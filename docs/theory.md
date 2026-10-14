---
icon: lucide/book-marked
---

# Theory

## Low-complexity sequences

DNA and protein sequences can contain contiguous stretches with simple repetitive patterns or strong compositional biases. These low-complexity sequences range from short homopolymeric runs (e.g., AAAAAAAAA) and short tandem repeats (e.g., GATGATGAT) to large-scale structural elements such as centromeric satellite DNA and tandem repeats with longer repeat units.

Low-complexity sequences can complicate many analyses. In sequence searches, for example, they generate spurious matches that obscure biologically meaningful similarities. To avoid this, low-complexity regions are commonly identified and masked before downstream analyses.

## Symmetric DUST (SDUST) { #sdust }

The SDUST algorithm[^1], implemented in the DustMasker tool included with NCBI's BLAST, measures the complexity of a sequence interval $x$ by how often its 3-mers repeat. The score $S_{\text{SDUST}}(x)$ is:

$$
S_{\text{SDUST}}(x) = \frac{\displaystyle\sum\nolimits_{t \in R} c_t(x)(c_t(x)-1)}{2(\ell(x)-1)}
$$

Here, $R$ is the set of all 64 possible 3-mers, $c_t(x)$ is the number of times 3-mer $t$ occurs in $x$, and $\ell(x) = |x| - 2$ is the total number of 3-mers in $x$. The numerator counts pairs of identical 3-mers, so intervals dominated by a few repeated 3-mers score high. A 3-base interval contains a single 3-mer and makes the denominator zero, so its score is defined as zero.

An interval is low-complexity when its score exceeds a threshold $T$. SDUST reports perfect intervals: substrings of length at most $W$ (64 bases by default) whose score exceeds $T$ and is not lower than the score of any of their own substrings. Every base covered by a perfect interval is masked. As in NCBI's tool, the `score_threshold` parameter of `DustMasker` is $10T$, so its default of 20 corresponds to $T = 2$.

To find perfect intervals, SDUST slides a window of $W$ bases along the sequence. Each time the window advances by one base, SDUST evaluates the intervals that end at the new base and keeps those that are perfect. Once the start of an interval leaves the window, the interval is final and its bases are masked. Because perfectness depends only on an interval's own sequence, a low-complexity interval is masked regardless of its surroundings. Reverse complementation only relabels the 3-mers, so the score and the masked intervals are strand symmetric.

## Longdust

SDUST finds short low-complexity regions but misses satellites and tandem repeats with long repeat units. A larger window is not a simple fix, because the cost of finding perfect intervals grows quickly with the window size and the SDUST score is biased by length. Even in random sequence, longer intervals accumulate more repeated 3-mers and score higher, so a large window would mask random sequence.

Longdust[^2] scores intervals with a statistical model of k-mer counts that removes this length bias. The score of an interval $x$ is:

$$
S_{\text{Longdust}}(x) = \sum\nolimits_{t \in R} \log(c_t(x)!) - f(\ell(x)) - T \cdot \ell(x)
$$

Here, $R$ is the set of all $4^k$ possible k-mers ($k = 7$ by default), $c_t(x)$ is the number of times k-mer $t$ occurs in $x$, and $\ell(x)$ is the total number of k-mers in $x$. The first term grows as k-mers repeat. The scaling function $f(\ell(x))$ is the expected value of the first term for a random sequence with $\ell(x)$ k-mers, optionally adjusted for GC content. Subtracting it brings the first two terms close to zero for random sequence of any length. The threshold $T$ (0.6 by default) then sets the stringency. An interval has a positive score only when the first two terms, divided by $\ell(x)$, exceed $T$.

To make the search fast, Longdust replaces perfect intervals with good low-complexity intervals: intervals with a positive score that no prefix or suffix outscores. Checking only prefixes and suffixes, rather than every substring, is much faster. All bases covered by good intervals are masked.

For each possible end position, Longdust scans backward to find starts that satisfy the suffix condition, then scans forward from each one to check the prefix condition. Both scans stay within a 5,000-base window by default. Heuristics such as an X-drop rule, which stops the backward scan once the score falls far below its best, speed up the search but can miss some good intervals. Because these heuristics are not strand symmetric, Longdust runs on both strands by default and merges the results.

## tantan

SDUST and Longdust score an interval by counting its k-mers. The tantan algorithm[^3] instead measures how similar a sequence $x$ is to itself at different offsets. For an offset $p$, it compares each letter $x_i$ with the letter $p$ positions earlier, $x_{i-p}$. In a tandem repeat whose repeat unit has length $p$, these letters tend to match.

tantan models the sequence with a hidden Markov model (HMM). The model has one background state and one repeat state for each offset $p$, up to a maximum offset (100 for DNA and 50 for proteins by default). In repeat state $p$, the substitution score $\sigma(x_i, x_{i-p})$ is converted into repeat evidence:

$$
L_{i,p} = \exp\!\bigl(\lambda\, \sigma(x_i, x_{i-p})\bigr)
$$

Here, $\lambda$ is a constant derived from the scoring matrix (+1 for a match and −1 for a mismatch in DNA, BLOSUM62 for proteins). It makes $L_{i,p}$ a likelihood ratio. Values above 1 favor the repeat state over the background, and values below 1 favor the background. Because similar letters still count as evidence, tantan can detect imperfect repeats.

Transitions between states control when repeats begin and end. At each position, the model enters a repeat with a small probability. This probability is split among offsets so that each offset is a fixed fraction as likely as the one before it, which favors short repeat units. A repeat state returns to the background with a fixed probability at each position and otherwise continues at the same offset. Optional gap states (off by default) also let repeat copies differ by short insertions and deletions.

Forward-backward decoding combines the evidence $L_{i,p}$ with these transition probabilities, summing over all possible state paths through the sequence. The result is the posterior probability $q_i$ that position $i$ is not in the background state:

$$
q_i = 1 - P(\text{position } i \text{ is in the background state} \mid x)
$$

Because $q_i$ adds up the support from every offset, it can flag regions that are weakly self-similar at many offsets yet lack a single well-defined repeat unit. Positions with $q_i$ at or above a threshold (0.5 by default) are masked.

[^1]: Morgulis, Aleksandr, *et al*. **A Fast and Symmetric Dust Implementation to Mask Low-Complexity DNA Sequences**. *Journal of Computational Biology*, vol. 13, no. 5, June 2006, pp. 1028–40. <https://doi.org/10.1089/cmb.2006.13.1028>.

[^2]: Li, Heng, and Brian Li. **Finding Low-Complexity DNA Sequences with Longdust**. *Bioinformatics*, vol. 42, no. 3, Feb. 2026, p. btag112. <https://doi.org/10.1093/bioinformatics/btag112>.

[^3]: Frith, Martin C. **A New Repeat-Masking Method Enables Specific Detection of Homologous Sequences**. *Nucleic Acids Research*, vol. 39, no. 4, Mar. 2011, pp. e23–e23. <https://doi.org/10.1093/nar/gkq1212>.
