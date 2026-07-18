---
icon: lucide/book-marked
---

# Theory

## Low-complexity sequences

Some genomic regions contain contiguous stretches of nucleotides with simple repetitive patterns or strong compositional biases. These low-complexity sequences range from short homopolymeric runs (e.g., AAAAAAAAA) and short tandem repeats (e.g., GATGATGAT) to large-scale structural elements such as centromeric satellite DNA and tandem repeats with longer repeat units.

Low-complexity sequences can complicate analyses in various ways. For example, they can generate spurious matches in sequence searches, obscuring biologically meaningful similarities. To mitigate this, it is common practice to identify and mask low-complexity regions prior to downstream analyses using dedicated algorithms that detect these regions based on characteristic features, such as increased repetition of sequence motifs relative to random sequences.

## Symmetric DUST (SDUST) { #sdust }

The SDUST algorithm[^1], implemented in the DustMasker tool included with NCBI’s BLAST, measures sequence complexity by assessing how frequently nucleotide 3-mers are repeated within a given sequence interval $x$. The complexity score, $S_{\text{SDUST}}(x)$, is calculated as:

$$
S_{\text{SDUST}}(x) = \frac{\displaystyle\sum\nolimits_{t \in R} c_t(x)(c_t(x)-1)}{2(\ell(x)-1)}
$$

In this formula, $R$ is the set of all 64 possible 3-mers, $c_t(x)$ is the frequency of 3-mer $t$ in $x$, and $\ell(x)$ is the total number of 3-mers in the candidate interval, where $\ell(x)=|x|-2$. For a 3-base interval, SDUST defines the score as zero.

SDUST identifies perfect low-complexity intervals, which are subsequences of length at most $W$ whose score $S_{\text{SDUST}}(x)$ exceeds a threshold $T$ and is not lower than that of any of their subsequences. During execution, SDUST moves a window of at most $W$ bases (64 by default) along the sequence and maintains the perfect intervals within the current window. When an interval leaves the window, its bases are finalized for masking; newly formed suffixes are then evaluated for additional perfect intervals. Because perfectness depends only on an interval and its substrings, the result is independent of flanking sequence. The score is also invariant under reverse complementation, so the masked intervals are strand symmetric.

## Longdust

While SDUST effectively identifies short-range low-complexity regions within DNA sequences, it is unsuited to identify satellite or tandem repeats with long repeat units. The algorithm becomes computationally prohibitive with large window sizes, and its scoring function exhibits a length bias that disproportionately classifies longer sequences as low-complexity. Moreover, because SDUST uses a fixed 3-mer size, it cannot adequately characterize repeats with longer motifs.

To overcome these limitations, Longdust[^2] employs a statistical model of k-mer count distributions, enabling efficient analysis within long genomic windows. The algorithm computes the score $S_{\text{Longdust}}(x)$ as follows:

$$
S_{\text{Longdust}}(x) = \sum\nolimits_{t \in R} \log(c_t(x)!) - f(\ell(x)) - T \cdot \ell(x)
$$

Here, $R$ is the set of all $4^k$ possible k-mers ($k=7$ by default), $c_t(x)$ is the count of k-mer $t$, and $\ell(x)$ is the total number of k-mers in the string. The scaling function $f(\ell(x))$ is derived from the expected k-mer-count distribution under a random-sequence background, optionally adjusted for GC content. It keeps random sequences near zero across lengths, avoiding the length bias of the SDUST score. The threshold $T$ (0.6 by default) controls the stringency of low-complexity classification.

Rather than reporting perfect intervals, Longdust identifies good low-complexity intervals, which are regions with a positive score for which no prefix or suffix has a higher score. For each end position ($j$), it scans backward to collect candidate starts, then scans forward from each candidate to determine whether $j$ is the best-scoring endpoint. Both passes are restricted to a fixed window ($w$; 5,000 bp by default), which keeps the search local and efficient. Candidate pruning and other heuristics make the algorithm inexact; by default, Longdust runs on both strands and merges the intervals to produce strand-symmetric output. Its X-drop heuristic also limits extensions that fall too far below their best score.

## tantan

The tantan algorithm[^3] identifies tandem repeats and other low-complexity regions in nucleotide or protein sequences with a hidden Markov model (HMM). Unlike SDUST and Longdust, which score candidate intervals from aggregate k-mer counts, tantan measures a sequence's self-similarity at several offsets. For an offset $p$, it compares $x_i$ with $x_{i-p}$, allowing it to detect inexact tandem patterns.

The HMM has a background state and one repeat state for each offset $p=1,\ldots,w$, where $w$ is the maximum offset. In repeat state $p$, the substitution score $S(x_i, x_{i-p})$ is converted to repeat evidence:

$$
L_{i,p}=\exp\!\bigl(\lambda S(x_i, x_{i-p})\bigr),
$$

Here, $\lambda$ is a scaling constant. Favorable substitutions provide stronger evidence than unfavorable ones, allowing imperfect repeats to be detected. A repeat starts with total probability $r$, divided among offsets with geometrically decaying weights ($s_{p+1}=d s_p$), favoring shorter offsets. Once entered, a repeat state can persist along a region or return to the background with probability $e$.

Forward-backward decoding combines this evidence with the transitions at all offsets and neighbouring positions. It calculates the posterior probability that position $i$ belongs to a non-background state:

$$
q_i=1-P(B_i\mid x),
$$

This represents the total posterior support across various offsets, rather than the score for one selected offset or a direct sum of substitution scores. Positions with $q_i$ at or above the chosen threshold are masked.

[^1]: Morgulis, Aleksandr, *et al*. **A Fast and Symmetric Dust Implementation to Mask Low-Complexity DNA Sequences**. *Journal of Computational Biology*, vol. 13, no. 5, June 2006, pp. 1028–40. <https://doi.org/10.1089/cmb.2006.13.1028>.

[^2]: Li, Heng, and Brian Li. **Finding Low-Complexity DNA Sequences with Longdust**. *Bioinformatics*, vol. 42, no. 3, Feb. 2026, p. btag112. <https://doi.org/10.1093/bioinformatics/btag112>.

[^3]: Frith, Martin C. **A New Repeat-Masking Method Enables Specific Detection of Homologous Sequences**. *Nucleic Acids Research*, vol. 39, no. 4, Mar. 2011, pp. e23–e23. <https://doi.org/10.1093/nar/gkq1212>.
