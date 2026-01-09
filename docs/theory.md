---
icon: lucide/book-marked
---

# Theory

## Low-complexity sequences

Some genomic regions contain contiguous stretches of nucleotides with simple repetitive patterns and strong compositional biases. These *low-complexity sequences* range from short homopolymeric runs (e.g., `AAAAAAAAA`) and short tandem repeats (e.g., `GATGATGAT`), which typically arise through DNA polymerase slippage, to large-scale structural elements such as centromeric satellite DNA and tandem repeats with longer repeat units.

Low-complexity sequences can interfere with computational analyses in various ways. For example, they can generate spurious matches in sequence searches, obscuring biologically meaningful similarities. To mitigate this, it is common practice to identify and mask low-complexity regions prior to downstream analyses using dedicated algorithms that detect these regions based on characteristic features, such as reduced information content relative to random sequences.

## Symmetric DUST (SDUST) { #sdust }

The SDUST algorithm[^1], implemented in the DustMasker tool included with NCBI’s BLAST, measures sequence complexity by assessing how frequently nucleotide 3-mers are repeated within a given sequence interval $x$. The complexity score, $S_{\text{SDUST}}(x)$, is calculated as:

$$
S_{\text{SDUST}}(x) = \frac{\displaystyle\sum\nolimits_{t \in R} c_t(x)(c_t(x)-1)}{2(\ell(x)-1)}
$$

In this formula, $R$ is the set of all 64 possible 3-mers, $c_t(x)$ is the frequency of 3-mer $t$ in $x$, and $\ell(x)$ is the total number of 3-mers in the candidate interval, where $\ell(x)=|x|-2$.

Low-complexity sequences are identified as **perfect low-complexity intervals**, defined as subsequences of length at most $W$ whose score $S_{\text{SDUST}}(x)$ is greater than a threshold $T$, and for which no subsequence of $x$ has a higher score than $x$ itself. During execution, SDUST applies a sliding window of fixed size $W$ (64 bp by default) that moves along the sequence while maintaining a dynamic list of all perfect intervals found entirely within its boundaries. At each step, intervals that are no longer fully contained are finalized for masking, while newly created suffixes are evaluated to identify any new perfect intervals. This method guarantees context-independence, as the decision to mask a region is based only on comparisons within that local window, isolated from flanking sequences.

## Longdust

While SDUST effectively identifies short-range low-complexity regions within DNA sequences, it is unsuited to identify satellite or tandem repeats with long repeat units. The algorithm becomes computationally prohibitive with large window sizes, and its scoring function exhibits a length bias that disproportionately classifies longer sequences as low-complexity. Moreover, because SDUST uses a fixed 3-mer size, it cannot adequately characterize repeats with longer motifs.

To overcome these limitations, Longdust[^2] employs a statistical model of k-mer count distributions, enabling efficient analysis within long genomic windows. The algorithm computes the score $S_{\text{Longdust}}(x)$ as follows:

$$
S_{\text{Longdust}}(x) = \sum\nolimits_{t \in R} \log(c_t(x)!) - f(\ell(x)) - T \cdot \ell(x)
$$

Where $R$ is the set of all $4^k$ possible k-mers ($k=7$ by default), $c_t(x)$ is the count of k-mer $t$, $\ell(x)$ is the total number of k-mers in the string, and $f(\ell(x))$ is a scaling function derived from the expected distribution of k-mer counts in random sequences. This scaling function ensures that random sequences score near zero regardless of their length, thereby avoiding the length bias inherent in the SDUST scoring function. The threshold $T$ (0.6 by default) controls the stringency of low-complexity classification.

Rather than reporting perfect intervals, Longdust identifies **good low-complexity intervals**, defined as regions with a positive score for which no sub-interval sharing the same start or end position has a higher score. To locate such intervals, for each end position ($j$) the algorithm first scans backward to collect candidate start positions, that is, potential left endpoints that could pair with ($j$) to form an interval. It then scans forward from each candidate to verify whether the resulting interval is a valid good low-complexity interval. For efficiency and to keep the assessment local, both the backward and forward scans are restricted to a fixed window ($w$) (5,000 bp by default), which limits how far upstream candidate start positions are considered. By combining the bounded backward collection of candidates with forward verification, the algorithm delimits low-complexity regions precisely while controlling runtime.

[^1]: Morgulis, Aleksandr, et al. **A Fast and Symmetric DUST Implementation to Mask Low-Complexity DNA Sequences**. *Journal of Computational Biology*, vol. 13, no. 5, June 2006, pp. 1028–40. <https://doi.org/10.1089/cmb.2006.13.1028>.

[^2]: Li, Heng, and Brian Li. **Finding Low-Complexity DNA Sequences with Longdust**. *arXiv*, 2025. <https://doi.org/10.48550/arxiv.2509.07357>.
