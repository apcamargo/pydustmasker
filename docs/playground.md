---
icon: lucide/gamepad-2
---

# Playground

Use the code editor below to try `pydustmasker` directly in your browser. You can edit the input sequence or the parameters and run it again to see how the results change.

/// pyodide-code-block

import pydustmasker

seq = (
    "GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGC"
    "AGCAGCAGCAGCAGCAGCAGCAGCTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACT"
    "TGCTGTGTCCACCCCATCGGACTTTAATTTTTATTAAAATAATAAAAATAAATATTTTAAT"
    "TTTTAAATTTAAATAATATTAGCGCCCTCCTGAAGTGCGTGGACACTCGCTATGAATCTCT"
)

maskers = {
    "SDUST": pydustmasker.DustMasker(seq),
    "Longdust": pydustmasker.LongdustMasker(seq),
    "tantan": pydustmasker.TantanMasker(seq),
}
for name, masker in maskers.items():
    print(f"{name}: {masker.intervals}")
    print(masker.mask())
    print()
///
