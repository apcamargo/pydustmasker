import math

import pytest

from pydustmasker import DustMasker, LongdustMasker, TantanMasker


def test_reject_non_ascii_sequences():
    for masker_type in (DustMasker, LongdustMasker, TantanMasker):
        with pytest.raises(ValueError):
            masker_type("ACGT🙂ACGT")


def test_dustmasker_creation():
    seq = "ATGC" * 10
    masker = DustMasker(seq)
    assert masker.sequence == seq
    assert masker.window_size == 64
    assert masker.score_threshold == 20
    assert isinstance(masker.intervals, tuple)


def test_dustmasker_repr():
    seq = "TACCCCCCCGCGTTTTTTT"
    masker = DustMasker(seq, window_size=64, score_threshold=20)
    assert (
        repr(masker)
        == "DustMasker(sequence: 'TACCCCCC…', intervals: ((2, 9), (12, 19)))"
    )


def test_dustmasker_slicing():
    seq = "TACCCCCCCGCGTTTTTTT"
    masker = DustMasker(seq, window_size=64, score_threshold=20)
    assert masker[0:1] == ((2, 9),)


def test_iterable():
    masker = DustMasker("TACCCCCCCGCGTTTTTTT", window_size=64, score_threshold=20)
    assert tuple(
        DustMasker("TACCCCCCCGCGTTTTTTT", window_size=64, score_threshold=20)
    ) == ((2, 9), (12, 19))
    assert tuple(iter(masker)) == masker.intervals


def test_dustmasker_masking():
    seq = "TACCCCCCCGCGTTTTTTT"
    masker = DustMasker(seq, window_size=64, score_threshold=20)
    assert masker.mask() == "TAcccccccGCGttttttt"
    assert masker.mask(hard=True) == "TANNNNNNNGCGNNNNNNN"


def test_dustmasker_window_size():
    seq = "TACCCCCCCGCGTTTTTTT"
    m1 = DustMasker(seq, window_size=64)
    m2 = DustMasker(seq, window_size=4)
    assert m1.window_size == 64
    assert m1.intervals == ((2, 9), (12, 19))
    assert m2.window_size == 4
    assert m2.intervals == ()


def test_dustmasker_score_threshold():
    seq = "TACCCCCCCGCGTTTTTTT"
    m1 = DustMasker(seq, score_threshold=20)
    m2 = DustMasker(seq, score_threshold=128)
    assert m1.score_threshold == 20
    assert m1.intervals == ((2, 9), (12, 19))
    assert m2.score_threshold == 128
    assert m2.intervals == ()


def test_ambigious():
    # No ambiguous
    seq1 = "GCCAGGCTGGCCAAGGAGATCttttttttttttttttttttttttAAGAGACCATGGCATGCACTGGCCAAGGAGATCttttttttttttttttttttttttAAGA"
    # With ambiguous
    seq2 = "GCCAGGCTGGCCAAGGAGATTCttttttttttttttttttttttttAAGAGCCARYCTGGCCAAGGAGANTCttttttttttttttttttttttttAAGA"
    # With ambiguous and masks
    seq3 = "GCCAGGCTGGCCAAGGAGATTCttttttttttttttttttttttttAFGAGCCAGGCTGGCCAAGGAGANTCtttttttttnNnttttttttAAGA"
    assert DustMasker(seq1, window_size=64).intervals == ((21, 45), (78, 102))
    assert DustMasker(seq2, window_size=64).intervals == ((22, 46), (72, 96))
    assert DustMasker(seq3, window_size=64).intervals == ((22, 46), (72, 81), (84, 92))


def test_longdustmasker_creation():
    seq = "ATGC" * 10
    masker = LongdustMasker(seq)
    assert masker.sequence == seq
    assert masker.window_size == 5000
    assert masker.score_threshold == 0.6
    assert masker.kmer == 7
    assert masker.gc is None
    assert masker.xdrop == 50
    assert masker.min_start_cnt == 3
    assert masker.approx is False
    assert masker.forward_only is False
    assert isinstance(masker.intervals, tuple)


def test_longdustmasker_repr():
    seq = "TACCCCCCCGCGTTTTTTT"
    masker = LongdustMasker(seq, window_size=64, score_threshold=0.1, kmer=3)
    assert repr(masker) == "LongdustMasker(sequence: 'TACCCCCC…', intervals: ((2, 19)))"


def test_longdustmasker_slicing():
    seq = "TACCCCCCCGCGTTTTTTT"
    masker = LongdustMasker(seq, window_size=64, score_threshold=0.1, kmer=3)
    assert masker[0:1] == ((2, 19),)


def test_longdustmasker_masking():
    seq = "TACCCCCCCGCGTTTTTTT"
    masker = LongdustMasker(seq, window_size=64, score_threshold=0.1, kmer=3)
    assert masker.mask() == "TAcccccccgcgttttttt"
    assert masker.mask(hard=True) == "TANNNNNNNNNNNNNNNNN"


def test_longdustmasker_window_size():
    seq = "GCTAGCAGTTCGAT" + "A" * 15 + "GCTAGCAGTTCGAT"
    m1 = LongdustMasker(seq)
    m2 = LongdustMasker(seq, window_size=8)
    assert m1.intervals == ((14, 29),)
    assert m2.intervals == ((12, 31),)


def test_longdustmasker_score_threshold():
    seq = "GCTAGCAGTTCGAT" + "A" * 20 + "GCTAGCAGTTCGAT"
    m1 = LongdustMasker(seq, score_threshold=0.6)
    m2 = LongdustMasker(seq, score_threshold=5.0)
    assert m1.intervals == ((14, 34),)
    assert m2.intervals == ()


def test_longdustmasker_validation_score_threshold():
    for score_threshold in (float("nan"), float("inf"), 0.0, -0.1):
        with pytest.raises(ValueError, match="invalid score threshold"):
            LongdustMasker("A" * 8, score_threshold=score_threshold)


def test_longdustmasker_kmer():
    seq = "TTTTGCGCACGTGTCGCTTGAATATATTTTTTTTT"
    m1 = LongdustMasker(seq, kmer=7, score_threshold=0.1, window_size=64)
    m2 = LongdustMasker(seq, kmer=3, score_threshold=0.1, window_size=64)
    assert m1.intervals == ((26, 35),)
    assert m2.intervals == ((0, 4), (21, 35))


def test_longdustmasker_validation_kmer_range():
    for kmer in (0, 13):
        with pytest.raises(ValueError, match="invalid k-mer size"):
            LongdustMasker("A" * 14, window_size=14, kmer=kmer)


def test_longdustmasker_validation_window_size_limit():
    assert (
        LongdustMasker("AA", window_size=65535, kmer=1, forward_only=True).intervals
        == ()
    )
    with pytest.raises(ValueError, match="invalid window size"):
        LongdustMasker("AA", window_size=65536, kmer=1)


def test_longdustmasker_xdrop():
    seq = "GTAGCGAT" + "A" * 20 + "GCTAGCAGTTCGATATAAGCT" + "A" * 20 + "GTAGCGAT"
    m1 = LongdustMasker(seq, xdrop=50)
    m2 = LongdustMasker(seq, xdrop=1)
    assert m1.intervals == ((7, 70),)
    assert m2.intervals == (
        (8, 28),
        (49, 69),
    )


def test_longdustmasker_min_start_cnt():
    seq = "GTCTTCTTCGTCTTCTTCGTCTTCTTCGTCTTCTTCATCTTCTTCGTCTTCTTCATCTTCTTCTTCTTCTTCTTCGTCTT"
    m1 = LongdustMasker(seq, min_start_cnt=2, score_threshold=0.1, window_size=64)
    m2 = LongdustMasker(seq, min_start_cnt=6, score_threshold=0.1, window_size=64)
    assert m1.intervals == ((0, 80),)
    assert m2.intervals == ((1, 80),)


def test_longdustmasker_approx():
    seq = "TAGTCATTTTTTTTGTCTTTGCGTGTTTGCAATTAATTGATCTTTTTTTTTAACCCGCCCCCCTTTATTTTGTCA"
    m1 = LongdustMasker(seq, approx=False, score_threshold=0.1, window_size=64)
    m2 = LongdustMasker(seq, approx=True, score_threshold=0.1, window_size=64)
    assert m1.intervals == (
        (6, 14),
        (42, 51),
    )
    assert m2.intervals == ((42, 51),)


def test_longdustmasker_gc_none_vs_float():
    seq = "CAAATATGGCTGACCGTCAGCAGATGATGGAAAAACGTATGGACATGATGCAATCCATGATGCAGATGATGATGGAC"
    m1 = LongdustMasker(seq, gc=None, score_threshold=0.1, window_size=64)
    m2 = LongdustMasker(seq, gc=0.1, score_threshold=0.1, window_size=64)
    assert m1.intervals == ((19, 76),)
    assert m2.intervals == ((64, 74),)


def test_longdustmasker_gc_none_vs_auto():
    seq = "CCATTGGATATAAATTCTCACTTCTGTTTTAGACATAAAATTATAATCAAAAGATTAATTATATTACTCAGTTCTTAAGAAGCAAAAGA"
    m1 = LongdustMasker(seq, gc=None, score_threshold=0.1, window_size=64, kmer=3)
    m2 = LongdustMasker(seq, gc="auto", score_threshold=0.1, window_size=64, kmer=3)
    assert m1.intervals == ((6, 89),)
    assert m2.intervals == ((7, 12), (34, 66), (76, 89))


def test_longdustmasker_forward_only():
    seq = "ACAGAAAAATGCGTACCCATCCACCTTTCAGTGCGTACCCACCCATCCACCTTTCAGTGCGTACCCATCCACCTTTCATTT"
    m1 = LongdustMasker(seq, forward_only=False, score_threshold=0.1, window_size=64)
    m2 = LongdustMasker(seq, forward_only=True, score_threshold=0.1, window_size=64)
    assert m1.intervals == ((4, 81),)
    assert m2.intervals == ((4, 78),)


def test_longdustmasker_ambigious():
    seq = "TACCNNNNCGCGTTTTTTT"
    masker = LongdustMasker(seq, window_size=64, score_threshold=0.1, kmer=3)
    assert masker.intervals == ((12, 19),)


def test_tantanmasker_creation():
    seq = "ATGCTAGCCGTAATGCGTACX"
    masker = TantanMasker(seq)
    assert masker.sequence == seq
    assert masker.protein is False
    assert masker.repeat_start == 0.005
    assert masker.repeat_end == 0.05
    assert masker.decay == 0.9
    assert masker.max_period == 100
    assert masker.gap_open == 0
    assert masker.gap_extend is None
    assert masker.score_threshold == 0.5
    assert masker.min_copy_number == 2.0
    assert isinstance(masker.intervals, tuple)
    assert isinstance(masker.probabilities, tuple)
    assert len(masker.probabilities) == len(seq)
    assert masker.intervals == ()
    assert masker.n_masked_bases == 0


def test_tantanmasker_protein_defaults():
    seq = "ACDEFGHIKLMNPQRSTVWY" * 5
    m1 = TantanMasker(seq)
    m2 = TantanMasker(seq, protein=True)
    assert m1.intervals == ()
    assert m2.protein is True
    assert m2.max_period == 50
    assert m2.intervals == ((20, 100),)


def test_tantanmasker_masking():
    seq = "ATTATTATTATTATT"
    masker = TantanMasker(seq)
    assert masker.intervals == ((3, 15),)
    assert masker.mask() == "ATTattattattatt"
    assert masker.mask(hard=True) == "ATTNNNNNNNNNNNN"


def test_tantanmasker_masking_protein():
    seq = "AC" * 7
    m1 = TantanMasker(seq)
    m2 = TantanMasker(seq, protein=True)
    assert m1.mask(hard=True) == "ACNNNNNNNNNNNN"
    assert m2.intervals == ((2, 14),)
    assert m2.mask() == "ACacacacacacac"
    assert m2.mask(hard=True) == "ACXXXXXXXXXXXX"


def test_tantanmasker_repeat_units():
    assert TantanMasker("ACACACACACACAC").repeat_units() == (("AC", 0, 14, 7.0),)
    assert TantanMasker("ACGTACGTACGTACGTACGT").repeat_units() == (("ACGT", 0, 20, 5.0),)
    assert TantanMasker("ATTATTATTATTATT").repeat_units() == (("ATT", 0, 15, 5.0),)


def test_tantanmasker_repeat_units_max_period_boundaries():
    # Period 1 is the valid minimum DP state space
    assert TantanMasker("A" * 20, max_period=1).repeat_units() == (("A", 0, 20, 20.0),)
    # Period 5 reaches scalar tails, one bulk group, and a bulk-plus-tail group
    assert TantanMasker("ACGTA" * 6, max_period=5).repeat_units() == (
        ("ACGTA", 0, 30, 6.0),
    )


def test_tantanmasker_repeat_units_with_ambiguous_letters():
    assert TantanMasker("ACGTNACGTNACGTNACGTNACGTNACGTN").repeat_units() == (
        ("ACGTN", 0, 29, 5.8),
    )
    assert TantanMasker("CATNCATNCATNCATNCATNCATNCATN").repeat_units() == (
        ("CATN", 0, 27, 6.75),
    )
    assert TantanMasker("ACGTRACGTRACGTRACGTRACGTRACGTR").repeat_units() == (
        ("ACGTN", 0, 29, 5.8),
    )
    assert TantanMasker(
        "ACDEFGHIKBACDEFGHIKBACDEFGHIKB", protein=True
    ).repeat_units() == (("ACDEFGHIKX", 0, 29, 2.9),)


def test_tantanmasker_probabilities_are_in_unit_range():
    dna = "ACGTACGTAGCTNNNNACACACACAC"
    assert all(0.0 <= p <= 1.0 for p in TantanMasker(dna).probabilities)
    assert all(
        0.0 <= p <= 1.0
        for p in TantanMasker("ACDEFGHIK" * 4, protein=True).probabilities
    )


def test_tantanmasker_score_threshold():
    seq = "ACGTACGTACGTACGTAAGT"
    m1 = TantanMasker(seq)
    m2 = TantanMasker(seq, score_threshold=0.0)
    assert m1.intervals == ((4, 19),)
    assert m2.intervals == ((0, 20),)


def test_tantanmasker_validation_max_period():
    # `repeat_offset_prob` passes max_period to `f64::powi`.
    for bad in (0, 2**31):
        with pytest.raises(ValueError, match="invalid max_period"):
            TantanMasker("ACGTACGTACGTACGT", max_period=bad)


def test_tantanmasker_validation_decay_non_normal():
    for bad in (float("nan"), float("inf"), 1e-310, 0.0, -1.0, 1.5):
        with pytest.raises(ValueError, match="invalid decay"):
            TantanMasker("ACACACACACACAC", decay=bad)


def test_tantanmasker_repeat_units_gapped():
    seq = "ACACACACACACCATCATCATCATCAT"
    m1 = TantanMasker(seq, min_copy_number=0)
    m2 = TantanMasker(seq, gap_open=7, gap_extend=1, min_copy_number=0)
    assert m1.repeat_units() == (("CAT", 9, 27, 6.0),)
    assert m2.repeat_units() == (("CAT", 0, 27, 10.666666666666666),)

    seq = "CATCATCATCATACACACACACACAC"
    m1 = TantanMasker(seq, min_copy_number=0)
    m2 = TantanMasker(seq, gap_open=7, gap_extend=1, min_copy_number=0)
    assert m1.repeat_units() == (("AC", 12, 26, 7.0),)
    assert m2.repeat_units() == (("AC", 0, 26, 11.0),)


def test_tantanmasker_repeat_units_gapped_min_copy_number():
    # Filtering uses the gapped copy-number calculation.
    seq = "CATCATCATCATACACACACACACAC"
    m1 = TantanMasker(seq, gap_open=7, gap_extend=1)
    m2 = TantanMasker(seq, gap_open=7, gap_extend=1, min_copy_number=12)
    assert m1.repeat_units() == (("AC", 0, 26, 11.0),)
    assert m2.repeat_units() == ()


def test_tantanmasker_repeat_units_min_copy_number():
    seq = "ACACACACACACAC"
    m1 = TantanMasker(seq)
    m2 = TantanMasker(seq, min_copy_number=10.0)
    assert m1.repeat_units() == (("AC", 0, 14, 7.0),)
    assert m2.repeat_units() == ()


def test_tantanmasker_repr():
    seq = "TACCCCCCCGCGTTTTTTT"
    masker = TantanMasker(seq)
    assert repr(masker) == "TantanMasker(sequence: 'TACCCCCC…', intervals: ())"


def test_tantanmasker_scalar_validation():
    for repeat_start in (1.0, -0.1):
        with pytest.raises(ValueError):
            TantanMasker("ACGTACGTACGTACGT", repeat_start=repeat_start)

    for repeat_end in (1.5, -0.1):
        with pytest.raises(ValueError):
            TantanMasker("ACGTACGTACGTACGT", repeat_end=repeat_end)

    for score_threshold in (2.0, -0.1):
        with pytest.raises(ValueError):
            TantanMasker("ACGTACGTACGTACGT", score_threshold=score_threshold)

    with pytest.raises(ValueError):
        TantanMasker("ACGTACGTACGTACGT", gap_extend=0)

    with pytest.raises(ValueError):
        TantanMasker("ACGTACGTACGTACGT", min_copy_number=-1.0)


def test_tantanmasker_validation_gap_probability():
    with pytest.raises(ValueError):
        TantanMasker("ACGTACGTACGTACGT", gap_open=0, gap_extend=1)


def test_tantanmasker_validation_empty_sequence():
    with pytest.raises(ValueError):
        TantanMasker("")


def test_tantanmasker_gapped_probabilities():
    seq = "ACGCGCGCGCGCAGCGCGCGCGCACGT"
    m1 = TantanMasker(seq)
    m2 = TantanMasker(seq, gap_open=0, gap_extend=2)
    assert m1.intervals == ((3, 24),)
    assert m2.intervals == ((4, 23),)
    oracle = [
        0,
        0.014,
        0.0897,
        0.394,
        0.559,
        0.724,
        0.806,
        0.863,
        0.891,
        0.908,
        0.916,
        0.921,
        0.923,
        0.938,
        0.946,
        0.949,
        0.949,
        0.946,
        0.94,
        0.926,
        0.898,
        0.838,
        0.712,
        0.443,
        0.358,
        0.197,
        0.0478,
    ]
    assert list(m2.probabilities) == pytest.approx(oracle, abs=1e-3)


def test_tantanmasker_gapped_small_decay():
    seq = "AC" * 10
    m1 = TantanMasker(seq, gap_open=7, gap_extend=1)
    m2 = TantanMasker(seq, gap_open=7, gap_extend=1, decay=1e-5)
    assert m1.intervals == ((2, 20),)
    assert all(math.isfinite(p) and 0.0 <= p <= 1.0 for p in m2.probabilities)
    assert max(m2.probabilities) > 0.5
    assert m2.intervals == ((2, 19),)
