from pydustmasker import DustMasker, LongdustMasker


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
    assert hasattr(iter(masker), "__next__")
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
    # no ambiguous
    seq1 = "GCCAGGCTGGCCAAGGAGATCttttttttttttttttttttttttAAGAGACCATGGCATGCACTGGCCAAGGAGATCttttttttttttttttttttttttAAGA"
    # with ambiguous
    seq2 = "GCCAGGCTGGCCAAGGAGATTCttttttttttttttttttttttttAAGAGCCARYCTGGCCAAGGAGANTCttttttttttttttttttttttttAAGA"
    # with ambiguous and masks
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


def test_longdustmasker_kmer():
    seq = "TTTTGCGCACGTGTCGCTTGAATATATTTTTTTTT"
    m1 = LongdustMasker(seq, kmer=7, score_threshold=0.1, window_size=64)
    m2 = LongdustMasker(seq, kmer=3, score_threshold=0.1, window_size=64)
    assert m1.intervals == ((26, 35),)
    assert m2.intervals == ((0, 4), (21, 35))


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
