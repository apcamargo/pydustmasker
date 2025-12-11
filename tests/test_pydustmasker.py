from pydustmasker import DustMasker, LongdustMasker
from pydustmasker._pydustmasker import _BaseMasker


def test_dust_masker_creation():
    masker = DustMasker("TACCCCCCCGCGTTTTTTT", window_size=64, score_threshold=20)
    assert masker.sequence == "TACCCCCCCGCGTTTTTTT"
    assert masker.window_size == 64
    assert masker.score_threshold == 20
    assert masker.intervals == [(2, 9), (12, 19)]
    assert isinstance(masker.intervals, list)


def test_n_masked_bases():
    masker = DustMasker("TACCCCCCCGCGTTTTTTT", window_size=64, score_threshold=20)
    assert masker.n_masked_bases == 14


def test_mask_soft():
    masker = DustMasker("TACCCCCCCGCGTTTTTTT", window_size=64, score_threshold=20)
    assert masker.mask(hard=False) == "TAcccccccGCGttttttt"


def test_mask_hard():
    masker = DustMasker("TACCCCCCCGCGTTTTTTT", window_size=64, score_threshold=20)
    assert masker.mask(hard=True) == "TANNNNNNNGCGNNNNNNN"


def test_repr():
    masker = DustMasker("GTACCCCCCCGTAACGTTTTT", window_size=64, score_threshold=20)
    assert "DustMasker(sequence: 'GTACCCCC…', intervals: [(3, 10)])" == repr(masker)


def test_window_size():
    masker = DustMasker("TACCCCCCCGCGTTTTTTT", window_size=4, score_threshold=20)
    assert masker.window_size == 4
    assert masker.intervals == []


def test_score_threshold():
    masker = DustMasker("TACCCCCCCGCGTTTTTTT", window_size=64, score_threshold=128)
    assert masker.score_threshold == 128
    assert masker.intervals == []


def test_ambigious():
    # no ambiguous
    seq1 = "GCCAGGCTGGCCAAGGAGATCttttttttttttttttttttttttAAGAGACCATGGCATGCACTGGCCAAGGAGATCttttttttttttttttttttttttAAGA"
    # with ambiguous
    seq2 = "GCCAGGCTGGCCAAGGAGATTCttttttttttttttttttttttttAAGAGCCARYCTGGCCAAGGAGANTCttttttttttttttttttttttttAAGA"
    # with ambiguous and masks
    seq3 = "GCCAGGCTGGCCAAGGAGATTCttttttttttttttttttttttttAFGAGCCAGGCTGGCCAAGGAGANTCtttttttttnNnttttttttAAGA"

    assert DustMasker(seq3, window_size=64).intervals == [(22, 46), (72, 81), (84, 92)]


def test_longdust_masker_basic():
    # Longdust defaults: window=5000, threshold=0.6, kmer=7
    # For a short sequence with simple repeats, we might need smaller window/kmer to trigger?
    # Actually let's just test API mechanics first.
    seq = "TACCCCCCCGCGTTTTTTT"
    masker = LongdustMasker(
        seq,
        window_size=10,
        score_threshold=20.0,
        kmer=3,
        xdrop=15,
        min_start_cnt=2,
        approx=True,
        gc=0.4,
        forward_only=True,
    )
    assert masker.sequence == seq
    assert masker.window_size == 10
    assert masker.score_threshold == 20.0
    assert masker.kmer == 3
    assert masker.xdrop == 15
    assert masker.min_start_cnt == 2
    assert masker.approx is True
    assert masker.forward_only is True
    assert isinstance(masker.intervals, list)


def test_longdust_gc_auto():
    # Test that gc="auto" works without errors
    seq = "GGCCGGCC"
    masker = LongdustMasker(seq, gc="auto")
    assert isinstance(masker.intervals, list)

    seq2 = "AAAAAAGGGGGG"
    masker2 = LongdustMasker(seq2, gc="auto")
    assert isinstance(masker2.intervals, list)


def test_longdust_gc_none():
    # Test that gc=None (uniform) works
    seq = "ACGTACGT"
    masker = LongdustMasker(seq)
    assert isinstance(masker.intervals, list)


def test_longdust_mask_hard():
    seq = "AAAAAAAA"
    # Create valid masker (params don't matter as much as checking mask() runs)
    masker = LongdustMasker(seq)
    masked = masker.mask(hard=True)
    assert len(masked) == len(seq)
    # If intervals found, should have Ns. If not, should be identical.
    if masker.intervals:
        assert "N" in masked
    else:
        assert masked == seq


def test_longdust_validation():
    seq = "ACGTACGT"

    # Test Invalid GC
    try:
        LongdustMasker(seq, gc=1.1)
        assert False, "Should raise ValueError for GC > 1.0"
    except ValueError:
        pass

    try:
        LongdustMasker(seq, gc=-0.1)
        assert False, "Should raise ValueError for GC < 0.0"
    except ValueError:
        pass

    # Test Invalid K-mer
    try:
        LongdustMasker(seq, kmer=0)
        assert False, "Should raise ValueError for kmer=0"
    except ValueError:
        pass

    # Test Invalid Score
    try:
        LongdustMasker(seq, score_threshold=0.0)
        assert False, "Should raise ValueError for score_threshold=0.0"
    except ValueError:
        pass

    # Test Invalid Auto GC (should be valid string)
    try:
        LongdustMasker(seq, gc="invalid")
        assert False, "Should raise ValueError for gc='invalid'"
    except ValueError:
        pass

    # Test Invalid min_start_cnt
    try:
        LongdustMasker(seq, min_start_cnt=0)
        assert False, "Should raise ValueError for min_start_cnt=0"
    except ValueError:
        pass

    try:
        LongdustMasker(seq, min_start_cnt=1)
        assert False, "Should raise ValueError for min_start_cnt=1"
    except ValueError:
        pass


def test_inheritance():
    d_masker = DustMasker("ACGT", window_size=64, score_threshold=20)
    ld_masker = LongdustMasker("ACGTACGT", window_size=10, score_threshold=0.6)

    assert isinstance(d_masker, _BaseMasker)
    assert isinstance(ld_masker, _BaseMasker)
    assert issubclass(DustMasker, _BaseMasker)
    assert issubclass(LongdustMasker, _BaseMasker)
