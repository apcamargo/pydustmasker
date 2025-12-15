import pytest

from pydustmasker import DustMasker


def test_dust_masker_creation():
    masker = DustMasker("TACCCCCCCGCGTTTTTTT", window_size=64, score_threshold=20)
    assert masker.sequence == "TACCCCCCCGCGTTTTTTT"
    assert masker.window_size == 64
    assert masker.score_threshold == 20
    assert masker.intervals == ((2, 9), (12, 19))
    assert isinstance(masker.intervals, tuple)
    assert len(masker) == 2


def test_slicing():
    masker = DustMasker("TACCCCCCCGCGTTTTTTT", window_size=64, score_threshold=20)
    assert masker[0] == (2, 9)
    assert masker[0:1] == ((2, 9),)
    assert masker[:2] == ((2, 9), (12, 19))


def test_iterable():
    masker = DustMasker("TACCCCCCCGCGTTTTTTT", window_size=64, score_threshold=20)
    assert hasattr(iter(masker), "__next__")
    assert tuple(
        DustMasker("TACCCCCCCGCGTTTTTTT", window_size=64, score_threshold=20)
    ) == ((2, 9), (12, 19))
    assert tuple(iter(masker)) == masker.intervals


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
    assert "DustMasker(sequence: 'GTACCCCC…', intervals: ((3, 10)))" == repr(masker)


def test_window_size():
    masker = DustMasker("TACCCCCCCGCGTTTTTTT", window_size=4, score_threshold=20)
    assert masker.window_size == 4
    assert masker.intervals == ()


def test_score_threshold():
    masker = DustMasker("TACCCCCCCGCGTTTTTTT", window_size=64, score_threshold=128)
    assert masker.score_threshold == 128
    assert masker.intervals == ()


def test_ambigious():
    # no ambiguous
    seq1 = "GCCAGGCTGGCCAAGGAGATCttttttttttttttttttttttttAAGAGACCATGGCATGCACTGGCCAAGGAGATCttttttttttttttttttttttttAAGA"
    assert DustMasker(seq1, window_size=64).intervals == ((21, 45), (78, 102))
    # with ambiguous
    seq2 = "GCCAGGCTGGCCAAGGAGATTCttttttttttttttttttttttttAAGAGCCARYCTGGCCAAGGAGANTCttttttttttttttttttttttttAAGA"
    assert DustMasker(seq2, window_size=64).intervals == ((22, 46), (72, 96))
    # with ambiguous and masks
    seq3 = "GCCAGGCTGGCCAAGGAGATTCttttttttttttttttttttttttAFGAGCCAGGCTGGCCAAGGAGANTCtttttttttnNnttttttttAAGA"
    assert DustMasker(seq3, window_size=64).intervals == ((22, 46), (72, 81), (84, 92))


def test_errors_creation():
    # sequence too short -> ValueError
    with pytest.raises(ValueError):
        DustMasker("AAA", window_size=64)
    # window_size too small -> ValueError
    with pytest.raises(ValueError):
        DustMasker("ACGTACGT", window_size=3)
    # negative integers for size/threshold should raise OverflowError when converting to unsigned types
    with pytest.raises(OverflowError):
        DustMasker("ACGTACGT", window_size=-1)
    with pytest.raises(OverflowError):
        DustMasker("ACGTACGT", window_size=64, score_threshold=-5)
    # wrong types should raise TypeError
    with pytest.raises(TypeError):
        DustMasker(12345, window_size=64)
    with pytest.raises(TypeError):
        DustMasker("ACGTACGT", window_size="foo")
