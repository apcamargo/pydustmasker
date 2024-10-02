from pydustmasker import DustMasker


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
