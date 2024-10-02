from __future__ import annotations

from collections.abc import Sequence

class DustMasker:
    sequence: str
    window_size: int
    score_threshold: int
    intervals: Sequence[tuple[int, int]]
    def __init__(
        self, sequence: str, window_size: int, score_threshold: int
    ) -> None: ...
    @property
    def n_masked_bases(self) -> int: ...
    def mask(self, hard: bool) -> str: ...
    def __repr__(self) -> str: ...
