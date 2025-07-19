# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.3] - 2025-07-19
### Fix
- Fix index overflow caused by non-ATCG characters.

### Added
- Include docstrings in the stub file.

## [1.0.2] - 2025-06-04
### Changed
- Use `saturating_sub` to compute `window_start`.
- Remove useless `PyResult` conversion from the `mask` method signature.
- Bump `pyo3` to `0.25.0` and `thiserror` to `2.0.12`.
- Set the `__version__` attribute in `lib.rs`.

## [1.0.1] - 2025-06-04
### Fix
- Treat `N` and `n` bases as `A` to prevent out of range masks.

### Added
- Added a home URL to `pyproject.toml`.

### Changed
- Added `Cargo.lock` to `.gitignore`.

## [1.0.0] - 2024-10-02
### Added
- First release.
