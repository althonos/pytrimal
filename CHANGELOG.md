# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/althonos/pytrimal/compare/v0.6.0...HEAD


## [v0.6.0] - 2022-10-17
[v0.6.0]: https://github.com/althonos/pytrimal/compare/v0.5.5...v0.6.0

### Fixed
- Overflow in `calculateSpuriousVector` for certain sequence block sizes.
- MMX backend never being used unless explicitly required even when being the best supported backend.
### Changed
- Use `archspec` Python package instead of `cpu-features` library for CPU feature detection at runtime.
- Use `importlib.resources.files` to load package data in `pytrimal.tests`.


## [v0.5.5] - 2022-10-17
[v0.5.5]: https://github.com/althonos/pytrimal/compare/v0.5.4...v0.5.5

### Fixed
- `calculateSpuriousVector` method of SIMD implementations not being declared `override`.
### Changed
- Replaced `aligned_alloc` with `posix_memalign` for compatibility with MacOS.


## [v0.5.4] - 2022-10-15
[v0.5.4]: https://github.com/althonos/pytrimal/compare/v0.5.3...v0.5.4

### Added
- AVX2 and MMX implementations of the SIMD statistics computation.
- Tests for all SIMD implementations supported on the local machine.
### Changed
- Refactor SIMD code using C++ templates and generic implementation.
### Fixed
- Broken rendering of function signatures in Sphinx documentation.
- `residues_mask` and `sequences_mask` attributes of `TrimmedAlignment` not being documented ([#1](https://github.com/althonos/pytrimal/issues/1)).


## [v0.5.3] - 2022-10-04
[v0.5.3]: https://github.com/althonos/pytrimal/compare/v0.5.2...v0.5.3

### Fixed
- `SimilarityMatrix.nt` inverting the `degenerated` argument value.


## [v0.5.2] - 2022-09-30
[v0.5.2]: https://github.com/althonos/pytrimal/compare/v0.5.1...v0.5.2

### Changed
- Replace NEON horizontal sums with implementations using `vaddvq` on Aarch64 and `vpaddl` on Armv7.
- Remove one layer of table lookup in all `Similarity::CalculateVectors` implementations.
- Make all SIMD code use local buffers and deallocate early.

### Fixed
- Invalid operator being used in Cython code to deallocate C++ arrays.

### Added
- SSE2 and NEON implementations for the `Gaps` statistic.


## [v0.5.1] - 2022-09-05
[v0.5.1]: https://github.com/althonos/pytrimal/compare/v0.5.0...v0.5.1

### Fixed
- Build of `cpu_features` for platforms without hardware detection support.


## [v0.5.0] - 2022-09-05
[v0.5.0]: https://github.com/althonos/pytrimal/compare/v0.4.0...v0.5.0

### Added
- `pytrimal.RepresentativeTrimmer` class to trim by maximum identity or fixed number of cluster representatives.
- `pickle` protocol support for all trimmer classes.
- Conversion methods to convert an `Alignment` from and to Biopython or PyHMMER objects.
- Arm NEON implementation of the statistics computation algorithm, with speed-up similar to that of the SSE implementation.

### Fixed
- `std::streambuf` implementation based on the `readinto` Python method not working on Arm because of `char` being used to read ASCII.

### Removed
- Support for Python 3.5, due to Cython compatibility issues.


## [v0.4.0] - 2022-08-14
[v0.4.0]: https://github.com/althonos/pytrimal/compare/v0.3.0...v0.4.0

### Added
- `BaseTrimmer.backend` property to get the backend used by a trimmer object.
- Zero-copy slicing for `AlignmentSequences` and `AlignmentResidues` objects.
- `noduplicateseqs` method for `AutomaticTrimmer` objects.
- `OverlapTrimmer` class to perform overlap trimming with SSE-accelerated implementation.
- `AutomaticTrimmer.METHODS` attribute to expose all supported automatic trimming methods.
- `__repr__` implementation to all trimmer classes.

### Fixed
- Missing deallocation code for standalone `AlignmentResidues` objects.
- `Alignment.load` not working properly in PyPy environments.
- `Alignment` constructor sometimes crashing when not given any sequence.

### Changed
- Use aligned memory for some temporary buffers used in SIMD code.
- Enable loop unrolling when supported by the compiler.
- Skip letter validation when creating an `Alignment` object with sequences from an `AlignmentSequences` object.

### Removed
- `consistency_threshold` and `consistency_window` arguments of `ManualTrimmer`.


## [v0.3.0] - 2022-06-26
[v0.3.0]: https://github.com/althonos/pytrimal/compare/v0.2.2...v0.3.0

### Added
- Support for loading an `Alignment` from a file-like object for certain formats.
- Generic optimized backend using caching optimizations from [inab/trimal#66](https://github.com/inab/trimal/pull/66).

### Fixed
- Compilation of code for OSX platforms in Python 3.10.
- File not being closed on error when loading a FASTA alignment.

### Changed
- Add tests for loading an `Alignment` without requiring `importlib.resources`.


## [v0.2.2] - 2022-06-08
[v0.2.2]: https://github.com/althonos/pytrimal/compare/v0.2.1...v0.2.2

### Added
- Keyword arguments to specify the half-window sizes in manual trimmer.
- `Alignment.dump` and `Alignment.dumps` function to write an alignment to a file, file-like object, or string.
- Optimized implementation of `Similarity::calculateVectors`.

### Changed
- Use faster implementation of SSE2 horizontal sum based on `_mm_sad_epu8`.


## [v0.2.1] - 2022-06-06
[v0.2.1]: https://github.com/althonos/pytrimal/compare/v0.2.0...v0.2.1

### Fixed
- Missing SSE2 files in source distribution.


## [v0.2.0] - 2022-06-06
[v0.2.0]: https://github.com/althonos/pytrimal/compare/v0.1.2...v0.2.0

### Added
- Vendored `cpu_features` library to perform runtime detection of CPU features.
- SIMD implementation of the *similarity* statistic code with SSE2 instructions.

### Fixed
- Compilation on platforms without OpenMP by adding an empty `omp.h` header file.


## [v0.1.2] - 2022-06-04
[v0.1.2]: https://github.com/althonos/pytrimal/compare/v0.1.1...v0.1.2

### Added
- Python constructor and buffer protocol support for `SimilarityMatrix`.
- `SimilarityMatrix.similarity` method to get the similarity between two characters instead.

### Fixed
- Source compilation failing because of source files in the `pytrimal` folder.


## [v0.1.1] - 2022-06-03
[v0.1.1]: https://github.com/althonos/pytrimal/compare/v0.1.0...v0.1.1

### Added
- Type annotations for all classes of the `pytrimal` extension module.

### Fixed
- Cython header files not being included in source distribution.


## [v0.1.0] - 2022-06-02
[v0.1.0]: https://github.com/althonos/pytrimal/compare/d774f73...v0.1.0

Initial release.
