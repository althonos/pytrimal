# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/althonos/pytrimal/compare/v0.3.0...HEAD


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
