# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/althonos/pytrimal/compare/v0.1.2...HEAD


## [v0.1.2] - 2022-06-04
[v0.1.1]: https://github.com/althonos/pytrimal/compare/v0.1.1...v0.1.2

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
