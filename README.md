# 🐍✂️ PytrimAl [![Stars](https://img.shields.io/github/stars/althonos/pytrimal.svg?style=social&maxAge=3600&label=Star)](https://github.com/althonos/pytrimal/stargazers)

*[Cython](https://cython.org/) bindings and Python interface to [trimAl](http://trimal.cgenomics.org/), a tool for automated alignment trimming. **Now with SIMD!***

[![Actions](https://img.shields.io/github/actions/workflow/status/althonos/pytrimal/test.yml?branch=main&logo=github&style=flat-square&maxAge=300)](https://github.com/althonos/pytrimal/actions)
[![Coverage](https://img.shields.io/codecov/c/gh/althonos/pytrimal?style=flat-square&maxAge=3600&logo=codecov)](https://codecov.io/gh/althonos/pytrimal/)
[![License](https://img.shields.io/badge/license-GPLv3-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/gpl-3.0/)
[![PyPI](https://img.shields.io/pypi/v/pytrimal.svg?style=flat-square&maxAge=3600&logo=PyPI)](https://pypi.org/project/pytrimal)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/pytrimal?style=flat-square&maxAge=3600&logo=anaconda)](https://anaconda.org/bioconda/pytrimal)
[![AUR](https://img.shields.io/aur/version/python-pytrimal?logo=archlinux&style=flat-square&maxAge=3600)](https://aur.archlinux.org/packages/python-pytrimal)
[![Wheel](https://img.shields.io/pypi/wheel/pytrimal.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pytrimal/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/pytrimal.svg?style=flat-square&maxAge=600&logo=python)](https://pypi.org/project/pytrimal/#files)
[![Python Implementations](https://img.shields.io/pypi/implementation/pytrimal.svg?style=flat-square&maxAge=600&label=impl)](https://pypi.org/project/pytrimal/#files)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pytrimal/)
[![Mirror](https://img.shields.io/badge/mirror-EMBL-009f4d?style=flat-square&maxAge=2678400)](https://git.embl.de/larralde/pytrimal/)
[![Issues](https://img.shields.io/github/issues/althonos/pytrimal.svg?style=flat-square&maxAge=600)](https://github.com/althonos/pytrimal/issues)
[![Docs](https://img.shields.io/readthedocs/pytrimal/latest?style=flat-square&maxAge=600)](https://pytrimal.readthedocs.io)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pytrimal/blob/main/CHANGELOG.md)
[![Downloads](https://img.shields.io/badge/dynamic/json?style=flat-square&color=303f9f&maxAge=86400&label=downloads&query=%24.total_downloads&url=https%3A%2F%2Fapi.pepy.tech%2Fapi%2Fprojects%2Fpytrimal)](https://pepy.tech/project/pytrimal)

***⚠️ This package is based on the release candidate of trimAl 2.0, and results
may not be consistent across versions or with the trimAl 1.4 results.***

## 🗺️ Overview

PytrimAl is a Python module that provides bindings to [trimAl](http://trimal.cgenomics.org/)
using [Cython](https://cython.org/). It implements a user-friendly, Pythonic
interface to use one of the different trimming methods from trimAl and
access results directly. It interacts with the trimAl internals, which has
the following advantages:

- **single dependency**: PytrimAl is distributed as a Python package, so you
  can add it as a dependency to your project, and stop worrying about the
  trimAl binary being present on the end-user machine.
- **no intermediate files**: Everything happens in memory, in a Python object
  you control, so you don't have to invoke the trimAl CLI using a
  sub-process and temporary files.
  [`Alignment`](https://pytrimal.readthedocs.io/en/latest/api/alignment.html#pytrimal.Alignment)
  objects can be created directly from Python code.
- **friendly interface**: The different trimming methods are implement as
  Python classes that can be configured independently.
- **error management**: Errors occuring in trimAl are converted
  transparently into Python exceptions, including an informative
  error message.
- **better performance**: PytrimAl uses *SIMD* instructions to compute
  statistics like pairwise sequence similarity. This makes the whole
  trimming process much faster for alignment with a large number of
  sequences, at the expense of slightly higher memory consumption.

## 📋 Roadmap

The following features are available or considered for implementation:

- [x] **automatic trimming**: Support for trimming alignments using one of the
  automatic heuristics implemented in trimAl.
- [x] **manual trimming**: Support for trimming alignments using manually
  defined conservation and gap thresholds for each residue position.
- [x] **overlap trimming**: Trimming sequences using residue and sequence
  overlaps to exclude regions with minimal conservation.
- [x] **representative trimming**: Select only representative sequences
  from the alignment, either using a fixed number, or a maximum identity
  threshold.
- [x] **alignment loading from disk**: Load an alignment from disk given
  a filename.
- [x] **alignment loading from a file-like object**: Load an alignment from
  a Python [file object](https://docs.python.org/3/glossary.html#term-file-object)
  instead of a file on the local filesystem.
- [x] **aligment creation from Python**: Create an alignment from a collection
  of sequences stored in Python strings.
- [x] **alignment formatting to disk**: Write an alignment to a file given
  a filename in one of the supported file formats.
- [x] **alignment formatting to a file-like object**: Write an alignment to
  a file-like object in one of the supported file formats.
- [ ] **reverse-translation**: Back-translate a protein alignment to align
  the sequences in genomic space.
- [x] **alternative similarity matrix**: Specify an alternative similarity
  matrix for the alignment (instead of BLOSUM62).
- [x] **similarity matrix creation**: Create a similarity matrix from scratch
  from Python code.
- [x] **windows for manual methods**: Use a sliding window for computing
  statistics in manual methods.

## 🔧 Installing

PytrimAl is available for all modern versions (3.6+), with no external dependencies.

It can be installed directly from [PyPI](https://pypi.org/project/pytrimal/),
which hosts some pre-built wheels for the x86-64 architecture (Linux/OSX)
and the Aarch64 architecture (Linux only), as well as the code required to compile
from source with Cython:
```console
$ pip install pytrimal
```

Otherwise, pytrimal is also available as a [Bioconda](https://bioconda.github.io/)
package:
```console
$ conda install -c bioconda pytrimal
```

## 💡 Example

Let's load an `Alignment` from a file on the disk, and use the *strictplus*
method to trim it, before printing the `TrimmedAlignment` as a Clustal block:
```python
from pytrimal import Alignment, AutomaticTrimmer

ali = Alignment.load("pytrimal/tests/data/example.001.AA.clw")
trimmer = AutomaticTrimmer(method="strictplus")

trimmed = trimmer.trim(ali)
for name, seq in zip(trimmed.names, trimmed.sequences):
    print(name.decode().rjust(6), seq)
```

This should output the following:
```
Sp8    GIVLVWLFPWNGLQIHMMGII
Sp10   VIMLEWFFAWLGLEINMMVII
Sp26   GLFLAAANAWLGLEINMMAQI
Sp6    GIYLSWYLAWLGLEINMMAII
Sp17   GFLLTWFQLWQGLDLNKMPVF
Sp33   GLHMAWFQAWGGLEINKQAIL
```

You can then use the
[`dump`](https://pytrimal.readthedocs.io/en/latest/api/alignment.html#pytrimal.Alignment.dump)
method to write the trimmed alignment to a file or file-like
object. For instance, save the results in
[PIR format](https://www.bioinformatics.nl/tools/crab_pir.html)
to a file named `example.trimmed.pir`:
```python
trimmed.dump("example.trimmed.pir", format="pir")
```

## 🧶 Thread-safety

Trimmer objects are thread-safe, and the `trim` method is re-entrant.
This means you can batch-process alignments in parallel using a
[`ThreadPool`](https://docs.python.org/3/library/multiprocessing.html#multiprocessing.pool.ThreadPool)
with a single trimmer object:
```python
import glob
import multiprocessing.pool
from pytrimal import Alignment, AutomaticTrimmer

trimmer = AutomaticTrimmer()
alignments = map(Alignment.load, glob.iglob("pytrimal/tests/data/*.fasta"))

with multiprocessing.pool.ThreadPool() as pool:
    trimmed_alignments = pool.map(trimmer.trim, alignments)
```

## ⏱️ Benchmarks

Benchmarks were run on a [i7-10710U CPU](https://ark.intel.com/content/www/us/en/ark/products/196448/intel-core-i710710u-processor-12m-cache-up-to-4-70-ghz.html)
@ 1.10GHz, using a single core to time the computation of several statistics,
on a variable number of sequences from
[`example.014.AA.EggNOG.COG0591.fasta`](https://github.com/inab/trimal/blob/trimAl/dataset/example.014.AA.EggNOG.COG0591.fasta),
an alignment of 3583 sequences and 7287 columns.

![Benchmarks](https://raw.githubusercontent.com/althonos/pytrimal/main/bench/v0.5.4.svg)

Each graph measures the computation time of a single trimAl statistic
(see the [Statistics page](https://pytrimal.readthedocs.io/en/stable/statistics.html)
of the [online documentation](https://pytrimal.readthedocs.io/) for more
information.)

The `None` curve shows the time using the internal trimAl 2.0 code,
the `Generic` curve shows a generic C implementation with some more
optimizations, and the `SSE` curve shows the time spent using a dedicated
class with [SIMD](https://en.wikipedia.org/wiki/Single_instruction,_multiple_data)
implementations of the statistic computation.

## 💭 Feedback

### ⚠️ Issue Tracker

Found a bug ? Have an enhancement request ? Head over to the [GitHub issue tracker](https://github.com/althonos/pytrimal/issues)
if you need to report or ask something. If you are filing in on a bug,
please include as much information as you can about the issue, and try to
recreate the same bug in a simple, easily reproducible situation.


### 🏗️ Contributing

Contributions are more than welcome! See
[`CONTRIBUTING.md`](https://github.com/althonos/pytrimal/blob/main/CONTRIBUTING.md)
for more details.


## 📋 Changelog

This project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html)
and provides a [changelog](https://github.com/althonos/pytrimal/blob/main/CHANGELOG.md)
in the [Keep a Changelog](http://keepachangelog.com/en/1.0.0/) format.


## ⚖️ License

This library is provided under the [GNU General Public License v3.0](https://choosealicense.com/licenses/gpl-3.0/).
trimAl is developed by the [trimAl team](http://trimal.cgenomics.org/trimal_team) and is distributed under the
terms of the GPLv3 as well. See `vendor/trimal/LICENSE` for more information.

*This project is in no way not affiliated, sponsored, or otherwise endorsed
by the [trimAl authors](http://trimal.cgenomics.org/trimal_team). It was developed
by [Martin Larralde](https://github.com/althonos/) during his PhD project
at the [European Molecular Biology Laboratory](https://www.embl.de/) in
the [Zeller team](https://github.com/zellerlab).*
