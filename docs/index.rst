PytrimAl |Stars|
================

.. |Stars| image:: https://img.shields.io/github/stars/althonos/pytrimal.svg?style=social&maxAge=3600&label=Star
   :target: https://github.com/althonos/pytrimal/stargazers

`Cython <https://cython.org/>`_ *bindings and Python interface to* `trimAl <http://trimal.cgenomics.org/>`_,
*a tool for automated alignment trimming.*

|Actions| |Coverage| |PyPI| |Bioconda| |AUR| |Wheel| |Versions| |Implementations| |License| |Source| |Mirror| |Issues| |Docs| |Changelog| |Downloads|

.. |Actions| image:: https://img.shields.io/github/workflow/status/althonos/pytrimal/Test/main?logo=github&style=flat-square&maxAge=300
   :target: https://github.com/althonos/pytrimal/actions

.. |Coverage| image:: https://img.shields.io/codecov/c/gh/althonos/pytrimal?style=flat-square&maxAge=600
   :target: https://codecov.io/gh/althonos/pytrimal/

.. |PyPI| image:: https://img.shields.io/pypi/v/pytrimal.svg?style=flat-square&maxAge=3600
   :target: https://pypi.python.org/pypi/pytrimal

.. |Bioconda| image:: https://img.shields.io/conda/vn/bioconda/pytrimal?style=flat-square&maxAge=3600
   :target: https://anaconda.org/bioconda/pytrimal

.. |AUR| image:: https://img.shields.io/aur/version/python-pytrimal?logo=archlinux&style=flat-square&maxAge=3600
   :target: https://aur.archlinux.org/packages/python-pytrimal

.. |Wheel| image:: https://img.shields.io/pypi/wheel/pytrimal?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pytrimal/#files

.. |Versions| image:: https://img.shields.io/pypi/pyversions/pytrimal.svg?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pytrimal/#files

.. |Implementations| image:: https://img.shields.io/pypi/implementation/pytrimal.svg?style=flat-square&maxAge=3600&label=impl
   :target: https://pypi.org/project/pytrimal/#files

.. |License| image:: https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square&maxAge=3600
   :target: https://choosealicense.com/licenses/mit/

.. |Source| image:: https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square
   :target: https://github.com/althonos/pytrimal/

.. |Mirror| image:: https://img.shields.io/badge/mirror-EMBL-009f4d?style=flat-square&maxAge=2678400
   :target: https://git.embl.de/larralde/pytrimal/

.. |Issues| image:: https://img.shields.io/github/issues/althonos/pytrimal.svg?style=flat-square&maxAge=600
   :target: https://github.com/althonos/pytrimal/issues

.. |Docs| image:: https://img.shields.io/readthedocs/pytrimal?style=flat-square&maxAge=3600
   :target: http://pytrimal.readthedocs.io/en/stable/?badge=stable

.. |Changelog| image:: https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square
   :target: https://github.com/althonos/pytrimal/blob/main/CHANGELOG.md

.. |Downloads| image:: https://img.shields.io/badge/dynamic/json?style=flat-square&color=303f9f&maxAge=86400&label=downloads&query=%24.total_downloads&url=https%3A%2F%2Fapi.pepy.tech%2Fapi%2Fprojects%2Fpytrimal
   :target: https://pepy.tech/project/pytrimal


Overview
--------

PytrimAl is a Python module that provides bindings to trimAl using
`Cython <https://cython.org/>`_. It directly interacts with the trimAl
internals, which has the following advantages:

- **single dependency**: PytrimAl is distributed as a Python package, so you
  can add it as a dependency to your project, and stop worrying about the
  trimAl binary being present on the end-user machine.
- **no intermediate files**: Everything happens in memory, in a Python object
  you control, so you don't have to invoke the trimAl CLI using a
  sub-process and temporary files. `Alignment` objects can be created
  directly from Python code.
- **friendly interface**: The different trimming methods are implement as
  Python classes that can be configured independently.
- **error management**: Errors occuring in trimAl are converted
  transparently into Python exceptions, including an informative
  error message.


Setup
-----

Run ``pip install pytrimal`` in a shell to download the latest release and all
its dependencies from PyPi, or have a look at the
:doc:`Installation page <install>` to find other ways to install ``pytrimal``.


Library
-------

.. toctree::
   :maxdepth: 2

   Installation <install>
   Statistics <statistics>
   Examples <examples/index>
   Contributing <contributing>
   API Reference <api/index>
   Changelog <changes>


License
-------

This library is provided under the `GNU General Public License v3.0 <https://choosealicense.com/licenses/gpl-3.0/>`_.
trimAl is developed by the `trimAl team <http://trimal.cgenomics.org/trimal_team>`_ and is distributed under the
terms of the GPLv3 as well. The ``cpu_features`` library was written by
`Guillaume Chatelet <https://github.com/gchatelet>`_ and is licensed under the terms of the
`Apache License 2.0 <https://choosealicense.com/licenses/apache-2.0/>`_.

*This project is in no way not affiliated, sponsored, or otherwise endorsed by
the original* `trimAl`_ *authors. It was developed by* `Martin Larralde <https://github.com/althonos>`_ *during his
PhD project at the* `European Molecular Biology Laboratory <https://www.embl.de/>`_
*in the* `Zeller team <https://github.com/zellerlab>`_.
