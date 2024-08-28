PytrimAl |Stars|
================

.. |Stars| image:: https://img.shields.io/github/stars/althonos/pytrimal.svg?style=social&maxAge=3600&label=Star
   :target: https://github.com/althonos/pytrimal/stargazers

`Cython <https://cython.org/>`_ *bindings and Python interface to* `trimAl <http://trimal.cgenomics.org/>`_,
*a tool for automated alignment trimming.*

|Actions| |Coverage| |PyPI| |Bioconda| |AUR| |Wheel| |Versions| |Implementations| |License| |Source| |Mirror| |Issues| |Docs| |Changelog| |Downloads|

.. |Actions| image:: https://img.shields.io/github/actions/workflow/status/althonos/pytrimal/test.yml?branch=main&logo=github&style=flat-square&maxAge=300
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

.. |Downloads| image:: https://img.shields.io/pypi/dm/pytrimal?style=flat-square&color=303f9f&maxAge=86400&label=downloads
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
:doc:`Installation page <guide/install>` to find other ways to install ``pytrimal``.


Library
-------

.. toctree::
   :maxdepth: 2

   User Guide <guide/index>
   API Reference <api/index>


Related Projects
----------------

The following Python libraries may be of interest for bioinformaticians.

.. grid:: 1 3 5 5
   :gutter: 1

   .. grid-item-card:: :fas:`diamond` PyHMMER
      :link: https://pyhmmer.readthedocs.io

      Profile Hidden Markov Models (with HMMER).

   .. grid-item-card:: :fas:`fire` Pyrodigal
      :link: https://pyrodigal.readthedocs.io

      Prokaryotic Gene Finding (with Prodigal).

   .. grid-item-card:: :fas:`virus-covid` Pyrodigal-gv
      :link: https://github.com/althonos/pyrodigal-gv

      Pyrodigal for Giant Viruses.

   .. grid-item-card:: :fas:`align-center` PyFAMSA
      :link: https://pyfamsa.readthedocs.io

      Multiple Sequence Alignment (with FAMSA).

   .. grid-item-card:: :fas:`scissors` PytrimAl
      :link: https://pytrimal.readthedocs.io

      Alignment Trimming (with trimAl).

   .. grid-item-card:: :fas:`music` LightMotif
      :link: https://lightmotif.readthedocs.io

      Platform-accelerated motif scoring.

   .. grid-item-card:: :fas:`knife;fa-custom` Diced
      :link: https://diced.readthedocs.io

      CRISPR Detection (with MinCED).

   .. grid-item-card:: :fas:`table-cells` Scoring Matrices
      :link: https://scoring-matrices.readthedocs.io

      Scoring matrices for Cython.

   .. grid-item-card:: :fas:`chain` Pyskani
      :link: https://pyskani.readthedocs.io

      Average Nucleotide Identity (with skani).

   .. grid-item-card:: :fas:`forward-fast` PyFastANI
      :link: https://pyfastani.readthedocs.io

      Average Nucleotide Identity (with FastANI).

   .. grid-item-card:: :fas:`magnifying-glass` PyJess
      :link: https://pyjess.readthedocs.io

      Geometric Template Matching (with Jess).

   .. grid-item-card:: :fas:`repeat` PyTantan
      :link: https://pytantan.readthedocs.io

      Tandem Repeat Masking (with Tantan).

   .. grid-item-card:: :fas:`gem` PyOpal
      :link: https://pyopal.readthedocs.io

      Query/Database Aligner (with Opal).

   .. grid-item-card:: :fas:`sword;fa-custom` PySWRD
      :link: https://pyswrd.readthedocs.io

      Database Heuristic Filtering (with SWORD).

   .. grid-item-card:: :fas:`rocket` Mini3di
      :link: https://github.com/althonos/mini3di

      Protein structure to 3di in pure Python.

   .. grid-item-card:: :fas:`calculator` ``peptides.py``
      :link: https://peptides.readthedocs.io

      Peptide descriptors for Python.

   .. grid-item-card:: :fas:`diagram-project` Pronto
      :link: https://pronto.readthedocs.io

      Open Biomedical Ontologies for Python.

   .. grid-item-card:: :fas:`box` NAFcodec
      :link: https://nafcodec.readthedocs.io

      Nucleotide Archival Format for Python.

   .. grid-item-card:: :fas:`bank` ``gb-io.py``
      :link: https://gb-io.readthedocs.io

      Fast GenBank parser for Python (with ``gb-io``).


License
-------

This library is provided under the `GNU General Public License v3.0 <https://choosealicense.com/licenses/gpl-3.0/>`_.
trimAl is developed by the `trimAl team <http://trimal.cgenomics.org/trimal_team>`_ and is distributed under the
terms of the GPLv3 as well. 

*This project is in no way not affiliated, sponsored, or otherwise endorsed by
the original* `trimAl`_ *authors. It was developed by* `Martin Larralde <https://github.com/althonos>`_ *during his
PhD project at the* `European Molecular Biology Laboratory <https://www.embl.de/>`_
*in the* `Zeller team <https://github.com/zellerlab>`_.
