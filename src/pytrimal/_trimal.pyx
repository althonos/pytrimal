# distutils: language = c++
# cython: language_level=3, linetrace=True, embedsignature=False, binding=False
"""Bindings to trimAl, a tool for automated alignment trimming.

Example:
    Create a multiple sequence alignment using the `~pytrimal.Alignment`
    class constructor::

        >>> import pytrimal
        >>> msa = pytrimal.Alignment(
        ...     names=[b"Sp8", b"Sp10", b"Sp26", b"Sp6", b"Sp17", b"Sp33"],
        ...     sequences=[
        ...         "-----GLGKVIV-YGIVLGTKSDQFSNWVVWLFPWNGLQIHMMGII",
        ...         "-------DPAVL-FVIMLGTIT-KFS--SEWFFAWLGLEINMMVII",
        ...         "AAAAAAAAALLTYLGLFLGTDYENFA--AAAANAWLGLEINMMAQI",
        ...         "-----ASGAILT-LGIYLFTLCAVIS--VSWYLAWLGLEINMMAII",
        ...         "--FAYTAPDLL-LIGFLLKTVA-TFG--DTWFQLWQGLDLNKMPVF",
        ...         "-------PTILNIAGLHMETDI-NFS--LAWFQAWGGLEINKQAIL",
        ...     ],
        ... )

    Then use a trimmer to trim the alignment, for instance with the
    *automated1* algorithm:

        >>> trimmer = AutomaticTrimmer("automated1")
        >>> trimmed = trimmer.trim(msa)

    The trimmed alignment supports the same methods as the original
    alignment object:

        >>> for sequence in trimmed.sequences:
        ...     print(sequence)
        VWLFPWNGLQIHMMGII
        EWFFAWLGLEINMMVII
        AAANAWLGLEINMMAQI
        SWYLAWLGLEINMMAII
        TWFQLWQGLDLNKMPVF
        AWFQAWGGLEINKQAIL

References:
    - Capella-Gutiérrez, Salvador, José M. Silla-Martínez, and Toni Gabaldón.
      *TrimAl: A Tool for Automated Alignment Trimming in Large-Scale
      Phylogenetic Analyses*. Bioinformatics 25, no. 15 (2009): 1972–73.
      :doi:`10.1093/bioinformatics/btp348`.

"""

# --- C imports --------------------------------------------------------------

cimport cython
from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_FORMAT, PyBUF_READ, PyBuffer_FillInfo
from cpython.bytes cimport PyBytes_FromStringAndSize, PyBytes_AsString
from cpython.list cimport PyList_New, PyList_SET_ITEM
from cpython.mem cimport PyMem_Free, PyMem_Malloc
from cpython.memoryview cimport PyMemoryView_FromMemory, PyMemoryView_GET_BUFFER
from cpython.ref cimport Py_INCREF
from cpython.unicode cimport (
    PyUnicode_New,
    PyUnicode_KIND,
    PyUnicode_DATA,
    PyUnicode_WRITE,
)

from libc.errno cimport errno
from libc.math cimport NAN, isnan, sqrt
from libc.stdio cimport printf
from libc.string cimport memset, memcpy
from libcpp cimport bool
from libcpp.string cimport string
from iostream cimport istream, ostream, stringbuf, filebuf, ios_base

cimport trimal
cimport trimal.alignment
cimport trimal.format_handling
cimport trimal.manager
cimport trimal.report_system
cimport trimal.similarity_matrix
cimport trimal.statistics
from trimal.statistics cimport ComputePlatform
from scoring_matrices.lib cimport ScoringMatrix

from pystreambuf cimport pyreadbuf, pyreadintobuf, pywritebuf

# --- Python imports ---------------------------------------------------------

import os
import threading
from scoring_matrices.lib import ScoringMatrix

__version__ = PROJECT_VERSION

# --- Constants --------------------------------------------------------------

_TARGET_CPU           = TARGET_CPU
_TARGET_SYSTEM        = TARGET_SYSTEM
_SSE2_BUILD_SUPPORT   = False
_SSE2_RUNTIME_SUPPORT = False
_AVX2_BUILD_SUPPORT   = False
_AVX2_RUNTIME_SUPPORT = False
_NEON_BUILD_SUPPORT   = False
_NEON_RUNTIME_SUPPORT = False

if TARGET_CPU == "x86":
    cimport cpu_features.x86
    _info = cpu_features.x86.GetX86Info()
    _SSE2_BUILD_SUPPORT   = SSE2_BUILD_SUPPORT
    _AVX2_BUILD_SUPPORT   = AVX2_BUILD_SUPPORT
    _SSE2_RUNTIME_SUPPORT = SSE2_BUILD_SUPPORT and _info["features"]["sse2"] != 0
    _AVX2_RUNTIME_SUPPORT = AVX2_BUILD_SUPPORT and _info["features"]["avx2"] != 0
elif TARGET_CPU == "x86_64" or TARGET_CPU == "amd64":
    cimport cpu_features.x86
    _info = cpu_features.x86.GetX86Info()
    _SSE2_BUILD_SUPPORT   = SSE2_BUILD_SUPPORT
    _AVX2_BUILD_SUPPORT   = AVX2_BUILD_SUPPORT
    _SSE2_RUNTIME_SUPPORT = SSE2_BUILD_SUPPORT  # always runtime support on x86-64
    _AVX2_RUNTIME_SUPPORT = AVX2_BUILD_SUPPORT and _info["features"]["avx2"] != 0
elif TARGET_CPU == "arm":
    cimport cpu_features.arm
    _info = cpu_features.arm.GetArmInfo()
    _NEON_BUILD_SUPPORT   = NEON_BUILD_SUPPORT
    _NEON_RUNTIME_SUPPORT = NEON_BUILD_SUPPORT and _info["features"]["neon"] != 0
elif TARGET_CPU == "aarch64" or TARGET_CPU == "arm64":
    _NEON_BUILD_SUPPORT   = NEON_BUILD_SUPPORT
    _NEON_RUNTIME_SUPPORT = NEON_BUILD_SUPPORT  # always runtime support on Aarch64

if _AVX2_RUNTIME_SUPPORT:
    _BEST_PLATFORM = ComputePlatform.AVX2
elif _SSE2_RUNTIME_SUPPORT:
    _BEST_PLATFORM = ComputePlatform.SSE2
elif _NEON_RUNTIME_SUPPORT:
    _BEST_PLATFORM = ComputePlatform.NEON
else:
    _BEST_PLATFORM = ComputePlatform.NONE

# --- Utilities --------------------------------------------------------------

ctypedef fused number:
    int
    float
    ssize_t

cdef number _check_range(number value, str name, number min_value, number max_value) except *:
    if value < min_value or value > max_value or isnan(value):
        raise ValueError(f"Invalid value for `{name}`: {value!r}")
    return value

cdef number _check_positive(number value, str name) except *:
    if value <= 0:
        raise ValueError(f"Invalid value for `{name}`: {value!r}")
    return value

cdef int _check_fileobj_read(object fileobj) except 1:
    cdef str ty = type(fileobj).__name__
    if not hasattr(fileobj, "seek") or not fileobj.seekable():
        raise TypeError(f"{ty!r} object is not seekable.")
    if not hasattr(fileobj, "readinto") and not hasattr(fileobj, "read"):
        raise TypeError(f"{ty!r} object has no attribute 'read'.")
    try:
        b = bytearray(0)
        fileobj.readinto(b)
    except Exception as err:
        raise TypeError(f"{ty!r} object is not open in binary mode.") from err
    return 0

cdef extern from *:
    """
    template <typename T>
    T* new_array(size_t n) {
        return new T[n];
    }
    template <typename T>
    void del_array(T* array) {
        delete[] array;
    }
    """
    T* new_array[T](size_t)
    void del_array[T](T*)

cdef extern from "<ios>":
    """
    std::ios_base::openmode READMODE = std::ios_base::in;
    std::ios_base::openmode WRITEMODE = std::ios_base::out | std::ios_base::trunc;
    """
    const ios_base.openmode READMODE
    const ios_base.openmode WRITEMODE


# --- Alignment classes ------------------------------------------------------

@cython.freelist(8)
cdef class AlignmentSequences:
    """A read-only view over the sequences of an alignment.

    Objects from this class are created in the `~Alignment.sequences`
    property of `~pytrimal.Alignment` objects. Use it to access the
    string data of individual rows from the alignment::

        >>> msa = Alignment.load("example.001.AA.clw")
        >>> len(msa.sequences)
        6
        >>> msa.sequences[0]
        '-----GLGKVIV-YGIVLGTKSDQFSNWVVWLFPWNGLQIHMMGII'
        >>> sum(seq.count('-') for seq in msa.sequences)
        43

    A slice over a subset of the sequences can be obtained as well without
    having to copy the internal data, allowing to create a new `Alignment`
    with only some sequences from the original one::

        >>> msa2 = Alignment(msa.names[:4:2], msa.sequences[:4:2])
        >>> len(msa2.sequences)
        2
        >>> msa2.sequences[1] == msa.sequences[2]
        True

    .. versionadded:: 0.4.0
       Support for zero-copy slicing.

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self, Alignment alignment):
        self._owner = alignment
        self._ali = alignment._ali
        self._index_mapping = alignment._sequences_mapping
        self._length = alignment._ali.numberOfSequences
        self._free_mapping = False

    def __dealloc__(self):
        if self._free_mapping:
            PyMem_Free(self._index_mapping)

    def __len__(self):
        return self._length

    def __getitem__(self, object index):
        assert self._ali is not NULL
        if isinstance(index, slice):
            start, stop, stride = index.indices(self._length)
            return self._slice(start, stop, stride)
        else:
            return self._sequence(index)

    # --- Utils --------------------------------------------------------------

    cdef AlignmentSequences _slice(self, int start, int stop, int stride):
        """Return a view of a subset of the alignment sequences, without copy.
        """
        cdef object             indices = range(start, stop, stride)
        cdef int                newlen  = len(indices)
        cdef AlignmentSequences view    = AlignmentSequences.__new__(AlignmentSequences, self._owner)

        view._length = newlen
        view._free_mapping = True
        view._index_mapping = <int*> PyMem_Malloc(newlen * sizeof(int))
        if view._index_mapping is NULL:
            raise MemoryError()

        for i, x in enumerate(indices):
            assert i < newlen
            if self._index_mapping is not NULL:
                assert x < self._length
                x = self._index_mapping[x]
            view._index_mapping[i] = x

        return view

    cdef str _sequence(self, int index):
        """Return a single sequence in the alignment, creating a new string.
        """
        cdef int    kind
        cdef object seq
        cdef char*  data
        cdef size_t x      = 0
        cdef int    index_ = index

        if index_ < 0:
            index_ += self._length
        if index_ < 0 or index_ >= self._length:
            raise IndexError(index)
        if self._index_mapping is not NULL:
            index_ = self._index_mapping[index_]

        assert index_ < self._ali.originalNumberOfSequences
        if SYS_VERSION_INFO_MAJOR <= 3 and SYS_VERSION_INFO_MAJOR <= 7 and SYS_IMPLEMENTATION_NAME == "pypy":
            seq  = PyBytes_FromStringAndSize(NULL, self._ali.numberOfResidues)
            data = PyBytes_AsString(seq)
            for i in range(self._ali.originalNumberOfResidues):
                if self._ali.saveResidues is NULL or self._ali.saveResidues[i] != -1:
                    data[x] = self._ali.sequences[index_][i]
                    x += 1
            return seq.decode('ascii')
        else:
            seq  = PyUnicode_New(self._ali.numberOfResidues, 0x7f)
            data = <char*> PyUnicode_DATA(seq)
            kind = PyUnicode_KIND(seq)
            for i in range(self._ali.originalNumberOfResidues):
                if self._ali.saveResidues is NULL or self._ali.saveResidues[i] != -1:
                    PyUnicode_WRITE(kind, data, x, self._ali.sequences[index_][i])
                    x += 1
            return seq


@cython.freelist(8)
cdef class AlignmentResidues:
    """A read-only view over the residues of an alignment.

    Objects from this class are created in the `~Alignment.residues`
    property of `~pytrimal.Alignment` objects. Use it to access the
    string data of individual columns from the alignment::

        >>> msa = Alignment.load("example.001.AA.clw")
        >>> len(msa.residues)
        46
        >>> msa.residues[0]
        '--A---'
        >>> msa.residues[-1]
        'IIIIFL'

    .. versionadded:: 0.4.0
       Support for zero-copy slicing.

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self, Alignment alignment):
        self._owner = alignment
        self._ali = alignment._ali
        self._index_mapping = alignment._residues_mapping
        self._length = alignment._ali.numberOfResidues
        self._free_mapping = False

    def __dealloc__(self):
        if self._free_mapping:
            PyMem_Free(self._index_mapping)

    def __len__(self):
        return self._length

    def __getitem__(self, object index):
        assert self._ali is not NULL
        if isinstance(index, slice):
            start, stop, stride = index.indices(self._length)
            return self._slice(start, stop, stride)
        else:
            return self._column(index)

    # --- Utils --------------------------------------------------------------

    cdef AlignmentResidues _slice(self, int start, int stop, int stride):
        """Return a view of a subset of the alignment residues, without copy.
        """
        cdef object            indices = range(start, stop, stride)
        cdef int               newlen  = len(indices)
        cdef AlignmentResidues view    = AlignmentResidues.__new__(AlignmentResidues, self._owner)

        view._length = newlen
        view._free_mapping = True
        view._index_mapping = <int*> PyMem_Malloc(newlen * sizeof(int))
        if view._index_mapping is NULL:
            raise MemoryError()

        for i, x in enumerate(indices):
            assert i < newlen
            if self._index_mapping is not NULL:
                assert x < self._length
                x = self._index_mapping[x]
            view._index_mapping[i] = x

        return view

    cdef str _column(self, int index):
        """Return a single residue column in the alignment, creating a new string.
        """
        cdef object col
        cdef int    kind
        cdef char*  data
        cdef size_t x      = 0
        cdef int    index_ = index
        cdef int    length = self._ali.numberOfResidues

        if index_ < 0:
            index_ += length
        if index_ < 0 or index_ >= self._length:
            raise IndexError(index)
        if self._index_mapping is not NULL:
            index_ = self._index_mapping[index_]

        assert index_ < self._ali.originalNumberOfResidues
        if SYS_VERSION_INFO_MAJOR <= 3 and SYS_VERSION_INFO_MAJOR <= 7 and SYS_IMPLEMENTATION_NAME == "pypy":
            col  = PyBytes_FromStringAndSize(NULL, self._ali.numberOfSequences)
            data = PyBytes_AsString(col)
            for i in range(self._ali.originalNumberOfSequences):
                if self._ali.saveSequences is NULL or self._ali.saveSequences[i] != -1:
                    data[x] = self._ali.sequences[i][index_]
                    x += 1
            return col.decode('ascii')
        else:
            col = PyUnicode_New(self._ali.numberOfSequences, 0x7f)
            data = <char*> PyUnicode_DATA(col)
            kind = PyUnicode_KIND(col)
            for i in range(self._ali.originalNumberOfSequences):
                if self._ali.saveSequences is NULL or self._ali.saveSequences[i] != -1:
                    PyUnicode_WRITE(kind, data, x, self._ali.sequences[i][index_])
                    x += 1
            return col


cdef class Alignment:
    """A multiple sequence alignment.
    """

    # --- Conversion methods -------------------------------------------------

    @classmethod
    def from_biopython(cls, object alignment not None):
        """Create a new `Alignment` from an iterable of Biopython records.

        Arguments:
            alignment (iterable of `~Bio.SeqRecord.SeqRecord`): An iterable
                of Biopython records objects to build the alignment from.
                Passing a `Bio.Align.MultipleSeqAlignment` object is also
                supported.

        Returns:
            `~pytrimal.Alignment`: A new alignment object ready for trimming.

        .. versionadded:: 0.5.0

        """
        names = []
        sequences = []
        for record in alignment:
            names.append(record.id.encode("utf-8"))
            try:
                sequences.append(bytes(record.seq))
            except TypeError:
                sequences.append(str(record.seq))
        return cls(names=names, sequences=sequences)

    def to_biopython(self):
        """Create a new `~Bio.Align.MultipleSeqAlignment` from this `Alignment`.

        Returns:
            `~Bio.Align.MultipleSeqAlignment`: A multiple sequence alignment
            object as implemented in Biopython.

        Raises:
            `ImportError`: When the `Bio` module cannot be imported.

        .. versionadded:: 0.5.0

        """
        import Bio.Align
        import Bio.SeqRecord
        import Bio.Seq

        if isinstance(Bio, ImportError):
            raise Bio
        records = []
        for name, seq in zip(self.names, self.sequences):
            record = Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(seq), name.decode("utf-8"))
            records.append(record)
        return Bio.Align.MultipleSeqAlignment(records)

    @classmethod
    def from_pyhmmer(cls, object alignment not None):
        """Create a new `Alignment` from a `pyhmmer.easel.TextMSA`.

        Arguments:
            alignment (`~pyhmmer.easel.TextMSA`): A PyHMMER object storing
                a multiple sequence alignment in text format.

        Returns:
            `~pytrimal.Alignment`: A new alignment object ready for trimming.

        .. versionadded:: 0.5.0

        """
        return cls(names=alignment.names, sequences=alignment.alignment)

    def to_pyhmmer(self):
        """Create a new `~pyhmmer.easel.TextMSA` from this `Alignment`.

        Returns:
            `~pyhmmer.easel.TextMSA`: A PyHMMER multiple sequence alignment
            in text mode.

        Raises:
            `ImportError`: When the `pyhmmer` module cannot be imported.

        .. versionadded:: 0.5.0

        """
        import pyhmmer.easel

        return pyhmmer.easel.TextMSA(
            sequences=[
                pyhmmer.easel.TextSequence(name=name, sequence=seq)
                for name, seq in zip(self.names, self.sequences)
            ]
        )

    # --- Parser / Loader ----------------------------------------------------

    @classmethod
    def load(cls, object file not None, str format = None):
        """Load a multiple sequence alignment from a file.

        Arguments:
            path (`str`, `bytes`, `os.PathLike` or file-like object): The
                file from which to read the alignment. If a file-like object
                is given, it must be open in *binary* mode and support random
                access with the ``seek`` method. Otherwise, ``file`` is
                treated as a path.
            format (`str`, *optional*): The file-format the alignment is
                stored in. Must be given when loading from a file-like
                object, will be autodetected when reading from a file.

        Returns:
            `~pytrimal.Alignment`: The deserialized alignment.

        Example:
            >>> msa = Alignment.load("example.001.AA.clw")
            >>> msa.names
            [b'Sp8', b'Sp10', b'Sp26', b'Sp6', b'Sp17', b'Sp33']

        .. versionchanged:: 0.3.0
           Add support for reading code from a file-like object.

        """
        cdef trimal.format_handling.FormatManager      manager
        cdef trimal.format_handling.BaseFormatHandler* handler

        cdef string         path_
        cdef char           cbuffer[DEFAULT_BUFFER_SIZE]
        cdef filebuf        fbuffer
        cdef pyreadbuf*     rbuffer   = NULL
        cdef pyreadintobuf* r2buffer  = NULL
        cdef istream*       stream    = NULL
        cdef Alignment      alignment = Alignment.__new__(Alignment)

        if SYS_VERSION_INFO_MAJOR == 3 and SYS_VERSION_INFO_MINOR < 6:
            TYPES = (str, bytes)
        else:
            TYPES = (str, bytes, os.PathLike)
        if isinstance(file, TYPES):
            # check that file exists and is not a directory
            if not os.path.exists(file):
                raise FileNotFoundError(file)
            elif os.path.isdir(file):
                raise IsADirectoryError(file)
            # load the alignment from the given path
            path_ = os.fsencode(file)
            alignment._ali = manager.loadAlignment(path_)
        else:
            # check the file-like object has all the required features
            _check_fileobj_read(file)
            # clear buffer (crash on PyPI otherwise)
            memset(cbuffer, 0, DEFAULT_BUFFER_SIZE*sizeof(char))
            # make sure a format was given
            if format is None:
                raise ValueError("Format must be specified when loading from a file-like object")
            # get the right format handler
            handler = manager.getFormatFromToken(format.lower().encode('ascii'))
            if handler is NULL:
                raise ValueError(f"Unknown alignment format: {format!r}")
            # create a file-like object wrapper
            # attempt to use `readinto` if available and not on PyPy
            if SYS_IMPLEMENTATION_NAME == "cpython" and hasattr(file, "readinto"):
                r2buffer = new pyreadintobuf(file)
                r2buffer.pubsetbuf(cbuffer, 512)
                stream = new istream(r2buffer)
            else:
                rbuffer = new pyreadbuf(file)
                rbuffer.pubsetbuf(cbuffer, 512)
                stream = new istream(rbuffer)
            # load the alignment from the istream
            try:
                if handler.CheckAlignment(stream) == 0:
                    raise RuntimeError(f"Failed to recognize format {format!r} in {file!r}")
                stream.seekg(0)
                alignment._ali = handler.LoadAlignment(stream[0])
            finally:
                del stream
                del rbuffer
                del r2buffer

        if alignment._ali is NULL:
            raise RuntimeError(f"Failed to load alignment from {file!r}.")
        return alignment

    cpdef void dump(self, object file, str format="fasta") except *:
        """Dump the alignment to a file or a file-like object.

        Arguments:
            file (`str`, `bytes`, `os.PathLike` or file-like object): The
                file to which to write the alignment. If a file-like object
                is given, it must be open in *binary* mode. Otherwise,
                ``file`` is treated as a path.
            format (`str`): The name of the alignment format to write. See
                below for a list of supported formats.

        Raises:
            `ValueError`: When ``format`` is not a recognized file format.
            `OSError`: When the path given as ``file`` could not be opened.

        Hint:
            The alignment can be written in one of the following formats:

            ``clustal``
              The alignment format produced by the Clustal and Clustal Omega
              alignment softwares.

            ``fasta``
              The aligned FASTA format, which outputs all sequences
              in the alignment as FASTA records with gap characters
              (see :wiki:`FASTA format`).

            ``html``
              An HTML report showing alignment in pseudo-Clustal format with
              colored residues.

            ``mega``
              The alignment format produced by the
              `MEGA <https://www.megasoftware.net>`_ software for evolutionary
              analysis of alignments.

            ``nexus``
              The NEXUS alignment format (see :wiki:`Nexus file`).

            ``phylip`` (or ``phylip40``):
              The PHYLIP 4.0 alignment format.

            ``phylip32``
              The PHYLIP 3.2 alignment format.

            ``phylippaml``
              A variant of PHYLIP 4.0 compatible with the
              `PAML <http://abacus.gene.ucl.ac.uk/software/paml.html>`_ tool
              for phylogenetic analysis.

            ``nbrf`` or ``pir``
              The format of Protein Information Resource database files,
              provided by the National Biomedical Research Foundation.

            Additionally, the ``fasta``, ``nexus``, ``phylippaml``, ``phylip32``,
            and  ``phylip40`` formats support an ``_m10`` variant, which limits
            the sequence names to 10 characters.

        .. versionadded:: 0.2.2

        """
        assert self._ali != NULL

        cdef bool                                      is_fileobj
        cdef bytes                                     path_
        cdef trimal.format_handling.FormatManager      manager
        cdef trimal.format_handling.BaseFormatHandler* handler
        cdef filebuf                                   fbuffer
        cdef pywritebuf*                               pbuffer    = NULL
        cdef ostream*                                  stream     = NULL

        handler = manager.getFormatFromToken(format.lower().encode('ascii'))
        if handler is NULL:
            raise ValueError(f"Could not recognize alignment format: {format!r}")

        if SYS_VERSION_INFO_MAJOR == 3 and SYS_VERSION_INFO_MINOR < 6:
            TYPES = (str, bytes)
        else:
            TYPES = (str, bytes, os.PathLike)
        if isinstance(file, TYPES):
            path_ = os.fsencode(file)
            if fbuffer.open(<const char*> path_, WRITEMODE) is NULL:
                raise OSError(errno, f"Failed to open {file!r}")
            stream = new ostream(&fbuffer)
        else:
            pbuffer = new pywritebuf(file)
            stream = new ostream(pbuffer)

        try:
            handler.SaveAlignment(self._ali[0], stream)
        finally:
            del stream
            if pbuffer is not NULL:
                del pbuffer
            fbuffer.close()

    cpdef str dumps(self, str format="fasta", str encoding="utf-8"):
        """Dump the alignment to a string in the provided format.

        Arguments:
            format (`str`): The format of the alignment. See
                the `~Alignment.dump` method for a list of supported
                formats.
            encoding (`str`): The encoding to use to decode sequence names.

        Raises:
            `ValueError`: When ``format`` is not a recognized file format.

        .. versionadded:: 0.2.2

        """
        assert self._ali != NULL

        cdef stringbuf                                 buffer
        cdef ostream*                                  stream
        cdef trimal.format_handling.FormatManager      manager
        cdef trimal.format_handling.BaseFormatHandler* handler

        handler = manager.getFormatFromToken(format.lower().encode('ascii'))
        if handler is NULL:
            raise ValueError(f"Could not recognize alignment format: {format!r}")

        buffer = stringbuf()
        stream = new ostream(&buffer)

        try:
            handler.SaveAlignment(self._ali[0], stream)
            return buffer.str().decode(encoding)
        finally:
            del stream

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._ali = NULL
        self._sequences_mapping = NULL
        self._residues_mapping = NULL

    def __dealloc__(self):
        if self._ali is not NULL:
            del self._ali
        if self._sequences_mapping is not NULL:
            PyMem_Free(self._sequences_mapping)
        if self._residues_mapping is not NULL:
            PyMem_Free(self._residues_mapping)

    def __init__(self, object names not None, object sequences not None):
        """__init__(self, names, sequences)\n--\n

        Create a new alignment with the given names and sequences.

        Arguments:
            names (`~collections.abc.Sequence` of `bytes`): The names of
                the sequences in the alignment.
            sequences (`~collections.abc.Sequence` of `bytes` or `str`): The
                actual sequences in the alignment.

        Examples:
            Create a new alignment with a list of sequences and a list of
            names::

                >>> alignment = Alignment(
                ...     names=[b"Sp8", b"Sp10", b"Sp26"],
                ...     sequences=[
                ...         "-----GLGKVIV-YGIVLGTKSDQFSNWVVWLFPWNGLQIHMMGII",
                ...         "-------DPAVL-FVIMLGTIT-KFS--SEWFFAWLGLEINMMVII",
                ...         "AAAAAAAAALLTYLGLFLGTDYENFA--AAAANAWLGLEINMMAQI",
                ...     ]
                ... )

            There should be as many sequences as there are names, otherwise
            a `ValueError` will be raised::

                >>> Alignment(
                ...     names=[b"Sp8", b"Sp10", b"Sp26"],
                ...     sequences=["GLQIHMMGII", "GLEINMMVII"]
                ... )
                Traceback (most recent call last):
                ...
                ValueError: `Alignment` given 3 names but 2 sequences

            Sequence characters will be checked, and an error will be
            raised if they are not one of the characters from a biological
            alphabet::

                >>> Alignment(
                ...     names=[b"Sp8", b"Sp10"],
                ...     sequences=["GLQIHMMGII", "GLEINMM123"]
                ... )
                Traceback (most recent call last):
                ...
                ValueError: The sequence "Sp10" has an unknown (49) character

        """
        cdef bytes  name
        cdef object sequence
        cdef int    nresidues     = -1
        cdef bool   validate_seqs = not isinstance(sequences, AlignmentSequences)

        if len(names) != len(sequences):
            raise ValueError(f"`Alignment` given {len(names)!r} names but {len(sequences)!r} sequences")

        self._ali = new trimal.alignment.Alignment()
        self._ali.numberOfSequences = len(sequences)
        self._ali.seqsName  = new_array[string](self._ali.numberOfSequences)
        self._ali.sequences = new_array[string](self._ali.numberOfSequences)

        for i, (name, sequence) in enumerate(zip(names, sequences)):

            if not self._ali.numberOfResidues:
                self._ali.numberOfResidues = len(sequence)
            if len(sequence) != self._ali.numberOfResidues:
                raise ValueError(f"Sequence length mismatch in sequence {i}: {len(sequence)} != {self._ali.numberOfResidues)}")

            self._ali.seqsName[i]  = name
            if isinstance(sequence, bytes):
                self._ali.sequences[i] = sequence
            else:
                self._ali.sequences[i] = sequence.encode('ascii')

        if self._ali.numberOfSequences == 0:
            self._ali.numberOfResidues = 0
        if self._ali.numberOfResidues > 0:
            self._ali.fillMatrices(self._ali.numberOfSequences > 1, validate_seqs)

        self._ali.originalNumberOfSequences = self._ali.numberOfSequences
        self._ali.originalNumberOfResidues = self._ali.numberOfResidues

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}(names={self.names!r}, sequences={list(self.sequences)!r})"

    def __copy__(self):
        return self.copy()

    # --- Properties ---------------------------------------------------------

    @property
    def names(self):
        """sequence of `bytes`: The names of the sequences in the alignment.
        """
        assert self._ali is not NULL
        assert self._ali.seqsName is not NULL

        cdef int    i
        cdef bytes  name
        cdef object names = []

        for i in range(self._ali.originalNumberOfSequences):
            if self._ali.saveSequences is NULL or self._ali.saveSequences[i] != -1:
                name = PyBytes_FromStringAndSize( self._ali.seqsName[i].data(), self._ali.seqsName[i].size() )
                names.append(name)

        return names

    @property
    def sequences(self):
        """`~pytrimal.AlignmentSequences`: The sequences in the alignment.
        """
        assert self._ali is not NULL
        return AlignmentSequences(self)

    @property
    def residues(self):
        """`~pytrimal.AlignmentResidues`: The residues in the alignment.
        """
        assert self._ali is not NULL
        return AlignmentResidues(self)

    # --- Functions ----------------------------------------------------------

    cpdef Alignment copy(self):
        """Create a copy of this alignment.
        """
        assert self._ali is not NULL
        cdef Alignment copy = (type(self)).__new__(type(self))
        copy._ali = new trimal.alignment.Alignment(self._ali[0])
        return copy


cdef class TrimmedAlignment(Alignment):
    """A multiple sequence alignment that has been trimmed.

    Internally, the trimming process produces a mask of sequences and
    a mask of residues. This class only exposes the filtered sequences
    and residues.

    Example:
        Create a trimmed alignment using two lists to filter out some
        residues and sequences::

            >>> trimmed = TrimmedAlignment(
            ...    names=[b"Sp8", b"Sp10", b"Sp26"],
            ...    sequences=["QFSNWV", "KFS--S", "NFA--A"],
            ...    sequences_mask=[True, True, False],
            ...    residues_mask=[True, True, True, False, False, True],
            ... )

        The `~TrimmedAlignment.names` and `~TrimmedAlignment.sequences`
        properties will only contain the retained sequences and residues::

            >>> list(trimmed.names)
            [b'Sp8', b'Sp10']
            >>> list(trimmed.sequences)
            ['QFSV', 'KFSS']

        Use the `~TrimmedAlignment.original_alignment` method to build
        the original unfiltered alignment containing all sequences and
        residues:

            >>> ali = trimmed.original_alignment()
            >>> list(ali.names)
            [b'Sp8', b'Sp10', b'Sp26']
            >>> list(ali.sequences)
            ['QFSNWV', 'KFS--S', 'NFA--A']

    """

    # --- Parser / Loader ----------------------------------------------------

    @classmethod
    def load(cls, object file not None, str format = None):
        # For compatibility, allow loading a trimmed alignment from a file
        # even though it makes it effectively not trimmed
        cdef Alignment alignment = Alignment.load(file, format)
        cdef TrimmedAlignment trimmed = TrimmedAlignment.__new__(TrimmedAlignment)
        trimmed._ali = alignment._ali
        alignment._ali = NULL
        trimmed._build_index_mapping()
        return trimmed

    # --- Magic methods ------------------------------------------------------

    def __init__(
        self,
        object names,
        object sequences not None,
        object sequences_mask = None,
        object residues_mask = None,
    ):
        """__init__(self, names, sequences, sequences_mask=None, residues_mask=None)\n--\n

        Create a new alignment with the given names, sequences and masks.

        Arguments:
            names (`~collections.abc.Sequence` of `bytes`): The names of
                the sequences in the alignment.
            sequences (`~collections.abc.Sequence` of `str`): The actual
                sequences in the alignment.
            sequences_mask (`~collections.abc.Sequence` of `bool`): A mask
                for which sequences to keep in the trimmed alignment. If
                given, must be as long as the ``sequences`` and ``names``
                list.
            residues_mask (`~collections.abc.Sequence` of `bool`): A mask
                for which residues to keep in the trimmed alignment. If
                given, must be as long as every element in the ``sequences``
                argument.

        """
        super().__init__(names, sequences)
        assert self._ali is not NULL

        cdef bool mask
        cdef int  i

        # mask sequences
        if sequences_mask is not None:
            if len(sequences_mask) != self._ali.originalNumberOfSequences:
                raise ValueError("Sequences mask must have the same length as the sequences list")
            self._ali.saveSequences = new_array[int](self._ali.originalNumberOfSequences)
            for i, mask in enumerate(sequences_mask):
                if mask:
                    self._ali.saveSequences[i] = i
                else:
                    self._ali.saveSequences[i] = -1
                    self._ali.numberOfSequences -= 1

        # mask residues
        if residues_mask is not None:
            if len(residues_mask) != self._ali.originalNumberOfResidues:
                raise ValueError("Sequences mask must have the same length as the sequences list")
            self._ali.saveResidues = new_array[int](self._ali.originalNumberOfResidues)
            for i, mask in enumerate(residues_mask):
                if mask:
                    self._ali.saveResidues[i] = i
                else:
                    self._ali.saveResidues[i] = -1
                    self._ali.numberOfResidues -= 1

        # build index mapping
        self._build_index_mapping()

    # --- Utils --------------------------------------------------------------

    cdef void _build_index_mapping(self) except *:
        assert self._ali is not NULL

        cdef ssize_t i
        cdef ssize_t x

        # create a mapping from new sequence index to old sequence index
        self._sequences_mapping = <int*> PyMem_Malloc(self._ali.numberOfSequences * sizeof(int))
        if self._sequences_mapping is NULL:
            raise MemoryError()
        x = 0
        for i in range(self._ali.originalNumberOfSequences):
            if self._ali.saveSequences is NULL or self._ali.saveSequences[i] != -1:
                self._sequences_mapping[x] = i
                x += 1

        # create a mapping from new residue index to old residue index
        self._residues_mapping = <int*> PyMem_Malloc(self._ali.numberOfResidues * sizeof(int))
        if self._residues_mapping is NULL:
            raise MemoryError()
        x = 0
        for i in range(self._ali.originalNumberOfResidues):
            if self._ali.saveResidues is NULL or self._ali.saveResidues[i] != -1:
                self._residues_mapping[x] = i
                x += 1

    # --- Properties ---------------------------------------------------------

    @property
    def residues_mask(self):
        """sequence of `bool`: Which residues are kept in the alignment.
        """
        assert self._ali is not NULL

        cdef int    i
        cdef object mask = PyList_New(self._ali.originalNumberOfResidues)

        for i in range(self._ali.originalNumberOfResidues):
            if self._ali.saveResidues is NULL or self._ali.saveResidues[i] != -1:
                Py_INCREF(True)
                PyList_SET_ITEM(mask, i, True)
            else:
                Py_INCREF(False)
                PyList_SET_ITEM(mask, i, False)

        return mask

    @property
    def sequences_mask(self):
        """sequence of `bool`: Which sequences are kept in the alignment.
        """
        assert self._ali is not NULL

        cdef int    i
        cdef object mask = PyList_New(self._ali.originalNumberOfSequences)

        for i in range(self._ali.originalNumberOfSequences):
            if self._ali.saveSequences is NULL or self._ali.saveSequences[i] != -1:
                Py_INCREF(True)
                PyList_SET_ITEM(mask, i, True)
            else:
                Py_INCREF(False)
                PyList_SET_ITEM(mask, i, False)

        return mask

    # --- Functions ----------------------------------------------------------

    cpdef Alignment original_alignment(self):
        """Rebuild the original alignment from which this object was obtained.

        Returns:
            `~pytrimal.Alignment`: The untrimmed alignment that produced
            this trimmed alignment.

        """
        assert self._ali is not NULL
        cdef Alignment orig = Alignment.__new__(Alignment)
        orig._ali = new trimal.alignment.Alignment(self._ali[0])
        del_array[int](orig._ali.saveSequences)
        del_array[int](orig._ali.saveResidues)
        orig._ali.saveSequences = NULL
        orig._ali.saveResidues = NULL
        orig._ali.numberOfSequences = orig._ali.originalNumberOfSequences
        orig._ali.numberOfResidues = orig._ali.originalNumberOfResidues
        return orig

    cpdef TrimmedAlignment terminal_only(self):
        """Get a trimmed alignment where only the terminal residues are removed.

        Returns:
            `~pytrimal.TrimmedAlignment`: The alignment where only terminal
            residues have been trimmed.

        """
        assert self._ali is not NULL
        cdef TrimmedAlignment term_only = TrimmedAlignment.__new__(TrimmedAlignment)
        term_only._ali = new trimal.alignment.Alignment(self._ali[0])
        term_only._ali.Cleaning.removeOnlyTerminal()
        term_only._build_index_mapping()
        return term_only

    cpdef TrimmedAlignment copy(self):
        """Create a copy of this trimmed alignment.
        """
        cdef TrimmedAlignment copy = TrimmedAlignment.__new__(TrimmedAlignment)
        copy._ali = new trimal.alignment.Alignment(self._ali[0])
        copy._build_index_mapping()
        return copy


# -- Trimmer classes ---------------------------------------------------------

cdef class BaseTrimmer:
    """A sequence alignment trimmer.

    All subclasses provide the same `trim` method, and are configured
    through their constructor.

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._platform = ComputePlatform.NONE

    def __init__(self, *, str platform = "detect"):
        """__init__(self, *, platform="detect")\n--\n

        Create a new base trimmer.

        Keyword Arguments:
            platform (`str`, *optional*): The compute platform to use
                to accelerate computation of pairwise similarity statistics.
                If `None` given, use the original code from trimAl. By default,
                attempt to auto-detect the best available platform for the host
                CPU.

        .. versionchanged:: 0.8.0
           Renamed ``backend`` keyword argument to ``platform``.

        """
        if TARGET_CPU == "x86" or TARGET_CPU == "x86_64":
            if platform =="detect":
                if SSE2_BUILD_SUPPORT and _SSE2_RUNTIME_SUPPORT:
                    self._platform = ComputePlatform.SSE2
                if AVX2_BUILD_SUPPORT and _AVX2_RUNTIME_SUPPORT:
                    self._platform = ComputePlatform.AVX2
            elif platform == "avx2":
                if not AVX2_BUILD_SUPPORT:
                    raise RuntimeError("Extension was compiled without AVX2 support")
                elif not _AVX2_RUNTIME_SUPPORT:
                    raise RuntimeError("Cannot run AVX2 instructions on this machine")
                else:
                    self._platform = ComputePlatform.AVX2
            elif platform == "sse2":
                if not SSE2_BUILD_SUPPORT:
                    raise RuntimeError("Extension was compiled without SSE2 support")
                elif not _SSE2_RUNTIME_SUPPORT:
                    raise RuntimeError("Cannot run SSE2 instructions on this machine")
                else:
                    self._platform = ComputePlatform.SSE2
            elif platform is None:
                self._platform = ComputePlatform.NONE
            else:
                raise ValueError(f"Unsupported platform on this architecture: {platform!r}")
        elif TARGET_CPU == "arm" or TARGET_CPU == "aarch64":
            if platform == "detect":
                self._platform = ComputePlatform.NONE
                if NEON_BUILD_SUPPORT and _NEON_RUNTIME_SUPPORT:
                    self._platform = ComputePlatform.NEON
            elif platform == "neon":
                if not NEON_BUILD_SUPPORT:
                    raise RuntimeError("Extension was compiled without NEON support")
                elif not _NEON_RUNTIME_SUPPORT:
                    raise RuntimeError("Cannot run NEON instructions on this machine")
                else:
                    self._platform = ComputePlatform.NEON
            elif platform is None:
                self._platform = ComputePlatform.NONE
            else:
                raise ValueError(f"Unsupported platform on this architecture: {platform!r}")
        else:
            if platform == "detect" or platform == "generic":
                self._platform = ComputePlatform.NONE
            elif platform is None:
                self._platform = ComputePlatform.NONE
            else:
                raise ValueError(f"Unsupported platform on this architecture: {platform!r}")

    def __repr__(self):
        cdef str ty  = type(self).__name__
        cdef str arg = ""
        if self._platform != _BEST_PLATFORM:
            arg = f"platform={self.platform!r}"
        return f"{ty}({arg})"

    def __getstate__(self):
        return {
            "platform": self.platform,
        }

    def __setstate__(self, dict state):
        try:
            self.__init__(platform=state["platform"])
        except (ValueError, RuntimeError):
            self.__init__(platform="detect")

    # --- Properties ---------------------------------------------------------

    @property
    def platform(self):
        """`str` or `None`: The compute platform for this trimmer.

        .. versionchanged:: 0.8.0
           Renamed ``backend`` property to ``platform``.

        """
        if self._platform == ComputePlatform.SSE2:
            return "sse2"
        elif self._platform == ComputePlatform.AVX2:
            return "avx2"
        elif self._platform == ComputePlatform.NEON:
            return "neon"
        elif self._platform == ComputePlatform.NONE:
            return None

    # --- Utils --------------------------------------------------------------

    cdef void _configure_manager(self, trimal.manager.trimAlManager* manager):
        pass

    # --- Functions ----------------------------------------------------------

    cpdef TrimmedAlignment trim(self, Alignment alignment, SimilarityMatrix matrix = None):
        """Trim the provided alignment.

        Arguments:
            alignment (`~pytrimal.Alignment`): A multiple sequence
                alignment to trim.
            matrix (`~pytrimal.SimilarityMatrix`, optional): An alternative
                similarity matrix to use for computing the similarity
                statistic. If `None`, a default matrix will be used based
                on the type of the alignment.

        Returns:
            `~pytrimal.TrimmedAlignment`: The trimmed alignment.

        Hint:
            This method is re-entrant, and can be called safely accross
            different threads. Most of the computations will be done after
            releasing the GIL.

        .. versionchanged:: 0.1.2
           Added the ``matrix`` optional argument.

        """
        # use a local manager object so that this method is re-entrant
        cdef Alignment                    copy
        cdef trimal.manager.trimAlManager manager

        # copy the alignment to the manager object so that the original
        # alignment is left untouched; for trimmed alignments, we must
        # first extract the saved sequences and residues, otherwise the
        # trimming will occur on the orignal alignment instead of the
        # trimmed one!
        if isinstance(alignment, TrimmedAlignment):
            copy = Alignment(alignment.names, alignment.sequences)
            manager.origAlig = copy._ali
            copy._ali = NULL
        else:
            manager.origAlig = new trimal.alignment.Alignment(alignment._ali[0])

        # configure the manager (to be implemented by the different subclasses)
        self._configure_manager(&manager)

        with nogil:
            # setup computation of optimized statistics with SIMD
            alignment._ali.Statistics.platform = self._platform
            # set flags
            manager.set_window_size()
            if manager.blockSize != -1:
                manager.origAlig.setBlockSize(manager.blockSize)
            # set similarity matrix from argument or load a default one
            if matrix is not None:
                manager.origAlig.Statistics.setSimilarityMatrix(&matrix._smx)
            else:
                manager.create_or_use_similarity_matrix()
            # clean alignment
            manager.clean_alignment()
            # use original alignment as single alignment if needed
            if manager.singleAlig == NULL:
                manager.singleAlig = manager.origAlig
                manager.origAlig = NULL

        # trim alignment and create a TrimmedAlignment object
        cdef TrimmedAlignment trimmed = TrimmedAlignment.__new__(TrimmedAlignment)
        trimmed._ali = new trimal.alignment.Alignment(manager.singleAlig[0])
        trimmed._build_index_mapping()
        return trimmed


cdef class AutomaticTrimmer(BaseTrimmer):
    """A sequence alignment trimmer with automatic parameter detection.

    trimAl provides several heuristic methods for automated trimming of
    multiple sequence algorithms:

    - ``strict``: A statistical method that combines *gaps* and *similarity*
      statistics to clean the alignment.
    - ``strictplus``: A statistical method that combines *gaps* and
      *similarity* statistics, optimized for Neighbour-Joining tree
      reconstruction.
    - ``gappyout``: A statistical method that only uses *gaps* statistic
      to clean the alignment.
    - ``automated1``: A meta-method that chooses between ``strict`` and
      ``gappyout``, optimized for Maximum Likelihood phylogenetic tree
      reconstruction.
    - ``automated2``: A meta-method that use automated selection on
      ``gappyout`` mode and keeps a minimum amount of columns.
    - ``nogaps``: A naive method that removes every column containing at
      least one gap.
    - ``noallgaps``: A naive method that removes every column containing
      only gaps.
    - ``noduplicateseqs``: A naive method that removes sequences that are
      equal on the alignment, keeping the latest occurence.

    Hint:
        A Python `frozenset` containing all valid automatic trimming methods
        can be obtained with the `AutomaticTrimmer.METHODS` attribute. This
        can be useful for listing or validating methods beforehand, e.g. to
        build a CLI with `argparse`.

    .. versionadded:: 0.4.0
       The `AutomaticTrimmer.METHODS` class attribute.

    .. versionadded:: 0.5.0
       Support for `pickle` protocol.

    .. versionadded:: 0.9.0
       The ``automated2`` method.

    """

    METHODS = frozenset({
        "strict",
        "strictplus",
        "gappyout",
        "nogaps",
        "noallgaps",
        "automated1",
        "automated2",
        "noduplicateseqs",
    })

    # --- Magic methods ------------------------------------------------------

    def __init__(self, str method="strict", *, str platform="detect"):
        """__init__(self, method="strict", *, platform="detect")\n--\n

        Create a new automatic alignment trimmer using the given method.

        Arguments:
            method (`str`): The automatic aligment trimming method. See
                the documentation for `AutomaticTrimmer` for a list of
                supported values.

        Keyword Arguments:
            platform (`str`, *optional*): The compute platform to use
                to accelerate computation of pairwise similarity statistics.
                If `None` given, use the original code from trimAl. By default,
                attempt to auto-detect the best available platform for the host
                CPU.

        Raises:
            `ValueError`: When ``method`` is not one of the automatic
                alignment trimming methods supported by trimAl.

        .. versionadded:: 0.4.0
           The ``noduplicateseqs`` method.

        .. versionchanged:: 0.8.0
           Renamed ``backend`` keyword argument to ``platform``.

        """
        super().__init__(platform=platform)

        if method not in self.METHODS:
            raise ValueError(f"Invalid value for `method`: {method!r}")
        self.method = method

    def __repr__(self):
        cdef str ty    = type(self).__name__
        cdef list args = [repr(self.method)]
        if self._platform != _BEST_PLATFORM:
            args.append(f"platform={self.platform!r}")
        return f"{ty}({', '.join(args)})"

    def __getstate__(self):
        return {
            "method":  self.method,
            "platform": self.platform,
        }

    def __setstate__(self, dict state):
        try:
            BaseTrimmer.__init__(self, platform=state["platform"])
        except (ValueError, RuntimeError):
            BaseTrimmer.__init__(self, platform="detect")
        self.method = state["method"]

    # --- Utils --------------------------------------------------------------

    cdef void _configure_manager(self, trimal.manager.trimAlManager* manager):
        BaseTrimmer._configure_manager(self, manager)
        manager.automatedMethodCount = 1
        if self.method == "strict":
            manager.strict = True
        elif self.method == "strictplus":
            manager.strictplus = True
        elif self.method == "gappyout":
            manager.gappyout = True
        elif self.method == "nogaps":
            manager.nogaps = True
        elif self.method == "noallgaps":
            manager.noallgaps = True
        elif self.method == "automated1":
            manager.automated1 = True
        elif self.method == "automated2":
            manager.automated2 = True
        elif self.method == "noduplicateseqs":
            manager.removeDuplicates = True


cdef class ManualTrimmer(BaseTrimmer):
    """A sequence alignment trimmer with manually defined thresholds.

    Manual trimming allows the user to specify independent thresholds for
    two different statistics:

    - *Gap threshold*: Remove columns where the gap ratio (or the absolute
      gap count) is higher than the provided threshold.
    - *Similarity threshold*: Remove columns with a similarity ratio lower
      than the provided threshold.

    In addition, the trimming can be restricted so that at least a
    configurable fraction of the original alignment is retained, in order
    to avoid stripping an alignment of distance sequences by aggressive
    trimming.

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._gap_threshold           = -1
        self._gap_absolute_threshold  = -1
        self._similarity_threshold    = -1
        self._conservation_percentage = -1
        self._window                  = -1
        self._gap_window              = -1
        self._similarity_window       = -1

    def __init__(
        self,
        *,
        object gap_threshold           = None,
        object gap_absolute_threshold  = None,
        object similarity_threshold    = None,
        object conservation_percentage = None,
        object window                  = None,
        object gap_window              = None,
        object similarity_window       = None,
        str    platform                 = "detect",
    ):
        """__init__(self, *, gap_threshold=None, gap_absolute_threshold=None, similarity_threshold=None, conservation_percentage=None, window=None, gap_window=None, similarity_window=None, platform="detect")\n--\n

        Create a new manual alignment trimmer with the given parameters.

        Keyword Arguments:
            gap_threshold (`float`, *optional*): The minimum fraction of
                 non-gap characters that must be present in a column to keep
                 the column.
            gap_absolute_threshold (`int`, *optional*): The absolute number
                of gaps allowed on a column to keep it in the alignment.
                Incompatible with ``gap_threshold``.
            similarity_threshold (`float`, *optional*): The minimum average
                similarity required.
            conservation_percentage (`float`, *optional*): The minimum
                percentage of positions in the original alignment to
                conserve.
            window (`int`, *optional*): The size of the half-window to use
                when computing statistics for an alignment.
            gap_window (`int`, *optional*): The size of the half-window to
                use when computing the *gap* statistic for an alignment.
                Incompatible with ``window``.
            similarity_window (`int`, *optional*): The size of the
                half-window to use when computing the *similarity* statistic
                for an alignment. Incompatible with ``window``.
            platform (`str`, *optional*): The compute platform to use
                to accelerate computation of pairwise similarity statistics.
                If `None` given, use the original code from trimAl. By default,
                attempt to auto-detect the best available platform for the host
                CPU.

        .. versionadded:: 0.2.2
           The keyword arguments for controling the half-window sizes.

        .. versionchanged:: 0.4.0
           Removed ``consistency_threshold`` and ``consistency_window``.

        .. versionchanged:: 0.8.0
           Renamed ``backend`` keyword argument to ``platform``.

        """
        super().__init__(platform=platform)

        if gap_threshold is not None and gap_absolute_threshold is not None:
            raise ValueError("Cannot specify both `gap_threshold` and `gap_absolute_threshold`")
        if window is not None and (gap_window is not None or similarity_window is not None):
            raise ValueError("Cannot specify both `window` and a specific window argument")

        if gap_threshold is not None:
            self._gap_threshold = 1 - _check_range[float](gap_threshold, "gap_threshold", 0, 1)
        if gap_absolute_threshold is not None:
            self._gap_absolute_threshold = _check_positive[ssize_t](gap_absolute_threshold, "gap_absolute_threshold")
        if similarity_threshold is not None:
            self._similarity_threshold = _check_range[float](similarity_threshold, "similarity_threshold", 0, 1)
        if conservation_percentage is not None:
            self._conservation_percentage = _check_range[float](conservation_percentage, "conservation_percentage", 0, 100)
        if window is not None:
            self._window = _check_positive[int](window, "window")
        if gap_window is not None:
            self._gap_window = _check_positive[int](gap_window, "gap_window")
        if similarity_window is not None:
            self._similarity_window = _check_positive[int](similarity_window, "similarity_window")

    def __repr__(self):
        cdef str ty    = type(self).__name__
        cdef list args = []
        if self._gap_threshold != -1:
            args.append(f"gap_threshold={1-self._gap_threshold!r}")
        if self._gap_absolute_threshold != -1:
            args.append(f"gap_absolute_threshold={self._gap_absolute_threshold!r}")
        if self._similarity_threshold != -1:
            args.append(f"similarity_threshold={self._similarity_threshold!r}")
        if self._conservation_percentage != -1:
            args.append(f"conservation_percentage={self._conservation_percentage!r}")
        if self._window != -1:
            args.append(f"window={self._window!r}")
        if self._gap_window != -1:
            args.append(f"gap_window={self._gap_window!r}")
        if self._similarity_window != -1:
            args.append(f"similarity_window={self._similarity_window!r}")
        if self._platform != _BEST_PLATFORM:
            args.append(f"platform={self.platform!r}")
        return f"{ty}({', '.join(args)})"

    def __getstate__(self):
        return {
            "platform":                 self.platform,
            "gap_threshold":           self._gap_threshold,
            "gap_absolute_threshold":  self._gap_absolute_threshold,
            "similarity_threshold":    self._similarity_threshold,
            "conservation_percentage": self._conservation_percentage,
            "window":                  self._window,
            "gap_window":              self._gap_window,
            "similarity_window":       self._similarity_window,
        }

    def __setstate__(self, dict state):
        try:
            BaseTrimmer.__init__(self, platform=state["platform"])
        except (ValueError, RuntimeError):
            BaseTrimmer.__init__(self, platform="detect")
        self._gap_threshold           = state["gap_threshold"]
        self._gap_absolute_threshold  = state["gap_absolute_threshold"]
        self._similarity_threshold    = state["similarity_threshold"]
        self._conservation_percentage = state["conservation_percentage"]
        self._window                  = state["window"]
        self._gap_window              = state["gap_window"]
        self._similarity_window       = state["similarity_window"]

    # --- Utils --------------------------------------------------------------

    cdef void _configure_manager(self, trimal.manager.trimAlManager* manager):
        manager.automatedMethodCount  = 0
        manager.gapThreshold          = self._gap_threshold
        manager.gapAbsoluteThreshold  = self._gap_absolute_threshold
        manager.similarityThreshold   = self._similarity_threshold
        manager.conservationThreshold = self._conservation_percentage
        manager.windowSize            = self._window
        manager.gapWindow             = self._gap_window
        manager.similarityWindow      = self._similarity_window


cdef class OverlapTrimmer(BaseTrimmer):
    """A sequence alignment trimmer for overlap blocks.

    Overlap trimming works by defining "good" positions, i.e. a position
    where most sequences agree (given a certain threshold) that the
    alignment contains a gap or a residue (independently of their agreement
    on that given residue). Sequences not containing enough "good" positions
    are then removed.

    Example:
        Consider the following alignment, where the last three sequences
        align decently, while the first sequence doesn't. In particular, it
        creates a large gap in the rest of the alignment::

            >>> ali = Alignment(
            ...     names=[b"Sp8", b"Sp17", b"Sp10", b"Sp26"],
            ...     sequences=[
            ...         "LG-----------TKSD---NNNNNNNNNNNNNNNNWV----------",
            ...         "APDLLL-IGFLLKTV-ATFG-----------------DTWFQLWQGLD",
            ...         "DPAVL--FVIMLGTI-TKFS-----------------SEWFFAWLGLE",
            ...         "AAALLTYLGLFLGTDYENFA-----------------AAAANAWLGLE",
            ...     ]
            ... )

        Let's create an overlap trimmer so that "good" positions correspond
        to an agreement between at least half of the sequences, and make it
        remove sequences with less than 40% of good positions::

            >>> trimmer = OverlapTrimmer(40.0, 0.5)

        Trimming will remove the first sequence because it doesn't contain
        enough good positions; then, the block containing only gaps will be
        removed (this is the default behaviour of all trimmer objects)::

            >>> trimmed = trimmer.trim(ali)
            >>> for name, seq in zip(trimmed.names, trimmed.sequences):
            ...     print(name.decode().ljust(5), seq)
            Sp17  APDLLL-IGFLLKTV-ATFGDTWFQLWQGLD
            Sp10  DPAVL--FVIMLGTI-TKFSSEWFFAWLGLE
            Sp26  AAALLTYLGLFLGTDYENFAAAAANAWLGLE

    .. versionadded:: 0.4.0

    .. versionadded:: 0.5.0
       Support for the `pickle` protocol.

    """

    # --- Magic methods ------------------------------------------------------

    def __init__(
        self,
        float sequence_overlap,
        float residue_overlap,
        *,
        str platform="detect"
    ):
        """__init__(self, sequence_overlap, residue_overlap, *, platform="detect")\n--\n

        Create a new overlap trimmer with the given thresholds.

        Arguments:
            sequence_overlap (`float`): The minimum percentage of "good"
                positions a sequence must contain to be kept in the
                alignment.
            residue_overlap (`float`): The fraction of matching residues
                a column must contain to be considered a "good" position.

        Keyword Arguments:
            platform (`str`, *optional*): The compute platform to use
                to accelerate computation of pairwise similarity statistics.
                If `None` given, use the original code from trimAl. By default,
                attempt to auto-detect the best available platform for the host
                CPU.

        """
        super().__init__(platform=platform)
        self._sequence_overlap = _check_range[float](sequence_overlap, "sequence_overlap", 0, 100)
        self._residue_overlap = _check_range[float](residue_overlap, "residue_overlap", 0, 1)

    def __repr__(self):
        cdef str ty    = type(self).__name__
        cdef list args = [repr(self._sequence_overlap), repr(self._residue_overlap)]
        if self._platform != _BEST_PLATFORM:
            args.append(f"platform={self.platform!r}")
        return f"{ty}({', '.join(args)})"

    def __getstate__(self):
        return {
            "platform":          self.platform,
            "sequence_overlap": self._sequence_overlap,
            "residue_overlap":  self._residue_overlap,
        }

    def __setstate__(self, dict state):
        try:
            BaseTrimmer.__init__(self, platform=state["platform"])
        except (ValueError, RuntimeError):
            BaseTrimmer.__init__(self, platform="detect")
        self._sequence_overlap = state["sequence_overlap"]
        self._residue_overlap  = state["residue_overlap"]

    # --- Utils --------------------------------------------------------------

    cdef void _configure_manager(self, trimal.manager.trimAlManager* manager):
        manager.automatedMethodCount  = 0
        manager.residuesOverlap       = self._residue_overlap
        manager.sequenceOverlap       = self._sequence_overlap


cdef class RepresentativeTrimmer(BaseTrimmer):
    """A sequence alignment trimmer for selecting representative sequences.

    Representative sequences on an alignment can be selected using a specific
    identity threshold, or a fixed number of representative sequences to keep.
    Representative trimming can be useful to reduce the weight of certain
    very similar sequences in an alignment, for instance to build a less
    conservative HMM.

    .. versionadded:: 0.5.0

    """

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._clusters = -1
        self._identity_threshold = -1

    def __init__(
        self,
        object clusters = None,
        object identity_threshold = None,
        *,
        str platform="detect"
    ):
        """__init__(self, clusters=None, identity_threshold=None, *, platform="detect")\n--\n

        Create a new representative alignment trimmer.

        Arguments:
            clusters (`int`, *optional*): The number of cluster
                representatives to keep in the trimmed alignment. Must be
                strictly positive. If the trimmer receives an alignment with
                less sequences than this, it will not perform any trimming.
            identity_threshold (`float`, *optional*): The identity threshold
                for which to get representative sequences.

        Keyword Arguments:
            platform (`str`, *optional*): The compute platform to use
                to accelerate computation of pairwise similarity statistics.
                If `None` given, use the original code from trimAl. By default,
                attempt to auto-detect the best available platform for the host
                CPU.

        Raises:
            `ValueError`: When both ``clusters`` and ``identity_threshold``
                are provided at the same time, or when they don't fall in
                a valid range.

        """
        super().__init__(platform=platform)
        if clusters is not None and identity_threshold is not None:
            raise ValueError("Cannot specify both `clusters` and `identity_threshold`")
        if clusters is not None:
            self._clusters = _check_positive[int](clusters, "clusters")
        if identity_threshold is not None:
            self._identity_threshold = _check_range[float](identity_threshold, "identity_threshold", 0, 1)

    def __repr__(self):
        cdef str ty    = type(self).__name__
        cdef list args = []
        if self._clusters != -1:
            args.append(f"clusters={self._clusters!r}")
        elif self._identity_threshold != -1:
            args.append(f"identity_threshold={self._identity_threshold!r}")
        if self._platform != _BEST_PLATFORM:
            args.append(f"platform={self.platform!r}")
        return f"{ty}({', '.join(args)})"

    def __getstate__(self):
        return {
            "platform":            self.platform,
            "clusters":           self._clusters,
            "identity_threshold": self._identity_threshold,
        }

    def __setstate__(self, dict state):
        try:
            BaseTrimmer.__init__(self, platform=state["platform"])
        except (ValueError, RuntimeError):
            BaseTrimmer.__init__(self, platform="detect")
        self._clusters           = state["clusters"]
        self._identity_threshold = state["identity_threshold"]

    # --- Utils --------------------------------------------------------------

    cdef void _configure_manager(self, trimal.manager.trimAlManager* manager):
        manager.automatedMethodCount  = 0
        manager.clusters              = self._clusters
        manager.maxIdentity           = self._identity_threshold


# -- Misc classes ------------------------------------------------------------

cdef class SimilarityMatrix(ScoringMatrix):
    """A similarity matrix for biological sequence characters.

    .. versionchanged:: 0.8.0
       Inherit from the `~scoring_matrices.ScoringMatrix` class.

    """

    DEFAULT_ALPHABET = trimal.aminoAcidResidues.decode('ascii')

    # --- Class methods ------------------------------------------------------

    @classmethod
    def aa(cls):
        """Create a default amino-acid similarity matrix (BLOSUM62).
        """
        blosum62 = ScoringMatrix.from_name("BLOSUM62")
        matrix = list(blosum62.shuffle(SimilarityMatrix.DEFAULT_ALPHABET))
        return cls(matrix, alphabet=SimilarityMatrix.DEFAULT_ALPHABET, name="BLOSUM62")

    @classmethod
    def nt(cls, bool degenerated=False):
        """Create a default nucleotide similarity matrix.

        Arguments:
            degenerated (`bool`): Set to `True` to create a similarity
                matrix for degenerated nucleotides.

        """
        cdef str                                       alphabet
        cdef trimal.similarity_matrix.similarityMatrix sm

        if degenerated:
            alphabet = trimal.degenerateNucleotideResidues.decode('ascii')
            sm.defaultNTDegeneratedSimMatrix()
        else:
            alphabet = trimal.nucleotideResidues.decode('ascii')
            sm.defaultNTSimMatrix()

        matrix = [[0.0 for _ in alphabet] for _ in alphabet]
        for i in range(len(alphabet)):
            for j in range(len(alphabet)):
                matrix[i][j] = sm.simMat[i][j]

        return cls(matrix, alphabet=alphabet)

    # --- Magic methods ------------------------------------------------------

    def __init__(
        self,
        object matrix not None,
        str alphabet not None = SimilarityMatrix.DEFAULT_ALPHABET,
        str name = None,
    ):
        """__init__(self, matrix, alphabet="ARNDCQEGHILKMFPSTWYVBZX*", name=None)\n--\n

        Create a new similarity matrix from the given alphabet and data.

        Arguments:
            matrix (`~numpy.typing.ArrayLike`): The similarity matrix,
                as a square matrix indexed by the alphabet characters.
            alphabet (`str`): The alphabet used for indexing the rows
                and columns of the similarity matrix.
            name (`str` or `None`): The name of the scoring matrix, if any.

        Example:
            Create a new similarity matrix using the HOXD70 scores by
            Chiaromonte, Yap and Miller (:pmid:`11928468`)::

                >>> matrix = SimilarityMatrix(
                ...     [[  91, -114,  -31, -123],
                ...      [-114,  100, -125,  -31],
                ...      [ -31, -125,  100, -114],
                ...      [-123,  -31, -114,   91]],
                ...     alphabet="ATCG",
                ...     name="HOXD70",
                ... )

            Create a new similarity matrix using one of the matrices from
            the `Bio.Align.substitution_matrices` module::

                >>> jones = Bio.Align.substitution_matrices.load('JONES')
                >>> matrix = SimilarityMatrix(jones, jones.alphabet, 'JONES')

        .. versionadded:: 0.1.2

        """
        cdef ssize_t i
        cdef ssize_t j
        cdef ssize_t k
        cdef float  total

        super().__init__(matrix, alphabet=alphabet, name=name)

        # check alphabet constraints
        if not self.alphabet.isupper():
            raise ValueError("Alphabet must only contain uppercase letters")
        if len(self.alphabet) > 28:
            raise ValueError(f"Cannot use alphabet of more than 28 symbols: {alphabet!r}")

        # allocate memory
        self._smx.memoryAllocation(self._size)

        # create the hashing vector with support for all ASCII codes
        for i, letter in enumerate(self.alphabet):
            j = ord(letter) - ord('A')
            if j < 0:
                raise ValueError(f"Invalid symbol in alphabet: {letter!r}")
            self._smx.vhash[ord(letter) - ord('A')] = i

        # copy the similarity matrix
        for i in range(self._size):
            memcpy(self._smx.simMat[i], self._matrix[i], sizeof(float) * self._size)

        # calculate Euclidean distance
        with nogil:
            for j in range(self._size):
                for i in range(j+1, self._size):
                    total = 0
                    for k in range(self._size):
                        total += (
                            (self._smx.simMat[k][j] - self._smx.simMat[k][i])
                          * (self._smx.simMat[k][j] - self._smx.simMat[k][i])
                        )
                    self._smx.distMat[i][j] = self._smx.distMat[j][i] = sqrt(total)

    # --- Functions ----------------------------------------------------------

    cpdef float similarity(self, str a, str b) except -1:
        """Return the similarity between two sequence characters.

        Example:
            >>> mx = SimilarityMatrix.nt()
            >>> mx.similarity('A', 'A')
            1.0
            >>> mx.similarity('A', 'T')
            0.0

        Raises:
            `ValueError`: When ``a`` or ``b`` is an invalid character or a
                character that was not defined in the matrix alphabet.
            `TypeError`: When ``a`` or ``b`` is a string containing more
                than one character.

        .. versionadded:: 0.1.2

        """
        cdef int ia   = ord(a)
        cdef int ib   = ord(b)

        if ia < ord('A') or ia > ord('Z'):
            raise ValueError(f"the symbol {a!r} is incorrect")
        if ib < ord('A') or ib > ord('Z'):
            raise ValueError(f"the symbol {b!r} is incorrect")

        cdef int numa = self._smx.vhash[ia - ord('A')]
        cdef int numb = self._smx.vhash[ib - ord('A')]

        if numa == -1:
            raise ValueError(f"the symbol {a!r} accesing the matrix is not defined in this object")
        if numb == -1:
            raise ValueError(f"the symbol {b!r} accesing the matrix is not defined in this object")

        return self._smx.simMat[numa][numb]

    cpdef float distance(self, str a, str b) except -1:
        """Return the distance between two sequence characters.

        Example:
            >>> mx = SimilarityMatrix.nt(degenerated=True)
            >>> mx.distance('A', 'A')
            0.0
            >>> mx.distance('A', 'T')
            1.5184...

        Raises:
            `ValueError`: When ``a`` or ``b`` is an invalid character or a
                character that was not defined in the matrix alphabet.
            `TypeError`: When ``a`` or ``b`` is a string containing more
                than one character.

        """
        cdef float distance = 0.0
        cdef char  a_code   = ord(a)
        cdef char  b_code   = ord(b)
        with nogil:
            distance  = self._smx.getDistance(a_code, b_code)
        return distance
