# distutils: language = c++
# cython: language_level=3, linetrace=True
"""Bindings to trimAl, a tool for automated alignment trimming.

References:
    - Capella-Gutiérrez, Salvador, José M. Silla-Martínez, and Toni Gabaldón.
      *TrimAl: A Tool for Automated Alignment Trimming in Large-Scale
      Phylogenetic Analyses*. Bioinformatics 25, no. 15 (2009): 1972–73.
      :doi:`10.1093/bioinformatics/btp348`.

"""

# --- C imports --------------------------------------------------------------

cimport cython
from _unicode cimport (
    PyUnicode_New,
    PyUnicode_READY,
    PyUnicode_KIND,
    PyUnicode_DATA,
    PyUnicode_WRITE,
)
from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_FORMAT, PyBUF_READ
from cpython.bytes cimport PyBytes_FromStringAndSize, PyBytes_AsString
from cpython.list cimport PyList_New, PyList_SET_ITEM
from cpython.mem cimport PyMem_Free, PyMem_Malloc
from cpython.ref cimport Py_INCREF

from libc.math cimport NAN, isnan, sqrt
from libc.stdio cimport printf
from libcpp cimport bool
from libcpp.string cimport string

cimport trimal
cimport trimal.alignment
cimport trimal.format_manager
cimport trimal.manager
cimport trimal.report_system
cimport trimal.similarity_matrix

IF SSE2_BUILD_SUPPORT:
    from pytrimal.impl.sse cimport SSESimilarity, SSECleaner

# --- Python imports ---------------------------------------------------------

import os
import threading

# --- Constants --------------------------------------------------------------

cdef set AUTOMATED_TRIMMER_METHODS = {
    "strict",
    "strictplus",
    "gappyout",
    "nogaps",
    "noallgaps",
    "automated1",
}

_TARGET_CPU           = TARGET_CPU
_SSE2_RUNTIME_SUPPORT = False
_SSE2_BUILD_SUPPORT   = False

IF TARGET_CPU == "x86" and TARGET_SYSTEM in ("freebsd", "linux_or_android", "macos", "windows"):
    from cpu_features.x86 cimport GetX86Info, X86Info
    cdef X86Info cpu_info = GetX86Info()
    _SSE2_BUILD_SUPPORT   = SSE2_BUILD_SUPPORT
    _AVX2_BUILD_SUPPORT   = AVX2_BUILD_SUPPORT
    _SSE2_RUNTIME_SUPPORT = cpu_info.features.sse2 != 0
    _AVX2_RUNTIME_SUPPORT = cpu_info.features.avx2 != 0


# --- Utilities --------------------------------------------------------------

cdef float _check_range(object value, str name, float min_value, float max_value) except *:
    if value < min_value or value > max_value or isnan(value):
        raise ValueError(f"Invalid value for `{name}`: {value!r}")
    return value

cdef extern from *:
    """
    template <typename T>
    T* new_array(size_t n) {
        return new T[n];
    }
    template <typename T>
    void del_array(T* array) {
        delete array;
    }
    """
    T* new_array[T](size_t)
    void del_array[T](T*)


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
        >>> sum(letter == '-' for seq in msa.sequences for letter in seq)
        43

    """

    def __cinit__(self, Alignment alignment):
        self._owner = alignment
        self._ali = alignment._ali
        self._index_mapping = alignment._sequences_mapping

    def __len__(self):
        assert self._ali is not NULL
        return self._ali.numberOfSequences

    def __getitem__(self, int index):
        assert self._ali is not NULL

        cdef int index_ = index
        cdef int length = self._ali.numberOfSequences

        if index_ < 0:
            index_ += length
        if index_ < 0 or index_ >= length:
            raise IndexError(index)
        if self._index_mapping is not NULL:
            index_ = self._index_mapping[index_]

        IF SYS_VERSION_INFO_MAJOR <= 3 and SYS_VERSION_INFO_MAJOR <= 7 and SYS_IMPLEMENTATION_NAME == "pypy":
            cdef bytes    seq  = PyBytes_FromStringAndSize(NULL, self._ali.numberOfResidues)
            cdef char*    data = PyBytes_AsString(seq)
            cdef size_t   x    = 0
            for i in range(self._ali.originalNumberOfResidues):
                if self._ali.saveResidues is NULL or self._ali.saveResidues[i] != -1:
                    data[x] = self._ali.sequences[index_][i]
                    x += 1
            return seq.decode('ascii')
        ELSE:
            cdef str      seq  = PyUnicode_New(self._ali.numberOfResidues, 0x7f)
            IF SYS_VERSION_INFO_MAJOR <= 3 and SYS_VERSION_INFO_MINOR < 12:
                PyUnicode_READY(seq)
            cdef void*    data = PyUnicode_DATA(seq)
            cdef int      kind = PyUnicode_KIND(seq)
            cdef size_t   x    = 0
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

    """

    def __cinit__(self, Alignment alignment):
        self._owner = alignment
        self._ali = alignment._ali
        self._index_mapping = alignment._residues_mapping

    def __len__(self):
        assert self._ali is not NULL
        return self._ali.numberOfResidues

    def __getitem__(self, int index):
        assert self._ali is not NULL

        cdef int index_ = index
        cdef int length = self._ali.numberOfResidues

        if index_ < 0:
            index_ += length
        if index_ < 0 or index_ >= length:
            raise IndexError(index)
        if self._index_mapping is not NULL:
            index_ = self._index_mapping[index_]

        IF SYS_VERSION_INFO_MAJOR <= 3 and SYS_VERSION_INFO_MAJOR <= 7 and SYS_IMPLEMENTATION_NAME == "pypy":
            cdef bytes    col  = PyBytes_FromStringAndSize(NULL, self._ali.numberOfSequences)
            cdef char*    data = PyBytes_AsString(col)
            cdef size_t   x    = 0
            for i in range(self._ali.originalNumberOfSequences):
                if self._ali.saveSequences is NULL or self._ali.saveSequences[i] != -1:
                    data[x] = self._ali.sequences[i][index_]
                    x += 1
            return col.decode('ascii')
        ELSE:
            cdef str      col = PyUnicode_New(self._ali.numberOfSequences, 0x7f)
            IF SYS_VERSION_INFO_MAJOR <= 3 and SYS_VERSION_INFO_MINOR < 12:
                PyUnicode_READY(col)
            cdef void*    data = PyUnicode_DATA(col)
            cdef int      kind = PyUnicode_KIND(col)
            cdef size_t   x    = 0
            for i in range(self._ali.originalNumberOfSequences):
                if self._ali.saveSequences is NULL or self._ali.saveSequences[i] != -1:
                    PyUnicode_WRITE(kind, data, x, self._ali.sequences[i][index_])
                    x += 1
            return col


cdef class Alignment:
    """A multiple sequence alignment.
    """

    @classmethod
    def load(cls, object path not None):
        """load(cls, path)\n--

        Load a multiple sequence alignment from a file.

        Arguments:
            path (`str`, `bytes` or `os.PathLike`): The path to the file
              containing the serialized alignment to load.

        Returns:
            `~pytrimal.Alignment`: The deserialized alignment.

        Example:
            >>> msa = Alignment.load("example.001.AA.clw")
            >>> msa.names
            [b'Sp8', b'Sp10', b'Sp26', b'Sp6', b'Sp17', b'Sp33']

        """
        cdef trimal.format_manager.FormatManager manager
        cdef Alignment alignment = Alignment.__new__(Alignment)
        cdef string    path_ =  os.fsencode(path)

        if not os.path.exists(path):
            raise FileNotFoundError(path)
        elif os.path.isdir(path):
            raise IsADirectoryError(path)

        alignment._ali = manager.loadAlignment(path_)
        return alignment

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

    def __init__(self, object names, object sequences):
        """__init__(self, names, sequences)\n--

        Create a new alignment with the given names and sequences.

        Arguments:
            names (`~collections.abc.Sequence` of `bytes`): The names of
                the sequences in the alignment.
            sequences (`~collections.abc.Sequence` of `str`): The actual
                sequences in the alignment.

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
        cdef bytes name
        cdef str   sequence
        cdef int   nresidues = -1

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
            self._ali.sequences[i] = sequence.encode('ascii') # FIXME: no decoding

        self._ali.fillMatrices(True, True)
        self._ali.originalNumberOfSequences = self._ali.numberOfSequences
        self._ali.originalNumberOfResidues = self._ali.numberOfResidues

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}(names={self.names!r}, sequences={list(self.sequences)!r})"

    def __copy__(self):
        return self.copy()

    @property
    def names(self):
        """sequence of `bytes`: The names of the sequences in the alignment.
        """
        assert self._ali is not NULL

        cdef int    i
        cdef int    x     = 0
        cdef bytes  name
        cdef object names = PyList_New(self._ali.numberOfSequences)

        for i in range(self._ali.originalNumberOfSequences):
            if self._ali.saveSequences is NULL or self._ali.saveSequences[i] != -1:
                name = PyBytes_FromStringAndSize( self._ali.seqsName[i].data(), self._ali.seqsName[i].size() )
                PyList_SET_ITEM(names, x, name)
                Py_INCREF(name)  # manually increase reference count because `PyList_SET_ITEM` doesn't
                x += 1

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

    cpdef Alignment copy(self):
        """copy(self)\n--

        Create a copy of this alignment.

        """
        assert self._ali is not NULL
        cdef Alignment copy = (type(self)).__new__(type(self))
        copy._ali = new trimal.alignment.Alignment(self._ali[0])
        return copy


cdef class TrimmedAlignment(Alignment):
    """A multiple sequence alignment that has been trimmed.

    Internally, the trimming process only produces a mask of sequences and
    a mask of residues. This class exposes the filtered sequences and
    residues.

    """

    @classmethod
    def load(cls, object path not None):
        # For compatibility, allow loading a trimmed alignment from a file
        # even though it makes it effectively not trimmed
        cdef Alignment alignment = Alignment.load(path)
        cdef TrimmedAlignment trimmed = TrimmedAlignment.__new__(TrimmedAlignment)
        trimmed._ali = alignment._ali
        alignment._ali = NULL
        trimmed._build_index_mapping()
        return trimmed

    def __init__(
        self,
        object names,
        object sequences,
        object sequences_mask = None,
        object residues_mask = None,
    ):
        """__init__(self, names, sequences, sequences_mask=None, residues_mask=None)\n--

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

    cpdef Alignment original_alignment(self):
        """original_alignment(self)\n--

        Rebuild the original alignment from which this object was obtained.

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
        """terminal_only(self)\n--

        Get a trimmed alignment where only the terminal residues are removed.

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
        """copy(self)\n--

        Create a copy of this trimmed alignment.

        """
        cdef TrimmedAlignment copy = TrimmedAlignment.__new__(TrimmedAlignment)
        copy._ali = new trimal.alignment.Alignment(self._ali[0])
        copy._build_index_mapping()
        return copy

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


# -- Trimmer classes ---------------------------------------------------------

cdef enum simd_backend:
    NONE = 0
    SSE2 = 1


cdef class BaseTrimmer:
    """A sequence alignment trimmer.

    All subclasses provide the same `trim` method, and are configured
    through their constructor.

    """

    def __init__(self, *, backend="detect"):
        """__init__(self, *, backend="detect")\n--

        Create a new base trimmer.

        Arguments:
            backend (`str`): The SIMD extension backend to use to accelerate
                computation of pairwise similarity statistics.

        .. versionadded:: 0.2.0
           The ``backend`` keyword argument.

        """
        IF TARGET_CPU == "x86":
            if backend =="detect":
                self._backend = simd_backend.NONE
                IF SSE2_BUILD_SUPPORT:
                    if _SSE2_RUNTIME_SUPPORT:
                        self._backend = simd_backend.SSE2
            elif backend == "sse":
                IF not SSE2_BUILD_SUPPORT:
                    raise RuntimeError("Extension was compiled without SSE2 support")
                if not _SSE2_RUNTIME_SUPPORT:
                    raise RuntimeError("Cannot run SSE2 instructions on this machine")
                self._backend = simd_backend.SSE2
            elif backend is None:
                self._backend = simd_backend.NONE
            else:
                raise ValueError(f"Unsupported backend on this architecture: {backend}")
        ELSE:
            if backend == "detect" or backend is None:
                self._backend = simd_backend.NONE
            else:
                raise ValueError(f"Unsupported backend on this architecture: {backend}")

    cdef void _setup_simd_code(self, trimal.manager.trimAlManager* manager):
        IF SSE2_BUILD_SUPPORT:
            if self._backend == simd_backend.SSE2:
                manager.origAlig.Statistics.similarity = new SSESimilarity(manager.origAlig)
                del manager.origAlig.Cleaning
                manager.origAlig.Cleaning = new SSECleaner(manager.origAlig)

    cdef void _configure_manager(self, trimal.manager.trimAlManager* manager):
        pass

    cpdef TrimmedAlignment trim(self, Alignment alignment, SimilarityMatrix matrix = None):
        """trim(self, alignment, matrix=None)\n--

        Trim the provided alignment.

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
        self._setup_simd_code(&manager)

        with nogil:
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
    - ``nogaps``: A naive method that removes every column containing at
      least one gap.
    - ``noallgaps``: A naive method that removes every column containing
      only gaps.

    """

    def __init__(self, str method="strict", *, str backend="detect"):
        """__init__(self, method="strict", *, backend="detect")\n--

        Create a new automatic alignment trimmer using the given method.

        Arguments:
            method (`str`): The automatic aligment trimming method. See
                the documentation for `AutomaticTrimmer` for a list of
                supported values.

        Keyword Arguments:
            backend (`str`): The SIMD extension backend to use to accelerate
                computation of pairwise similarity statistics.

        Raises:
            `ValueError`: When ``method`` is not one of the automatic
                alignment trimming methods supported by trimAl.

        .. versionadded:: 0.2.0
           The ``backend`` keyword argument.

        """
        super().__init__(backend=backend)

        if method not in AUTOMATED_TRIMMER_METHODS:
            raise ValueError(f"Invalid value for `method`: {method!r}")
        self.method = method

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


cdef class ManualTrimmer(BaseTrimmer):
    """A sequence alignment trimmer with manually defined thresholds.

    Manual trimming allows the user to specify independent thresholds for
    four different statistics:

    - *Consistency threshold*: Remove columns with a consistency ratio
      lower than the provided threshold.
    - *Gap threshold*: Remove columns where the gap ratio (or the absolute
      gap count) is higher than the provided threshold.
    - *Similarity threshold*: Remove columns with a similarity ratio lower
      than the provided threshold.

    In addition, the trimming can be restricted so that at least a
    configurable fraction of the original alignment is retained, in order
    to avoid stripping an alignment of distance sequences by aggressive
    trimming.

    """

    def __cinit__(self):
        self._gap_threshold           = -1
        self._gap_absolute_threshold  = -1
        self._similarity_threshold    = -1
        self._consistency_threshold   = -1
        self._conservation_percentage = -1

    def __init__(
        self,
        *,
        object gap_threshold           = None,
        object gap_absolute_threshold  = None,
        object similarity_threshold    = None,
        object consistency_threshold   = None,
        object conservation_percentage = None,
        str backend="detect",
    ):
        """__init__(self, *, gap_threshold=None, gap_absolute_threshold=None, similarity_threshold=None, consistency_threshold=None, conservation_percentage=None, backend="detect")\n--

        Create a new manual alignment trimmer with the given parameters.

        Keyword Arguments:
            gap_threshold (`float`): The minimum fraction of non-gap
                 characters that must be present in a column to keep
                 the column.
            gap_absolute_threshold (`int`, *optional*): The absolute number
                of gaps allowed on a column to keep it in the alignment.
                Incompatible with ``gap_threshold``.
            similarity_threshold (`float`, *optional*): The minimum average
                similarity required.
            consistency_threshold (`float`, *optional*): The minimum
                consistency value required.
            conservation_percentage (`float`, *optional*): The minimum
                percentage of positions in the original alignment to
                conserve.
            backend (`str`): The SIMD extension backend to use to accelerate
                computation of pairwise similarity statistics.

        .. versionadded:: 0.2.0
           The ``backend`` keyword argument.

        """
        super().__init__(backend=backend)

        if gap_threshold is not None and gap_absolute_threshold is not None:
            raise ValueError("Cannot specify both `gap_threshold` and `gap_absolute_threshold`")

        if gap_threshold is not None:
            self._gap_threshold = 1 - _check_range(gap_threshold, "gap_threshold", 0, 1)
        if gap_absolute_threshold is not None:
            if gap_absolute_threshold < 0:
                raise ValueError(f"Invalid value for `gap_absolute_threshold`: {gap_absolute_threshold!r}")
            self._gap_absolute_threshold = gap_absolute_threshold
        if similarity_threshold is not None:
            self._similarity_threshold = _check_range(similarity_threshold, "similarity_threshold", 0, 1)
        if consistency_threshold is not None:
            self._consistency_threshold = _check_range(consistency_threshold, "consistency_threshold", 0, 1)
        if conservation_percentage is not None:
            self._conservation_percentage = _check_range(conservation_percentage, "conservation_percentage", 0, 100)

    cdef void _configure_manager(self, trimal.manager.trimAlManager* manager):
        manager.automatedMethodCount = 0
        manager.gapThreshold = self._gap_threshold
        manager.gapAbsoluteThreshold = self._gap_absolute_threshold
        manager.similarityThreshold = self._similarity_threshold
        manager.consistencyThreshold = self._consistency_threshold
        manager.conservationThreshold = self._conservation_percentage


# -- Misc classes ------------------------------------------------------------


cdef class SimilarityMatrix:
    """A similarity matrix for biological sequence characters.
    """

    @classmethod
    def aa(cls):
        """aa(cls)\n--

        Create a default amino-acid similarity matrix (BLOSUM62).

        """
        cdef SimilarityMatrix matrix = cls.__new__(cls)
        with nogil:
            matrix._smx.defaultAASimMatrix()
        return matrix

    @classmethod
    def nt(cls, bool degenerated=False):
        """nt(cls, degenerated=False)\n--

        Create a default nucleotide similarity matrix.

        Arguments:
            degenerated (`bool`): Set to `True` to create a similarity
                matrix for degenerated nucleotides.

        """
        cdef SimilarityMatrix matrix = cls.__new__(cls)
        with nogil:
            if degenerated:
                matrix._smx.defaultNTSimMatrix()
            else:
                matrix._smx.defaultNTDegeneratedSimMatrix()
        return matrix

    def __init__(self, str alphabet not None, object matrix not None):
        """__init__(alphabet, matrix)\n--

        Create a new similarity matrix from the given alphabet and data.

        Arguments:
            alphabet (`str`): The alphabet of the similarity matrix.
            matrix (`~numpy.typing.ArrayLike`): The similarity matrix,
                as a square matrix indexed by the alphabet characters.

        Example:
            Create a new similarity matrix using the HOXD70 scores by
            Chiaromonte, Yap and Miller (:pmid:`11928468`)::

                >>> matrix = SimilarityMatrix(
                ...     "ATCG",
                ...     [[  91, -114,  -31, -123],
                ...      [-114,  100, -125,  -31],
                ...      [ -31, -125,  100, -114],
                ...      [-123,  -31, -114,   91]]
                ... )

            Create a new similarity matrix using one of the matrices from
            the `Bio.Align.substitution_matrices` module::

                >>> jones = Bio.Align.substitution_matrices.load('JONES')
                >>> matrix = SimilarityMatrix(jones.alphabet, jones)

        .. versionadded:: 0.1.2

        """
        cdef int    i
        cdef int    j
        cdef object row
        cdef float  value
        cdef float  sum
        cdef str    letter
        cdef int    length = len(matrix)

        if len(alphabet) != length:
            raise ValueError("Alphabet must have the same length as matrix")
        if not alphabet.isupper():
            raise ValueError("Alphabet must only contain uppercase letters")
        if not all(len(row) == length for row in matrix):
            raise ValueError("`matrix` must be a square matrix")
        if length > 28:
            raise ValueError(f"Cannot use alphabet of more than 28 symbols: {alphabet!r}")

        # allocate memory
        self._smx.memoryAllocation(length)

        # create the hashing vector with support for all ASCII codes
        for i, letter in enumerate(alphabet):
            j = ord(letter) - ord('A')
            if j < 0:
                raise ValueError(f"Invalid symbol in alphabet: {letter!r}")
            self._smx.vhash[ord(letter) - ord('A')] = i

        # create the similarity matrix
        for i, row in enumerate(matrix):
            for j, value in enumerate(row):
                self._smx.simMat[i][j] = value

        # calculate Euclidean distance
        with nogil:
            for j in range(length):
                for i in range(j+1, length):
                    sum = 0
                    for k in range(length):
                        sum += (
                            (self._smx.simMat[k][j] - self._smx.simMat[k][i])
                          * (self._smx.simMat[k][j] - self._smx.simMat[k][i])
                        )
                    self._smx.distMat[i][j] = self._smx.distMat[j][i] = sqrt(sum)

    def __len__(self):
        return self._smx.numPositions

    IF SYS_IMPLEMENTATION_NAME == "cpython":

        def __getbuffer__(self, Py_buffer* buffer, int flags):
            # setup indexing information
            self._shape[0] = self._smx.numPositions
            self._shape[1] = self._smx.numPositions
            self._suboffsets[0] = 0
            self._suboffsets[1] = -1
            self._strides[0] = sizeof(float*)
            self._strides[1] = sizeof(float)
            # update buffer information
            if flags & PyBUF_FORMAT:
                buffer.format = b"f"
            else:
                buffer.format = NULL
            buffer.buf = self._smx.simMat
            buffer.internal = NULL
            buffer.itemsize = sizeof(float)
            buffer.len = self._shape[0] * self._shape[1] * sizeof(float)
            buffer.ndim = 2
            buffer.obj = self
            buffer.readonly = 1
            buffer.shape = &self._shape[0]
            buffer.suboffsets = &self._suboffsets[0]
            buffer.strides = &self._strides[0]

    cpdef float similarity(self, str a, str b) except -1:
        """similarity(self, a, b)\n--

        Return the similarity between two sequence characters.

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
        """distance(self, a, b)\n--

        Return the distance between two sequence characters.

        Example:
            >>> mx = SimilarityMatrix.nt()
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
