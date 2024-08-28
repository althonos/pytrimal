Statistics
==========

*This section is partially adapted from the supplementary material of the* `trimAl manual <http://trimal.cgenomics.org/_media/manual.b.pdf>`_.

trimAl provides several methods for trimming an alignment, using one or more
statistics computed from the alignment.

Introduction
------------

Consider an alphabet :math:`V`, which is usually the protein alphabet.
We define the gapped alphabet :math:`V^\gamma = V \cup \{ \Gamma \}` where
:math:`\Gamma` is the gap symbol, and the indeterminate alphabet
:math:`V_\xi = V \cup \{ X \}` where :math:`X` is the indetermination symbol.

A multiple sequence alignment composed of :math:`m` sequences and :math:`n`
residues can be defined as the matrix:

.. math::

    A = (a_{i,j}) \in  (V_\xi^\gamma) ^{m \times n}

Gap
---

Definition
^^^^^^^^^^

The *Gap* statistic is defined for each column $k$ as the fraction of sequences
without a gap character in column $k$:

.. math::

    Gap_k(A) = \frac{1}{m} \sum_{i=1}^m \mathbb{1}_{V_\xi}(a_{i,k}) =
     \frac{1}{m} \sum_{i=1}^{m} \begin{cases} 0 & \text{if $a_{i,k} = \Gamma$}  \\ 1 & \text{otherwise} \end{cases}

Performance
^^^^^^^^^^^

The *Gap* statistic is the easiest and fastest one to compute for a given
alignment, with a runtime complexity of :math:`O(mn)`, and a memory complexity
of :math:`O(1)`. It is the only statistic that scales linearly with the number
of sequences, since it does not require any pairwise comparisons to be made.

Similarity
----------

Definition
^^^^^^^^^^

.. rubric:: Scoring Matrix

The similarity statistic uses a scoring matrix (such as BLOSUM62) to compute the
similarity between every sequence pair in the alignment. A scoring matrix can
be defined as a function mapping two symbols of :math:`V` to an arbitrary score
in :math:`\mathbb{R}`. Based on this scoring matrix, the distance between
two symbols :math:`x` and :math:`y` is computed as the euclidean distance:

.. math::

    D_S(x,y) = \sqrt { \sum_{z \in V}{(S(z, y) - S(z, x))^2} }

.. rubric:: Mismatch

Given two sequences :math:`i` and :math:`j`, the alignment mismatch :math:`W_{i,j}`
is defined as the fraction of non-identical residues in the pairwise alignment
between the two sequences:

.. math::

    W_{i,j}(A) = 1 - \frac{
      \sum_{k=1}^{n}
        \begin{cases}
        1 & \text{if $a_{i,k} = a_{j,k}$, $a_{i,k} \in V$, $a_{j,k} \in V$} \\
        0 & \text{otherwise}
        \end{cases}
    }{
      \sum_{k=1}^n
        \begin{cases}
        1 & \text{if $a_{i,k} \in V$ or $a_{j,k} \in V$} \\
        0 & \text{otherwise}
      \end{cases}
    }

.. rubric:: Similarity

The mean of the pairwise sequence distances weighted by the pairwise sequence
mismatch is then calculated, only taking into account positions where
both sequences have residue:

.. math::

    Q_{S,k}(A) = \frac
      { \sum_{i=1}^{m}{ \sum_{j=i+1}^{m}{ \mathbb{1}_{V}(a_{i,k}) \mathbb{1}_{V}(a_{j,k}) W_{i,j}(A) D_S(a_{i,k}, a_{j,k})}}}
      { \sum_{i=1}^{m}{ \sum_{j=i+1}^{m}{ \mathbb{1}_{V}(a_{i,k}) \mathbb{1}_{V}(a_{j,k}) W_{i,j}(A) }}}

Finally, the similarity score for a column is only computed if a column contains
gaps in 80% of the sequences or more. This can be obtained from the Gap statistic:

.. math::

    Sim_{S,k}(A) = \begin{cases}
    0 & \text{if $Gap_k(A) < 0.2$} \\
    e^{-Q_{S,k}(A)} & \text{otherwise}
    \end{cases}


Performance
^^^^^^^^^^^

The *Similarity* statistic has a runtime complexity of :math:`O(nm^2)` and a
memory complexity of :math:`O(m^2)`, as it requires an external matrix to
store the pairwise sequence identity for efficient computation.

An optimized implementation using SIMD is available for the SSE backend. It uses
additional column vectors to store intermediate results for a single column before
when computing the identity statistic, bringing the memory complexity to
:math:`O(m^2 + n)`.


Overlap
-------

Definition
^^^^^^^^^^

.. rubric:: Hit Ratio

For each sequence :math:`i` of the alignment and every position $k$, we can compute
the hit ratio as the fraction of sequences that agree with $i$ at the position $k$:

.. math::

    Hit_{i,k}(A) = \frac{1}{m-1} \sum_{\begin{align}j = 1 \\ j \ne i\end{align}}^{m}{
        \begin{cases}
        1 & \text{if $a_{i,k} = a_{j,k}$} \\
        1 & \text{if $a_{i,k} \in V$, $a_{j,k} \in V$} \\
        0 & \text{otherwise}
        \end{cases}
    }

.. rubric:: Sequence Overlap

Given a residue threshold :math:`r \in [0, 1]`, the overlap for a
sequence :math:`i` is the fraction of positions with a hit ratio above
:math:`r`:

.. math::

    Ovl_i(A) = \frac{1}{n} \sum_{k=1}^{n}{
      \begin{cases}
        1 & \text{if $Hit_{i,k}(A) > r$} \\
        0 & \text{otherwise}
      \end{cases}
    }

Performance
^^^^^^^^^^^

The *Overlap* statistic has a runtime complexity of :math:`O(nm^2)` and a
memory complexity of :math:`O(1)`, however the trimAl implementation uses
an additional vector to store the overlap value for each sequence, bringing
the memory complexity to :math:`O(m)`.

An optimized implementation using SIMD is available for the SSE backend. It
computes the hit ratio of all positions for every sequence, allowing some
optimizations in data access and processing, but requiring additional memory
to store intermediate column data, bringing the memory complexity to :math:`O(m+n)`.
