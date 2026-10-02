Quantum operators
=================

``qOperator`` represents quantum-number-resolved operator blocks in a compact, symmetry-aware form.
It is used throughout the NRG implementation to manage block structure, basis transformations, and
sector-dependent operator data.

How operator blocks are represented
-----------------------------------

Internally, a ``qOperator`` is a map from a pair of sector indices ``(i, j)`` to a ``qmatrix``. Each stored
matrix contains the operator elements connecting states in sector ``j`` to states in sector ``i``. Blocks
that are absent from the map are treated as structurally zero; ``get(i, j)`` reports an optional pointer so
callers can distinguish a missing block without allocating a dense matrix for it.

The ``set`` overloads insert a block by copying or moving its matrix and reject duplicate sector pairs.
``unitaryTransform`` rotates each stored block with the corresponding sector eigenvector matrices, applying
the left-sector adjoint and right-sector rotation. This lets NRG code update an operator after diagonalizing
the Hamiltonian without expanding it into the full Hilbert space. ``clear`` removes all blocks, and
``display`` prints them for inspection.

For the generated class reference, see the ``qOperator`` `class documentation <api/classqOperator.html>`_.


