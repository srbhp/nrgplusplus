Fermion basis
=============

The ``fermionBasis`` class defines the fermionic basis and the quantum-number sectors used to build and
block-diagonalize impurity Hamiltonians. It provides the symmetry structure used throughout the NRG
solver setup.

How the basis is constructed
----------------------------

The constructor starts from ``dof`` fermionic orbitals and builds their many-body occupation basis. Each
orbital's creation matrix is assembled with Kronecker products; the ``sigz`` factors supply the signs needed
for fermionic anticommutation between orbitals. Applying each creation operator and its adjoint gives the
occupation values stored for every basis state.

The selected ``modelSymmetry`` then determines the quantum-number vector assigned to each state:

* ``chargeOnly`` groups states by total particle number.
* ``spinOnly`` groups by the difference between the configured up- and down-orbital occupations.
* ``chargeAndSpin`` keeps the up- and down-orbital occupation counts as separate quantum numbers.

States with identical quantum-number vectors are collected into blocks. The class uses those block indices to
convert full-space operators into ``qOperator`` objects containing only nonzero sector-to-sector matrices.
``set_f_dag_operators`` applies this conversion to the fermion creation operators. The same machinery is
available through ``get_block_operators`` for other operators. ``get_block_Hamiltonian`` extracts diagonal
sector blocks and checks that the input Hamiltonian is Hermitian and has no couplings between distinct
quantum-number sectors.

This class constructs basis and operator structure; the impurity or bath model builds its Hamiltonian and
passes it through these block operations.

For the generated class reference, see the ``fermionBasis`` `class documentation <api/classfermionBasis.html>`_.
