Spin-half impurity model
========================

The ``spinhalf`` model is the canonical single-impurity fermionic model used for spinful Anderson-type
and Kondo-type calculations. It defines the local Hilbert space, quantum-number sectors, and the
fermionic operators needed by the NRG solver.

How the local model is built
----------------------------

The constructor creates a two-orbital ``fermionBasis`` for spin up and spin down, giving the four local
occupation states. It forms the number operators from each creation operator and its adjoint, then assembles
the local Hamiltonian

.. math::

   H = \epsilon_d(n_\uparrow + n_\downarrow) + h(n_\uparrow - n_\downarrow)
	   + U n_\uparrow n_\downarrow.

The Hamiltonian is converted to quantum-number blocks and diagonalized sector by sector. The model exposes
the resulting sector labels and energies, the spin-resolved creation operators, and the fermion-parity factor
``chi_Q`` used when operators from neighboring sites are combined. It also constructs the double-occupancy
operator ``n_up n_down`` for correlation and response calculations.

``spinhalf`` describes a single local site only: it does not add Wilson bath sites or perform NRG truncation.
Those responsibilities belong to ``nrgcore``, which consumes the model's basis, energies, and operators.

For the detailed API and member list, see the ``spinhalf`` `class documentation <api/classspinhalf.html>`_.
