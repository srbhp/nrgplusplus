Backward iteration
==================

The ``fdmBackwardIteration`` class performs the backward-iteration steps used in full-density-matrix
NRG calculations. It propagates reduced density matrices and computes the thermal weights needed for
spectral and dynamical observables.

How backward iteration works
-----------------------------

The helper is attached to an existing ``nrgcore`` object and processes one saved Wilson shell at a time,
starting from the last shell and moving toward the impurity. The NRG object supplies the shell's sector
energies, eigenbasis, kept-state indices, and the mapping between system and bath basis states.

For each shell, the helper updates its kept-state index sets and assembles a density matrix for every quantum-
number sector. On the final shell, the no-argument density-matrix path assigns normalized weight to states
within ``energyErrorBar`` of the shell ground state. On earlier shells, it places the supplied Boltzmann
weights on discarded states and carries the reduced density matrix from the next shell on the kept states.
This is how contributions from discarded shells are combined with the states retained for the next backward
step.

To obtain that next-shell density matrix, the helper rotates each shell matrix into the coupled system-bath
basis, then traces over bath states. The result is stored as ``reducedRho`` by system sector. The class also
provides ``rhoDotStaticOperators`` to contract the shell density matrix with sector-diagonal operators and
return their traces. ``clearKeptIndex`` resets the cached matrices, weights, and kept-state history before a
new backward pass.

The helper retains a pointer to the NRG object rather than copying its state, so that object must outlive the
backward-iteration calculation. Temperature and degeneracy tolerance are supplied in the same energy units as
the NRG spectrum.

For the generated class reference and member list, see the ``fdmBackwardIteration`` `class documentation <api/classfdmBackwardIteration.html>`_.

