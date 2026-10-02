Core NRG solver
===============

The ``nrgcore`` class is the iterative solver at the center of an NRG calculation. It combines an impurity model
with a sequence of bath sites, diagonalizes the growing system in quantum-number sectors, and retains the
low-energy states needed for the next iteration. The impurity and bath model classes provide their own basis,
energies, quantum numbers, and creation operators; ``nrgcore`` uses these to construct and solve each enlarged
system Hamiltonian.

The solver does not choose the bath discretization or calculate the hopping schedule. The caller supplies the
hopping amplitudes and energy rescaling factor for each site through :cpp:func:`add_bath_site`.

Iteration overview
------------------

An iteration grows the current system by one bath site. Internally, the calculation follows these steps:

1. **Initialize the system.** The constructor stores references to the impurity and bath models and checks that
	they have matching quantum-number dimensions and the same number of creation-operator channels. The first
	call to ``add_bath_site`` takes the impurity's energies, basis, and creation operators as the starting system.
	The model objects must remain alive for as long as the solver uses them.

2. **Build the enlarged basis.** The solver combines each retained system sector with each bath sector. Their
	quantum-number vectors are added, and combinations with the same resulting quantum numbers are grouped into
	a single sector. The grouped indices record which system-bath basis states belong to each block.

3. **Assemble one Hamiltonian per sector.** The diagonal entries combine the rescaled system energies with the
	bath-site energies. Off-diagonal entries couple system and bath states using their creation operators, the
	bath ``chi`` factors, and the corresponding hopping amplitude. The hopping vector must have one entry for each
	operator channel. The ``rescale`` argument multiplies the existing system energies when the new Hamiltonian
	is assembled.

4. **Diagonalize the blocks.** Each sector Hamiltonian is diagonalized independently. The resulting energies
	are stored in ``eigenvaluesQ``, with one energy vector per quantum-number sector.

5. **Update operators for the next iteration.** The system creation operators are transformed into the newly
	diagonalized basis and combined with the new bath-site operators. These transformed operators are then used
	to construct the hopping terms when the following site is added. This update happens before truncation
	because it needs the current basis and kept-state indices.

6. **Truncate and advance.** ``update_internal_state`` advances the iteration counter, retains low-energy states
	across all sectors, and makes the current quantum-number basis the previous system basis for the next call.
	``max_kept_states`` is a target for the total number of states across sectors, not a separate limit per
	sector. An energy tolerance is applied at the cutoff, so near-degenerate states can also be retained and the
	actual count can exceed the target. When states are discarded, the retained energies are shifted relative to
	the ground-state energy.

The first site follows the same sequence, with the impurity as the initial system. The impurity is initially
truncated if necessary before it is combined with that first bath site.

Caller workflow
---------------

The solver separates building a new iteration from committing it as the previous system. Call
``update_internal_state`` after every ``add_bath_site``; it is not called automatically by ``add_bath_site``.
The SIAM example follows this pattern:

.. code-block:: cpp

	nrgcore<spinhalf, spinhalf> siam(impurity, bathModel);
	siam.set_parameters(1024);

	siam.add_bath_site({V, V}, 1.0);
	siam.update_internal_state();

	for (int site = 0; site < nMax; ++site) {
	  siam.add_bath_site({hopping(site, Lambda), hopping(site, Lambda)}, rescale);
	  siam.update_internal_state();
	}

``set_parameters`` sets the total kept-state target and resets the iteration counter. Its default target is
1024 states. Set it before beginning the iteration sequence. The hopping vector has one value per creation
operator channel; for the two-channel spinful example above, it contains two values.

State and results
-----------------

The solver exposes energies and basis information by sector as well as several public state vectors used by
other library components. The most useful read accessors are:

* ``get_basis_nQ()`` returns the current sector quantum numbers.
* ``get_eigenvaluesQ()`` returns the energies grouped by sector.
* ``get_f_dag_operator()`` returns the transformed creation operators used for the next iteration.
* ``all_eigenvalue`` contains the combined spectrum from all sectors, sorted in ascending order before the
	cutoff is applied. It therefore also includes energies of states that are subsequently discarded.

The returned data describes the current iteration. In particular, call ``update_internal_state`` before adding
the next site so the truncation and basis transition have taken place.

For generated implementation details and the member reference, see the ``nrgcore`` `class documentation <api/classnrgcore.html>`_.


