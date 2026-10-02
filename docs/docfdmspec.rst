FDM spectrum
============

The ``fdmSpectrum`` class computes spectral functions from the reduced density matrix and the kept NRG
states. It is the main entry point for dynamical correlation functions in the full density matrix
framework.

How the spectrum is assembled
-----------------------------

An ``fdmSpectrum`` instance refers to an ``nrgcore`` object and accumulates contributions as its Wilson
shells are visited. Before processing shells, ``setOperator`` registers a required operator set ``B`` and an
optional corresponding set ``A``. If ``A`` is omitted, the adjoint of each ``B`` block is used. The operator
vectors are held by pointer, so they must remain alive while the spectrum is being calculated.

For each shell, ``calcSpectrum`` updates the kept-state indices and constructs a sector density matrix. On the
final shell, weights are initialized from the shell ground-state degeneracy; on later shells, the reduced
density matrix from the next shell is carried on kept states. The operator contraction evaluates the two
operator-density-matrix orderings for each pair of sectors and bins their energy differences, scaled by the
provided ``energyScale``, into separate positive- and negative-frequency arrays. Small contributions below the
internal weight tolerance are skipped.

After the contraction, the full density matrix is rotated into the coupled system-bath basis and traced over
the bath to produce the reduced matrix used by the next shell. Once all shells have contributed,
``saveFinalData`` writes the logarithmic energy grid and the positive and negative weight arrays through the
provided file object, then clears the accumulated weights. The internal grid defaults to 100,000 points over
the energy range from ``1e-10`` to ``10``.

For the detailed generated API reference, see the ``fdmSpectrum`` `class documentation <api/classfdmSpectrum.html>`_.


