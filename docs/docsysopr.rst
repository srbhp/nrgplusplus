System-operator updates
=======================

The ``update_system_operator`` helper rebuilds the system operators used by the NRG solver after each
iteration, keeping the operator basis consistent with the current retained symmetry sectors.

How an operator is carried to the next iteration
-------------------------------------------------

The helper takes the current ``nrgcore`` state and a vector of operators expressed in the previous system
basis. For each pair of new quantum-number sectors, it walks the grouped system-bath basis indices and copies
the old operator matrix elements only between basis states with the same bath-sector index. This is the
identity action on the newly added bath site, while retaining the operator's action on the old system.

It then rotates each assembled block into the new eigenbasis using the sector eigenvector matrices, applying
the left-sector adjoint and right-sector transformation. The input vector is replaced by these transformed
``qOperator`` blocks, ready for use in later iterations or observable calculations. The function relies on
the current coupled-sector map, kept-state indices, bath dimensions, and diagonalization matrices from
``nrgcore``; it does not construct or evolve the operators independently of that state.

For the generated function reference, see ``update_system_operator`` in the `API documentation <api/function_sysOperator_8hpp_1a83ad7b6117e9ef905244c84064f1f616.html>`_.

