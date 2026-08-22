# Code Review TODOs

## High priority

- [ ] Restrict `.github/workflows/cmake.yml:6-7,17-20,66-72` so pull-request builds cannot write to or push `main`; use read-only permissions for PR builds.
- [ ] Fix `nrgcore/include/utils/qmatrix.hpp:361-369` to pass complex scalar values to `cblas_zgemm` instead of `double*`.
- [ ] Initialize `aMatrix` correctly in `nrgcore/include/dynamics/fdmSpectrum.hpp:221-229` when an annihilation operator is supplied.
- [ ] Correct `nrgcore/include/utils/sparseSolver.hpp:75,181` so ground-state solvers select the smallest algebraic eigenvalue rather than the eigenvalue closest to zero.

## Medium priority

- [ ] Fix the square-matrix constructor in `nrgcore/include/utils/qmatrix.hpp:159`; it currently constructs and discards a temporary.
- [ ] Correct the truncation boundary in `nrgcore/include/nrgcore/nrgcore.hpp:285-299` so non-degenerate spectra retain exactly `max_kept_states`.
- [ ] Update block-diagonal validation in `nrgcore/include/models/fermionBasis.hpp:357-367` to detect cancellation-hidden off-block elements.
- [ ] Fix ground-state initialization for all-positive spectra in `nrgcore/include/dynamics/fdmSpectrum.hpp:122-149` and `nrgcore/include/dynamics/fdmBackwardIteration.hpp:126-153`.
- [ ] Make `setTemperature()` affect spectral calculations in `nrgcore/include/dynamics/fdmSpectrum.hpp:27,284`, or remove the ineffective public API.
- [ ] Implement RAII cleanup for temporary files in `nrgcore/include/nrgcore/nrgData.hpp:99-128`, including cleanup after `close()`.
- [ ] Make HDF5 open failures explicit in `nrgcore/include/utils/h5stream.hpp:276-293` by validating modes and propagating descriptive errors.
- [ ] Fix the malformed `wget` pipeline in `scripts/pre_install.sh:6` and enable fail-fast behavior.
- [ ] Make sanitizers, debug symbols, and `-O0` optional in `CMakeLists.txt:14` instead of forcing them for every build.
- [ ] Make the test target buildable from a normal checkout by fixing paths/compiler assumptions in `test/spinBasis.cpp:1-2` and `test/makefile:2`; integrate it with CMake/CTest.
- [ ] Align Doxygen configuration and CMake outputs in `docs/CMakeLists.txt:8-28` and `docs/Doxyfile.in:994,2301`; enable XML output and use source-tree input paths.
