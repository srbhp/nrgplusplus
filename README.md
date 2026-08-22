# nrgplusplus

[![.github/workflows/cmake.yml](https://github.com/srbhp/nrgplusplus/actions/workflows/cmake.yml/badge.svg?branch=main&event=push)](https://github.com/srbhp/nrgplusplus/actions/workflows/cmake.yml)

**A modern C++ library for high-performance Numerical Renormalization Group (NRG) calculations**

`nrgplusplus` is a implementation of the Numerical Renormalization Group method for solving quantum impurity problems. It provides a flexible, modular architecture for studying strongly correlated electron systems using modern C++ with optimized performance through BLAS/LAPACK libraries.

## ✨ Features

- **Multiple impurity models**: Single Impurity Anderson Model (SIAM), Kondo effect, multi-channel systems, multi-impurity configurations, and superconducting proximity effects
- **Dynamic properties**: Spectral functions, correlation functions, and real-time dynamics via FDM (Full Density Matrix)
- **Thermodynamic calculations**: Entropy, magnetic susceptibility, specific heat for arbitrary temperatures
- **Flexible symmetries**: Support for spin conservation, charge conservation, and arbitrary quantum number blocks
- **High performance**: Template-based C++20 implementation with BLAS/LAPACK optimization
- **Comprehensive examples**: 10+ working examples covering different physical systems
- **HDF5 output**: Native support for saving and analyzing large datasets

## 🚀 Quick Start

### 1. Build the Project

```bash
git clone https://github.com/srbhp/nrgplusplus.git
cd nrgplusplus
mkdir build && cd build
cmake ..
make
```

For detailed build instructions and dependencies, see the [documentation](https://srbhp.github.io/nrgplusplus/build.html).

### 2. Run Your First Example

```bash
cd examples/rgflowSIAM
./rgflowSIAM
python3 plot.py  # Plot the RG flow
```

### 3. Example: Single Impurity Anderson Model (SIAM)

Define the impurity with onsite energy and Coulomb repulsion:

```cpp
#include "nrgcore.hpp"
#include "spinhalf.hpp"

// Create impurity and bath models
spinhalf impurity(eps=-1.0, U_int=2.0);
spinhalf bathModel(0, 0);

// Initialize NRG solver
nrgcore<spinhalf, spinhalf> siam(impurity, bathModel);
siam.set_parameters(1024);  // Keep up to 1024 states per iteration
```

Iteratively add bath sites and save results:

```cpp
h5stream::h5stream output("siam_results.h5");

double Lambda = 2.0;  // RG flow parameter
for (int iteration = 0; iteration < nMax; iteration++) {
  double V = 0.5;  // Impurity-bath coupling
  double rescale = (iteration > 0) ? std::sqrt(Lambda) : 1.0;
  
  siam.add_bath_site({V, V}, rescale);
  siam.update_internal_state();
  
  output.write(siam.all_eigenvalue, "eigenvalues_" + std::to_string(iteration));
}
output.close();
```

See the [documentation](https://srbhp.github.io/nrgplusplus/) for complete API reference and more detailed examples.

## 📚 Available Models & Examples

| Example | Description | Location |
|---------|-------------|----------|
| **RG Flow SIAM** | Renormalization group flow for the Single Impurity Anderson Model | `examples/rgflowSIAM/` |
| **RG Flow for Kondo Model** | RG analysis of the Kondo effect with spin-flip scattering | `examples/rgflowKondo/` |
| **Entropy of SIAM** | Temperature-dependent entropy for Anderson impurity | `examples/entropySIAM/` |
| **Entropy Kondo Model** | Temperature evolution of the Kondo entropy | `examples/entropyKondo/` |
| **Spectrum  $A(\omega)$ for SIAM** | Spectral function calculation via Frequency Domain Mode | `examples/fdmSpectrumSiam/` |
| **Spectrum for Two-Channel SIAM** | Multi-channel system spectral properties | `examples/fdmSpectrumTwoChannel/` |
| **Resonant Level** | Non-interacting resonant level model | `examples/resonantLevel/` |
| **Two-Channel SIAM** | Multi-channel Anderson impurity system | `examples/twoChannelSiam/` |
| **Rabi & Anderson** | Rabi oscillations coupled to Anderson impurity | `examples/rabiAnderson/` |
| **Free Model** | Non-interacting model for testing | `examples/freeModel/` |

Each example includes:
- C++ source code (`main.cpp`)
- CMake build configuration
- Python analysis scripts (where applicable)
- Sample output files and detailed README

## 🏗️ Project Architecture

```
nrgcore/
├── nrgcore/          # Main NRG solver and data structures
│   ├── nrgcore.hpp   # Core NRG iteration logic
│   ├── nrgData.hpp   # Data containers for NRG results
│   ├── qOperator.hpp # Quantum operator definitions
│   └── qsymmetry.hpp # Symmetry and quantum number blocks
├── models/           # Impurity models
│   ├── spinhalf.hpp  # Spin-1/2 models
│   ├── spinnless.hpp # Spinless fermion models
│   └── ...           # Additional models (Kondo, superconducting, etc.)
├── dynamics/         # Dynamic properties
│   ├── fdmSpectrum.hpp      # Frequency domain spectral functions
│   └── fdmBackwardIteration.hpp  # Backward iteration for dynamics
└── utils/            # Utility functions
    ├── h5stream.hpp  # HDF5 I/O utilities
    ├── qmatrix.hpp   # Quantum matrix operations
    └── functions.hpp # Helper functions
```

## 📖 Documentation

Complete documentation available at: https://srbhp.github.io/nrgplusplus/

Key documentation sections:
- [Build Instructions](https://srbhp.github.io/nrgplusplus/build.html)
- [API Reference](https://srbhp.github.io/nrgplusplus/api/library_root.html)
- [NRG Core Documentation](https://srbhp.github.io/nrgplusplus/docnrgcore.html)
- [Model Classes](https://srbhp.github.io/nrgplusplus/docspinhalf.html)

## 🔧 Dependencies

- **C++20 compiler** (GCC 10+, Clang 12+)
- **BLAS/LAPACK** (OpenBLAS, Intel MKL, or equivalent)
- **LAPACKE** (C interface to LAPACK)
- **HDF5** (for data I/O)
- **CMake 3.11+** (build system)
- **Python 3.11+** (optional, for visualization and analysis)

## 💡 Theoretical Foundation

This implementation follows the conventions and methods described in:

> Bulla, R., Costi, T. A., & Pruschke, T. (2008).  
> "The Numerical Renormalization Group Method for Quantum Impurity Systems."  
> *Reviews of Modern Physics*, 80(2), 395–450.  
> https://doi.org/10.1103/RevModPhys.80.395

## 🤝 Contributing

Contributions are welcome! To contribute:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/your-feature`)
3. Make your changes and add tests if applicable
4. Submit a pull request with a description of your changes

For major changes, please open an issue first to discuss your ideas.

## 📋 Similar Projects

Explore other NRG implementations:
- [DMNRG](http://www.phy.bme.hu/~dmnrg/) - Original implementation
- [NRG Ljubljana](https://github.com/rokzitko/nrgljubljana) - Alternative modern implementation

## 📝 License

This project is licensed under the terms specified in the [LICENSE](LICENSE) file.

## 🙏 Acknowledgments

This project was developed with support from:
- **Prof. Jonas Fransson**, Uppsala University
- **Prof. Frithjof Anders**, Technical University Dortmund

Their guidance and support were instrumental in making this project possible.

---

**Questions or issues?** Please open a GitHub issue or consult the [documentation](https://srbhp.github.io/nrgplusplus/).
