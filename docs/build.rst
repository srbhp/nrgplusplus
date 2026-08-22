===========================
Build Guide
===========================

.. contents:: 

Core Dependencies
====================================

**Required for compilation:**

1. **BLAS & LAPACK** — Linear algebra libraries (required)
   
   - **Option A (Recommended)**: OpenBLAS or system BLAS/LAPACK
   - **Option B**: Intel MKL (high performance, optional)
   
2. **LAPACKE** — C interface to LAPACK (required)

3. **HDF5** — For data I/O (required)
   - Save and load NRG iteration states
   - Store large datasets efficiently

4. **C++20 Compiler** — GCC 10+, Clang 12+, or Intel compiler

**Optional:**

- Intel OneAPI MKL (for higher performance)
- Python 3.11+ (for visualization scripts)
- Sphinx, Doxygen, Graphviz (for documentation generation)


Which BLAS/LAPACK to Use?
====================================

**OpenBLAS (Recommended for most users)**
   - Open-source, good performance
   - Easy to install on any Linux distribution
   - Compatible with all compilers (GCC, Clang, Intel)
   - Install package: ``libopenblas-dev``

**Intel MKL (Highest performance, requires Intel account)**
   - Optimized for Intel processors
   - Requires registration at https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html
   - Install via ``intel-oneapi-mkl-devel``

**System BLAS/LAPACK (Minimal)**
   - May have lower performance
   - Good for testing and development




Ubuntu/Debian Based Linux
-----------------------------------

**Option 1: Using OpenBLAS (Recommended)**

This is the easiest and most compatible option:

.. code-block:: bash

    sudo apt-get update
    sudo apt-get install -y \
        libopenblas-dev \
        liblapack-dev \
        liblapacke-dev \
        libhdf5-dev \
        cmake \
        build-essential

    # Optional: C++20 compiler (if GCC 10+ not already installed)
    sudo apt-get install -y gcc-11 g++-11

    # Optional: for documentation generation
    sudo apt-get install -y doxygen graphviz sphinx-doc
    pip3 install sphinx-rtd-theme breathe sphinx-sitemap exhale


**Option 2: Using Intel MKL (High Performance)**

For maximum performance on Intel processors:

.. code-block:: bash

    # Add Intel package repository
    wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
        | gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
    
    echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" \
        | sudo tee /etc/apt/sources.list.d/oneAPI.list
    
    # Install Intel MKL and compiler
    sudo apt-get update
    sudo apt-get install -y \
        intel-oneapi-mkl-devel \
        intel-oneapi-compiler-dpcpp-cpp \
        libhdf5-dev

    # Set up Intel environment variables
    source /opt/intel/oneapi/setvars.sh

    # Optional: Add to ~/.bashrc to persist across sessions
    echo "source /opt/intel/oneapi/setvars.sh" >> ~/.bashrc


Fedora / RHEL / CentOS
-----------------------------------

Using OpenBLAS:

.. code-block:: bash

    sudo dnf install -y \
        openblas-devel \
        lapack-devel \
        lapacke-devel \
        hdf5-devel \
        cmake \
        gcc-c++

Alternatively, using Intel MKL via the Intel repository:

.. code-block:: bash

    # Follow Intel's official repository setup for your distribution
    # https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html
    
    sudo dnf install -y intel-oneapi-mkl-devel intel-oneapi-compiler-dpcpp-cpp


Arch Linux
-----------------------------------

.. code-block:: bash

    sudo pacman -S \
        openblas \
        lapack \
        hdf5 \
        cmake \
        base-devel

    # Optional: Intel MKL
    yay -S intel-oneapi-mkl



Building the Project
====================================

**1. Clone the Repository**

.. code-block:: bash

    git clone https://github.com/srbhp/nrgplusplus.git
    cd nrgplusplus


**2. Create Build Directory**

.. code-block:: bash

    mkdir build
    cd build


**3. Configure with CMake**

.. code-block:: bash

    cmake ..

If CMake fails to find BLAS, see **Troubleshooting** below.


**4. Build**

.. code-block:: bash

    make -j$(nproc)

The compiled executables will be in ``build/examples/*/`` directories.


**5. Run an Example**

.. code-block:: bash

    cd ../examples/rgflowSIAM
    ../../build/examples/rgflowSIAM/rgflowSIAM
    python3 plot.py  # Visualize results


Building a Specific Example Only
====================================

If you only want to build one example (faster for testing):

.. code-block:: bash

    mkdir build
    cd build
    cmake ..
    make rgflowSIAM -j$(nproc)
    ../examples/rgflowSIAM/rgflowSIAM


Troubleshooting
====================================

**CMake Error: "Could NOT find BLAS"**

This means BLAS/LAPACK libraries are not installed or CMake cannot find them.

**Solution:**

1. **Install BLAS/LAPACK first:**

   On Ubuntu/Debian:
   
   .. code-block:: bash
   
       sudo apt-get install -y libopenblas-dev liblapack-dev liblapacke-dev
   
   On Fedora/RHEL:
   
   .. code-block:: bash
   
       sudo dnf install -y openblas-devel lapack-devel lapacke-devel
   
   On macOS:
   
   .. code-block:: bash
   
       brew install openblas lapack

2. **Clean and reconfigure CMake:**

   .. code-block:: bash
   
       rm -rf build
       mkdir build
       cd build
       cmake ..
       make -j$(nproc)

**CMake Error: "Could NOT find HDF5"**

Install HDF5 development files:

.. code-block:: bash

    # Ubuntu/Debian
    sudo apt-get install -y libhdf5-dev
    
    # Fedora/RHEL
    sudo dnf install -y hdf5-devel
    
    # macOS
    brew install hdf5


**CMake Finds Old BLAS Version**

If you have both OpenBLAS and MKL installed, you can specify which one to use:

.. code-block:: bash

    cd build
    rm CMakeCache.txt
    cmake -DBLA_VENDOR=OpenBLAS ..
    make -j$(nproc)

Alternative vendors: ``OpenBLAS``, ``Intel10_64lp_seq``, ``ATLAS``, ``PhiPACK``, ``ACML``, ``Apple``, ``NAS``, ``Generic``


**Compiler Not Found**

If you get C++20 compiler errors, install a newer compiler:

.. code-block:: bash

    # Ubuntu/Debian - Install GCC 11
    sudo apt-get install -y gcc-11 g++-11
    
    # Set as default
    sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-11 100
    sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-11 100


Building Documentation
====================================

To generate HTML documentation locally:

.. code-block:: bash

    # Install documentation tools
    sudo apt-get install -y doxygen graphviz sphinx-doc
    pip3 install sphinx-rtd-theme breathe exhale
    
    # Build documentation
    cd build
    cmake ..
    make docs
    
    # Open documentation in browser
    open ../docs/build/html/index.html  # macOS
    xdg-open ../docs/build/html/index.html  # Linux

