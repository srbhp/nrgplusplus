===========================
QuickStart -
===========================

.. contents:: 

Install dependencies
====================================

-   We use Intel Compiler and MKL library for our matrix operations.
  ``intel-oneapi-mkl-devel`` and ``intel-oneapi-compiler-dpcpp-cpp`` are these two packages for this purpose. GCC and Clang also works fine but we anyway need MKL 
  for lapack. If you prefer not to use MKL then comment MKL #include part  
  in the qmatrix.hpp file and uncomment the lapacke part.

- We use HDF5 (``libhdf5-dev``) for our to read and write data. We can also fully
  save the NRG iteration state into a hdf5 file and load the NRG state to continue
  the caculation. This usefull for during the calculation of the static and 
  dynamic quantities. 

- To generate the documentation we need few more python packages. 
  Install these packages if you want to generate the documentation. 
  
  - ``sphinx-doc doxygen graphviz``
  - ``sphinx-rtd-theme breathe sphinx-sitemap sphinx exhale``




Ubuntu/Debian based Linux OS
-----------------------------------
Here is a list commands for to install the dependencies.

  .. code-block:: python

    sudo apt-get update 
    sudo apt-get install sphinx-doc doxygen graphviz libhdf5-dev --yes
    pip3 install sphinx-rtd-theme breathe sphinx-sitemap sphinx exhale
    wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \ | gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
    echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list
    sudo apt update
    sudo apt install intel-oneapi-mkl-devel intel-oneapi-compiler-dpcpp-cpp --yes
    source  /opt/intel/oneapi/setvars.sh # We should add this line to the end of our ~/.bashrc line  


Other Linux OS
-----------------------------------
We can use the currect package manager for to configure
the repository for ``intel-mkl``. 
https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html

We can use ``pip3`` to  install python Packages.



How to Build a Model
====================================

Clone the `nrgplusplus` repo.
-----------------------------------

- [optional] Create a new folder. ``mkdir test_nrg; cd test_nrg``

- Clone the repo. ``git clone https://github.com/srbhp/nrgplusplus.git``


Examples
-----------------------------------

``examples`` directory contains examples of different Models. 

- Copy a folder here from ``examples`` to start woking on a specific model.
-  Edit ``CMakeLists.txt`` to add the model folder.
    for example ``add_subdirectory(ExampleModel)``  to the end of the ``CMakeLists.txt`` file.
- `build` the project by invoking `make`.
   - create a build directory and invoke cmake from there
    - ``mkdir build ; cd build ; cmake .. ; make``
    - if we manage to compile the project, then the executable will be in the ``build/ExampleModel/`` directory.



================

