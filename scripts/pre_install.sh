#!/bin/bash
#
sudo apt-get update 
sudo apt-get install sphinx-doc doxygen graphviz libhdf5-dev --yes
pip3 install sphinx-rtd-theme breathe sphinx-sitemap sphinx exhale
wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \ | gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list
sudo apt update
sudo apt install intel-oneapi-mkl-devel intel-oneapi-compiler-dpcpp-cpp --yes
source  /opt/intel/oneapi/setvars.sh 


