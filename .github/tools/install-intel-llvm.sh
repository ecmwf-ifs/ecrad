#!/usr/bin/env bash
# Installation instructions from https://www.intel.com/content/www/us/en/docs/oneapi/installation-guide-linux/2025-0/apt-005.html

version=2023.2.0
KEY=GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB

# download the key to system keyring
wget -nv -O- https://apt.repos.intel.com/intel-gpg-keys/$KEY \
  | gpg --dearmor \
  | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null

# add signed entry to apt sources and configure the APT client to use Intel repository:
echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" \
  | sudo tee /etc/apt/sources.list.d/oneAPI.list

sudo apt-get update
sudo apt-get install \
    intel-oneapi-compiler-fortran-$version \
    intel-oneapi-compiler-dpcpp-cpp-$version \
    intel-oneapi-mpi-devel-2021.10.0 \
    intel-oneapi-mkl-$version
