#!/usr/bin/env bash
#
# Build HDF5, NetCDF-C and NetCDF-Fortran (with the compilers of the current
# environment) for the pyecrad wheels; called by build_wheel.sh inside the
# manylinux container. The libraries are bundled in the wheel by auditwheel.
#
# The installation directory is cached by the pyecrad_wheels workflow: its
# cache key depends on this file and on the .github/tools/install-*.sh scripts.
#
# Usage: build_wheel_deps.sh PREFIX
#   HDF5 is installed in PREFIX/hdf5, NetCDF-C and NetCDF-Fortran in
#   PREFIX/netcdf. Nothing is done if PREFIX/complete exists.

set -e

PREFIX="$1"
TOOLS="$(cd "$(dirname "${BASH_SOURCE[0]}")/../tools" && pwd)"
SRC="${TMPDIR:-/tmp}/wheel-deps-src"

if [ -f "${PREFIX}/complete" ]; then
    echo "HDF5 and NetCDF already built in ${PREFIX}"
    exit 0
fi
rm -rf "${PREFIX:?}/hdf5" "${PREFIX:?}/netcdf"

# Libraries in lib (not lib64), as expected by the scripts and FindNetCDF.cmake
# NetCDF-C: no remote access (DAP, byte-range, NCZarr), which would bring
# libcurl and OpenSSL in the wheel
export HDF5_EXTRA_CMAKE_ARGS="-DCMAKE_INSTALL_LIBDIR=lib -DHDF5_ENABLE_Z_LIB_SUPPORT=ON"
export NETCDF_C_EXTRA_CMAKE_ARGS="-DCMAKE_INSTALL_LIBDIR=lib -DENABLE_DAP=OFF -DENABLE_DAP4=OFF -DENABLE_BYTERANGE=OFF -DENABLE_NCZARR=OFF"
export NETCDF_FORTRAN_EXTRA_CMAKE_ARGS="-DCMAKE_INSTALL_LIBDIR=lib"

"${TOOLS}/install-hdf5.sh" --prefix "${PREFIX}/hdf5" --tmpdir "${SRC}"
"${TOOLS}/install-netcdf-c.sh" --prefix "${PREFIX}/netcdf" --tmpdir "${SRC}" \
    --hdf5-root "${PREFIX}/hdf5"
"${TOOLS}/install-netcdf-fortran.sh" --prefix "${PREFIX}/netcdf" --tmpdir "${SRC}" \
    --netcdf-root "${PREFIX}/netcdf"

touch "${PREFIX}/complete"
