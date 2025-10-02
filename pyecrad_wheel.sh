#!/usr/bin/env bash

# This script builds the python wheel using a container with low version of GLIBC
# allowing execution on a wide range of linux systems.
# Script can be called with up to two arguments: target and pypi
# target van be:
# - get: to get the container and dependencies source code
# - deps: to build the dependencies
# - lib: to build the library
# - wheel: to build the wheels and the source distribution
# - pypi: to push the source distributions and the wheels on PyPI
#         in this case, the pypi argument must be provided
# - all: to perform all these actions
# pypi is the repository on which wheels are pushed ('pypi' and 'all'
# targets.

# versions of dependencies
netcdf_version="v4.9.2"
netcdf_uri=http://github.com/Unidata/netcdf-c
netcdffortran_version="v4.6.2"
netcdffortran_uri=https://github.com/Unidata/netcdf-fortran
hdf5=https://resources.unidata.ucar.edu/netcdf/netcdf-4/hdf5-1.8.15.tar.bz2
curl=https://curl.haxx.se/download/curl-8.5.0.tar.gz

# container image to be used
container_uri=docker://quay.io/pypa/manylinux2014_x86_64

# local certificate (for pip)
local_certificate=/etc/ssl/certs/ca-certificates.crt

# python versions for which to build a wheel
#export python_versions="cp310-cp310 cp311-cp311 cp312-cp312 cp313-cp313 cp313-cp313t"
# This module does not depend on the exact version used, we pick one
export python_versions="cp312-cp312"

# end of user customisation
# =============================================================================

# working directories
export ECRAD_ROOT=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
export TMPWORKDIR=$ECRAD_ROOT/tmp
export SRC=$TMPWORKDIR/src
mkdir -p $TMPWORKDIR
mkdir -p $SRC

# container stuff
export CONTAINER_SIF_PATH="$TMPWORKDIR/$(echo $container_uri | awk -F '/' '{print $NF}').sif"
export CONTAINER_ROOT=/work  # bind path from within the container
export SINGULARITY_BINDPATH=$ECRAD_ROOT:$CONTAINER_ROOT  # binding out:in
export SINGULARITY_TMPDIR=$TMPWORKDIR/singularity
export REQUESTS_CA_BUNDLE=$TMPWORKDIR/ca-certificates.crt  # accessible copy within the container
mkdir -p $SINGULARITY_TMPDIR

# cmake/make stuff
export BUILD=$CONTAINER_ROOT/tmp/build
export INSTALL_DIR=$CONTAINER_ROOT/tmp/install  # path within the container for the sake of wheel link edition
export CMAKE_INSTALL_PREFIX=$INSTALL_DIR

# end of definitions
# =============================================================================

# target action
target=$1
if [ "$target" != "all" ] && [ "$target" != "get" ] && [ "$target" != "deps" ] && \
   [ "$target" != "lib" ] && [ "$target" != "wheel" ]  && [ "$target" != "pypi" ]; then
  echo "First argument must be one of 'all', 'get', 'deps', 'lib', 'wheel', 'pypi'"
  exit 1
fi

# ecrad version
version=$(cat VERSION)

# get container and dependencies
if [ "$target" == "get" ] || [ "$target" == "all" ]; then
  # clone dependencies packages (out of container)
  [ -d $SRC/netcdf ] && rm -rf $SRC/netcdf
  git clone $netcdf_uri -b $netcdf_version $SRC/netcdf

  [ -d $SRC/netcdffortran ] && rm -rf $SRC/netcdffortran
  git clone $netcdffortran_uri -b $netcdffortran_version $SRC/netcdffortran

  rm -rf $SRC/hdf5-* $SRC/hdf5.tgz
  wget $hdf5 -O $SRC/hdf5.tgz
  (cd $SRC; tar xf hdf5.tgz)

  rm -rf $SRC/curl-* $SRC/curl.tgz
  wget $curl -O $SRC/curl.tgz
  (cd $SRC; tar xf curl.tgz; rm -f curl.tgz)

  # get container
  singularity pull $CONTAINER_SIF_PATH $container_uri

  # copy certificate for pip install from within the container
  cp $local_certificate $REQUESTS_CA_BUNDLE
fi

# build dependencies
if [ "$target" == "deps" ] || [ "$target" == "all" ]; then
  cat - <<..EOF > ${TMPWORKDIR}/build_deps.sh
  set -e

  cd $SRC/curl-*
  ./configure --prefix=$INSTALL_DIR --without-ssl
  make -j 4
  make install

  cd $SRC/hdf5-*
  ./configure --prefix=$INSTALL_DIR
  make -j 4
  make install

  rm -rf $BUILD
  mkdir -p $BUILD
  cd $BUILD
  cmake -S $SRC/netcdf -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR
  cmake --build $BUILD
  cmake --install $BUILD

  rm -rf $BUILD
  mkdir -p $BUILD
  cd $BUILD
  cmake -S $SRC/netcdffortran -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR
  cmake --build $BUILD
  cmake --install $BUILD
..EOF
  chmod +x ${TMPWORKDIR}/build_deps.sh
  singularity run -i $CONTAINER_SIF_PATH ${TMPWORKDIR}/build_deps.sh
fi

# build ecrad lib
if [ "$target" == "lib" ] || [ "$target" == "all" ]; then
  cat - <<..EOF > ${TMPWORKDIR}/build_lib.sh
  export PATH=\$PATH:$INSTALL_DIR/bin # to find nf-config
  make python
..EOF
  chmod +x ${TMPWORKDIR}/build_lib.sh
  singularity run -i $CONTAINER_SIF_PATH ${TMPWORKDIR}/build_lib.sh
fi

# build wheels
if [ "$target" == "wheel" ] || [ "$target" == "all" ]; then
  cat - <<..EOF > ${TMPWORKDIR}/build_wheel.sh
  export PATH=\$PATH:$INSTALL_DIR/bin # to find nf-config
  export LD_LIBRARY_PATH=$INSTALL_DIR/lib:$INSTALL_DIR/lib64 # for auditwheel repair
  for p in $python_versions; do
    /opt/python/\$p/bin/pip install build
    /opt/python/\$p/bin/pip install auditwheel
    /opt/python/\$p/bin/python -m build --wheel
    wheel=dist/pyecrad-${version}-\$p-linux_x86_64.whl
    if [ -f \$wheel ]; then
      /opt/python/\$p/bin/python -m auditwheel repair \$wheel
    fi
  done
  wheel=dist/pyecrad-${version}-py3-none-any.whl
  if [ -f \$wheel ]; then
    /opt/python/\$p/bin/python -m auditwheel repair \$wheel
  fi
..EOF
  chmod +x ${TMPWORKDIR}/build_wheel.sh
  singularity run -i $CONTAINER_SIF_PATH ${TMPWORKDIR}/build_wheel.sh
fi

# distribute wheels
if [ "$target" == "pypi" ] || [ "$target" == "all" ]; then
  pypi=$2
  if [ "$pypi" == "" ]; then
    echo "Second argument must be the PyPI server to use"
    exit 2
  fi
  python3 -m twine upload --repository $pypi wheelhouse/pyecrad-${version}-py3-none-manylinux2014_x86_64.manylinux_2_17_x86_64.whl
fi
