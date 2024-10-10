#!/usr/bin/env bash

version=4.6.1

TEMPORARY_FILES="${TMPDIR:-/tmp}"
export NETCDF_INSTALL_DIR=$(pwd)/netcdf-install
while [ $# != 0 ]; do
    case "$1" in
    "--prefix")
        export NETCDF_INSTALL_DIR="$2"; shift
        ;;
    "--tmpdir")
        TEMPORARY_FILES="$2"; shift
        ;;
    "--version")
        version="$2"; shift
        ;;
    "--netcdf-root")
        export NETCDF_ROOT="$2"; shift
        ;;
    *)
        echo "Unrecognized argument '$1'"
        exit 1
        ;;
    esac
    shift
done

NETCDF_MIRROR=https://downloads.unidata.ucar.edu/netcdf-fortran/
NETCDF_VERSION=${version}

URL=${NETCDF_MIRROR}/${NETCDF_VERSION}/netcdf-fortran-${NETCDF_VERSION}.tar.gz
FOLDER=/netcdf-fortran-${NETCDF_VERSION}

if [ ! -d "${TEMPORARY_FILES}/${FOLDER}" ]; then
  echo "Downloading ${TEMPORARY_FILES}/${FOLDER} from URL [${URL}]"
  mkdir -p ${TEMPORARY_FILES}
  curl --location \
       "${URL}" | tar zx -C "${TEMPORARY_FILES}"
else
   echo "Download already present in ${TEMPORARY_FILES}/${FOLDER}"
fi


mkdir -p ${TEMPORARY_FILES}/build-${FOLDER} && cd ${TEMPORARY_FILES}/build-${FOLDER}
rm -rf ./*
NETCDF_LIB_DIR="$(${NETCDF_ROOT}/bin/nc-config --libdir)"
cmake -G Ninja ${TEMPORARY_FILES}/${FOLDER} \
    -DnetCDF_LIBRARIES="${NETCDF_LIB_DIR}/libnetcdf.so" -DnetCDF_INCLUDE_DIR="${NETCDF_ROOT}/include" \
    -DENABLE_TESTS=OFF -DCMAKE_INSTALL_PREFIX="${NETCDF_INSTALL_DIR}"
cmake --build . --config Release
cmake --install .
