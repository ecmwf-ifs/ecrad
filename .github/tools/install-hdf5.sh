#!/usr/bin/env bash

version=1.14.4-3

TEMPORARY_FILES="${TMPDIR:-/tmp}"
export HDF5_INSTALL_DIR=$(pwd)/hdf5-install
while [ $# != 0 ]; do
    case "$1" in
    "--prefix")
        export HDF5_INSTALL_DIR="$2"; shift
        ;;
    "--tmpdir")
        TEMPORARY_FILES="$2"; shift
        ;;
    "--version")
        version="$2"; shift
        ;;
    *)
        echo "Unrecognized argument '$1'"
        exit 1
        ;;
    esac
    shift
done

HDF5_MIRROR=https://support.hdfgroup.org/ftp/HDF5/releases
HDF5_VERSION=${version}

# Pick out version parts separated by '.'
VERSION_PARTS=($(echo ${HDF5_VERSION} | tr "." "\n"))
# Major version, e.g., 1.14
MAJOR_VERSION=${VERSION_PARTS[0]}.${VERSION_PARTS[1]}

# Minor version parts, including patch level (if any), e.g., 3 or 4-3
MINOR_VERSION_PARTS=($(echo ${VERSION_PARTS[2]} | tr "-" "\n"))

# Minor version without patch level
MINOR_VERSION=${MINOR_VERSION_PARTS[0]}

URL=${HDF5_MIRROR}/hdf5-${MAJOR_VERSION}/hdf5-${MAJOR_VERSION}.${MINOR_VERSION}/src/hdf5-${HDF5_VERSION}.tar.gz
FOLDER=hdf5-${HDF5_VERSION}

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
cmake -G Ninja ${TEMPORARY_FILES}/${FOLDER} -DHDF5_BUILD_FORTRAN=ON -DHDF5_BUILD_HL_LIB=ON -DBUILD_TESTING=OFF
cmake --build . --config Release
cmake --install . --prefix "${HDF5_INSTALL_DIR}"
