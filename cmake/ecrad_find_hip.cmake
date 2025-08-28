# (C) Copyright 2020- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

macro( ecrad_find_hip )
  # This macro finds all HIP related libraries, if found, HAVE_HIP=1

  cmake_minimum_required( VERSION 3.24 FATAL_ERROR )

  set( options "" )
  set( single_value_args REQUIRED )
  set( multi_value_args "" )
  cmake_parse_arguments( _PAR "${options}" "${single_value_args}" "${multi_value_args}"  ${_FIRST_ARG} ${ARGN} )

  set(HIP_REQUIRED "")
  if( _PAR_REQUIRED )
    set(HIP_REQUIRED "REQUIRED" )
  endif()

  set(HAVE_HIP 1)

  # Setup ROCM_PATH
  if (NOT DEFINED ROCM_PATH )
    find_path(ROCM_PATH
      hip
      ENV{ROCM_DIR}
      ENV{ROCM_PATH}
      ENV{HIP_PATH}
      ${HIP_PATH}/..
      ${HIP_ROOT_DIR}/../
      ${ROCM_ROOT_DIR}
      /opt/rocm)
  endif()
  ecbuild_info("ROCM path: ${ROCM_PATH}")
  # Update CMAKE_PREFIX_PATH to make sure all the configs that hip depends on are found.
  set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${ROCM_PATH}")

  set(HAVE_HIP 1)

  if( HAVE_HIP )
    find_package(rocprofiler-sdk-roctx CONFIG ${HIP_REQUIRED} PATHS ${ROCM_PATH}/lib)
    if( NOT rocprofiler-sdk-roctx_FOUND )
      ecbuild_info("rocprofiler-sdk-roctx libraries not found: HAVE_HIP=0")
      set( HAVE_HIP 0 )
    endif()

    if( HAVE_HIP )
      list( APPEND ECRAD_GPU_ROCPROFILER_SDK_LIBRARIES ${rocprofiler-sdk-roctx_LIBRARIES})
    endif()

  endif()
  ecbuild_info("ROCProfiler SDK libraries: ${ECRAD_GPU_ROCPROFILER_SDK_LIBRARIES}")
endmacro()
