# (C) Copyright 2014- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Capture ecbuild defaults and/or flags set by a toolchain
set( ${PNAME}_Fortran_FLAGS "${${PNAME}_Fortran_FLAGS} ${ECBUILD_Fortran_FLAGS}" )
set( ${PNAME}_Fortran_FLAGS_BIT "${${PNAME}_Fortran_FLAGS_BIT} ${ECBUILD_Fortran_FLAGS_BIT}" )
set( ${PNAME}_Fortran_FLAGS_DEBUG "${${PNAME}_Fortran_FLAGS_DEBUG} ${ECBUILD_Fortran_FLAGS_DEBUG}" )

if(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
  set(autopromote_flags   "-sreal64")
  set(checkbounds_flags   "-Rb")
  set(fpe_flags           "-Ktrap=fp")
  set(initsnan_flags      "-ei")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(autopromote_flags   "-fdefault-real-8 -fdefault-double-8")
  set(checkbounds_flags   "-fcheck=bounds")
  set(fpe_flags           "-ffpe-trap=invalid,zero,overflow")
  set(initsnan_flags      "-finit-real=snan")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(autopromote_flags   "-real-size 64")
  set(checkbounds_flags   "-check bounds")
  set(initsnan_flags      "-init=snan")
  set(fpe_flags           "-fpe0")
  set(vectorization_flags "-march=core-avx2 -no-fma")
  set(fpmodel_flags       "-fp-model precise -fp-speculation=safe")
  set(transcendentals_flags "-fast-transcendentals")
  set(heap_arrays_flags   "-heap-arrays 32")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI|NVHPC")
  set(autopromote_flags   "-r8")
  set(fpe_flags           "-Ktrap=fp")
  set(vectorization_flags "-O3 -fast")
  string(REPLACE "-O2" "" ${PNAME}_Fortran_FLAGS_BIT ${${PNAME}_Fortran_FLAGS_BIT})
  set(checkbounds_flags   "-Mbounds")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Flang")
  set(autopromote_flags   "-fdefault-real-8")
  set(fpe_flags           "-ffp-exception-behavior=strict")

endif()

ecbuild_add_fortran_flags( "-g -O0"   NAME base_debug BUILD DEBUG)
# if( NOT HAVE_SINGLE_PRECISION )
#   ecbuild_add_fortran_flags( "${autopromote_flags}"   NAME autopromote )
# endif()
if( DEFINED vectorization_flags )
  # vectorization flags must be per-sourcefile overrideable, so are set via ${PNAME}_Fortran_FLAGS
  # set( ${PNAME}_Fortran_FLAGS_BIT "${${PNAME}_Fortran_FLAGS_BIT} ${vectorization_flags}" )
  string( CONCAT ${PNAME}_Fortran_FLAGS_BIT "${${PNAME}_Fortran_FLAGS_BIT} -fpe0 -convert big_endian -align array64byte "
      "-assume byterecl,noold_maxminloc,realloc_lhs "
      "-traceback -fpic -qopenmp -qopenmp-threadprivate compat -fp-model precise "
      "-fp-speculation=safe -qopt-report=2 -qopt-report-phase=all -fast-transcendentals -ftz "
      "-finline-functions -finline-limit=1500 -Winline "
      "-diag-disable=7713 "   # disable unused statement functions warning
      "-diag-disable=11021 "  # disable unresolved libraries in ipo warning
      "-diag-disable=10397 "  # disable message reporting location of opt report files
      "-diag-disable=8291" )  # disable remark: Recommended relationship between field width 'W' and the number of fractional digits 'D' in this edit descriptor is 'W>=D+7'
endif()
# if( DEFINED fpmodel_flags )
#   ecbuild_add_fortran_flags( "${fpmodel_flags}"   NAME fpmodel BUILD BIT)
# endif()
# if( DEFINED transcendentals_flags )
#   ecbuild_add_fortran_flags( "${transcendentals_flags}"   NAME transcendentals BUILD BIT)
# endif()
# if( DEFINED heap_arrays_flags )
#   ecbuild_add_fortran_flags( "${heap_arrays_flags}"   NAME heap_arrays )
# endif()

if( CMAKE_BUILD_TYPE MATCHES "Debug" )
  foreach( debug_flag    fpe initsnan checkbounds )
    if( ${debug_flag}_flags )
      set( ${PNAME}_Fortran_FLAGS_DEBUG "${${PNAME}_Fortran_FLAGS_DEBUG} ${${debug_flag}_flags}" )
    endif()
  endforeach()
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    # In case '-check all' has been added, we need to remove the '-check arg_temp_created' warnings
    set( ${PNAME}_Fortran_FLAGS_DEBUG "${${PNAME}_Fortran_FLAGS_DEBUG} -check noarg_temp_created" )
  endif()
endif()

# if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
#   if( NOT CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 10 )
#     ecbuild_add_fortran_flags( "-fallow-argument-mismatch" NAME argument_mismatch )
#   endif()
#   if( NOT CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "aarch64" )
#     ecbuild_add_fortran_flags( "-m64" NAME gnu_arch )
#   endif()
#   if( LOKI_MODE MATCHES "idem-stack|scc-stack" AND HAVE_LOKI )
#     ecbuild_add_fortran_flags( "-fcray-pointer" NAME cray_pointer )
#   endif()
# endif()

# if(CMAKE_Fortran_COMPILER_ID MATCHES "Flang")
#   # Linker complains of unknown arguments:
#   #    warning: argument unused during compilation: '-fdefault-real-8' [-Wunused-command-line-argument]
#   foreach( LINKER_FLAGS CMAKE_EXE_LINKER_FLAGS CMAKE_SHARED_LINKER_FLAGS CMAKE_STATIC_LINKER_FLAGS )
#     set( ${LINKER_FLAGS} "${${LINKER_FLAGS}} -Wno-unused-command-line-argument")
#   endforeach()
# endif()
