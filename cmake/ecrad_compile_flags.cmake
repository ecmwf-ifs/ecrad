# (C) Copyright 2014- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

if(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
  set(checkbounds_flags   "-Rb")
  set(fpe_flags           "-Ktrap=fp")
  set(initsnan_flags      "-ei")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(checkbounds_flags   "-fcheck=bounds")
  set(fpe_flags           "-ffpe-trap=invalid,zero,overflow")
  set(initsnan_flags      "-finit-real=snan")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
  set(align_flags         "-align array64byte")
  set(alloc_flags         "-heap-arrays 32")
  set(checkbounds_flags   "-check bounds")
  set(initsnan_flags      "-init=snan")
  set(inline_flags        "-finline-functions")
  set(vectorization_flags "-assume byterecl,realloc_lhs -march=core-avx2 -no-fma")
  set(fpmodel_flags       "-fpe0 -fp-model=precise -fp-speculation=safe")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(align_flags         "-align array64byte")
  set(alloc_flags         "-heap-arrays 32")
  set(checkbounds_flags   "-check bounds")
  set(initsnan_flags      "-init=snan")
  set(inline_flags        "-finline-functions -finline-limit=1500 -Winline")
  set(vectorization_flags "-assume byterecl,realloc_lhs -march=core-avx2 -no-fma")
  set(fpmodel_flags       "-fpe0 -fp-model precise -fp-speculation=safe -ftz -fast-transcendentals")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI|NVHPC")
  set(fpe_flags           "-Ktrap=fp")
  set(vectorization_flags "-O3 -fast")
  string(REPLACE "-O2" "" CMAKE_Fortran_FLAGS_BIT ${CMAKE_Fortran_FLAGS_BIT})
  set(checkbounds_flags   "-Mbounds")

endif()

if( DEFINED align_flags )
  ecbuild_add_fortran_flags( "${align_flags}"   NAME align )
endif()
if( DEFINED alloc_flags )
  ecbuild_add_fortran_flags( "${alloc_flags}"   NAME alloc )
endif()
if( DEFINED convert_flags )
  ecbuild_add_fortran_flags( "${convert_flags}"   NAME convert )
endif()
if( DEFINED vectorization_flags )
  ecbuild_add_fortran_flags( "${vectorization_flags}" NAME vectorization BUILD BIT )
endif()
if( DEFINED fpmodel_flags )
  ecbuild_add_fortran_flags( "${fpmodel_flags}"   NAME fpmodel BUILD BIT )
endif()
if( DEFINED inline_flags )
  ecbuild_add_fortran_flags( "${inline_flags}"   NAME inline BUILD BIT )
endif()

foreach( debug_flag fpe initsnan checkbounds )
  if( ${debug_flag}_flags )
    ecbuild_add_fortran_flags( "${${debug_flag}_flags}" NAME ${debug_flag} BUILD DEBUG )
  endif()
endforeach()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  # In case '-check all' has been added, we need to remove the '-check arg_temp_created' warnings
  ecbuild_add_fortran_flags( "-check noarg_temp_created" NAME noarg BUILD DEBUG)
endif()
