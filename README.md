# ECRAD - ECMWF atmospheric radiation scheme

This document last updated 9 June 2022

Robin Hogan <r.j.hogan@ecmwf.int>

For more complete information about compilation and usage of ecRad,
please see the documentation on the
[ecRad web site](https://confluence.ecmwf.int/display/ECRAD).


## INTRODUCTION

This package contains the offline version of a radiation scheme
suitable for use in atmospheric weather and climate models.  The code
is designed to be extensible and flexible.  For example, the gas
optics, cloud optics and solver are completely separated (see
`radiation/radiation_interface.F90` where they are called in sequence),
thereby facilitating future changes where different gas models or
solvers may be switched in and out independently. The offline code is
parallelized using OpenMP.

Five solvers are currently available:

1. The Monte Carlo Independent Column Approximation (McICA) of Pincus
et al. (2003). This is is a now widely used method for treating cloud
structure efficiently. The implementation in this package is more
efficient than the one currently operational in the ECMWF model, and
produces less noise in partially cloudy situations. Note that since
McICA is stocastic, individual flux profiles using McICA may differ
simply due to random variations in the sampling of the cloud field.

2. The Tripleclouds scheme of Shonk and Hogan (2008). This represents
cloud structure by dividing each layer into three regions, one clear
and two cloudy with different optical depth. It is somewhat slower
than McICA but does not generate noise.

3. The Speedy Algorithm for Radiative Transfer through Cloud Sides
(SPARTACUS) of Hogan et al. (JGR 2016). This is a method for
efficiently treating 3D radiative effects associated with clouds. It
uses the same differential equations proposed by Hogan and Shonk (JAS
2013), but solves them using a matrix exponential method that is much
more elegant than their method, and is also here extended to the
longwave (see Schaefer et al., JGR 2016).  It also incorporates the
Tripleclouds methodology of Shonk and Hogan (2008) to represent cloud
inhomogeneity.

4. A homogeneous (plane parallel) solver in which clouds are assumed
to fill the gridbox horizontally.  This is useful for computing
Independent Column Approximation benchmarks.

5. A "cloudless" solver if your focus is on clear skies.

Two gas optics models are available:

1. The Rapid Radiative Transfer Model for GCMs (RRTMG), the
implementation being that from the ECMWF Integrated Forecasting System
(IFS).

2. The ECMWF Correlated k-Distribution (ecCKD) scheme (since ecRad
1.5), which uses a flexible discretization of the spectrum that is
read from a file at run-time.


## PACKAGE OVERVIEW

The subdirectories are as follows:

- `radiation` - the ecRad souce code

- `ifsaux` - source code providing a (sometimes dummy) IFS environment

- `ifsrrtm` - the IFS implementation of the RRTMG gas optics scheme

- `utilities` - source code for useful utilities, such as reading netCDF
       files

- `drhook` - dummy version of the Dr Hook profiling system

- `driver` - the source code for the offline driver program

- `ifs` - slightly modified source files from the IFS that are used to provide inputs to
        ecRad, but not used in this offline version except if you compile the ecrad_ifs_driver executable

- `mod` - where Fortran module files are written

- `lib` - where the static libraries are written

- `bin` - where the executable ecrad is written

- `data` - contains configuration data read at run-time

- `test` - test cases including Matlab code to plot the outputs

- `include` - automatically generated interface blocks for non-module routines

- `practical` - exercises to get started with ecRad 


## TO COMPILE

1. Ensure you have a reasonably recent Fortran compiler - it needs to
support modules with `contains` and `procedure` statements for
example.  Ensure you have the Fortran netCDF library installed
(versions 3 or 4) and that the module file is compatible with your
Fortran compiler.

2. You can compile the code using 

       make PROFILE=<prof>

   where `<prof>` is one of `gfortran`, `pgi`, `cray` or `intel`.
   This will read the system-specific configurations from the file
   `Makefile_include.<prof>`.  If you omit the `PROFILE=` option then
   `gfortran` will be assumed. If you have a compiler other than these
   then create such a file for your compiler following the example in
   `Makefile_include.gfortran`. Two additional profiles are provided,
   `ecmwf` which builds on the `gfortran` profile and `uor`
   (University of Reading) which is built on the `pgi` profile.
   
   If the compile is successful then static libraries should appear in
   the `lib` directory, and then the executable `bin/ecrad`.

3. To clean-up, type `make clean`.  To build an unoptimized version
   for debugging, you can do
   
       make PROFILE=<prof> DEBUG=1
   
   or you can specifically override the variables in `Makefile_include.<prof>`
   using, for example
   
       make PROFILE=<prof> OPTFLAGS=-O0 DEBUGFLAGS="-g -pg"
   
   To compile in single precision add `SINGLE_PRECISION=1` to the
   `make` command line.  To compile with the Dr Hook profiling system,
   first install ECMWF's [fiat library]([ecRad web
   site](https://github.com/ecmwf-ifs/fiat), then add
   `FIATDIR=/path/to/fiat` to the `make` command line, such that the
   files `$FIATDIR/lib/libfiat.so` and
   `$FIATDIR/module/fiat/yomhook.mod` can be found at build time.
   

## TO TEST

The offline driver is run via

    ecrad <namelist.nam> <input_file.nc> <output_file.nc>

where the radiation scheme is configured using the Fortran namelist
file `<namelist.nam>`, and the inputs and outputs are in netCDF
format.

The `practical` directory contains a set of practical exercises to
help new users become familiar with the capabilities of ecRad. Start
by reading the instructions in `practical/ecrad_practical.pdf`.

The `test/ifs` directory contains a pole-to-pole slice of
low-resolution IFS model data in a form to use as input to the offline
version of ecRad. It includes aerosols extracted from the CAMS
climatology used operationally since IFS Cycle 43R3. Typing `make
test` in this directory runs a number of configurations of ecRad
described in the Makefile. The Matlab script `plot_ifs.m` can be used
to visualize the results. The file
`ecrad_meridian_default_out_REFERENCE.nc` contains a reference version
of the output file `ecrad_meridian_default_out.nc` (case "a"), which
you can compare to be sure your compilation is working as
expected. This case has essentially been superceded by the slice in the
`practical` directory.

The `test/i3rc` directory contains the 1D profile of the I3RC cumulus
test case used by Hogan et al. (2016). Typing `make test` in this
directory runs the various 1D and 3D configurations of ecRad. The
Matlab script `plot_i3rc.m` can then be used to visualize the results,
reproducing three of the figures from Hogan et al. (2016). Note that
you will need to ensure that a reasonably up-to-date version of the
`nco` tools are available and in your path.  This test involves
running the duplicate_profiles.sh script, which duplicates the single
profile in `i3rc_mls_cumulus.nc`, each with a different solar zenith
angle.

The `test/surface` directory contains tests of the surface tile types,
although this is under development and so nothing here is guaranteed
to work.

Alternatively, type `make test` in the top-level directory to run all
cases.

In addition to writing the output file, a file containing the
intermediate radiative properties of the atmosphere for each g-point
can be stored in `radiative_properties.nc` (edit the config namelist to
enable this), but note that the g-points have been reordered in
approximate order of optical depth if the SPARTACUS solver is chosen.


## LICENCE

(C) Copyright 2014- ECMWF.

This software is licensed under the terms of the Apache Licence Version 2.0
which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

In applying this licence, ECMWF does not waive the privileges and immunities
granted to it by virtue of its status as an intergovernmental organisation
nor does it submit to any jurisdiction.
Copyright statements are given in the file NOTICE.

The ifsrrtm directory of this package includes a modified version of
the gas optics part of the Rapid Radiative Transfer Model for GCMS
(RRTMG).  RRTMG was developed at Atmospheric & Environmental Research
(AER), Inc., Lexington, Massachusetts and is available under the
"3-clause BSD" license; for details, see ifsrrtm/AER-BSD3-LICENSE.


## PUBLICATIONS

The ecRad radiation scheme itself is described here:

 - Hogan, R. J., and A. Bozzo, 2018: A flexible and efficient radiation
scheme for the ECMWF model.  J. Adv. Modeling Earth Syst., 10, 1990-2008,
doi:10.1029/2018MS001364.

 - Hogan, R. J., and A. Bozzo, 2016: ECRAD: A new radiation scheme for
the IFS. ECMWF Technical Memorandum number 787, 35pp:
http://www.ecmwf.int/en/elibrary/16901-ecrad-new-radiation-scheme-ifs

A two-part paper is published in Journal of Geophysics Research
describing the SPARTACUS technique:

 - Schäfer, S. A. K., R. J. Hogan, C. Klinger, J.-C. Chiu and B. Mayer,
2016: Representing 3D cloud-radiation effects in two-stream schemes: 1. Longwave considerations and effective cloud edge length.
J. Geophys. Res., 121, 8567-8582.
http://www.met.reading.ac.uk/~swrhgnrj/publications/spartacus_part1.pdf

 - Hogan, R. J., S. A. K. Schäfer, C. Klinger, J.-C. Chiu and B. Mayer,
2016: Representing 3D cloud-radiation effects in two-stream schemes: 2. Matrix formulation and broadband evaluation. J. Geophys. Res., 121,
8583-8599.
http://www.met.reading.ac.uk/~swrhgnrj/publications/spartacus_part2.pdf

More recent developments on the shortwave SPARTACUS solver, available
since ecRad 1.1.10, are described here:

 - Hogan, R. J., M. D. Fielding, H. W. Barker, N. Villefranque and
S. A. K. Schäfer, 2019: Entrapment: An important mechanism to explain
the shortwave 3D radiative effect of clouds. J. Atmos. Sci., 76,
2123–2141.

The ecCKD gas optics scheme is described here:

 - Hogan, R. J., and M. Matricardi, 2022: a tool for generating fast
k-distribution gas-optics models for weather and climate
applications. J. Adv. Modeling Earth Sys., in review.


## CONTACT

Please email Robin Hogan <r.j.hogan@ecmwf.int> with any queries or bug
fixes, but note that ECMWF does not commit to providing support to
users of this software.
