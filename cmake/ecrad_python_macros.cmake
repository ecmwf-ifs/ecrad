# (C) Copyright 2024- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

##############################################################################
#.rst:
#
# ecrad_find_python_venv
# ======================
#
# Call ``find_package( Python3 )``, making sure to discover a specific
# virtual environment at the given location ``VENV_PATH``::
#
#   ecrad_find_python_venv( VENV_PATH <path> [ PYTHON_VERSION <version str> ] )
#
# Options
# -------
# :VENV_PATH: The path to the virtual environment
# :PYTHON_VERSION: Optional specification of permissible Python versions for find_package
#
# Output variables
# ----------------
# :Python3_FOUND:       Exported into parent scope from FindPython3
# :Python3_EXECUTABLE:  Exported into parent scope from FindPython3
# :Python3_VENV_BIN:    The path to the virtual environment's `bin` directory
# :ENV{VIRTUAL_ENV}:    Environment variable with the virtual environment directory,
#                       emulating the activate script
#
##############################################################################

function( ecrad_find_python_venv )

    set( options "" )
    set( oneValueArgs VENV_PATH PYTHON_VERSION )
    set( multiValueArgs "" )

    cmake_parse_arguments( _PAR "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    if( _PAR_UNPARSED_ARGUMENTS )
        message( FATAL_ERROR "Unknown keywords given to ecrad_find_python_venv(): \"${_PAR_UNPARSED_ARGUMENTS}\"" )
    endif()

    if( NOT _PAR_VENV_PATH )
        message( FATAL_ERROR "No VENV_PATH provided to ecrad_find_python_venv()" )
    endif()

    # Update the environment with VIRTUAL_ENV variable (mimic the activate script)
    set( ENV{VIRTUAL_ENV} "${_PAR_VENV_PATH}" )

    # Change the context of the search to only find the venv
    set( Python3_FIND_VIRTUALENV ONLY )

    # Unset Python3_EXECUTABLE because it is also an input variable
    # see https://cmake.org/cmake/help/latest/module/FindPython.html#artifacts-specification
    unset( Python3_EXECUTABLE )
    # To allow cmake to discover the newly created venv if Python3_ROOT_DIR
    # was passed as an argument at build-time
    set( Python3_ROOT_DIR "${_PAR_VENV_PATH}" )

    # Launch a new search
    find_package( Python3 ${_PAR_PYTHON_VERSION} COMPONENTS Interpreter REQUIRED )

    # Find the binary directory of the virtual environment
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -c "import sys; import os.path; print(os.path.dirname(sys.executable), end='')"
        OUTPUT_VARIABLE Python3_VENV_BIN
    )

    # Forward variables to parent scope
    foreach ( _VAR_NAME Python3_FOUND Python3_EXECUTABLE Python3_VENV_BIN )
        set( ${_VAR_NAME} ${${_VAR_NAME}} PARENT_SCOPE )
    endforeach()

endfunction()

##############################################################################
#.rst:
#
# ecrad_create_python_venv
# ========================
#
# Discover a Python 3 interpreter and create a virtual environment at the
# specified location ``VENV_PATH``. ::
#
#   ecrad_create_python_venv( VENV_PATH <path> [ PYTHON_VERSION <version str> ] [ INSTALL_VENV ] )
#
# Installation procedure
# ----------------------
#
# Create a virtual environment at the given location (`VENV_PATH`)
#
# Options
# -------
#
# :VENV_PATH: The path to use for the virtual environment
# :PYTHON_VERSION: Optional specification of permissible Python versions for find_package
# :INSTALL_VENV: If provided, an equivalent virtual environment will also be created in
#                ``${CMAKE_INSTALL_PREFIX}/var/${VENV_NAME}`` upon installation
#
##############################################################################

function( ecrad_create_python_venv )

    set( options INSTALL_VENV )
    set( oneValueArgs VENV_NAME PYTHON_VERSION )
    set( multiValueArgs "" )

    cmake_parse_arguments( _PAR "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    if( _PAR_UNPARSED_ARGUMENTS )
        message( FATAL_ERROR "Unknown keywords given to ecrad_create_python_venv(): \"${_PAR_UNPARSED_ARGUMENTS}\"" )
    endif()

    if( NOT _PAR_VENV_NAME )
        message( FATAL_ERROR "No VENV_NAME provided to ecrad_create_python_venv()" )
    endif()

    set( VENV_PATH "${CMAKE_CURRENT_BINARY_DIR}/${_PAR_VENV_NAME}" )

    # Discover only system install Python 3
    set( Python3_FIND_VIRTUALENV STANDARD )
    find_package( Python3 ${_PAR_PYTHON_VERSION} COMPONENTS Interpreter REQUIRED )

    # Ensure the activate script is writable in case the venv exists already
    foreach( activate_script activate activate.csh activate.fish Activate.ps1 )
        if( EXISTS "${VENV_PATH}/bin/${activate_script}" )
            file( CHMOD "${VENV_PATH}/bin/${activate_script}" FILE_PERMISSIONS OWNER_READ OWNER_WRITE )
        endif()
    endforeach()

    # Create a virtualenv
    message( STATUS "Create Python virtual environment ${VENV_PATH}" )
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -m venv "${VENV_PATH}"
        OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY
    )

    # Upon installation, we create an equivalent Python venv in the installation directory
    if( _PAR_INSTALL_VENV )
        install(
            CODE
                "execute_process( COMMAND ${Python3_EXECUTABLE} -m venv \${CMAKE_INSTALL_PREFIX}/var/${_PAR_VENV_NAME} RESULT_VARIABLE _RET OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY )"
        )
        set( Python3_INSTALL_VENV "\${CMAKE_INSTALL_PREFIX}/var/${_PAR_VENV_NAME}" PARENT_SCOPE )
    endif()

endfunction()

##############################################################################
#.rst:
#
# ecrad_setup_python_venv
# =======================
#
# Find Python 3, create a virtual environment and make it available. ::
#
#   ecrad_setup_python_venv( VENV_PATH <path> [ PYTHON_VERSION <version str> ] [ INSTALL_VENV ] )
#
# It combines calls to ``ecrad_create_python_venv`` and ``ecrad_find_python_venv``
#
# Options
# -------
#
# :VENV_PATH: The path to use for the virtual environment
# :PYTHON_VERSION: Optional specification of permissible Python versions for find_package
# :INSTALL_VENV: If provided, an equivalent virtual environment will also be created in
#                ``${CMAKE_INSTALL_PREFIX}/var/${VENV_NAME}`` upon installation
#
# Output variables
# ----------------
# :Python3_FOUND:        Exported into parent scope from FindPython3
# :Python3_EXECUTABLE:   Exported into parent scope from FindPython3
# :Python3_VENV_BIN:     The path to the virtual environment's `bin` directory
# :Python3_INSTALL_VENV: The path to the virtual environment in the install directory,
#                        if INSTALL_VENV has been provided.
# :ENV{VIRTUAL_ENV}:     Environment variable with the virtual environment directory,
#                        emulating the activate script
#
##############################################################################

function( ecrad_setup_python_venv )

    set( options INSTALL_VENV )
    set( oneValueArgs VENV_NAME PYTHON_VERSION )
    set( multiValueArgs "" )

    cmake_parse_arguments( _PAR "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    if( _PAR_UNPARSED_ARGUMENTS )
        message( FATAL_ERROR "Unknown keywords given to ecrad_setup_python_venv(): \"${_PAR_UNPARSED_ARGUMENTS}\"" )
    endif()

    if( NOT _PAR_VENV_NAME )
        message( FATAL_ERROR "No VENV_NAME provided to ecrad_setup_python_venv()" )
    endif()

    # Create the virtual environment
    set( _ARGS VENV_NAME "${_PAR_VENV_NAME}" )
    if( _PAR_PYTHON_VERSION )
        list( APPEND _ARGS PYTHON_VERSION "${_PAR_PYTHON_VERSION}" )
    endif()
    if( _PAR_INSTALL_VENV )
        list( APPEND _ARGS INSTALL_VENV )
    endif()

    ecrad_create_python_venv( ${_ARGS} )

    if( Python3_INSTALL_VENV )
        set( Python3_INSTALL_VENV "${Python3_INSTALL_VENV}" PARENT_SCOPE )
    endif()

    # Discover Python in the virtual environment and set-up variables
    set( _ARGS VENV_PATH "${CMAKE_CURRENT_BINARY_DIR}/${_PAR_VENV_NAME}" )
    if( _PAR_PYTHON_VERSION )
        list( APPEND _ARGS PYTHON_VERSION "${_PAR_PYTHON_VERSION}" )
    endif()
    ecrad_find_python_venv( ${_ARGS} )

    foreach ( _VAR_NAME Python3_FOUND Python3_EXECUTABLE Python3_VENV_BIN )
        set( ${_VAR_NAME} ${${_VAR_NAME}} PARENT_SCOPE )
    endforeach()

endfunction()

##############################################################################
#.rst:
#
# ecrad_download_python_wheels
# ============================
#
# Download all dependencies for the given ``REQUIREMENT_SPEC`` and cache them in a
# wheelhouse at ``WHEELS_DIR``
#
#   ecrad_download_python_wheels( REQUIREMENT_SPEC <spec> [ WHEELS_DIR <path> ] [ PYTHON_VERSION <version str> ] )
#
# Implementation note
# -------------------
#
# This function does intentionally not expose all PIP options directly because the PIP command line
# interface allows to specify option values via environment variables. These can therefore be used
# to further control the PIP behaviour, see https://pip.pypa.io/en/stable/cli/pip_download/
#
# Because PIP does not provide a mechanism for downloading PEP 518 build dependencies,
# this function builds the wheel also for the provided REQUIREMENT_SPEC instead of only downloading
# the required dependencies. See https://github.com/pypa/pip/issues/7863 for details.
# To provide a sane minimum, setuptools and wheel packages are always downloaded.
#
# The provided PYTHON_VERSION is used to discover a Python interpreter matching the version
# specification when calling pip. To download wheels for specific platforms or Python versions,
# use the PIP_PLATFORM, PIP_PYTHON_VERSION, PIP_IMPLEMENTATION, or PIP_ABI environment variables.
#
# It is safe to call this function during an offline build, as long as all wheels are already
# available in the wheelhouse. A dry-run call to ``pip install`` is used to determine the need
# for any wheel downloads before executing the ``pip download`` command.
#
# Options
# -------
#
# :REQUIREMENT_SPEC: The requirement spec as given to ``pip download`` and ``pip wheel``
# :WHEELS_DIR: The path of the wheelhouse directory to cache the wheels. Defaults to
#              ``${CMAKE_CURRENT_BINARY_DIR}/wheelhouse``
# :PYTHON_VERSION: Optional specification of permissible Python versions for find_package
#
##############################################################################

function( ecrad_download_python_wheels )

    set( options "" )
    set( oneValueArgs WHEELS_DIR PYTHON_VERSION )
    set( multiValueArgs REQUIREMENT_SPEC )

    cmake_parse_arguments( _PAR "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    if( _PAR_UNPARSED_ARGUMENTS )
        message( FATAL_ERROR "Unknown keywords given to ecrad_download_python_wheels(): \"${_PAR_UNPARSED_ARGUMENTS}\"" )
    endif()

    if( NOT _PAR_REQUIREMENT_SPEC )
        message( FATAL_ERROR "No REQUIREMENT_SPEC provided to ecrad_download_python_wheels()" )
    endif()

    message( STATUS "Checking for cached wheels in ${WHEELS_DIR}" )

    # Check for a suitable python interpreter
    find_package( Python3 ${_PAR_PYTHON_VERSION} COMPONENTS Interpreter REQUIRED QUIET )

    # If no wheelhouse dir is given, create one in the current binary directory
    if( _PAR_WHEELS_DIR )
        set( WHEELS_DIR "${_PAR_WHEELS_DIR}" )
    else()
        set( WHEELS_DIR "${CMAKE_CURRENT_BINARY_DIR}/wheelhouse" )
    endif()
    file( MAKE_DIRECTORY "${WHEELS_DIR}" )

    # Ensure PIP is available
    execute_process( COMMAND ${Python3_EXECUTABLE} -m ensurepip -U )

    # We use a dry-run installation to check if all dependencies have already been downloaded
    execute_process(
        COMMAND
            ${Python3_EXECUTABLE} -m pip install
                --dry-run --break-system-packages
                --no-index --find-links "${WHEELS_DIR}" --only-binary :all:
                ${_PAR_REQUIREMENT_SPEC}
        OUTPUT_QUIET ERROR_QUIET
        RESULT_VARIABLE _RET_VAL
    )

    if( "${_RET_VAL}" EQUAL "0" )

        message( STATUS "All dependency wheels for ${_PAR_REQUIREMENT_SPEC} found in cache" )

    else()

        message( STATUS "Downloading dependency wheels for ${_PAR_REQUIREMENT_SPEC} to ${WHEELS_DIR}" )

        # Download typical build dependencies for wheels: setuptools and wheel
        execute_process(
            COMMAND
                ${Python3_EXECUTABLE} -m pip download
                --disable-pip-version-check --dest "${WHEELS_DIR}"
                setuptools wheel
            OUTPUT_QUIET
            COMMAND_ERROR_IS_FATAL ANY
        )

        # Download dependencies for the specified REQUIREMENT_SPEC
        execute_process(
            COMMAND
                ${Python3_EXECUTABLE} -m pip download
                --disable-pip-version-check --dest "${WHEELS_DIR}"
                ${_PAR_REQUIREMENT_SPEC}
            OUTPUT_QUIET
            COMMAND_ERROR_IS_FATAL ANY
        )

        # Here we _build_ the package instead of just downloading its build dependencies. Sadly, this is necessary because
        # PIP does not yet provide a mechanism to download the build dependencies for PEP 518 packages.
        # See https://github.com/pypa/pip/issues/7863 for details.
        # When this is resolved, we should instead download only build dependencies here, which will defer the actual wheel
        # building to the `build_python_wheel` function
        execute_process(
            COMMAND
                ${Python3_EXECUTABLE} -m pip wheel
                    --disable-pip-version-check --wheel-dir "${WHEELS_DIR}"
                    ${_PAR_REQUIREMENT_SPEC}
            OUTPUT_QUIET
            COMMAND_ERROR_IS_FATAL ANY
        )

    endif()

endfunction()

##############################################################################
#.rst:
#
# build_python_wheel
# ==================
#
# Build a Python wheel for the given ``REQUIREMENT_SPEC`` and store it in the
# specified ``BUILD_DIR``. This uses no online sources to download packages,
# any required dependencies must be available in ``WHEELS_DIR``.
# Use ``download_python_wheels`` to make them available if necessary.
#
#   build_python_wheels( REQUIREMENT_SPEC <spec> [ BUILD_DIR <path> ] [ WHEELS_DIR <path> ] [ PYTHON_VERSION <version str> ] )
#
# Implementation note
# -------------------
#
# This function does intentionally not expose all PIP options directly because the PIP command line
# interface allows to specify option values via environment variables. These can therefore be used
# to further control the PIP behaviour, see https://pip.pypa.io/en/stable/cli/pip_download/
#
# The provided PYTHON_VERSION is used to discover a Python interpreter matching the version
# specification when calling pip. To build wheels for specific platforms or Python versions,
# use the PIP_PLATFORM, PIP_PYTHON_VERSION, PIP_IMPLEMENTATION, or PIP_ABI environment variables.
#
# Options
# -------
#
# :REQUIREMENT_SPEC: The requirement spec as given to ``pip download`` and ``pip wheel``
# :BUILD_DIR: The path to store the built wheel. Defaults to ``${CMAKE_CURRENT_BINARY_DIR}/wheelhouse``
# :WHEELS_DIR: The path of the wheelhouse directory where to look for cached wheels. Defaults to ``BUILD_DIR``
# :PYTHON_VERSION: Optional specification of permissible Python versions for find_package
#
##############################################################################

function( build_python_wheel )

    set( options "" )
    set( oneValueArgs WHEELS_DIR BUILD_DIR )
    set( multiValueArgs REQUIREMENT_SPEC )

    cmake_parse_arguments( _PAR "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    if( _PAR_UNPARSED_ARGUMENTS )
        message( FATAL_ERROR "Unknown keywords given to build_python_wheels(): \"${_PAR_UNPARSED_ARGUMENTS}\"" )
    endif()

    if( NOT _PAR_REQUIREMENT_SPEC )
        message( FATAL_ERROR "No REQUIREMENT_SPEC provided to build_python_wheels()" )
    endif()

    message( STATUS "Building wheel for ${REQUIREMENT_SPEC}" )

    # Check for a suitable python interpreter
    find_package( Python3 ${_PAR_PYTHON_VERSION} COMPONENTS Interpreter REQUIRED QUIET )

    # If no build dir is given, create one in the current binary directory
    if( _PAR_BUILD_DIR )
        set( BUILD_DIR "${_PAR_BUILD_DIR}" )
    else()
        set( BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/wheelhouse" )
    endif()
    file( MAKE_DIRECTORY "${BUILD_DIR}" )

    # If no wheelhouse is given, use the build directory
    if( _PAR_WHEELS_DIR )
        set( WHEELS_DIR "${_PAR_WHEELS_DIR}" )
    else()
        set( WHEELS_DIR "${BUILD_DIR}" )
    endif()
    file( MAKE_DIRECTORY "${WHEELS_DIR}" )

    # Ensure PIP is available
    execute_process( COMMAND ${Python3_EXECUTABLE} -m ensurepip -U )

    execute_process(
        COMMAND
            ${Python3_EXECUTABLE} -m pip wheel
                --no-index --find-links "${WHEELS_DIR}" --wheel-dir "${BUILD_DIR}"
                ${_PAR_REQUIREMENT_SPEC}
            OUTPUT_QUIET
            COMMAND_ERROR_IS_FATAL ANY
    )

endfunction()

##############################################################################
#.rst:
#
# install_python_package
# ======================
#
# Install a Python package using the provided ``REQUIREMENT_SPEC``.
#
#   install_python_package( REQUIREMENT_SPEC <spec> [ WHEELS_DIR <path> ] [ INSTALL ] )
#
# This assumes that the ``Python3_EXECUTABLE`` has been made available to use, e.g.,
# via a ``find_package( Python3 )``, ``find_python_venv()`` or ``setup_python_venv()``.
#
# By default this will search for the package and its dependencies in the
# standard package index.
# Providing a wheelhouse ``WHEELS_DIR`` ensures that this installation is
# an offline operation, taking wheels only from the provided path.
# If required, these can be fetched explicitly via ``download_python_wheels``.
#
# If Python3_INSTALL_VENV variable is set, the package will also be installed
# into the virtual environment at installation time.
#
# Implementation note
# -------------------
#
# This function does intentionally not expose all PIP options directly because the PIP command line
# interface allows to specify option values via environment variables. These can therefore be used
# to further control the PIP behaviour, see https://pip.pypa.io/en/stable/cli/pip_download/
#
# The provided PYTHON_VERSION is used to discover a Python interpreter matching the version
# specification when calling pip. To build wheels for specific platforms or Python versions,
# use the PIP_PLATFORM, PIP_PYTHON_VERSION, PIP_IMPLEMENTATION, or PIP_ABI environment variables.
#
# Options
# -------
#
# :REQUIREMENT_SPEC: The requirement spec as given to ``pip download`` and ``pip wheel``
# :WHEELS_DIR: The path of the wheelhouse directory where to look for cached wheels.
#
##############################################################################
function( ecrad_install_python_package )

    set( options "" )
    set( oneValueArgs WHEELS_DIR )
    set( multiValueArgs REQUIREMENT_SPEC )

    cmake_parse_arguments( _PAR "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    if( _PAR_UNPARSED_ARGUMENTS )
        message( FATAL_ERROR "Unknown keywords given to ecrad_install_python_package(): \"${_PAR_UNPARSED_ARGUMENTS}\"" )
    endif()

    if( NOT _PAR_REQUIREMENT_SPEC )
        message( FATAL_ERROR "No REQUIREMENT_SPEC provided to install_python_package()" )
    endif()

    # Check for a suitable python interpreter
    find_package( Python3 ${_PAR_PYTHON_VERSION} COMPONENTS Interpreter REQUIRED QUIET )

    if( _PAR_WHEELS_DIR )
        # Force installation from provided wheelhouse
        set( INSTALL_OPTS "--no-index --find-links=${_PAR_WHEELS_DIR}" )
    else()
        # Default pip install
        set( INSTALL_OPTS "--disable-pip-version-check" )
    endif()

    # Ensure PIP is available
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -m ensurepip -U
        OUTPUT_QUIET
        COMMAND_ERROR_IS_FATAL ANY
    )

    message( STATUS "Installing Python package ${_PAR_REQUIREMENT_SPEC}" )

    # Install package
    execute_process( COMMAND
        ${Python3_EXECUTABLE} -m pip install ${INSTALL_OPTS} ${_PAR_REQUIREMENT_SPEC}
        OUTPUT_QUIET
        COMMAND_ERROR_IS_FATAL ANY
    )

    # Upon installation, repeat the installation
    if( Python3_INSTALL_VENV )
        install(
            CODE
                "execute_process( COMMAND ${Python3_INSTALL_VENV}/bin/python -m ensurepip -U OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY )"
        )
        install(
            CODE
                "execute_process( COMMAND ${Python3_INSTALL_VENV}/bin/python -m pip install ${INSTALL_OPTS} ${_PAR_REQUIREMENT_SPEC} OUTPUT_QUIET COMMAND_ERROR_IS_FATAL ANY )"
        )
    endif()

endfunction()
