if( HAVE_PYECRAD AND DEFINED PYECRAD_PYTHON )
    set( _PYECRAD_VENV_BIN "${CMAKE_BINARY_DIR}/pyecrad_env/bin" )
    add_test(NAME pyecrad_test
        COMMAND ${CMAKE_COMMAND} -E env
            "VIRTUAL_ENV=${CMAKE_BINARY_DIR}/pyecrad_env"
            "PATH=${_PYECRAD_VENV_BIN}:${CMAKE_BINARY_DIR}/bin:$ENV{PATH}"
        ${PYECRAD_PYTHON} -m pytest -v
            ${CMAKE_CURRENT_SOURCE_DIR}/test_pyecrad.py
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
    set_tests_properties(pyecrad_test PROPERTIES
        DEPENDS "pyecrad;ecrad_dp"
    )
endif()
