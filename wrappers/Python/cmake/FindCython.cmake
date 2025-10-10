# Find Cython
# This module finds the Cython compiler and creates an imported target

find_package(Python COMPONENTS Interpreter REQUIRED)

execute_process(
    COMMAND ${Python_EXECUTABLE} -c "import cython; print(cython.__version__)"
    OUTPUT_VARIABLE CYTHON_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
    RESULT_VARIABLE CYTHON_NOT_FOUND
)

if(CYTHON_NOT_FOUND)
    set(Cython_FOUND FALSE)
else()
    set(Cython_FOUND TRUE)

    # Get path to cython executable
    execute_process(
        COMMAND ${Python_EXECUTABLE} -c "import sys; import cython; print(sys.executable)"
        OUTPUT_VARIABLE CYTHON_EXECUTABLE
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    if(NOT TARGET Cython::Cython)
        add_executable(Cython::Cython IMPORTED)
        set_target_properties(Cython::Cython PROPERTIES
            IMPORTED_LOCATION "${Python_EXECUTABLE}"
        )
        # Set command to run cython via python -m
        set(CYTHON_COMMAND "${Python_EXECUTABLE}" -m cython)
    endif()

    message(STATUS "Found Cython: ${CYTHON_VERSION}")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Cython
    REQUIRED_VARS Cython_FOUND
    VERSION_VAR CYTHON_VERSION
)
