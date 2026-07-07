# Create destination directory
file(MAKE_DIRECTORY "${DST_DIR}")

# Find all .java files in the source directory
file(GLOB JAVA_FILES "${SRC_DIR}/*.java")

foreach(f ${JAVA_FILES})

    # Skip Example.java
    if(NOT f STREQUAL "${SRC_DIR}/Example.java")
	
        # Extract just the filename (e.g., AbstractState.java)
        get_filename_component(fname "${f}" NAME)
        file(RENAME "${f}" "${DST_DIR}/${fname}")
    endif()
endforeach()
