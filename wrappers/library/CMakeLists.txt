cmake_minimum_required(VERSION 3.0.0)

#######################################
#         PROJECT INFORMATION         #
#-------------------------------------#
# This CMakeLists.txt file is for the #
# shared and static library of the    #
# CoolProp thermodynamic library.     #
#######################################

project("CoolProp Library" VERSION ${COOLPROP_VERSION})


#######################################
#        ADD CUSTOM SOURCES           #
#-------------------------------------#
# The main folder should contain most #
# of the source files. Only custom    #
# sources should be added in this     # 
# CMakeLists.txt                      #
#######################################

# This file is only needed for the library, normal builds do not need it.
set(COOLPROP_LIBRARY_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/CoolPropLib.cpp")
set(COOLPROP_LIBRARY_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/CoolPropLib.h")
set(COOLPROP_LIBRARY_EXPORTS "${CMAKE_CURRENT_SOURCE_DIR}/CoolPropLib.def")

# Create source groups, mainly for Visual Studio
source_group("Library Sources" FILES ${COOLPROP_LIBRARY_SOURCES} ${COOLPROP_LIBRARY_HEADERS} ${COOLPROP_LIBRARY_EXPORTS})


#######################################
#          SET DEPENDENCIES           #
#-------------------------------------#
# We make use of a number of other    #
# libraries for various purposes.     #
#######################################

if(CMAKE_DL_LIBS)
  find_package (${CMAKE_DL_LIBS} REQUIRED)
endif()


#######################################
#       CREATE DEFAULT TARGET         #
#-------------------------------------#
# All targets should be placed in     #
# subdirectories and they should also #
# have unique names. However, we can  #
# define one single default target to #
# simplify integration.               #
#######################################

# Note that switching between static 
# and shared library should be done by 
# turning on the CMake option 
# `BUILD_SHARED_LIBS` only. 
# Use `CMAKE_INSTALL_PREFIX` to install 
# shared and static library to different locations.

# [library-examples]> cmake -Hright-way -B_builds/shared -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX="`pwd`/_install/configuration-A"
# [library-examples]> cmake --build _builds/shared --target install

add_library(COOLPROP_LIBRARY $<TARGET_OBJECTS:COOLPROP_OBJECT_LIBRARY> ${COOLPROP_LIBRARY_SOURCES} ${COOLPROP_LIBRARY_HEADERS} ${COOLPROP_LIBRARY_EXPORTS})
set_target_properties(COOLPROP_LIBRARY PROPERTIES OUTPUT_NAME "CoolProp")
target_link_libraries(COOLPROP_LIBRARY ${CMAKE_DL_LIBS})

install(
    TARGETS COOLPROP_LIBRARY
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
)

if(BUILD_SHARED_LIBS)
  target_compile_definitions(COOLPROP_LIBRARY PUBLIC "-DCOOLPROP_LIB")
endif()

if(MSVC)
  #set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS YES CACHE BOOL "Export all symbols")
  # Add postfix for debugging
  set_property(TARGET COOLPROP_LIBRARY PROPERTY DEBUG_POSTFIX d)
  set_property(TARGET COOLPROP_LIBRARY PROPERTY RELEASE_POSTFIX )
  # No lib prefix for the shared library
  set_property(TARGET COOLPROP_LIBRARY PROPERTY PREFIX "")
  #add_custom_command(TARGET COOLPROP_LIBRARY POST_BUILD
  #  COMMAND dumpbin /EXPORTS $<TARGET_FILE:COOLPROP_LIBRARY> > $<TARGET_DIR:COOLPROP_LIBRARY>/exports.txt)
  #install(FILES ${CMAKE_CURRENT_BINARY_DIR}/exports.txt
  #  DESTINATION shared_library/${CMAKE_SYSTEM_NAME}/${BITNESS}bit${CONVENTION})
endif()
# For Linux
if("${CMAKE_SYSTEM_NAME}" MATCHES "Linux")
  set_property(TARGET COOLPROP_LIBRARY PROPERTY VERSION ${COOLPROP_VERSION})
  set_property(TARGET COOLPROP_LIBRARY PROPERTY SOVERSION ${COOLPROP_VERSION_MAJOR})
endif()
# For Mac systems
if("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
  # see https://support.enthought.com/hc/en-us/articles/204469410-OS-X-GCC-Clang-and-Cython-in-10-9-Mavericks
  set_property(TARGET COOLPROP_LIBRARY APPEND_STRING PROPERTY COMPILE_FLAGS " -stdlib=libstdc++ -mmacosx-version-min=10.6")
  set_property(TARGET COOLPROP_LIBRARY APPEND_STRING PROPERTY LINK_FLAGS " -stdlib=libstdc++ -mmacosx-version-min=10.6")
endif()