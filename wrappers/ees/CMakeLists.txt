cmake_minimum_required(VERSION 3.0.0)

#######################################
#         PROJECT INFORMATION         #
#-------------------------------------#
# This CMakeLists.txt file is for the #
# EES user library of the CoolProp    #
# thermodynamic library.              #
#######################################

project("CoolProp EES Wrapper" VERSION ${COOLPROP_VERSION})


#######################################
#        ADD CUSTOM SOURCES           #
#-------------------------------------#
# The main folder should contain most #
# of the source files. Only custom    #
# sources should be added in this     # 
# CMakeLists.txt                      #
#######################################

# This file is only needed for EES, normal builds do not need it.
set(COOLPROP_EES_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp")

# Create source groups, mainly for Visual Studio
source_group("EES Sources" FILES ${COOLPROP_EES_SOURCES})


#######################################
#          SET DEPENDENCIES           #
#-------------------------------------#
# We make use of a number of other    #
# libraries for various purposes.     #
#######################################

#if(CMAKE_DL_LIBS)
#  find_package (${CMAKE_DL_LIBS} REQUIRED)
#endif()


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

add_library(coolprop_ees_library ${COOLPROP_OBJECT_SOURCES} ${COOLPROP_EES_SOURCES})
set_target_properties(coolprop_ees_library PROPERTIES OUTPUT_NAME "CoolProp")
target_link_libraries(coolprop_library ${CMAKE_DL_LIBS})
add_dependencies(coolprop_library generate_headers)

if(BUILD_SHARED_LIBS)
  set_property(TARGET coolprop_library PROPERTY POSITION_INDEPENDENT_CODE ON)
  target_compile_definitions(coolprop_library PUBLIC "-DCOOLPROP_LIB")
  target_include_directories(coolprop_library PRIVATE ${COOLPROP_OBJECT_INCLUDES}) # PRIVATE is only respected for shared libraries!
  target_include_directories(coolprop_library PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})  # for CoolPropLib.h
else()
  target_include_directories(coolprop_library PUBLIC ${COOLPROP_OBJECT_INCLUDES})
  target_include_directories(coolprop_library PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
endif()

if(MSVC)
  #set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS YES CACHE BOOL "Export all symbols")
  # Add postfix for debugging
  set_property(TARGET coolprop_library PROPERTY DEBUG_POSTFIX d)
  set_property(TARGET coolprop_library PROPERTY RELEASE_POSTFIX )
  # No lib prefix for the shared library
  set_property(TARGET coolprop_library PROPERTY PREFIX "")
  get_filename_component(VCBINPATH ${CMAKE_CXX_COMPILER} DIRECTORY)  
  add_custom_command(TARGET coolprop_library POST_BUILD
      COMMAND ${VCBINPATH}/dumpbin.exe /EXPORTS $<TARGET_FILE:coolprop_library> > $<TARGET_FILE_DIR:coolprop_library>/exports.txt
      COMMAND ${VCBINPATH}/dumpbin.exe /EXPORTS $<TARGET_FILE:coolprop_library>
  )
  #install(FILES ${CMAKE_CURRENT_BINARY_DIR}/exports.txt
  #  DESTINATION shared_library/${CMAKE_SYSTEM_NAME}/${BITNESS}bit${CONVENTION})
endif()
# For Linux
if("${CMAKE_SYSTEM_NAME}" MATCHES "Linux")
  set_property(TARGET coolprop_library PROPERTY VERSION ${COOLPROP_VERSION})
  set_property(TARGET coolprop_library PROPERTY SOVERSION ${COOLPROP_VERSION_MAJOR})
endif()
# For Mac systems
if("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
  # see https://support.enthought.com/hc/en-us/articles/204469410-OS-X-GCC-Clang-and-Cython-in-10-9-Mavericks
  set_property(TARGET coolprop_library APPEND_STRING PROPERTY COMPILE_FLAGS " -stdlib=libstdc++ -mmacosx-version-min=10.6")
  set_property(TARGET coolprop_library APPEND_STRING PROPERTY LINK_FLAGS " -stdlib=libstdc++ -mmacosx-version-min=10.6")
endif()




# install(
    # TARGETS coolprop_library
    # LIBRARY DESTINATION lib
    # ARCHIVE DESTINATION lib
    # RUNTIME DESTINATION bin
# )


# #
# # Installation section
# target_include_directories(IF97 INTERFACE  
    # $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    # $<INSTALL_INTERFACE:include>
# )
# set(include_install_dir "include")
# set(config_install_dir "lib/cmake/${PROJECT_NAME}")
# set(generated_dir "${CMAKE_CURRENT_BINARY_DIR}/generated")
# set(version_config "${generated_dir}/${PROJECT_NAME}ConfigVersion.cmake")
# set(project_config "${generated_dir}/${PROJECT_NAME}Config.cmake")
# set(TARGETS_EXPORT_NAME "${PROJECT_NAME}Targets")
# set(namespace "${PROJECT_NAME}::")
# include(CMakePackageConfigHelpers)

# # Use:
# #   * PROJECT_VERSION
# write_basic_package_version_file(
    # "${version_config}" COMPATIBILITY SameMajorVersion
# )

# # Use:
# #   * TARGETS_EXPORT_NAME
# #   * PROJECT_NAME
# configure_package_config_file(
    # "cmake/Config.cmake.in"
    # "${project_config}"
    # INSTALL_DESTINATION "${config_install_dir}"
# )

# install(
    # TARGETS IF97
    # EXPORT "${TARGETS_EXPORT_NAME}"
    # LIBRARY DESTINATION "lib"
    # ARCHIVE DESTINATION "lib"
    # RUNTIME DESTINATION "bin"
    # INCLUDES DESTINATION "${include_install_dir}"
# )

# install(
    # FILES IF97.h
    # DESTINATION "${include_install_dir}"
# )

# install(
    # FILES "${project_config}" "${version_config}"
    # DESTINATION "${config_install_dir}"
# )

# install(
    # EXPORT "${TARGETS_EXPORT_NAME}"
    # NAMESPACE "${namespace}"
    # DESTINATION "${config_install_dir}"
# )




# EES is only compiled for 32bit Windows -> Really?

  IF (NOT "${BITNESS}" STREQUAL "32")
    MESSAGE(FATAL_ERROR "You cannot build the EES wrapper as a 64-bit library.")
  ENDIF()
  # Prepare the sources
  include_directories(${APP_INCLUDE_DIRS})
  list (APPEND APP_SOURCES "wrappers/EES/main.cpp")
  list (APPEND APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/${COOLPROP_LIBRARY_SOURCE}")
  add_library(COOLPROP_EES SHARED ${APP_SOURCES})
  # Modify the target and add dependencies
  add_dependencies (COOLPROP_EES generate_headers)
  set_target_properties (COOLPROP_EES PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DCOOLPROP_LIB -DCONVENTION=__cdecl")
  set_target_properties (COOLPROP_EES PROPERTIES SUFFIX ".dlf" PREFIX "")
  # Creates "COOLPROP_EES.dlf"
  if ( NOT MSVC )
    set_target_properties (COOLPROP_EES PROPERTIES COMPILE_FLAGS "-m32" LINK_FLAGS "-m32")
  elseif ( MSVC )
    set_target_properties (COOLPROP_EES PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} )
    set_target_properties (COOLPROP_EES PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_BINARY_DIR})
    set_target_properties (COOLPROP_EES PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_BINARY_DIR})
    # etc for the other available configuration types (MinSizeRel, RelWithDebInfo)
  endif ()
  # copy required files
  add_custom_command (TARGET COOLPROP_EES PRE_BUILD
    COMMAND ${CMAKE_COMMAND} ARGS "-E" "copy" "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/EES/CoolProp.htm" "${CMAKE_CURRENT_BINARY_DIR}/."
	COMMAND ${CMAKE_COMMAND} ARGS "-E" "copy" "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/EES/CoolProp.LIB" "${CMAKE_CURRENT_BINARY_DIR}/."
	COMMAND ${CMAKE_COMMAND} ARGS "-E" "copy" "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/EES/CoolProp_EES_Sample.EES" "${CMAKE_CURRENT_BINARY_DIR}/."
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Copying the EES files to the build directory"
    VERBATIM
  )
  # install the generated library and the other files
  install(TARGETS COOLPROP_EES DESTINATION "${CMAKE_INSTALL_PREFIX}/EES/${CMAKE_SYSTEM_NAME}")
  install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/EES/CoolProp.htm" DESTINATION "${CMAKE_INSTALL_PREFIX}/EES/${CMAKE_SYSTEM_NAME}")
  install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/EES/CoolProp.LIB" DESTINATION "${CMAKE_INSTALL_PREFIX}/EES/${CMAKE_SYSTEM_NAME}")
  install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/wrappers/EES/CoolProp_EES_Sample.EES" DESTINATION "${CMAKE_INSTALL_PREFIX}/EES/${CMAKE_SYSTEM_NAME}")
