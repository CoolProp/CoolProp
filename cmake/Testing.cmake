#######################################
#         TESTING MODULES             #
#-------------------------------------#
# Test executables and development    #
# tools for CoolProp                  #
#######################################

# Main executable module
if(COOLPROP_MAIN_MODULE)
  # Allow you to independently add back the testing CPP files
  if(COOLPROP_TEST)
    list(APPEND APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Tests/Tests.cpp")
    list(APPEND APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Tests/CoolProp-Tests.cpp")
  endif()

  list(APPEND APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cxx")
  add_executable(Main ${APP_SOURCES})
  add_dependencies(Main generate_headers)

  if(COOLPROP_TEST)
    set_target_properties(Main PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DENABLE_CATCH")
  endif()

  if(COOLPROP_IWYU)
    find_program(iwyu_path NAMES include-what-you-use iwyu)
    if(NOT iwyu_path)
      message(FATAL_ERROR "Could not find the program include-what-you-use")
    endif()
    set_property(TARGET Main PROPERTY CXX_INCLUDE_WHAT_YOU_USE ${iwyu_path})
  endif()

  if(UNIX)
    target_link_libraries(Main ${CMAKE_DL_LIBS})
  endif()
endif()

# Catch2 test runner
if(COOLPROP_CATCH_MODULE)
  add_subdirectory("externals/Catch2")
  enable_testing()

  list(APPEND APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Tests/CoolProp-Tests.cpp")

  # CATCH TEST, compile everything with catch and set test entry point
  add_executable(CatchTestRunner ${APP_SOURCES})
  add_dependencies(CatchTestRunner generate_headers)

  target_link_libraries(CatchTestRunner PRIVATE Catch2::Catch2WithMain)
  target_compile_definitions(CatchTestRunner PRIVATE ENABLE_CATCH)
  # Incbin is disabled for Catch2 because of a weird behavior where out-of-date zlib-compressed fluid information were being used.
  target_compile_definitions(CatchTestRunner PRIVATE COOLPROP_NO_INCBIN)

  if(UNIX)
    target_link_libraries(CatchTestRunner PRIVATE ${CMAKE_DL_LIBS})
  endif()

  target_include_directories(CatchTestRunner PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/externals/multicomplex/multicomplex/include")

  include(CTest)
  include(${CMAKE_CURRENT_SOURCE_DIR}/externals/Catch2/extras/Catch.cmake)

  if(NOT CMAKE_GENERATOR STREQUAL Xcode)
    # Test discovery doesn't work in Xcode, due to a signing bug in Xcode which causes discovery to fail: https://github.com/catchorg/Catch2/issues/2411
    catch_discover_tests(CatchTestRunner DISCOVERY_MODE PRE_TEST)
  endif()

  if(COOLPROP_IWYU)
    find_program(iwyu_path NAMES include-what-you-use iwyu)
    if(NOT iwyu_path)
      message(FATAL_ERROR "Could not find the program include-what-you-use")
    endif()
    set_property(TARGET CatchTestRunner PROPERTY CXX_INCLUDE_WHAT_YOU_USE ${iwyu_path})
  endif()

  if(COOLPROP_LAZY_LOAD_SUPERANCILLARIES)
    target_compile_definitions(CatchTestRunner PRIVATE LAZY_LOAD_SUPERANCILLARIES)
  else()
    # lazy load superancillaries by default in debug mode
    target_compile_definitions(CatchTestRunner PUBLIC $<$<CONFIG:Debug>:LAZY_LOAD_SUPERANCILLARIES>)
  endif()
endif()

# C++ documentation example test
if(COOLPROP_CPP_EXAMPLE_TEST)
  add_executable(docuTest.exe "Web/examples/C++/Example.cpp")
  add_dependencies(docuTest.exe CoolProp)
  target_link_libraries(docuTest.exe CoolProp)
  if(UNIX)
    target_link_libraries(docuTest.exe ${CMAKE_DL_LIBS})
  endif()
  add_test(DocumentationTest docuTest.exe)
endif()

# Code snippets builder
if(COOLPROP_SNIPPETS)
  list(APPEND APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/${COOLPROP_LIBRARY_SOURCE}")

  # Make the static library with which the snippets will be linked
  add_library(CoolProp STATIC ${APP_SOURCES})
  add_dependencies(CoolProp generate_headers)
  set_property(TARGET CoolProp APPEND_STRING PROPERTY COMPILE_FLAGS " -DEXTERNC")

  # Collect all the snippets
  file(GLOB_RECURSE snippets "${CMAKE_CURRENT_SOURCE_DIR}/Web/coolprop/snippets/*.cxx")

  message(STATUS "snippets found = ${snippets}")
  foreach(snippet ${snippets})
    get_filename_component(snippet_name ${snippet} NAME)
    get_filename_component(snippet_exe ${snippet} NAME_WE)
    message(STATUS "snippet_name = ${snippet_name}")

    add_executable(${snippet_exe} ${snippet})
    add_dependencies(${snippet_exe} CoolProp)
    target_link_libraries(${snippet_exe} CoolProp)
    if(UNIX)
      target_link_libraries(${snippet_exe} ${CMAKE_DL_LIBS})
    endif()

    if(MSVC)
      set_target_properties(${snippet_exe} PROPERTIES
                            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin
                            RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_BINARY_DIR}/bin
                            RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_BINARY_DIR}/bin)
      set(BIN_PATH "${CMAKE_CURRENT_BINARY_DIR}/bin")
    else()
      set(BIN_PATH "${CMAKE_CURRENT_BINARY_DIR}")
    endif()

    set_property(TARGET ${snippet_exe} APPEND_STRING PROPERTY COMPILE_FLAGS " -DEXTERNC")

    # Run it and save the output to a file with .output appended
    add_custom_command(TARGET ${snippet_exe} POST_BUILD
                       COMMAND ${BIN_PATH}/${snippet_exe} >
                               ${CMAKE_CURRENT_SOURCE_DIR}/Web/coolprop/snippets/${snippet_name}.output)
  endforeach()
endif()

# Clang address sanitizer build
if(COOLPROP_CLANG_ADDRESS_SANITIZER)
  set(CMAKE_CXX_FLAGS "-fsanitize=address -g")
  list(APPEND APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Tests/catch_always_return_success.cxx")

  # CATCH TEST, compile everything with catch and set test entry point
  add_executable(CatchTestRunner ${APP_SOURCES})
  add_dependencies(CatchTestRunner generate_headers)
  set_target_properties(CatchTestRunner PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -DENABLE_CATCH")

  set(CMAKE_CXX_FLAGS "-O1")
  set(CMAKE_EXE_LINKER_FLAGS "-fsanitize=address -fno-omit-frame-pointer -lstdc++")

  if(UNIX)
    target_link_libraries(CatchTestRunner ${CMAKE_DL_LIBS})
  endif()

  add_custom_command(TARGET CatchTestRunner POST_BUILD
                     COMMAND ${CMAKE_CURRENT_BINARY_DIR}/CatchTestRunner)
endif()

# Profiling support
if(COOLPROP_PROFILE)
  if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-g -O2")
    set(CMAKE_C_FLAGS "-g -O2")
  endif()
endif()

# Code coverage support
if(COOLPROP_COVERAGE)
  if(CMAKE_COMPILER_IS_GNUCXX)
    # See also http://stackoverflow.com/a/16536401 (detailed guide on using gcov with cmake)
    include(CodeCoverage)
    set(CMAKE_CXX_FLAGS "-g -O0 -fprofile-arcs -ftest-coverage")
    set(CMAKE_C_FLAGS "-g -O0 -fprofile-arcs -ftest-coverage")
    setup_target_for_coverage(CoolProp_coverage Main coverage)
  endif()
endif()
