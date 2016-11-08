# - Functions to handle compiler flags
# This module provides functions to modify the compiler flags. It simplifies
# working with these flags, but it should be used with caution.

# http://www.cmake.org/Wiki/CMake_FAQ#Dynamic_Replace.
# http://stackoverflow.com/questions/18233513/cmake-g-set-cxx-flags-02-but-03-is-still-there

#######################################
#        FUNCTION DEFINITIONS         #
#-------------------------------------#
# Define some macros that simplify    #
# working with the CMakeLists file.   #
# We really should use subfolders to  #
# slim this file, but for now macros  #
# can help us a little.               #
#######################################

set(debug_c_flags CMAKE_C_FLAGS_DEBUG)
set(debug_cxx_flags CMAKE_CXX_FLAGS_DEBUG)

set(release_c_flags CMAKE_C_FLAGS_RELEASE CMAKE_C_FLAGS_MINSIZEREL CMAKE_C_FLAGS_RELWITHDEBINFO)
set(release_cxx_flags CMAKE_CXX_FLAGS_RELEASE CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)

set(all_c_flags CMAKE_C_FLAGS ${debug_c_flags} ${release_c_flags})
set(all_c_flags CMAKE_CXX_FLAGS ${debug_cxx_flags} ${release_cxx_flags})

function(remove_compiler_flag_release flag)
  foreach (flag_var IN LISTS release_c_flags release_cxx_flags)
    string(REPLACE "${flag}" "" ${flag_var} "${${flag_var}}")
    set(${flag_var} "${${flag_var}}" PARENT_SCOPE)
  endforeach()
endfunction()

function(remove_compiler_flag_debug flag)
  foreach (flag_var IN LISTS debug_c_flags debug_cxx_flags)
    string(REPLACE "${flag}" "" ${flag_var} "${${flag_var}}")
    set(${flag_var} "${${flag_var}}" PARENT_SCOPE)
  endforeach()
endfunction()

function(add_compiler_flag_release flag)
  foreach (flag_var IN LISTS release_c_flags release_cxx_flags)
    set(${flag_var} "${${flag_var}} ${flag}" PARENT_SCOPE)
  endforeach()
endfunction()

function(add_compiler_flag_debug flag)
  foreach (flag_var IN LISTS debug_c_flags debug_cxx_flags)
    set(${flag_var} "${${flag_var}} ${flag}" PARENT_SCOPE)
  endforeach()
endfunction()

function(replace_compiler_flag_release flag_old flag_new)
  foreach (flag_var IN LISTS release_c_flags release_cxx_flags)
    if(${flag_var} MATCHES "${flag_old}")
       string(REGEX REPLACE "${flag_old}" "${flag_new}" ${flag_var} "${${flag_var}}")
    else(${flag_var} MATCHES "${flag_old}")
       set(${flag_var} "${${flag_var}} ${flag_new}")
    endif(${flag_var} MATCHES "${flag_old}")
    set(${flag_var} "${${flag_var}}" PARENT_SCOPE)
  endforeach()
endfunction()

function(replace_compiler_flag_debug flag_old flag_new)
  foreach (flag_var IN LISTS debug_c_flags debug_cxx_flags)
    if(${flag_var} MATCHES "${flag_old}")
       string(REGEX REPLACE "${flag_old}" "${flag_new}" ${flag_var} "${${flag_var}}")
    else(${flag_var} MATCHES "${flag_old}")
       set(${flag_var} "${${flag_var}} ${flag_new}")
    endif(${flag_var} MATCHES "${flag_old}")
    set(${flag_var} "${${flag_var}}" PARENT_SCOPE)
  endforeach()
endfunction()

#function(modify_compiler_flag_release flag_old flag_new flag_def)
#  if(${flag_new})
#    replace_compiler_flag(${flag_old} ${flag_new})
#  elseif(${flag_def})
#    replace_compiler_flag(${flag_old} ${flag_def})
#  else()
#    MESSAGE(FATAL_ERROR "You did not provide a proper replacement for ${flag_old}.")
#  endif()
#endfunction()
