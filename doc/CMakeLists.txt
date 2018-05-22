cmake_minimum_required(VERSION 3.0.0)


#######################################
#         PROJECT INFORMATION         #
#-------------------------------------#
# This CMakeLists.txt file is part of #
# the CoolProp thermodynamic library. #
#######################################

project("CoolProp Doxygen" VERSION ${COOLPROP_VERSION})

#add_custom_target(Name [ALL] [command1 [args1...]]
#                  [COMMAND command2 [args2...] ...]
#                  [DEPENDS depend depend depend ... ]
#                  [BYPRODUCTS [files...]]
#                  [WORKING_DIRECTORY dir]
#                  [COMMENT comment]
#                  [VERBATIM] [USES_TERMINAL]
#                  [COMMAND_EXPAND_LISTS]
#                  [SOURCES src1 [src2...]])
#add_custom_command(TARGET doxygen
#                   PRE_BUILD | PRE_LINK | POST_BUILD
#                   COMMAND command1 [ARGS] [args1...]
#                   [COMMAND command2 [ARGS] [args2...] ...]
#                   [BYPRODUCTS [files...]]
#                   [WORKING_DIRECTORY dir]
#                   [COMMENT comment]
#                   [VERBATIM] [USES_TERMINAL])

add_custom_target(docs-doxygen doxygen --version
                  COMMAND doxygen Doxyfile
                  [DEPENDS depend depend depend ... ]
                  [BYPRODUCTS [files...]]
                  [WORKING_DIRECTORY dir]
                  [COMMENT comment]
                  [VERBATIM] [USES_TERMINAL]
                  [COMMAND_EXPAND_LISTS]
                  [SOURCES src1 [src2...]])
                  
                  factory.addStep(ShellCommand(command=' '.join([docActivateCmd(), "doxygen", "--version", "&&", "doxygen", "Doxyfile"]), workdir= 'build', haltOnFailure = True))







#include(DoxygenTargets)
#add_doxygen(
#    "${COOLPROP_SOURCES_ROOT}/Doxyfile"
#    "${COOLPROP_SOURCES_ROOT}/Web/_static/doxygen"
#    "${COOLPROP_SOURCES_ROOT}/src;${COOLPROP_SOURCES_ROOT}/include"
#    ""
#    ""
#    "${COOLPROP_SOURCES_ROOT}/Web/_static/doxygen"
#    [INSTALL_COMPONENT <installcomponent>]
#    [INSTALL_PDF_NAME <installpdfname>] ]
#    [DOC_TARGET <targetname>]
#    [PROJECT_NUMBER <versionnumber>]
#    [NO_WARNINGS]
#    [NO_PDF])