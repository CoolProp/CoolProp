#######################################
#         PROJECT INFORMATION         #
#-------------------------------------#
# Version information and project     #
# metadata for CoolProp               #
#######################################

# Project version
set(COOLPROP_VERSION_MAJOR 7)
set(COOLPROP_VERSION_MINOR 1)
set(COOLPROP_VERSION_PATCH 0)
set(COOLPROP_VERSION_REVISION )
set(COOLPROP_VERSION
    "${COOLPROP_VERSION_MAJOR}.${COOLPROP_VERSION_MINOR}.${COOLPROP_VERSION_PATCH}${COOLPROP_VERSION_REVISION}"
)
message(STATUS "CoolProp version: ${COOLPROP_VERSION}")

# Project metadata
string(TIMESTAMP COOLPROP_YEAR 2010-%Y)
set(COOLPROP_PUBLISHER "The CoolProp developers")
