#######################################
#         PROJECT INFORMATION         #
#-------------------------------------#
# Version information and project     #
# metadata for CoolProp               #
#######################################

# Note: COOLPROP_VERSION_MAJOR/MINOR/PATCH/REVISION are defined in main CMakeLists.txt
# for backwards compatibility with dev/generate_headers.py which parses CMakeLists.txt

# Compose version string
set(COOLPROP_VERSION
    "${COOLPROP_VERSION_MAJOR}.${COOLPROP_VERSION_MINOR}.${COOLPROP_VERSION_PATCH}${COOLPROP_VERSION_REVISION}"
)
message(STATUS "CoolProp version: ${COOLPROP_VERSION}")

# Project metadata
string(TIMESTAMP COOLPROP_YEAR 2010-%Y)
set(COOLPROP_PUBLISHER "The CoolProp developers")
