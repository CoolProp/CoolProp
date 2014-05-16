#!/bin/bash
#
set -x 
# Specify the files
SOURCES=( "CoolPropTools" "MatrixMath" "Solvers" "PolyMath" )
ALLSRCS=""
#
CATCH="-I externals/Catch/include -D ENABLE_CATCH -D CATCH_CONFIG_MAIN "
#CATCH=""
#
LENGTH=${#SOURCES[@]}
#
# Loop through the sources and compile them
for (( i=0; i<${LENGTH}; i++ )); do
  g++ $CATCH-I include -c -o ${SOURCES[$i]}.o src/${SOURCES[$i]}.cpp
  ALLSRCS="$ALLSRCS ${SOURCES[$i]}.o"
done
g++ $ALLSRCS -o test
exit 0
