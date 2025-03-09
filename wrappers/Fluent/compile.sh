#!/bin/bash
# This Script build the CoolProp sources into a shared library, then it links it to the Fluent udf

# Warning: this script deletes directory libudf before creating a new one

# Warning: read README.rst before editing the variables SOLVER and FLUENT_BIN_FOLDER

SOLVER="$1"
FLUENT_BIN_FOLDER="$2"

if [ "$SOLVER" == "" ]
then
   #compiles 3d double precision by default -->  adapt the following line
   SOLVER="3ddp" 
fi

if [ "$FLUENT_BIN_FOLDER" == "" ]
then
   #call by default qualified Fluent version for studies -->  adapt the following line
   FLUENT_BIN_FOLDER="/users5/appli/produits/ansys_inc/v150/fluent/bin" 
fi
FLUENT_COMPILER=g++

HOME=`pwd`
#CoolProp root source file (where CMakelist.txt is) -->  adapt the following line
CoolProp_sources_repository="/usersCFD/epd/ah33666/CoolProp-6.1.0"
#Python 2.7 is mandatory to use cmake, by default on our machines the Python version is 2.4 -->  adapt the following line
PYTHON_273="/users5/appli/produits/PYTHON/2.7.3/bin/python"
#Specify cmake shortcut if default cmake<2.8 -->  adapt the following line
alias cmake='/users5/appli/tools/CMAKE/3.6.0/bin/cmake'
# beginning usefull Fluent environnement variables extraction
echo exit > exit.jou
if [ "$FLUENT_BIN_FOLDER" = "" ]
then
	fluent $SOLVER -g -env -i exit.jou > FluentEnvironment.dat 2>&1
else
	$FLUENT_BIN_FOLDER"/fluent" $SOLVER -g -env -i exit.jou > FluentEnvironment.dat 2>&1
fi
rm exit.jou
FLUENT_INC=`cat FluentEnvironment.dat | grep FLUENT_INC | sed 's/FLUENT[_A-Z]*=//'`
FLUENT_ARCH=`cat FluentEnvironment.dat | grep FLUENT_ARCH | sed 's/FLUENT[_A-Z]*=//'`
FLUENT_PROD_DIR=`cat FluentEnvironment.dat | grep FLUENT_PROD_DIR | sed 's/FLUENT[_A-Z]*=//'`
PATH_MAKEFILE1=$FLUENT_PROD_DIR"/src/makefile.udf"
PATH_MAKEFILE2=$FLUENT_PROD_DIR"/src/makefile.udf2"
PATH_LIBRARY="libudf/"$FLUENT_ARCH"/"$SOLVER
# extraction ending

#building the libudf directory
cd $HOME
rm -rf libudf/
mkdir libudf
cp $PATH_MAKEFILE2 libudf/makefile
mkdir libudf/src
cp $PATH_MAKEFILE1 libudf/src/makefile	
mkdir libudf/$FLUENT_ARCH
mkdir $PATH_LIBRARY
mkdir ${PATH_LIBRARY}"_host"
mkdir ${PATH_LIBRARY}"_node"
udfNames=""
echo -n CSOURCES= > ${PATH_LIBRARY}/user.udf 
for file in *.c
do
  echo -n $file" " >> ${PATH_LIBRARY}/user.udf 
  udfNames=$udfNames" "$file
done
echo -n -e "\n"HSOURCES= >> ${PATH_LIBRARY}/user.udf 
for file in *.h
do
  if [$file != "*.h"]
  then
  	echo -n $file" " >> ${PATH_LIBRARY}/user.udf 
  fi
done
echo -e "\n"FLUENT_INC=$FLUENT_INC >> ${PATH_LIBRARY}/user.udf
echo GPU_SUPPORT=off >> ${PATH_LIBRARY}/user.udf
cp ${PATH_LIBRARY}/user.udf ${PATH_LIBRARY}_host/user.udf
cp ${PATH_LIBRARY}/user.udf ${PATH_LIBRARY}_node/user.udf
#moving and compiling the source 
cp *.c libudf/src/
cp *.h libudf/src/
cd libudf
#make "FLUENT_ARCH=lnamd64" "SOURCES=$udfNames" "FLUENT_INC=$FLUENT_INC" "CC=$FLUENT_COMPILER" "CFLAGS_LNAMD64=-D_lnamd64 -D_GNU_SOURCE -fpic -shared -ansi -O -Wall -DPTR_RESTRICT= "
make "FLUENT_ARCH=$FLUENT_ARCH" 
cd $HOME
#remove the CoolProp build directory before compilation in the build directory
rm -rf CoolProp_Build
mkdir CoolProp_Build
cd CoolProp_Build
#uncomment following line for make file generation if python version <2.7
cmake -G "Unix Makefiles" -DCOOLPROP_EXTERNC_LIBRARY:BOOL="ON" -DCOOLPROP_SHARED_LIBRARY:BOOL="ON" -DPYTHON_EXECUTABLE:FILEPATH=$PYTHON_273 $CoolProp_sources_repository
#uncomment following line make file generation if python is 2.7
#cmake -G "Unix Makefiles" -DCOOLPROP_EXTERNC_LIBRARY:BOOL="ON" -DCOOLPROP_SHARED_LIBRARY:BOOL="ON" $CoolProp_sources_repository
make
cd $HOME
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/CoolProp_Build/

cp CoolProp_Build/*.so* ${PATH_LIBRARY}/.
cd  ${PATH_LIBRARY}/
g++ -Wall -shared *.o libCoolProp.so -o libudf.so

cd $HOME
cp CoolProp_Build/*.so* ${PATH_LIBRARY}_host/.
cd ${PATH_LIBRARY}_host
g++ -Wall -shared *.o libCoolProp.so -o libudf.so

cd $HOME
cp CoolProp_Build/*.so* ${PATH_LIBRARY}_node/.
cd ${PATH_LIBRARY}_node
g++ -Wall -shared *.o libCoolProp.so -o libudf.so
