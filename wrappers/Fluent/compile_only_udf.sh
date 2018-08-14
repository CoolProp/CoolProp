#!/bin/bash
# This Script only links the udf with the coolprop library, useful if the coolprop library is already compiled

# Warning: this script deletes directory libudf before creating a new one

# Warning: read README.rst before editing the variables SOLVER and FLUENT_BIN_FOLDER

SOLVER="$1"
FLUENT_BIN_FOLDER="$2"

if [ "$SOLVER" == "" ]
then
   #compiles 3d double precision by default
   SOLVER="3ddp" 
fi

if [ "$FLUENT_BIN_FOLDER" == "" ]
then
   #call by default qualified Fluent version for studies -->  adapt the following line
   FLUENT_BIN_FOLDER="/users5/appli/produits/ansys_inc/v150/fluent/bin" 
fi
FLUENT_COMPILER=g++

HOME=`pwd`

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
#make "FLUENT_ARCH=$FLUENT_ARCH" "SOURCES=$udfNames" "FLUENT_INC=$FLUENT_INC" "CC=$FLUENT_COMPILER" "CFLAGS_LNAMD64=-D_lnamd64 -D_GNU_SOURCE -fpic -shared -ansi -O -Wall -DPTR_RESTRICT= " # previous command to compile udf, it has been replaced with the following because gcc compiles udf with name mangling (to be sure the command nm -D libudf.so must return a function PropsSI and not _Z7PropsSIPcS_dS_dS_)
make "FLUENT_ARCH=$FLUENT_ARCH" 
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
