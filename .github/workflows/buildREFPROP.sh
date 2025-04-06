git clone https://github.com/CoolProp/gpgREFPROP
cd gpgREFPROP

# --batch to prevent interactive command
# --yes to assume "yes" for questions
gpg --quiet --batch --yes --decrypt --passphrase="$REFPROP_GPG_PASSPHRASE" --output REFPROP.zip REFPROP.zip.gpg
unzip REFPROP.zip -d REFPROP 

pip3 install numpy 

git clone --recursive https://github.com/usnistgov/REFPROP-cmake
cd REFPROP-cmake
cmake -B bld -S . -DCMAKE_BUILD_TYPE=Release -DREFPROP_FORTRAN_PATH=`pwd`/../REFPROP/FORTRAN
cmake --build bld --config Release
cp -r `pwd`/../REFPROP/FLUIDS `pwd`/bld
cp -r `pwd`/../REFPROP/MIXTURES `pwd`/bld

pip install ctREFPROP
COOLPROP_REFPROP_ROOT=`pwd`/bld python ../test_REFPROP.py
