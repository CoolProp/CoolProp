Installation
============

Requirements
------------
* CMake
* git
* Octave (and development headers)

For ubuntu and friends, you can install Octave using
```
# Install Octave
sudo apt-get install octave liboctave-dev
```

Once octave is installed, you can run the builder and tests using
```
# Check out the sources for CoolProp
git clone https://github.com/CoolProp/CoolProp
# Move into the folder you just created
cd CoolProp
# Make a build folder
mkdir -p build/Octave
# Move into that folder
cd build/Octave
# Build the makefile using CMake
cmake ../.. -DCOOLPROP_OCTAVE_MODULE=ON -DBUILD_TESTING=ON
# Make the C# files (by default files will be generated in folder install_root/Octave relative to CMakeLists.txt file)
make install
# Run the integration tests
ctest --extra-verbose
```
