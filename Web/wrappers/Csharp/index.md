Installation
============

Requirements
------------
* CMake
* git
* C#

For ubuntu and friends, you will need to install Mono C# using
```
# Install Mono version of C#
sudo apt-get install mono-mcs mono-runtime
```

For windows, download the Visual Studio 2010 version of C# (other versions should be fine too)

Once mono c# is installed, you can run the 
```
# Check out the sources for CoolProp
git clone https://github.com/CoolProp/CoolProp
# Move into that folder you just created
cd CoolProp
# Make a build folder
mkdir -p build/Csharp
# Move into that folder
cd build/Csharp
# Build the makefile using CMake
cmake ../.. -DCOOLPROP_JAVA_MODULE=ON -DBUILD_TESTING=ON
# Make the C# files
make install
# Run the integration tests
ctest --extra-verbose
```
