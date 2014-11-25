.. _cmake:

*****
CMake
*****

CMake is a very powerful cross-platform system for generating automated build
systems.  On Unix it builds makefiles which it then runs, on windows, it
normally builds projects for Visual Studio (though it can also build Makefiles
using MinGW on windows).

CMake is used to build each of the wrappers for CoolProp.

To install CMake, please visit `CMake.org <http://www.cmake.org/>`_ or install
it through your systems package manager. For Mac OSX, you can use ``homebrew``
and run ``brew install cmake``. If you do not like that, you can also install
it from the official binary release, you have to start it as
super user to be able to install the symlinks in the system path. Run::

    sudo /Applications/CMake.app/Contents/MacOS/CMake

once and choose Tools->Install For Command Line Use from the menu. Alternatively,
you can run these commands from a normal terminal session::

    sudo ln -s /Applications/CMake.app/Contents/bin/ccmake /usr/bin/
    sudo ln -s /Applications/CMake.app/Contents/bin/cmake /usr/bin/
    sudo ln -s /Applications/CMake.app/Contents/bin/cmake-gui /usr/bin/
    sudo ln -s /Applications/CMake.app/Contents/bin/cmakexbuild /usr/bin/
    sudo ln -s /Applications/CMake.app/Contents/bin/cpack /usr/bin/
    sudo ln -s /Applications/CMake.app/Contents/bin/ctest /usr/bin/

