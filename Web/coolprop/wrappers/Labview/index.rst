
.. _Labview:

***************
Labview Wrapper
***************

.. contents:: :depth: 2

Pre-compiled Binaries for windows
=================================

Shared library
--------------

You will need to download the release shared library for windows from :sfdownloads:`sourceforge <shared_library/Windows>` or the development version of the shared library from :sfnightly:`the buildbot server <shared_library/Windows>`.  If you are on 32-bit windows, download the __cdecl version.  For 64-bit, there is only one calling convention.

Place the downloaded shared library next to your library

Available libraries
-------------------

:download:`CoolProp.vi` :
Basic Library to get the properties from CoolProp.dll

:download:`CoolProp.llb` :
More advanced library, allowing to compute thermodynamic diagrams, real-time calculation of properties,
measurement processing, etc.  There is :download:`additional information <CoolProp & Labview.docx>` available for the use of the llb library

User-Compiled Binaries for windows
==================================

Check out :ref:`the instructions to build your own shared library <shared_library>`.  It should be a __cdecl DLL for 32-bit windows.

User-Compiled Binaries for VxWorks (cRIO)
=========================================

Once the toolchain is configured
--------------------------------

To compile::

    cmake ../.. -DCOOLPROP_VXWORKS_LIBRARY_MODULE=ON -DCMAKE_TOOLCHAIN_FILE=../../dev/cmake/Toolchains/powerpc-vxworks-crio.cmake -DCMAKE_CXX_FLAGS="-D__powerpc__"

Use pre-built toolchain from FirstForge
---------------------------------------

Instructions from http://firstforge.wpi.edu/sf/wiki/do/viewPage/projects.c--11_toochain/wiki/BinaryInstall

Debian GNU/Linux (Wheezy/Testing)

Note: Because of minimal dependencies this may work on other Debian based distributions (e.g. Ubuntu). This is, however, experimental.

This currently uses GCC 4.8.0 and Binutils 2.22-8 (default in debian testing).

Add this line to /etc/apt/sources.list to add the repository::

    deb http://debian.repo.frc.s3.amazonaws.com jessie main

Add the maintainer key for the repository::

    sudo wget http://debian.repo.frc.s3.amazonaws.com/rbmj.gpg.key
    sudo apt-key add rbmj.gpg.key

Run the following commands::

    sudo apt-get update
    sudo apt-get install gcc-powerpc-wrs-vxworks

Set the WIND_BASE environmental variable (or add to ~/.profile)::

    export WIND_BASE=/usr/powerpc-wrs-vxworks/wind_base

Building your own toolchain for VxWorks Target
----------------------------------------------

Instructions from http://firstforge.wpi.edu/sf/wiki/do/viewPage/projects.c--11_toochain/wiki/ManualInstall

Generic Linux/UNIX

These are the instructions to build from source manually. Note that some lines begin with $ and some begin with #. Lines that start with a $ can be run as a regular user. Lines that start with # must be run as root (on most distributions, just prefix the command with 'sudo').

1: Download all the components::

    $ wget http://ftp.gnu.org/gnu/gcc/gcc-4.8.2/gcc-4.8.2.tar.bz2
    $ wget http://ftp.gnu.org/gnu/binutils/binutils-2.23.1.tar.bz2
    $ wget ftp://ftp.ni.com/pub/devzone/tut/updated_vxworks63gccdist.zip

2: Set WIND_BASE::

    # echo 'export WIND_BASE=/usr/powerpc-wrs-vxworks/wind_base' >> /etc/profile
    $ source /etc/profile

3: Install the WindRiver headers and development resources::

    $ unzip updated_vxworks63gccdist.zip
    # mkdir -p /usr/powerpc-wrs-vxworks/wind_base/target
    # mkdir -p /usr/powerpc-wrs-vxworks/share/ldscripts
    # cp -R gccdist/WindRiver/vxworks-6.3/host /usr/powerpc-wrs-vxworks/wind_base
    # cp -R gccdist/WindRiver/vxworks-6.3/target/h/. /usr/powerpc-wrs-vxworks/sys-include
    # ln -fsT /usr/powerpc-wrs-vxworks/sys-include/wrn/coreip /usr/powerpc-wrs-vxworks/wind_base/target/h
    # sed '/ENTRY(_start)/d' < /usr/powerpc-wrs-vxworks/wind_base/target/h/tool/gnu/ldscripts/link.OUT > /usr/powerpc-wrs-vxworks/share/ldscripts/dkm.ld

4: extract binutils and gcc, and the dependency libraries::

    $ tar -jxf gcc-4.8.2.tar.bz2
    $ tar -jxf binutils-2.23.1.tar.bz2
    $ cd gcc-4.8.0
    $ ./contrib/download_prerequisites
    $ cd ..

5: Build & install binutils::

    $ mkdir binutils-build
    $ cd binutils-build
    $ ../binutils-2.23.1/configure --prefix=/usr --target=powerpc-wrs-vxworks --disable-nls
    $ make -j4
    # make install
    $ cd ..

6: Build & install gcc::

    $ mkdir gcc-build
    $ cd gcc-build
    $ ../gcc-4.8.2/configure \
          --prefix=/usr \
          --target=powerpc-wrs-vxworks \
          --with-gnu-as \
          --with-gnu-ld \
          --with-headers \
          --disable-shared \
          --disable-libssp \
          --disable-multilib \
          --with-float=hard \
          --enable-languages=c,c++ \
          --enable-libstdcxx \
          --enable-threads=vxworks \
          --without-gconv \
          --disable-libgomp \
          --disable-nls \
          --disable-libmudflap \
          --with-cpu-PPC603 \
          --disable-symvers \
          CFLAGS_FOR_TARGET='-mstrict-align -mlongcall -g -O2' \
          CXXFLAGS_FOR_TARGET='-mstrict-align -mlongcall -g -O2'

    $ make -j4
    # make install
