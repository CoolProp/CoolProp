Steps
-----

1. Download a WinXP 32-bit version virtual box image
2. Install Winxp 32-bit virtual box image, install 7zip, tortoisegit, git
3. Take the updated_vxworks63gccdist.ZIP from (http://www.ni.com/white-paper/5694/en/) and extract to c:\\gccdist so that you have folders c:\\gccdist\\docs, c:\\gccdist\\supplemental, c:\\gccdist\\WindRiver
4. Check out coolprop sources
5. Open a console and cd to CoolProp sources
6. cd to wrappers/Labview/vxWorks
7. run build.bat (It can take a really long time (several minutes) for some reason for Mixtures.cpp, be patient)
8. Upload the generated file CoolProp.out to ni-rt/system folder using the National Instruments Measurement & Automation Explorer.  Right click on the unit, then file transfer. Also see http://www.ni.com/white-paper/3365/en/

Debugging
---------
1. Go to the IP address of the cRIO in your browser (I found IE worked better than Chrome for this)
2. Create a new session
3. Go to the console (You will need to have this enabled on your cRIO)
4. Try to run the VI on the cRIO, it will try to load the module and print out the errors that occur (if any) when loading the module

Fixing Problems
---------------
Make variables static if at all possible
Do not enable Catch (it uses language features that are incompatible with vxWorks cRIO)

To be determined: What do do about _Dtest, _Nan, _Inf ??

See Also
--------
https://decibel.ni.com/content/docs/DOC-13537

Use pre-built binaries from FirstForge
--------------------------------------

Instructions from http://firstforge.wpi.edu/sf/wiki/do/viewPage/projects.c--11_toochain/wiki/BinaryInstall

Debian GNU/Linux (Wheezy/Testing)

Note: Because of minimal dependencies this may work on other Debian based distributions (e.g. Ubuntu). This is, however, experimental.

This currently uses GCC 4.8.0 and Binutils 2.22-8 (default in debian testing).

Add this line to /etc/apt/sources.list to add the repository:

```
deb http://debian.repo.frc.s3.amazonaws.com jessie main
```

Add the maintainer key for the repository

```
sudo wget -O - http://debian.repo.frc.s3.amazonaws.com/rbmj.gpg.key | sudo apt-key add -
```

Run the following commands:

```
sudo apt-get update
sudo apt-get install gcc-powerpc-wrs-vxworks
```

Building your own version of GCC for VxWorks Target
---------------------------------------------------

Instructions from http://firstforge.wpi.edu/sf/wiki/do/viewPage/projects.c--11_toochain/wiki/ManualInstall

Generic Linux/UNIX

These are the instructions to build from source manually. Note that some lines begin with $ and some begin with #. Lines that start with a $ can be run as a regular user. Lines that start with # must be run as root (on most distributions, just prefix the command with 'sudo').

1: Download all the components.

```
$ wget http://ftp.gnu.org/gnu/gcc/gcc-4.8.2/gcc-4.8.2.tar.bz2
$ wget http://ftp.gnu.org/gnu/binutils/binutils-2.23.1.tar.bz2
$ wget ftp://ftp.ni.com/pub/devzone/tut/updated_vxworks63gccdist.zip
```

2: Set WIND_BASE

```
# echo 'export WIND_BASE=/usr/powerpc-wrs-vxworks/wind_base' >> /etc/profile
$ source /etc/profile
```

3: Install the WindRiver headers and development resources.

```
$ unzip updated_vxworks63gccdist.zip
# mkdir -p /usr/powerpc-wrs-vxworks/wind_base/target
# mkdir -p /usr/powerpc-wrs-vxworks/share/ldscripts
# cp -R gccdist/WindRiver/vxworks-6.3/host /usr/powerpc-wrs-vxworks/wind_base
# cp -R gccdist/WindRiver/vxworks-6.3/target/h/. /usr/powerpc-wrs-vxworks/sys-include
# ln -fsT /usr/powerpc-wrs-vxworks/sys-include/wrn/coreip /usr/powerpc-wrs-vxworks/wind_base/target/h
# sed '/ENTRY(_start)/d' < /usr/powerpc-wrs-vxworks/wind_base/target/h/tool/gnu/ldscripts/link.OUT > /usr/powerpc-wrs-vxworks/share/ldscripts/dkm.ld
```

4: extract binutils and gcc, and the dependency libraries

```
$ tar -jxf gcc-4.8.2.tar.bz2
$ tar -jxf binutils-2.23.1.tar.bz2
$ cd gcc-4.8.0
$ ./contrib/download_prerequisites
$ cd ..
```

5: Build & install binutils

```
$ mkdir binutils-build
$ cd binutils-build
$ ../binutils-2.23.1/configure --prefix=/usr --target=powerpc-wrs-vxworks --disable-nls
$ make -j4
# make install
$ cd ..
```

6: Build & install gcc

```
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
```
