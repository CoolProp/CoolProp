# ============================================================================
# Name        : Makefile
# Author      : Jorrit Wronski (jowr@mek.dtu.dk) 
# Version     : 0.1
# Copyright   : Use and modify at your own risk.
# Description : Makefile for a CoolProp solver based on ExternalMedia.
# ============================================================================
# The installation procedure should be as follows:
# 1) make header library
# 2) sudo make install
# ============================================================================
# general commands:
RM    := rm -f
CP    := cp 
CH    := chmod 0644 
MK    := mkdir -p 
LD    := ldconfig 
LN    := ln -sf 
LT    := libtool 
AR    := ar rvs 

# used for the output
MAJORVERSION     =4
MINORVERSION     =1
THENAME          =CoolPropLib
LIBRARYEXTENSION =.a
THETEST          =coolPropLibTest

###########################################################
#  Setting the directories for library, header and 
#  binary files created in this makefile. 
###########################################################
SRCDIR     =./src
# For Dymola integration
DYMDIR     =/opt/dymola
LIBINST    =$(DYMDIR)/bin/lib
HEADINST   =$(DYMDIR)/source
# For OpenModelica integration
# LIBINST    =/usr/lib/omc
# HEADINST   =/usr/include/omc
# For system-wide installation
# LIBINST    =/usr/local/lib
# HEADINST   =/usr/local/include
BINDIR     =./bin/LinuxGCC

COOLPROPDIR=../../CoolProp

LIBS       =
OPTFLAGS   =-O3 -funroll-loops -ffast-math# -ffloat-store # optimisation, remove for debugging

# List of files to compile
OBJS             = \
	$(SRCDIR)/basesolver.o \
	$(SRCDIR)/coolpropsolver.o \
	$(SRCDIR)/errorhandling.o \
	$(SRCDIR)/solvermap.o \
	$(SRCDIR)/testsolver.o


###########################################################
#  Change these lines if you are using a different C++ 
#  compiler or if you would like to use other flags. 
###########################################################
CPPC         =g++
DEBUGFLAGS   =-g -O3
CPPFLAGS     =$(OPTFLAGS)# -Wall -pedantic -fbounds-check -ansi -Wpadded -Wpacked -malign-double -mpreferred-stack-boundary=8

.PHONY          : install
LIBFILE         =$(THENAME)
LIBRARY         =lib$(THENAME)
HEADERFILE      =$(THENAME)
HEADEREXTENSION =.h


###########################################################
#  Change these lines if you have other needs regarding
#  the generated shared library file.  
###########################################################
ifeq ($(LIBRARYEXTENSION),.so)
  install         : install_dynamic
  LIBFLAGS        =-rdynamic -fPIC -shared -Wl,-soname,$(LIBRARY)$(LIBRARYEXTENSION).$(MAJORVERSION)
endif


###########################################################
#  Change these lines if you have other needs regarding
#  the generated static library file.  
###########################################################
ifeq ($(LIBRARYEXTENSION),.a)
  install         : install_static
  LIBFLAGS        =-fPIC -static
endif


###########################################################
#  Copy files to places recognised by the system.
###########################################################
.PHONY     : install_dynamic
install_dynamic : header library
	$(MK) $(HEADINST) $(LIBINST)
	$(CP) $(BINDIR)/$(HEADERFILE)$(HEADEREXTENSION) $(HEADINST)/$(HEADERFILE)$(HEADEREXTENSION)
	$(CP) $(BINDIR)/$(LIBRARY).so $(LIBINST)/$(LIBRARY).so.$(MAJORVERSION).$(MINORVERSION)
	$(CH) $(HEADINST)/$(HEADERFILE)$(HEADEREXTENSION)
	$(CH) $(LIBINST)/$(LIBRARY).so.$(MAJORVERSION).$(MINORVERSION)
	$(LD) -l $(LIBINST)/$(LIBRARY).so.$(MAJORVERSION).$(MINORVERSION)
	$(LN) $(LIBINST)/$(LIBRARY).so.$(MAJORVERSION) $(LIBINST)/$(LIBRARY).so
	$(LD) 

.PHONY     : install_static
install_static : header library
	$(MK) $(HEADINST) $(LIBINST)
	$(CP) $(BINDIR)/$(HEADERFILE)$(HEADEREXTENSION) $(HEADINST)/$(HEADERFILE)$(HEADEREXTENSION)
	$(CP) $(BINDIR)/$(LIBRARY).a  $(LIBINST)/$(LIBRARY).a 
	$(CH) $(HEADINST)/$(HEADERFILE)$(HEADEREXTENSION)
	$(CH) $(LIBINST)/$(LIBRARY).a 
	$(LD) 

.PHONY     : uninstall
uninstall  : 
	$(RM) $(HEADINST)/$(HEADERFILE)$(HEADEREXTENSION)
	$(RM) $(LIBINST)/$(LIBRARY)$(LIBRARYEXTENSION)*
	$(RM) $(LIBINST)/$(LIBRARY).a

.PHONY     : all
all        : header library


############################################################
##  Compile coolprop sources to static files
############################################################	
#COOLCPP_FILES  := $(shell find $(COOLPROPDIR)/ -type f -name '*.cpp')
COOLCPP_FILES  := $(wildcard $(COOLPROPDIR)/*.cpp)
COOLOBJ_FILES  := $(addprefix $(SRCDIR)/,$(notdir $(COOLCPP_FILES:.cpp=.o)))

# Define search path for prerequisite files
vpath %.cpp $(COOLPROPDIR):$(COOLPROPDIR)/purefluids:$(COOLPROPDIR)/pseudopurefluids:$(SRCDIR)

.PHONY         : coolprop
coolprop       : $(COOLOBJ_FILES)


##########################################################
#  Compile the C++ sources into a library file that can 
#  be used as shared or static object.
###########################################################
.PHONY     : header
header     : $(BINDIR)/$(HEADERFILE)$(HEADEREXTENSION)
$(BINDIR)/$(HEADERFILE)$(HEADEREXTENSION): $(SRCDIR)/$(HEADERFILE)$(HEADEREXTENSION)
	$(MK) $(BINDIR)
	$(CP) $(SRCDIR)/$(HEADERFILE)$(HEADEREXTENSION) $(BINDIR)

.PHONY     : library
library    : $(BINDIR)/$(LIBRARY)$(LIBRARYEXTENSION)

$(BINDIR)/$(LIBRARY).so: $(OBJS) $(COOLOBJ_FILES) $(SRCDIR)/$(LIBFILE).o
	$(MK) $(BINDIR)
	$(CPPC) $(LIBFLAGS) $(CPPFLAGS) -o $(BINDIR)/$(LIBRARY).so $(SRCDIR)/$(LIBFILE).o $(OBJS) 

$(BINDIR)/$(LIBRARY).a: $(OBJS) $(COOLOBJ_FILES) $(SRCDIR)/$(LIBFILE).o
	$(MK) $(BINDIR)
	$(AR) $(BINDIR)/$(LIBRARY).a $^


###########################################################
#  General rulesets for compilation.
###########################################################
.PHONY: clean
clean:
	$(RM) **.o **.so **.a **.mod $(BINDIR)/* $(SRCDIR)/*.o 

$(SRCDIR)/%.o : %.cpp
	$(CPPC) $(CPPFLAGS) -o $(SRCDIR)/$*.o -I$(COOLPROPDIR) -I$(COOLPROPDIR)/purefluids -I$(COOLPROPDIR)/pseudopurefluids -I$(SRCDIR) -c $<


############################################################
##  Create the documentation from annotations in the source
##  files with DOXYGEN, a configuration file is needed.
############################################################
#.PHONY     : doc
#doc        : doc/Doxyfile
#	doxygen doc/Doxyfile
#	cd doc/latex ; \
#	make all 
	
