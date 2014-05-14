# ============================================================================
# Name        : Makefile
# Author      : Rod Persky (rodney.persky@connect.qut.edu.au) 
#               Jorrit Wronski (jowr@mek.dtu.dk)        
# Version     : 0.1
# Copyright   : Use and modify at your own risk.
# Description : Makefile for CoolProp shared library
# ============================================================================
# The installation procedure should be as follows:
# 1) make
# 2) sudo make install
 
# If you want a static library then add LIBTYPE=a to *ALL* of the
#  make and install commands. Also whilst you're adding commands may as well
#  include -j $(nproc) to increase build speed. So for a static library use:
#  make LIBTYPE=a -j $(nproc)
#  sudo make LIBTYPE=a install
#
# If you need a clean install use:
#  sudo make uninstall install
#
# You can simply uninstall by using:
#  sudo make uninstall
#
# When using the library in Eclipse in your source files add:
#  #include <CoolProp/CoolProp.h>
#
# And in the linker Librarys add (dl is only needed for the static library):
#  CoolProp
#  dl
#
# Header files are going to be installed into 
#  /usr/local/include
# while the libraries will reside in 
#  /usr/local/lib
# ============================================================================
# General commands:
# These commands are faily generic accross all platforms
 
CC    := g++
RM    := rm -f
CP    := cp 
CD    := cd 
CHF   := chmod 0644
CHD   := chmod 0655
MK    := mkdir -p -m 655
RS    := rsync -Rq
FN    := find
LD    := ldconfig 
LN    := ln -sf 
LT    := libtool 
AR    := ar rvs
MV    := mv
PUD   := pushd
POD   := popd
SHELL := /bin/sh
DEST  := /usr/local
UNAME := $(shell uname)
 
# What library type are we going to build by default (shared)
LIBTYPE = so
 
# We're going to install into /usr/local/include 
INCDIR      = $(DEST)/include/$(BASENAME)
LIBDIR      = $(DEST)/lib
COOLPROPDIR = ../../CoolProp/
BINDIR      = ./build
 
# Optimisation flags, change to your wants and desires
USER_OPTFLAGS = -O0 #-funroll-loops -ffast-math -ffloat-store

# Language specific flags use -D COOLPROP_LIB to build a C library
LANG_FLAGS = -D COOLPROP_LIB

# Debug flags
DBGFLAGS = -O0 -g3
WRNFLAGS = -Wall -Wpedantic -Wextra
 
# =============================
# ~ END OF CUSTOMISABLE CODE ~
# =============================
# Sort out exactly what we want to do,
# shared or static, debug or standard
# build. Get the final names of the execs
# =============================
 
BASENAME=CoolProp
VFILE = $(COOLPROPDIR)/version.h 
VERSION=$(shell cut -d'"' -f2 $(VFILE))
MAJOR=$(shell echo $(VERSION) | cut -d'.' -f1 - )
MINOR=$(shell echo $(VERSION) | cut -d'.' -f2 - )
 
ifeq ($(LIBTYPE),a)
  LIBTYPECLEAN=a
else
  ifeq ($(UNAME), Linux)
    LIBTYPECLEAN=so
  endif
  ifeq ($(UNAME), Darwin)
    LIBTYPECLEAN=dylib
  endif
endif
ifeq ($(DEBUG),debug)
  DEBUGCLEAN=_dbg
endif
 

EXECUTABLE=lib$(BASENAME)$(DEBUGCLEAN).$(LIBTYPECLEAN).$(VERSION)
LIBNAME=lib$(BASENAME)$(DEBUGCLEAN).$(LIBTYPECLEAN).$(MAJOR).$(MINOR)
EXECUTABLE_MAJOR=lib$(BASENAME)$(DEBUGCLEAN).$(LIBTYPECLEAN).$(MAJOR)
EXECUTABLE_BASE=lib$(BASENAME)$(DEBUGCLEAN).$(LIBTYPECLEAN)
 
 
# ======================
# Compile all the flags
# ======================
 
# The standard flags are meant to represent the safest compilation,
#  not fastest resulting library speed.
 
_LIBRARIES = -ldl
_INCLUDES = -I$(COOLPROPDIR)
 
 
STD_OPTFLAGS = -fPIC
_CPPFLAGS=$(_LIBRARIES) $(STD_OPTFLAGS) $(LANG_FLAGS) $(USER_OPTFLAGS)
 
ifeq ($(LIBTYPE),a)
  #AR doesn't take flags
else
  ifeq ($(UNAME), Linux)
    _LIBFLAGS=-shared -Wl,-soname,$(LIBNAME)
  endif
  ifeq ($(UNAME), Darwin)
    _LIBFLAGS=-dynamiclib -Wl,-headerpad_max_install_names,-undefined,dynamic_lookup,-compatibility_version,$(MAJOR).$(MINOR),-current_version,$(MAJOR).$(MINOR),-install_name,$(LIBNAME)
  endif
endif
 
ifeq ($(DEBUG),debug)
  _CPPFLAGS=$(STD_OPTFLAGS) $(DBGFLAGS) $(WRNFLAGS)
endif
 
#Get sources
COOLCPP_FILES=$(shell echo $(COOLPROPDIR)*.cpp)
COOLOBJ_FILES  := $(addprefix $(BINDIR)/,$(notdir $(COOLCPP_FILES:.cpp=.o)))
 
# Path for prerequisite files
vpath %.cpp $(COOLPROPDIR)

# ========================
# Compile coolprop sources
# ========================
 
.PHONY: all
all: $(COOLOBJ_FILES) $(BINDIR)/$(EXECUTABLE)

$(COOLOBJ_FILES): | $(BINDIR)

$(BINDIR):
	$(MK) -m 755 $(BINDIR)

# Make the .o files (for both static and shared)
$(BINDIR)/%.o : %.cpp
	$(CC) $(_CPPFLAGS) -c $(_INCLUDES) $< -o $@

# Link shared library
$(BINDIR)/lib$(BASENAME)$(DEBUGCLEAN).so.$(VERSION): $(COOLOBJ_FILES) 
	$(CC) $(_LIBFLAGS) $(COOLOBJ_FILES) -o $@

# Link shared library on Mac
$(BINDIR)/lib$(BASENAME)$(DEBUGCLEAN).dylib.$(VERSION): $(COOLOBJ_FILES) 
	$(CC) $(_LIBFLAGS) $(COOLOBJ_FILES) -o $@

# Link static library
$(BINDIR)/lib$(BASENAME)$(DEBUGCLEAN).a.$(VERSION): $(COOLOBJ_FILES)
	$(AR) $@ $^

# ========================
# Install and Uninstall
# ========================

.PHONY : install
install: $(BINDIR)/$(EXECUTABLE)
	$(MK)  $(INCDIR) $(LIBDIR)
	($(CD) $(COOLPROPDIR) && $(FN) . -name '*.h' -exec $(RS) {} $(INCDIR)/ \; )
	$(FN)  $(INCDIR) -type d -exec $(CHD)  {} \;
	$(FN)  $(INCDIR) -type f -exec $(CHF)  {} \;
ifeq ($(UNAME), Linux)
	$(CP)  $(BINDIR)/$(EXECUTABLE) $(LIBDIR)
endif
ifeq ($(UNAME), Darwin)
	install -m 644 -o root -g wheel $(BINDIR)/$(EXECUTABLE) $(LIBDIR)/
endif
	$(CHF) $(LIBDIR)/$(EXECUTABLE)
	($(CD) $(LIBDIR) && $(LN) $(EXECUTABLE) $(LIBNAME) )
	($(CD) $(LIBDIR) && $(LN) $(LIBNAME) $(EXECUTABLE_MAJOR) )
	($(CD) $(LIBDIR) && $(LN) $(EXECUTABLE_MAJOR) $(EXECUTABLE_BASE) )
ifeq ($(UNAME), Linux)
	$(LD) -n $(INCDIR)
	$(LD) -l $(LIBDIR)/$(EXECUTABLE)
	$(LD)
endif

.PHONY : uninstall
uninstall:
	$(RM)r $(INCDIR)
	$(RM)  $(LIBDIR)/$(EXECUTABLE)
	$(RM)  $(LIBDIR)/$(LIBNAME)
	$(RM)  $(LIBDIR)/$(EXECUTABLE_MAJOR)
	$(RM)  $(LIBDIR)/$(EXECUTABLE_BASE)
ifeq ($(UNAME), Linux)
	$(LD)
endif
	
# ========================
# Accessory functions
# ========================

.PHONY : clean
clean:
	$(RM) $(BINDIR)/*
 
.PHONY : check_flags
check_flags:
	@printf "\nBuild settings:\n\
	  VERSION=$(VERSION)\n\
	  MAJOR=$(MAJOR)\n\
	  MINOR=$(MINOR)\n\
	  EXECUTABLE=$(EXECUTABLE)\n\
	  LIBNAME=$(LIBNAME)\n\
	  LIBFLAGS=$(_LIBFLAGS)\n\
	  CPPFLAGS=$(_CPPFLAGS)\n\
	  INCLUDES=$(_INCLUDES)\n\n\
	  COOLCPP_FILES:\n$(COOLCPP_FILES)\n\n\
	  COOLOBJ_FILES:\n$(COOLOBJ_FILES)\n\n"
