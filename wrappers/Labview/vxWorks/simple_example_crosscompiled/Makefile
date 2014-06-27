# Name of the final output library (written to OUTFILE.out)
OUTFILE = plus_one

# hostname or IP address of the target, if automatically transferring binary via FTP.
# To automatically transfer via FTP, uncomment the FTP section below.
TARGET = 0.0.0.0

# directory to write the output library on the target, if automatically transferring via FTP.
TARGET_DIR = /ni-rt/system

# These variables control the compiler and linker flags. Change them as
# appropriate.

DEBUG_MODE = 0

ADDED_CFLAGS =

ifeq ($(DEBUG_MODE),1)
	OBJ_DIR := PPC603gnu_DEBUG
	CFLAGS = -g -mlongcall
else
	OBJ_DIR := PPC603gnu
	CFLAGS = -Os -fstrength-reduce -fno-builtin -fno-strict-aliasing -mlongcall
endif

# Expects NAME.c or NAME.cpp, and produces OBJDIR/NAME.o
SRCS = $(wildcard *.cpp)
OBJS = $(patsubst %.c,$(OBJ_DIR)/%.o,$(SRCS))

LINKFLAGS = $(CFLAGS)

# This is the name of the output shared library.
PROJECT_TARGETS := $(OBJ_DIR)/$(OUTFILE).out

# If you have other VxWorks .a files to reference, list them here.
LIBS = 
LIBPATH = 

CPP = powerpc-wrs-vxworks-g++
CPU = PPC603
TOOL_FAMILY = gnu
TOOL = gnu
CC_ARCH_SPEC = -mcpu=603 -mstrict-align -mno-implicit-fp

# This basic rule compiles a .c file into a .o file. It can be adapted to
# all other source files that gcc can compile, including assembly (.s) and
# C++ (.cpp, .cc, .C, .cxx) files. To enable support for those extensions,
# copy this rule and modify its extension and compile flags for the
# required source file type.

# Adapted rule for .cpp files
$(OBJ_DIR)/%.o : %.cpp
	 $(CPP) $(CFLAGS) $(CC_ARCH_SPEC) -ansi  -Wall  -MD -MP $(ADDED_CFLAGS) $(IDE_INCLUDES) $(ADDED_INCLUDES) -DCPU=$(CPU) -DTOOL_FAMILY=$(TOOL_FAMILY) -DTOOL=$(TOOL) -D_WRS_KERNEL   $(DEFINES) -o "$@" -c "$<"

all : check_objectdir $(PROJECT_TARGETS) 
	@echo Build of $(OUTFILE).out completed successfully.

# FTP the output file to the target (uncomment if desired)
#	@echo 
#	@echo Transferring $(OUTFILE).out to target $(TARGET) via FTP...
#	@echo user anonymous anonymous > $(OBJ_DIR)/ftpCommands.txt
#	@echo binary >> $(OBJ_DIR)/ftpCommands.txt
#	@echo send $(OBJ_DIR)/$(OUTFILE).out $(TARGET_DIR)/$(OUTFILE).out >> $(OBJ_DIR)/ftpCommands.txt
#	@echo quit >> $(OBJ_DIR)/ftpCommands.txt
#	@ftp -n -i -s:$(OBJ_DIR)/ftpCommands.txt $(TARGET)

$(PROJECT_TARGETS) : $(OBJS)
	$(CPP) -mlongcall -r -nostdlib -Wl,-X -o "$@" $(OBJS) $(LIBPATH) $(LIBS)  $(ADDED_LIBPATH) $(ADDED_LIBS)

check_objectdir :
	@if [ ! -d "$(OBJ_DIR)" ]; then\
		mkdir -p $(OBJ_DIR);\
	fi

clean :
	rm -f $(OBJS) $(PROJECT_TARGETS) $(wildcard $(OBJ_DIR)/*.unstripped)

.DUMMY: check_objectdir clean
