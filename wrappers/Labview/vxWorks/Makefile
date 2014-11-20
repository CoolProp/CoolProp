# These variables control the compiler and linker flags. Change them as
# appropriate.
DEBUG_MODE = 0

ADDED_CFLAGS = -D__powerpc__ -DDISABLE_CATCH
ADDED_INCLUDES = -I../../../include -I../../../externals/Eigen
ADDED_LIBS = -lstdc++ -lgcc -lsupc++

ifeq ($(DEBUG_MODE),1)
OBJ_DIR := PPC603gnu_DEBUG
CFLAGS = -g -mlongcall
else
OBJ_DIR := PPC603gnu
CFLAGS = -Os -fstrength-reduce -fno-builtin -fno-strict-aliasing -mlongcall
endif

LINKFLAGS = $(CFLAGS) -v
LIBPATH = 
SRC_DIR = ../../../src/

# List all the *compiled* object files here, under the OBJ_DIR
# directory. Make will automatically locate the source file and
# compile it.
CPP_FILES = $(shell find $(SRC_DIR) -type f -name '*.cpp')
OBJECTS := $(subst $(SRC_DIR),$(OBJ_DIR)/,$(subst .cpp,.o,$(CPP_FILES)))
DIRECTORIES_SRC := $(addprefix $(OBJDIR),$(shell find $(SRC_DIR) -type d | tr '\n' ' '))
DIRECTORIES_BIN := $(subst $(SRC_DIR),$(OBJ_DIR)/,$(DIRECTORIES_SRC))

# This is the name of the output shared library.
PROJECT_TARGETS := $(OBJ_DIR)/CoolProp.out

# If you have other VxWorks .a files to reference, list them here.
LIBS = 
LIBPATH = 

# Everything after this line should not need to be modified for
# basic compilation. However, significant changes to the build structure
# will probably involve modifying these lines.

WIND_BASE := $(subst \,/,$(WIND_BASE))

CPU = PPC603
TOOL_FAMILY = gnu
TOOL = gnu
CC_ARCH_SPEC = -mcpu=603 -mstrict-align 

IDE_INCLUDES = -I$(WIND_BASE)/target/h -I$(WIND_BASE)/target/h/wrn/coreip -I../../../CoolProp

# This basic rule compiles a .c file into a .o file. It can be adapted to
# all other source files that gcc can compile, including assembly (.s) and
# C++ (.cpp, .cc, .C, .cxx) files. To enable support for those extensions,
# copy this rule and modify its extension and compile flags for the
# required source file type.

# Adapted rule for .cpp files
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp
	mkdir -p $(subst $(SRC_DIR),$(OBJ_DIR)/,$<)
	powerpc-wrs-vxworks-g++ $(CFLAGS) $(CC_ARCH_SPEC) -ansi  -Wall  -MD -MP $(ADDED_CFLAGS) $(IDE_INCLUDES) $(ADDED_INCLUDES) -DCPU=$(CPU) -DTOOL_FAMILY=$(TOOL_FAMILY) -DTOOL=$(TOOL) -D_WRS_KERNEL   $(DEFINES) -o "$@" -c "$<"

all:
  @echo "Please execute next commands:"
all : check_objectdir $(PROJECT_TARGETS) 

$(PROJECT_TARGETS) : $(OBJECTS)
	rm -f "$@" ctdt.c;powerpc-wrs-vxworks-nm $(OBJECTS) | tclsh $(WIND_BASE)/host/resource/hutils/tcl/munch.tcl -c ppc > ctdt.c
	powerpc-wrs-vxworks-g++ $(LINKFLAGS) $(CC_ARCH_SPEC) -fdollars-in-identifiers -Wall $(ADDED_CFLAGS) $(IDE_INCLUDES) $(ADDED_INCLUDES) -DCPU=$(CPU) -DTOOL_FAMILY=$(TOOL_FAMILY) -DTOOL=$(TOOL) -D_WRS_KERNEL   $(DEFINES)  -o ctdt.o -c ctdt.c
	powerpc-wrs-vxworks-g++ -r -nostdlib -Wl,-X -T /usr/powerpc-wrs-vxworks/share/ldscripts/dkm.ld -o "$@" $(OBJECTS) $(LIBPATH) $(LIBS)  $(ADDED_LIBPATH) $(ADDED_LIBS) ctdt.o
	rm -f ctdt.c ctdt.o


check_objectdir :
	@if [ ! -d "$(OBJ_DIR)" ]; then\
		mkdir -p $(OBJ_DIR);\
	fi

clean :
	rm -f $(OBJECTS) $(PROJECT_TARGETS) $(wildcard $(OBJ_DIR)/*.unstripped)

.DUMMY: check_objectdir clean
