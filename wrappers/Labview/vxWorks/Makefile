# These variables control the compiler and linker flags. Change them as
# appropriate.
DEBUG_MODE = 0

ADDED_CFLAGS = -D__powerpc__ -DDISABLE_CATCH

ifeq ($(DEBUG_MODE),1)
OBJ_DIR := PPC603gnu_DEBUG
CFLAGS = -g -mlongcall
else
OBJ_DIR := PPC603gnu
CFLAGS = -Os -fstrength-reduce -fno-builtin -fno-strict-aliasing -mlongcall
endif

LINKFLAGS = $(CFLAGS) -v
LIBPATH = 

# List all the *compiled* object files here, under the OBJ_DIR
# directory. Make will automatically locate the source file and
# compile it.
OBJECTS := $(subst ../../../CoolProp/,$(OBJ_DIR)/,$(subst .cpp,.o,$(wildcard ../../../CoolProp/*.cpp)))

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
CC_ARCH_SPEC = -mcpu=603 -mstrict-align -mno-implicit-fp

IDE_INCLUDES = -I$(WIND_BASE)/target/h -I$(WIND_BASE)/target/h/wrn/coreip -I../../../CoolProp

# This basic rule compiles a .c file into a .o file. It can be adapted to
# all other source files that gcc can compile, including assembly (.s) and
# C++ (.cpp, .cc, .C, .cxx) files. To enable support for those extensions,
# copy this rule and modify its extension and compile flags for the
# required source file type.
$(OBJ_DIR)/%.o : ../../../CoolProp/%.c
	ccppc $(CFLAGS) $(CC_ARCH_SPEC) -ansi  -Wall  -MD -MP $(ADDED_CFLAGS) $(IDE_INCLUDES) $(ADDED_INCLUDES) -DCPU=$(CPU) -DTOOL_FAMILY=$(TOOL_FAMILY) -DTOOL=$(TOOL) -D_WRS_KERNEL   $(DEFINES) -o "$@" -c "$<"

# Adapted rule for .cpp files
$(OBJ_DIR)/%.o : ../../../CoolProp/%.cpp
	c++ppc $(CFLAGS) $(CC_ARCH_SPEC) -ansi  -Wall  -MD -MP $(ADDED_CFLAGS) $(IDE_INCLUDES) $(ADDED_INCLUDES) -DCPU=$(CPU) -DTOOL_FAMILY=$(TOOL_FAMILY) -DTOOL=$(TOOL) -D_WRS_KERNEL   $(DEFINES) -o "$@" -c "$<"

all : check_objectdir $(PROJECT_TARGETS) 

$(PROJECT_TARGETS) : $(OBJECTS)
	rm -f "$@" ctdt.c;nmppc $(OBJECTS) | tclsh $(WIND_BASE)/host/resource/hutils/tcl/munch.tcl -c ppc > ctdt.c
	ccppc $(LINKFLAGS) $(CC_ARCH_SPEC) -fdollars-in-identifiers -Wall $(ADDED_CFLAGS) $(IDE_INCLUDES) $(ADDED_INCLUDES) -DCPU=$(CPU) -DTOOL_FAMILY=$(TOOL_FAMILY) -DTOOL=$(TOOL) -D_WRS_KERNEL   $(DEFINES)  -o ctdt.o -c ctdt.c
	ccppc -r -nostdlib -Wl,-X -T $(WIND_BASE)/target/h/tool/gnu/ldscripts/link.OUT -o "$@" $(OBJECTS) $(LIBPATH) $(LIBS)  $(ADDED_LIBPATH) $(ADDED_LIBS) ctdt.o
	rm -f ctdt.c ctdt.o


check_objectdir :
	@if [ ! -d "$(OBJ_DIR)" ]; then\
		mkdir -p $(OBJ_DIR);\
	fi

clean :
	rm -f $(OBJECTS) $(PROJECT_TARGETS) $(wildcard $(OBJ_DIR)/*.unstripped)

.DUMMY: check_objectdir clean
