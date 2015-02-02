include lualib.mk

COOLPROP_CFLAGS  ?= $(shell $(PKGCONFIG) --cflags CoolProp)

CFLAGS       ?= -g -O2 -Wall -Wextra -Wswitch-enum -Wwrite-strings -Wshadow
XCFLAGS      += -std=c99 -pedantic-errors -fPIC
XCFLAGS      += $(LUA_CFLAGS) $(COOLPROP_CFLAGS)

all: coolprop/capi.so

coolprop/capi.o: coolprop/capi.c coolprop/compat.h

install: all
	$(MKDIR) '$(DESTDIR)$(LUA_CMOD_DIR)/coolprop'
	$(MKDIR) '$(DESTDIR)$(LUA_LMOD_DIR)/coolprop'
	$(INSTALLX) coolprop/capi.so '$(DESTDIR)$(LUA_CMOD_DIR)/coolprop/capi.so'
	$(INSTALL) coolprop.lua '$(DESTDIR)$(LUA_LMOD_DIR)/'
	$(INSTALL) coolprop/ffi.lua '$(DESTDIR)$(LUA_LMOD_DIR)/coolprop/'

uninstall:
	$(RM) '$(DESTDIR)$(LUA_LMOD_DIR)/coolprop.lua'
	$(RM) -r '$(DESTDIR)$(LUA_CMOD_DIR)/coolprop/'
	$(RM) -r '$(DESTDIR)$(LUA_LMOD_DIR)/coolprop/'

clean:
	$(RM) coolprop/capi.so coolprop/capi.o

.PHONY: all install uninstall clean

.DELETE_ON_ERROR:
