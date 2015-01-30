include lualib.mk

COOLPROP_CFLAGS  ?= $(shell $(PKGCONFIG) --cflags CoolProp)
COOLPROP_LDFLAGS ?= $(shell $(PKGCONFIG) --libs-only-L CoolProp)
COOLPROP_LDLIBS  ?= $(or $(shell $(PKGCONFIG) --libs-only-l CoolProp), -lCoolProp)
COOLPROP_INCDIR  ?= $(shell $(PKGCONFIG) --variable=includedir CoolProp)
COOLPROP_LIBDIR  ?= $(shell $(PKGCONFIG) --variable=libdir CoolProp)
COOLPROP_HEADER  ?= $(or $(COOLPROP_INCDIR)

CFLAGS       ?= -g -O2 -Wall -Wextra -Wswitch-enum -Wwrite-strings -Wshadow
XCFLAGS      += -std=c99 -pedantic-errors -fPIC
XCFLAGS      += $(LUA_CFLAGS) $(COOLPROP_CFLAGS)
XLDFLAGS     += $(GUMBO_LDFLAGS) $(COOLPROP_LDLIBS)

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
