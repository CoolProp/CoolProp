include lualib.mk

COOLPROP_CFLAGS  ?= $(shell $(PKGCONFIG) --cflags CoolProp)
COOLPROP_LDFLAGS ?= $(shell $(PKGCONFIG) --libs-only-L CoolProp)
COOLPROP_LDLIBS  ?= $(or $(shell $(PKGCONFIG) --libs-only-l CoolProp), -lCoolProp)
COOLPROP_INCDIR  ?= $(shell $(PKGCONFIG) --variable=includedir CoolProp)
COOLPROP_LIBDIR  ?= $(shell $(PKGCONFIG) --variable=libdir CoolProp)
COOLPROP_HEADER  ?= $(or $(COOLPROP_INCDIR)

CFLAGS       ?= -g -O2 -Wall -Wextra -Wswitch-enum -Wwrite-strings -Wshadow
XCFLAGS      += -std=c99 -pedantic-errors -fpic
XCFLAGS      += $(LUA_CFLAGS) $(COOLPROP_CFLAGS)
XLDFLAGS     += $(GUMBO_LDFLAGS) $(COOLPROP_LDLIBS)

all: coolprop/capi.so

coolprop/capi.o: coolprop/capi.c coolprop/compat.h

clean:
	$(RM) coolprop/capi.so coolprop/capi.o

.PHONY: all

.DELETE_ON_ERROR:
