APP_STL := c++_static
APP_CPPFLAGS += -fexceptions
LOCAL_LDLIBS += -latomic
LIBCXX_FORCE_REBUILD := true
