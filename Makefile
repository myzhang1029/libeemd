.PHONY: all clean install uninstall

version := 1.4.1
src_files := $(wildcard src/*.c) $(wildcard src/%.h)
obj_files := $(patsubst src/%.c,obj/%.o,$(filter %.c,$(src_files)))
gsl_flags := $(shell pkg-config --libs --cflags gsl)
ifeq ($(gsl_flags),)
$(error Failed to query GSL complilation flags from pkg-config)
endif
gsl_flags += -DHAVE_INLINE
commonflags := -Wall -Wextra -std=c99 -pedantic -Wno-unknown-pragmas -Wshadow -Wpointer-arith -Wconversion
commonflags += $(CFLAGS)
commonflags += -g -DEEMD_DEBUG=0
commonflags += -DEEMD_VERSION=\"$(version)\"
commonflags += -fopenmp
PREFIX ?= /usr

SONAME = -soname
ifeq ($(shell uname -s),Darwin)
    SONAME = -install_name
endif

define uninstall_msg
If you used $(PREFIX) as the prefix when running `make install`,
you can undo the install by removing these files:
$(PREFIX)/include/eemd.h
$(PREFIX)/lib/libeemd.a
$(PREFIX)/lib/libeemd.so
$(PREFIX)/lib/libeemd.so.$(version)
endef
export uninstall_msg

all: libeemd.so.$(version) libeemd.a eemd.h

clean:
	rm -f libeemd.so libeemd.so.$(version) libeemd.a eemd.h obj/eemd.o
	rm -rf obj

install:
	install -d $(PREFIX)/include
	install -d $(PREFIX)/lib
	install -m644 eemd.h $(PREFIX)/include
	install -m644 libeemd.a $(PREFIX)/lib
	install libeemd.so.$(version) $(PREFIX)/lib
	cp -Pf libeemd.so $(PREFIX)/lib

uninstall:
	@echo "$$uninstall_msg"

obj:
	mkdir -p obj

obj/%.o: src/%.c src/%.h | obj
	gcc $(commonflags) -c $< $(gsl_flags) -o $@

libeemd.a: $(obj_files)
	$(AR) rcs $@ $^

libeemd.so.$(version): $(src_files)
	gcc $(commonflags) $(filter %.c,$^) -fPIC -shared -Wl,$(SONAME),$@ $(gsl_flags) -o $@
	ln -sf $@ libeemd.so

eemd.h: src/eemd.h
	cp $< $@
