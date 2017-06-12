.PHONY: all
gsl_flags := $(shell pkg-config --libs --cflags gsl) -DHAVE_INLINE
commonflags := -Wall -Wextra -std=c99 -pedantic -Wno-unknown-pragmas -Wshadow -Wpointer-arith -Wconversion
commonflags += $(CFLAGS) -g

all: eemd_example ceemdan_example bemd_example

%_example: %_example.c
	gcc $(commonflags) -I../src $^ -L.. -Wl,-rpath,.. $(gsl_flags) -leemd -o $@
