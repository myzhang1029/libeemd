AUTOMAKE_OPTIONS = foreign subdir-objects
ACLOCAL_AMFLAGS = -I m4
SUBDIRS = . examples

lib_LTLIBRARIES = libeemd.la
include_HEADERS = src/eemd.h
libeemd_la_SOURCES = src/eemd.c src/eemd.h
libeemd_la_CFLAGS = @OPENMP_CFLAGS@
libeemd_la_LDFLAGS = @OPENMP_CFLAGS@ -no-undefined
libeemd_la_LIBADD = -lgsl
