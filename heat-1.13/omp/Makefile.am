noinst_PROGRAMS = heat fibo1 fibo2 fibo3

heat_SOURCES = $(COMMON)/main.c compute.c

fibo1_SOURCES = fibo1.c
fibo1_LDADD = 

fibo2_SOURCES = fibo2.c
fibo2_LDADD = 

fibo3_SOURCES = fibo3.c
fibo3_LDADD = 

## Uncomment the following:
# noinst_PROGRAMS += matrix
# matrix_SOURCES = ...

#
# Maintenance rules
#
CONFIGURE_DEPENDENCIES = $(srcdir)/../VERSION
ACLOCAL_AMFLAGS = --force -I .
include $(srcdir)/../src/libcommon.mk


