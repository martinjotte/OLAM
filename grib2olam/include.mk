F_COMP=h5fc
F_OPTS=-O3 -xAVX2

# C_COMP=gcc
C_COMP=icc
C_OPTS=-O3 -xAVX2 -Wno-unknown-pragmas

WGRIB_LIBS=
WGRIB_INCS=

JASPER_INCS=-I/usr/local/include
JASPER_LIBS=-L/usr/local/lib -ljasper -ljpeg

PNG_INCS=-I/usr/local/include
PNG_LIBS=-L/usr/local/lib -lpng -lz

LOADER_OPTS=

# For Apple OSX: the stack size may need to be increased at link time
#LOADER_OPTS=-Wl,-stack_size -Wl,0x10000000
