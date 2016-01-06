F_COMP=h5fc

F_OPTS=-O3

C_COMP=cc
C_OPTS=-O3

JASPER_INCS=-I/usr/local/include
JASPER_LIBS=-L/usr/local/lib -ljasper -ljpeg

PNG_INCS=-I/usr/local/include
PNG_LIBS=-L/usr/local/lib -lpng -lz

HDF5_LIBS=
HDF5_INCS=

LOADER=$(F_COMP)
LOADER_OPTS=$(F_OPTS)

# For Apple OSX: the stack size may need to be increased at link time
#LOADER_OPTS=$(F_OPTS) -Wl,-stack_size -Wl,0x10000000
