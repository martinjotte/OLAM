#-----------------  LINUX Intel Fortran ifort/gcc ---------------

F_COMP=h5pfc

# If the compiler supports (and the user wants to use)
# the module IEEE_ARITHMETIC, uncomment below
IEEE_ARITHMETIC=yes

# If using MPI libraries:
OLAM_MPI=yes

# If parallel hdf5 is supported, uncomment the next line
OLAM_PARALLEL_HDF5=yes

# VERY OPTIMIZED:
F_OPTS=-Ofast -xHost -fno-alias -fno-fnalias -ip -traceback -assume norealloc_lhs -vecabi=cmdtarget -qopt-multi-version-aggressive -fno-pic -fminshared -ffinite-math-only -falign-loops -qopt-zmm-usage=high -qopt-prefetch=5 -qopt-subscript-in-range -qopenmp

# SAFE OPTIMIZATIONS:
#F_OPTS=-Ofast -xHost -traceback -vecabi=cmdtarget -qopenmp

# FOR REPRODUCIBLE RESULTS
#F_OPTS=-O2 -xHost -fp-model precise -fp-speculation=safe -traceback -stand f08 -assume norealloc_lhs

# EXTENDED DEBUGGING:
#F_OPTS=-O2 -xHost -g -fp-model strict -fp-speculation=strict -check bounds -traceback -warn interfaces,unused -stand f08 -debug extended -check pointers -check uninit -init=arrays -init=snan -qopenmp

# EXTRA OPTIONS FOR FIXED-SOURCE CODE
FIXED_SRC_FLAGS=-fixed -132

C_COMP=icc
C_OPTS=-O3 -xHost

LIBNCARG=-L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c -lcairo \
	-lpng -lX11 -lXext -lXrender -lfreetype

LOADER=$(F_COMP)

LOADER_OPTS=$(F_OPTS)

# For Apple OSX: the stack size needs to be increased at link time
#LOADER_OPTS=$(F_OPTS) -Wl,-stack_size -Wl,0x80000000

# to allow ifort compiler to link with pg-compiled ncar graphics:
# LIBS=-z muldefs -L/opt/pgi/linux86-64/5.2/lib -lpgftnrtl -lpgc

#-------------------------------------------------------------
