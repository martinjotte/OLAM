# Activate appropriate parts below, comment out others:

#-----------------  LINUX Portland Group pgf77/gcc ---------------
#F_COMP=pgf90
#
## If the compiler supports (and the user wants to use)
## the module IEEE_ARITHMETIC, uncomment below
#IEEE_ARITHMETIC=yes
#
## If using MPI libraries:
#OLAM_MPI=yes
#MPI_PATH=/usr/local/mpich
#PAR_INCS=-I$(MPI_PATH)/include
#PAR_LIBS=-L$(MPI_PATH)/lib -lmpich -lpmpich
#
## If parallel hdf5 is supported, uncomment the next line
#OLAM_PARALLEL_HDF5=yes
#
## OPTIMIZED:
#F_OPTS=-Mvect=cachesize:524288 -Munroll -Mnoframe -O2 -pc 64 \
#       -Mfree  -Mbyteswapio
#
## DEBUG:
#F_OPTS= -g  -pc 64 -Mfree  -Mbyteswapio -Mbounds
#
## EXTRA OPTIONS FOR FIXED-SOURCE CODE
#FIXED_SRC_FLAGS=-ffixed-form -ffixed-line-length-132

## FORTRAN FLAGS FOR BIG FILES WHICH WOULD HAVE EXCESSIVE COMPILATION TIME
#SLOW_FFLAGS=-O1 -g
#
#C_COMP=gcc
#C_OPTS=-O3 -DUNDERSCORE -DLITTLE
#
#NCARG_DIR=/usr/local/ncarg/lib
#LIBNCARG=-L$(NCARG_DIR) -lncarg -lncarg_gks -lncarg_c \
#          -L/usr/X11/lib64 -lX11 -ldl -lpthread
#
#HDF5_LIBS=-L/usr/local/lib64 -lhdf5 -lz -lm
#HDF5_INCS=-I/usr/local/include
#
#LOADER=$(F_COMP)
#LOADER_OPTS=-v -Wl,-static $(F_OPTS)
#LIBS=


#-----------------  LINUX Intel Fortran ifort/gcc ---------------
F_COMP=h5pfc

# If the compiler supports (and the user wants to use)
# the module IEEE_ARITHMETIC, uncomment below
IEEE_ARITHMETIC=yes

# If using MPI libraries:
OLAM_MPI=yes

# If parallel hdf5 is supported, uncomment the next line
OLAM_PARALLEL_HDF5=yes

# If you don't use a MPI compiler wrapper script, specify MPI
# includes and libraries here if OLAM_MPI=yes
#PAR_INCS=-I/usr/local/mpich/include
#PAR_LIBS=-L/usr/local/mpich/lib -lmpich

# OPTIMIZED:
F_OPTS=-xHost -O3 -fno-alias -ip -openmp -traceback
#F_OPTS=-g -O3 -xHost -traceback

# DEBUG:
#F_OPTS=-g -fp-model source -check bounds -traceback -inline-level=0 \
#        -debug extended -check pointers -check uninit -warn interfaces,unused

# EXTRA OPTIONS FOR FIXED-SOURCE CODE
FIXED_SRC_FLAGS=-fixed -132

# FORTRAN FLAGS FOR BIG FILES WHICH WOULD HAVE EXCESSIVE COMPILATION TIME
SLOW_FFLAGS=-O1 -g -no-ip -traceback

C_COMP=gcc
C_OPTS=-O3 -DUNDERSCORE -DLITTLE

NCARG_DIR=/data/local/lib
LIBNCARG=-L$(NCARG_DIR) -lncarg -lncarg_gks -lncarg_c \
          -L/usr/lib64 -lX11 -ldl -lpng -lpthread -lgfortran

# IF you don't use a HDF5 compiler wrapper script, specify HDF5
# includes and libraries here
#HDF5_INCS=-I/usr/local/hdf5/include
#HDF5_LIBS=-L/usr/local/hdf5/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz -lm

LOADER=$(F_COMP)
LOADER_OPTS=-static-intel $(F_OPTS)

# For Apple OSX: the stack size needs to be increased at link time
# LOADER_OPTS=-static-intel $(F_OPTS) -Wl,-stack_size -Wl,0x10000000

# to allow ifort compiler to link with pg-compiled ncar graphics:
# LIBS=-z muldefs -L/opt/pgi/linux86-64/5.2/lib -lpgftnrtl -lpgc

#----------------- IBM xlf/xlc ---------------------------------
#F_COMP=xlf95_r       # without MPI
#
## F_COMP=mpxlf95_r   # with MPI
#
## If the compiler supports (and the user wants to use)
## the module IEEE_ARITHMETIC, uncomment below
#IEEE_ARITHMETIC=yes
#
## If using MPI libraries:
#OLAM_MPI=yes
#
## If parallel hdf5 is supported, uncomment the next line
#OLAM_PARALLEL_HDF5=yes
#
## IBM compilers use a non-standard flag for pre-processing symbols,
## so uncomment the next line:
#F_DEFINE_FLAG=-WF,-D
#
#C_COMP=xlc_r
#C_OPTS=-O3 -qarch=auto -qtune=auto -qcache=auto -q64
#
## OPTIMIZED:
#F_OPTS=-O3 -qhot -q64 -qsimd=auto -qunroll=yes -qinline=auto \
#       -qsuffix=cpp=F90:f=f90 -qalias=noaryovrlp,nopteovrlp \
#       -qarch=auto -qtune=auto -qcache=auto
#
## DEBUG:
##F_OPTS=-O2 -g -C -qflttrap -qfullpath -qsuffix=cpp=F90:f=f90 \
##       -q64 -qwarn64 -qfloat=fenv -qinitauto=7fbfffff -qsigtrap \
##       -qarch=auto -qtune=auto -qcache=auto -qstrict
#
## EXTRA OPTIONS FOR FIXED-SOURCE CODE
#FIXED_SRC_FLAGS=-qfixed=132
#
## FORTRAN FLAGS FOR BIG FILES WHICH WOULD HAVE EXCESSIVE COMPILATION TIME
#SLOW_FFLAGS=-O1 -g
#
#F_DEFINE_FLAG=-WF,-D
#C_COMP=xlc_r
#C_OPTS=-O3 -qarch=auto -qtune=auto -qcache=auto -q64
#
#NCARG_DIR=/usr/local/apps/ncl-5.0.0/lib64/r4i4
#LIBNCARG=-L$(NCARG_DIR) -lncarg -lncarg_gks -lncarg_c -lX11
#
#HDF5_LIBS=-L/contrib/hdf5/lib64/r4i4 -lhdf5 -lsz -lz
#HDF5_INCS=-I/contrib/hdf5/include/hdf5-64
#
#LOADER=$(F_COMP)
#LOADER_OPTS=$(F_OPTS)
#LIBS=
#
#-----------------------------------------------------------
