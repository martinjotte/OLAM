
# Activate appropriate parts below, comment out others:

#-----------------  LINUX Portland Group pgf77/gcc ---------------
#F_COMP=pgf90
# If the compiler supports (and the user wants to use)
# the module IEEE_ARITHMETIC, uncomment below
#IEEE_ARITHMETIC=yes

# If using MPI libraries:
#MPI_PATH=/usr/local/mpich
#PAR_INCS=-I$(MPI_PATH)/include
#PAR_LIBS=-L$(MPI_PATH)/lib -lmpich -lpmpich
#OLAM_MPI=yes
#
# OPTIMIZED:
#F_OPTS=-Mvect=cachesize:524288 -Munroll -Mnoframe -O2 -pc 64 \
#       -Mfree  -Mbyteswapio
#
# DEBUG:
#F_OPTS=-O0 -g #  -pc 64 -Mfree  -Mbyteswapio -Mbounds
#
#C_COMP=gcc
#C_OPTS=-O0 -DUNDERSCORE -DLITTLE
#
#NCARG_DIR=/opt/ncl/lib
#LIBNCARG=-L$(NCARG_DIR) -lncarg -lncarg_gks -lncarg_c \
#          -lX11 -ldl -lpthread -lgfortran
#
#HDF5_LIBS=-L/home/pedropl/stuff/hdf5-pgi-serial/lib -lhdf5_fortran -lhdf5 -lz -lm
#HDF5_INCS=-I/home/pedropl/stuff/hdf5-pgi-serial/include
# If using HDF5 fortran interface, uncomment below
# There are chances to need "-lhdf5_fortran" on HDF5_LIBS above
#OLAM_HDF5_FORTRAN=yes
#
#NETCDF_LIBS=-lnetcdf
#NETCDF_INCS=
#
#LOADER=$(F_COMP)
#LOADER_OPTS=$(F_OPTS)
#LIBS=


#-----------------  LINUX Intel Fortran ifort/gcc ---------------
F_COMP=h5pfc

# If using MPI libraries:
#MPI_PATH=/usr/local/mpich
#PAR_INCS=-I$(MPI_PATH)/include
#PAR_LIBS=-L$(MPI_PATH)/lib -lmpich -lpmpich
OLAM_MPI=yes
IEEE_ARITHMETIC=yes
OLAM_PARALLEL_HDF5=yes
#USE_ED2=yes

# OPTIMIZED:
F_OPTS=-O3 -xHost -fno-alias -ip -traceback
#F_OPTS=-O3 -fno-alias -ip -traceback -g -debug extended # -override-limits -guide-data-trans -guide-vec
#F_OPTS=-O2 -xHost -fno-alias -traceback -check all -g -debug extended

# DEBUG:
#F_OPTS=-g -fp-model precise -check bounds -traceback \
#        -debug extended -ftrapuv -check pointers -check uninit

# FORTRAN FLAGS FOR BIG FILES WHICH WOULD HAVE EXCESSIVE COMPILATION TIME
SLOW_FFLAGS=-O1 -g -no-ip -traceback

C_COMP=gcc
C_OPTS=-O3

NCARG_DIR=$(NCARG_ROOT)/lib
LIBNCARG=-L$(NCARG_DIR) -lncarg -lncarg_gks -lncarg_c

HDF5_LIBS= #-L$(HDF5_ROOT)/lib -lhdf5_fortran -lhdf5 -lhdf5_hl -lz -lm
HDF5_INCS= #-I$(HDF5_ROOT)/include

NETCDF_LIBS=-L$(NETCDF_ROOT)/lib -lnetcdf -lnetcdff
NETCDF_INCS=-I$(NETCDF_ROOT)/include

LOADER=$(F_COMP) $(LIBNCARG)
#LOADER_OPTS=$(F_OPTS)

# For Apple OSX: the stack size needs to be increased at link time
LOADER_OPTS=$(F_OPTS) -Wl,-stack_size -Wl,0x40000000

## to allow ifort compiler to link with pg-compiled ncar graphics:
## LIBS=-z muldefs -L/opt/pgi/linux86-64/5.2/lib -lpgftnrtl -lpgc
#
#----------------- IBM xlf/xlc ---------------------------------
#F_COMP=xlf95_r       # without MPI
#
## F_COMP=mpxlf95_r   # with MPI
## OLAM_MPI=yes       # with MPI
#
#F_DEFINE_FLAG=-WF,-D
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
#NETCDF_LIBS=-L/usr/local/apps/netcdf-3.6.2/lib -lnetcdf
#NETCDF_INCS=-I/usr/local/apps/netcdf-3.6.2/include
#
#LOADER=$(F_COMP)
#LOADER_OPTS=$(F_OPTS)
#LIBS=
#
#-----------------------------------------------------------
