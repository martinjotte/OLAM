.SUFFIXES:
.SUFFIXES: .F90 .f90 .c .o

# Include definitions

include ./paths.mk
include ./include.mk
include ./sources.mk

OBJECT_FILES = $(addsuffix .o, \
	       $(addprefix objects/, $(notdir $(basename $(SOURCE_FILES)))))

# Compiler commands

INCLUDES   = -I$(OUTILS) $(PAR_INCS) $(HDF5_INCS) $(NETCDF_INCS)
F_COMMAND  = $(F_COMP) -c $(F_OPTS) $(INCLUDES)
C_COMMAND  = $(C_COMP) -c $(C_OPTS) $(INCLUDES)

# F_DEFINE_FLAG is needed for the IBM XLF compiler, because the 
# -D switch is used for something other than preprocessor symbols

ifndef F_DEFINE_FLAG
F_DEFINE_FLAG=-D
endif

ifdef OLAM_MPI
MPI_DEF=$(F_DEFINE_FLAG)OLAM_MPI
OLAM_SUFFIX=-mpi
endif

ifdef IEEE_ARITHMETIC
IEEE_ARITHMETIC_DEF=$(F_DEFINE_FLAG)IEEE_ARITHMETIC
endif

ifdef OLAM_HDF5_FORTRAN
HDF5_FORTRAN_DEF=$(F_DEFINE_FLAG)OLAM_HDF5_FORTRAN
HDF5_C_DEF=$(F_DEFINE_FLAG)OLAM_HDF5_FORTRAN
endif


# Define executable name

BASE=olam
EXE=$(BASE)-$(OLAM_VERSION)$(OLAM_SUFFIX)

# Rules

%.o: %.F90
	(cd objects; $(F_COMMAND) $(MPI_DEF) $(HDF5_FORTRAN_DEF) $(IEEE_ARITHMETIC_DEF) $(<F))

%.o: %.f90
	(cd objects; $(F_COMMAND) $(<F))

%.o: %.c
	$(C_COMMAND) $(HDF5_C_DEF) $< -o $@

# Define targets

all: objects/Make.depend
	$(MAKE) $(EXE)

objects/Make.depend: $(SOURCE_FILES) sources.mk
	mkdir -p objects
	@echo "Generating dependency information"
	@./sfmakedepend.txt -I $(OUTILS) -f objects/Make.depend $(SOURCE_FILES)
	rm -f objects/Make.depend.old
	cat objects/Make.depend | sed s:objects/mpi.mod::g | sed s:objects/hdf5.mod::g > objects/Make.depend.new
	mv objects/Make.depend.new objects/Make.depend
	@echo "Copying source code to objects directory for compilation"
	@cp -p $(SOURCE_FILES) objects/

$(EXE): $(OBJECT_FILES)
	@echo ""
	$(LOADER) -o $(EXE) objects/*.o $(LOADER_OPTS)  \
	    $(LIBNCARG) $(PAR_LIBS) $(LIBS) $(HDF5_LIBS) $(NETCDF_LIBS)
	@echo ""
	@echo Finished building === $(EXE)
	@echo ""

clean:
	@echo ""
	rm -f objects/* $(EXE)
	@echo ""

objects/olam_mpi_atm.o:  include.mk
objects/olam_mpi_land.o: include.mk
objects/olam_mpi_sea.o:  include.mk

FORCE:

doc-html:
	@echo ""
	@echo "Building HTML documentation..."
	doxygen
	@echo ""
	@echo Finished HTML documentation. Point your browser to ./html/index.html
	@echo ""

clean-doc-html:
	@echo ""
	@echo "Cleaning HTML documentation..."
	rm -rf html
	@echo ""


sinclude objects/Make.depend
