.SUFFIXES:
.SUFFIXES: .F90 .f90 .c .o
.PHONY: all clean

# Include definitions

include ./paths.mk
include ./include.mk
include ./sources.mk

OBJECT_FILES = $(addsuffix .o, \
	       $(addprefix objects/, $(notdir $(basename $(SOURCE_FILES)))))

# Compiler commands

INCLUDES   = -I$(OUTILS) $(PAR_INCS) $(HDF5_INCS) $(NETCDF_INCS)
C_COMMAND  = $(C_COMP) -c $(C_OPTS) $(INCLUDES)

# F_DEFINE_FLAG is needed for the IBM XLF compiler, because the 
# -D switch is used for something other than preprocessor symbols

ifndef F_DEFINE_FLAG
F_DEFINE_FLAG=-D
endif

# For ED2
HDF5_FORTRAN_DEF=$(F_DEFINE_FLAG)USE_HDF5

ifdef OLAM_MPI
MPI_DEF=$(F_DEFINE_FLAG)OLAM_MPI
OLAM_SUFFIX=-mpi
endif

ifdef IEEE_ARITHMETIC
IEEE_ARITHMETIC_DEF=$(F_DEFINE_FLAG)IEEE_ARITHMETIC
endif

ifdef OLAM_PARALLEL_HDF5
HDF5_FORTRAN_DEF+=$(F_DEFINE_FLAG)OLAM_PARALLEL_HDF5
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
	$(C_COMMAND) $< -o $@

# Define targets

all: objects/FLIST objects/Make.depend objects/Make.targets
	$(MAKE) $(EXE)

objects/Make.depend: $(SOURCE_FILES) sources.mk
	@echo ""
	@echo "Generating dependency information"
	@mkdir -p objects
	@./sfmakedepend.txt -I $(OUTILS) -f objects/Make.depend $(SOURCE_FILES)
	@sed -i -e 's:objects/mpi.mod::g' -e 's:objects/hdf5.mod::g' objects/Make.depend

objects/Make.targets: $(SOURCE_FILES) sources.mk
	@echo ""
	@echo "Generating target dependency list"
	@mkdir -p objects
	@rm -f objects/Make.targets
	@for sfile in $(SOURCE_FILES); do \
	    dfile=objects/$$(basename $$sfile); \
	    echo $$dfile: $$sfile $$'\n\tcp' -p $$sfile $$dfile; \
	done > objects/Make.targets

$(EXE): $(OBJECT_FILES)
	@echo ""
	$(LOADER) -o $(EXE) objects/*.o $(LOADER_OPTS)  \
	    $(LIBNCARG) $(PAR_LIBS) $(LIBS) $(HDF5_LIBS) $(NETCDF_LIBS)
	@echo ""
	@echo Finished building === $(EXE)
	@echo ""

objects/FLIST: FORCE
	@for sfile in $(SOURCE_FILES); do \
	    if [ ! -f $$sfile ]; then \
	         echo $$'\nError: Source file "'$$sfile$$'" not found.\n'; \
	    fi \
	done
	@touch objects/FLIST

clean:
	@echo ""
	rm -f objects/* $(EXE)
	@echo ""

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

FORCE:

# Default fortan compiler command:
objects/%.o: F_COMMAND = $(F_COMP) -c $(F_OPTS) $(INCLUDES)

# For objects that take long to compile, use different fortran flags:
objects/ed_state_vars.o:        F_COMMAND = $(F_COMP) -c $(SLOW_FFLAGS) $(INCLUDES)
objects/ed_init_full_history.o: F_COMMAND = $(F_COMP) -c $(SLOW_FFLAGS) $(INCLUDES)
objects/edio.o:                 F_COMMAND = $(F_COMP) -c $(SLOW_FFLAGS) $(INCLUDES)

# Should we rebuild all of the .F90 files if include.mk or Makefile have changed?
$(addsuffix .o, $(addprefix objects/, $(notdir $(basename \
	$(filter %.F90,$(SOURCE_FILES)))))): include.mk Makefile

sinclude objects/Make.depend
sinclude objects/Make.targets
