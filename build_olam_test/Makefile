.SUFFIXES:
.SUFFIXES: .F90 .f90 .F .f .c .o
.PHONY: all clean

# Include definitions

include ./paths.mk
include ./include.mk
include ./sources.mk

OBJECT_FILES = $(addsuffix .o, \
	       $(addprefix objects/, $(notdir $(basename $(SOURCE_FILES)))))

# Compiler commands

INCLUDES  = -I../$(OUTILS) -I../$(CMAQ) -I../$(MEGAN) \
            $(PAR_INCS) $(HDF5_INCS) $(NETCDF_INCS)
C_COMMAND = $(C_COMP) -c $(C_OPTS) $(INCLUDES)

# F_DEFINE_FLAG is needed for the IBM XLF compiler, because the 
# -D switch is used for something other than preprocessor symbols

ifndef F_DEFINE_FLAG
F_DEFINE_FLAG=-D
endif

ifdef OLAM_MPI
ifeq ($(OLAM_MPI),$(filter $(OLAM_MPI),YES yes))
MPI_DEF=$(F_DEFINE_FLAG)OLAM_MPI
OLAM_SUFFIX=-mpi
endif
endif

ifdef IEEE_ARITHMETIC
ifeq ($(IEEE_ARITHMETIC),$(filter $(IEEE_ARITHMETIC),YES yes))
IEEE_ARITHMETIC_DEF=$(F_DEFINE_FLAG)IEEE_ARITHMETIC
endif
endif

ifdef OLAM_PARALLEL_HDF5
ifeq ($(OLAM_PARALLEL_HDF5),$(filter $(OLAM_PARALLEL_HDF5),YES yes))
HDF5_FORTRAN_DEF=$(F_DEFINE_FLAG)OLAM_PARALLEL_HDF5
endif
endif

DEFS=$(MPI_DEF) $(IEEE_ARITHMETIC_DEF) $(HDF5_FORTRAN_DEF)

# Define executable name

BASE=olam
EXE=$(BASE)-$(OLAM_VERSION)$(OLAM_SUFFIX)

# Rules

%.o: %.F90
	(cd objects; $(F_COMMAND) $(DEFS) $(<F))

%.o: %.f90
	(cd objects; $(F_COMMAND) $(<F))

%.o: %.F
	(cd objects; $(F_COMMAND) $(FIXED_SRC_FLAGS) $(DEFS) $(<F))

%.o: %.f
	(cd objects; $(F_COMMAND) $(FIXED_SRC_FLAGS) $(<F))

%.o: %.c
	(cd objects; $(C_COMMAND) $(DEFS) $(<F))

# Define targets

all: objects/FLIST objects/Make.depend objects/Make.targets
	$(MAKE) $(EXE)

objects/Make.depend: $(SOURCE_FILES) sources.mk include.mk Makefile
	@echo ""
	@echo "Generating dependency information"
	@mkdir -p objects
	@./sfmakedepend.txt -I $(OUTILS) -f objects/Make.depend $(SOURCE_FILES)
	@sed -i.tmp -e 's:objects/mpi.mod::g' -e 's:objects/hdf5.mod::g' -e 's:objects/h5ds.mod::g' objects/Make.depend

objects/Make.targets: $(SOURCE_FILES) sources.mk include.mk Makefile
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

# Should we rebuild all of the .F90 files if include.mk or Makefile have changed?
$(addsuffix .o, $(addprefix objects/, $(notdir $(basename \
	$(filter %.F90,$(SOURCE_FILES)))))): include.mk Makefile

sinclude objects/Make.depend
sinclude objects/Make.targets
