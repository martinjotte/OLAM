.SUFFIXES:
.SUFFIXES: .F90 .f90 .c .o
.PHONY: all clean

include ./include.mk

#DEFINES= -D$(CMACH)

#FORT_INCS= -I. -I$(HDF5_INCS)
F_COMMAND = $(F_COMP) -c $(F_OPTS) $(FORT_INCS)
C_COMMAND = $(C_COMP) -c $(C_OPTS) $(DEFINES)

#
#     Define run archive and executable file names.
#

EXE=grib2olam

# Define objects

OBJS = grib2olam.o  \
       ir_popen.o   \
       grib_get.o   \
       tokenize.o   \
       hdf5_utils.o \
       hdf5_f2f.o   \
       look_ssth5.o \
       relax.o

# Rules

%.o: %.F90
	$(F_COMMAND) $(DEFS) $< -o $@

%.o: %.f90
	$(F_COMMAND) $< -o $@

%.o: %.c
	$(C_COMMAND) $< -o $@

# Define targets

all:  $(EXE) wgrib

$(EXE): $(OBJS)
	@echo ""
	$(LOADER) -o $(EXE) $(OBJS)
	@echo ""
	@echo Finished building === $(EXE)
	@echo ""

wgrib: wgrib_patched.c
	@echo "Compiling patched wgrib executable"
	$(C_COMP) $(C_OPTS) $(DEFINES) -o wgrib wgrib_patched.c
	@echo ""
	@echo Finished building wgrib
	@echo ""

clean:
	@echo ""
	rm -f $(EXE) $(OBJS) *.mod wgrib wgrib_patched.o
	@echo ""

# Dependencies

grib2olam.o: grib_get.o hdf5_utils.o

hdf5_utils.o: hdf5_f2f.o
