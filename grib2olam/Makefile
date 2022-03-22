.SUFFIXES:
.SUFFIXES: .F90 .f90 .c .o
.PHONY: all clean wgrib2

include ./include.mk

#DEFINES= -D$(CMACH)

#FORT_INCS= -I. -I$(HDF5_INCS)
F_COMMAND = $(F_COMP) -c $(F_OPTS) $(FORT_INCS)
C_COMMAND = $(C_COMP) -c $(C_OPTS) $(DEFINES)

WGRIB2_FLAGS  = $(C_OPTS) $(PNG_INCS) $(JASPER_INCS) -DUSE_PNG -DUSE_JASPER
WGRIB2_LDFLGS = $(C_OPTS) $(JASPER_LIBS) $(PNG_LIBS)
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
       relax.o      \
       gaussian.o   \
       bit_shave.o

# Rules

%.o: %.F90
	$(F_COMMAND) $(DEFS) $< -o $@

%.o: %.f90
	$(F_COMMAND) $< -o $@

%.o: %.c
	$(C_COMMAND) $< -o $@

# Define targets

all:  $(EXE) wgrib wgrib2

$(EXE): $(OBJS)
	@echo ""
	$(LOADER) -o $(EXE) $(OBJS) $(LOADER_OPTS)
	@echo ""
	@echo Finished building === $(EXE)
	@echo ""

wgrib: wgrib_patched.c
	@echo "Compiling patched wgrib executable"
	$(C_COMP) $(C_OPTS) $(DEFINES) -o wgrib wgrib_patched.c
	@echo ""
	@echo Finished building wgrib
	@echo ""

wgrib2:
	@echo "Compiling wgrib2 executable"
	$(MAKE) -C wgrib2_src CC="$(C_COMP)" CFLAGS="$(WGRIB2_FLAGS)" LDFLAGS="$(WGRIB2_LDFLGS)"
	ln -sf wgrib2_src/wgrib2 ./
	@echo Finished building wgrib2
	@echo ""

clean:
	@echo ""
	rm -f $(EXE) $(OBJS) *.mod wgrib wgrib_patched.o wgrib2
	make -C wgrib2_src clean
	@echo ""

# Dependencies

grib2olam.o: grib_get.o hdf5_utils.o

hdf5_utils.o: hdf5_f2f.o

grib_get.o: gaussian.o tokenize.o
