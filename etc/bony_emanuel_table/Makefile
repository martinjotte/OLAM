.SUFFIXES:
.SUFFIXES: .F90 .f90 .o
.PHONY: all clean

F_COMP = h5fc
F_OPTS = -O3

F_COMMAND = $(F_COMP) -c $(F_OPTS)

#
#     Define run archive and executable file names.
#

EXE=bony_eman_mktab

# Define objects

OBJS = bony_eman_mktab.o \
       hdf5_utils.o       \
       hdf5_f2f.o

# Rules

%.o: %.F90
	$(F_COMMAND) $(DEFS) $< -o $@

%.o: %.f90
	$(F_COMMAND) $< -o $@

%.o: %.c
	$(C_COMMAND) $< -o $@

# Define targets

all:  $(EXE)

$(EXE): $(OBJS)
	@echo ""
	$(F_COMP) $(F_OPTS) -o $(EXE) $(OBJS)
	@echo ""
	@echo Finished building === $(EXE)
	@echo ""

clean:
	@echo ""
	rm -f $(EXE) $(OBJS) *.mod
	@echo ""

# Dependencies

bony_eman_mktab.o: hdf5_utils.o
hdf5_utils.o:       hdf5_f2f.o
