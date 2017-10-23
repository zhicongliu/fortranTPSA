# Find all source files, create a list of corresponding object files
SRCS=$(wildcard *.f90) $(notdir $(wildcard src/*.f90))
OBJS=$(patsubst %.f90,%.o,$(SRCS))

# Ditto for mods (They will be in both lists)
MODS=$(wildcard mod*.f90) $(notdir $(wildcard src/mod*.f90))
MOD_OBJS=$(patsubst %.f90,%.o,$(MODS))

VPATH  = %.90 src

# Compiler/Linker settings
FC = gfortran
FLFLAGS = -std=f2003 -g 
FCFLAGS = -std=f2003 -g -c -Wall -Wextra -Wconversion -Og -pedantic -fcheck=bounds -fmax-errors=5
PROGRAM = run
PRG_OBJ = $(PROGRAM).o

# Clean the suffixes
.SUFFIXES:

# Set the suffixes we are interested in
.SUFFIXES: .f90 .o

# make without parameters will make first target found.
default : $(PROGRAM)

# Compiler steps for all objects
$(OBJS) : %.o : %.f90
	$(FC) $(FCFLAGS) -o $@ $<

# Linker
$(PROGRAM) : $(OBJS)
	$(FC) $(FLFLAGS) -o $@ $^
	mv *.mod* *.o obj/
	@echo "make done"
debug:
	@echo "SRCS = $(SRCS)"
	@echo "OBJS = $(OBJS)"
	@echo "MODS = $(MODS)"
	@echo "MOD_OBJS = $(MOD_OBJS)"
	@echo "PROGRAM = $(PROGRAM)"
	@echo "PRG_OBJ = $(PRG_OBJ)"

clean:
	rm -rf $(OBJS) $(PROGRAM) $(patsubst %.o,%.mod,$(MOD_OBJS))
	rm -rf obj/*.o obj/*.mod
	rm -rf *.mod

.PHONY: debug default clean

# Dependencies

# Main program depends on all modules
$(PRG_OBJ) : $(MOD_OBJS)

# Blocks and allocations depends on shared
mod_polymap.o: mod_mathfunc.o
quadrupole.o: mod_tpsa.o
example.o: mod_tpsa.o
