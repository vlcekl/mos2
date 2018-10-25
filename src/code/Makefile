#FC = mpif90
FC = gfortran
#FFLAGS = -fopenmp -O2 -ffast-math -w
FFLAGS = -O2 -ffast-math -w

FILES = rndmod.f90 sysdef.f90 config.f90 ener.f90 move.f90 measure.f90 lgmain.f90

all: clean

lg: $(FILES)
	$(FC) $(FFLAGS) -o $@ $(FILES)

clean: lg
	rm *.mod
