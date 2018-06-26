FC = ifort
MPIFORT = ifort
#FFLAGS = -O0 -openmp -g -warn all -check all -debug all -traceback -fpe-all=0 -no-ipo -mkl=sequential
FFLAGS = -O2 -openmp -vec-report-0 -par-report-0 -openmp-report-0 -mkl=sequential
LIBS = -larpack -lblas
PREFIX = /home/users/0001/uk006500/Programming/DVR_Standalone
LDFLAGS = -L$(PREFIX)/arpack

all: arpack dvr_diag 

# How to make object files from Fortran 77 files.
%.o: %.f
	$(FC) $(FFLAGS) $(LDFLAGS) -c -o $@ $<

# How to make object files from Fortran 90 files.
%.o: %.f90
	$(FC) $(FFLAGS) $(LDFLAGS) -c -o $@ $<

# How to make object files from Fortran 90 preprocessor files.
%.o: %.F90
	$(MPIFORT) $(FFLAGS) $(LDFLAGS) -c -o $@ $<

# How to create Arpack library file 
%.a:
	bash -c "cd arpack && make clean all && cd .." 

OBJ = dvr_spline_mod.o dvr_diag_mod.o dvr_diag.o

dvr_diag: $(OBJ)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJ) $(LIBS)

arpack: $(PREFIX)/arpack/libarpack.a
#	bash -c "echo 'Compiling Arpack library'"

clean:
	@rm -f *.o *.g90 *.mod *.dat *.out dvr_diag

fullclean:
	@rm arpack/libarpack.a arpack/*.o -f *.o *.g90 *.mod *.dat *.out dvr_diag

.PHONY: all clean 
