FC = ifort
MPIFORT = ifort
#FFLAGS = -O0 -openmp -g -warn all -check all -debug all -traceback -fpe-all=0 -no-ipo -mkl=sequential
#FFLAGS = -O2 -openmp -vec-report-0 -par-report-0 -openmp-report-0 -mkl=sequential
FFLAGS = -O2 -qopenmp -mkl=sequential
LIBS = -larpack -lblas
#PREFIX = /home/users/0001/uk006500/Programming/DVR_Standalone
PREFIX = `pwd`
LDFLAGS = -L$(PREFIX)/arpack
VPATH = src
BPATH = bin

all: arpack dvr

# How to make object files from Fortran 77 files.
%.o: %.f
	$(FC) $(FFLAGS) $(LDFLAGS) -c -o $@ $<

# How to make object files from Fortran 90 files.
$(BPATH)/%.o: $(VPATH)/%.f90
	$(FC) $(FFLAGS) $(LDFLAGS) -c -o $@ $<

# How to make object files from Fortran 90 preprocessor files.
%.o: %.F90
	$(MPIFORT) $(FFLAGS) $(LDFLAGS) -c -o $@ $<

# How to create Arpack library file 
%.a:
	bash -c "cd arpack && make clean all && cd .." 

OBJ_DIAG = bin/*.o

dvr: $(OBJ_DIAG)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJ_DIAG) $(LIBS)

arpack: $(PREFIX)/arpack/libarpack.a
#	bash -c "echo 'Compiling Arpack library'"

clean:
	@rm -f *.o *.g90 *.mod *.dat *.out dvr_diag rad_elements ang_elements rad_check

fullclean:
	@rm arpack/libarpack.a arpack/*.o -f *.o *.g90 *.mod *.dat *.out dvr_diag

.PHONY: all clean 
