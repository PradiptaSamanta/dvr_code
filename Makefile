FC = ifort
MPIFORT = ifort
#FFLAGS = -O0 -openmp -g -warn all -check all -debug all -traceback -fpe-all=0 -no-ipo -mkl=sequential
FFLAGS = -O2 -openmp -vec-report-0 -par-report-0 -openmp-report-0 -mkl=sequential
LIBS = -larpack -lblas
PREFIX = /home/users/0001/uk006500/Programming/DVR_Standalone
LDFLAGS = -L$(PREFIX)/arpack

all: arpack dvr_diag ang_elements rad_elements

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

OBJ_DIAG = dvr_spline_mod.o dvr_diag_mod.o dvr_diag.o
OBJ_ANG = dvr_spline_mod.o angular_mod.o angular_elements.o
OBJ_RAD = dvr_spline_mod.o dvr_diag_mod.o radial_mod.o radial_elements.o

dvr_diag: $(OBJ_DIAG)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJ_DIAG) $(LIBS)

ang_elements: $(OBJ_ANG)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJ_ANG) $(LIBS)

rad_elements: $(OBJ_RAD)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJ_RAD) $(LIBS)

arpack: $(PREFIX)/arpack/libarpack.a
#	bash -c "echo 'Compiling Arpack library'"

clean:
	@rm -f *.o *.g90 *.mod *.dat *.out dvr_diag ang_elements rad_elements

fullclean:
	@rm arpack/libarpack.a arpack/*.o -f *.o *.g90 *.mod *.dat *.out dvr_diag

.PHONY: all clean 
