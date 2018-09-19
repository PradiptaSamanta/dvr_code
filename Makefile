FC = ifort
MPIFORT = ifort
#FFLAGS = -O0 -openmp -g -warn all -check all -debug all -traceback -fpe-all=0 -no-ipo -mkl=sequential
#FFLAGS = -O2 -openmp -vec-report-0 -par-report-0 -openmp-report-0 -mkl=sequential
FFLAGS = -O2 -qopenmp -mkl=sequential
LIBS = -larpack -lblas
#PREFIX = /home/users/0001/uk006500/Programming/DVR_Standalone
PREFIX = `pwd`
LDFLAGS = -L$(PREFIX)/arpack

all: arpack dvr_diag ang_elements rad_elements rad_check transf_ints comb_ints

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
OBJ_RADC = dvr_spline_mod.o dvr_diag_mod.o radial_mod.o radial_check.o
OBJ_TRANSINT= dvr_spline_mod.o gen_ints_mod.o transf_ints_2e_mod.o transf_ints.o
OBJ_COMBINE= dvr_spline_mod.o dvr_diag_mod.o radial_mod.o combine_mod.o combine_elements.o

dvr_diag: $(OBJ_DIAG)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJ_DIAG) $(LIBS)

ang_elements: $(OBJ_ANG)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJ_ANG) $(LIBS)

rad_elements: $(OBJ_RAD)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJ_RAD) $(LIBS)

rad_check: $(OBJ_RADC)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJ_RADC) $(LIBS)

transf_ints: $(OBJ_TRANSINT)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJ_TRANSINT) $(LIBS)

comb_ints: $(OBJ_COMBINE)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJ_COMBINE) $(LIBS)

arpack: $(PREFIX)/arpack/libarpack.a
#	bash -c "echo 'Compiling Arpack library'"

clean:
	@rm -f *.o *.g90 *.mod *.dat *.out dvr_diag rad_elements ang_elements rad_check

fullclean:
	@rm arpack/libarpack.a arpack/*.o -f *.o *.g90 *.mod *.dat *.out dvr_diag

.PHONY: all clean 
