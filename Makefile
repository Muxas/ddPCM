#
#   Makefile for ddCOSMO
#
#RunF77 = ifort
#FFLAGS = -O3 -xHost -qopenmp
RunF77 = gfortran
#FFLAGS = -O3 -march=native -llapack -lblas
FFLAGS = -O1 -march=native ${MKLROOT}/lib/libmkl_intel_lp64.a \
	 ${MKLROOT}/lib/libmkl_sequential.a ${MKLROOT}/lib/libmkl_core.a \
	 -lpthread -lm -ldl -ftree-vectorize -finline-functions -fcheck=all
#RunF77 = pgfortran
#FFLAGS = -O3 -mp

MODS   = ddcosmo.o pcm_fmm.o ddpcm_lib.o
OBJS   = ${MODS} mkrhs.o llgnew.o forces_dd.o efld.o\
	matvec.o cosmo.o jacobi_diis.o
#
all:    main.exe main_fmm.exe
	
main.exe: $(MODS) $(OBJS)
	$(RunF77) $(FFLAGS) main.f90 -o main.exe $(OBJS)

main_fmm.exe: $(MODS) $(OBJS)
	$(RunF77) $(FFLAGS) main_fmm.f90 -o main_fmm.exe $(OBJS)
#
%.o: %.f
	$(RunF77) $(FFLAGS) -c $*.f
%.o: %.f90
	$(RunF77) $(FFLAGS) -c $*.f90
#
clean:
	rm -fr $(OBJS) *.exe *.mod *.so
#
test_pcm_fmm: pcm_fmm.f90 llgnew.o test_pcm_fmm.f90
	f2py -m test_pcm_fmm -c pcm_fmm.f90 llgnew.o --opt="${FFLAGS}"
#
test:	test_pcm_fmm.f90 pcm_fmm.o llgnew.o
	$(RunF77) $(FFLAGS) test_pcm_fmm.f90 -o test_pcm_fmm.exe pcm_fmm.o\
	    llgnew.o
#
.PHONY:	all clean test_pcm_fmm test

