#
#   Makefile for ddCOSMO
#
#RunF77 = ifort
#FFLAGS = -O3 -xHost -qopenmp
RunF77 = gfortran
#FFLAGS = -O3 -march=native -llapack -lblas
FFLAGS = -O3 -march=native ${MKLROOT}/lib/libmkl_intel_lp64.a \
	 ${MKLROOT}/lib/libmkl_sequential.a ${MKLROOT}/lib/libmkl_core.a \
	 -lpthread -lm -ldl -ftree-vectorize
#RunF77 = pgfortran
#FFLAGS = -O3 -mp

MODS   = ddcosmo.o ddpcm_lib.o pcm_fmm.o
OBJS   = mkrhs.o llgnew.o main.o ddcosmo.o ddpcm_lib.o forces_dd.o efld.o\
	matvec.o cosmo.o jacobi_diis.o
#
all:    $(MODS) $(OBJS)
	$(RunF77) $(FFLAGS) -o main.exe $(OBJS)
#
%.o: %.f
	$(RunF77) $(FFLAGS) -c $*.f
%.o: %.f90
	$(RunF77) $(FFLAGS) -c $*.f90
#
clean:
	rm -fr $(OBJS) *.exe *.mod *.so pcm_fmm.o

mydx:	all mydx.f90
	f2py -m mydx -c mydx.f90 $(OBJS)

pcm_fmm:	pcm_fmm.f90 llgnew.f
	f2py -m pcm_fmm -c pcm_fmm.f90 llgnew.f --opt="${FFLAGS}"

.PHONY:	all clean mydx pcm_fmm

