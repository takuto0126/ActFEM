# Coded on September 06, 2017 by T. MINAMI

#FC=ifort  # commented out 2026.11.25 CC
# use ifx if ifx is available, otherwise ifort (ifort) 2026.11.25 CC
ifneq ($(shell command -v ifx 2>/dev/null),)
FC = ifx
else ifneq ($(shell command -v ifort 2>/dev/null),)
FC = ifort
else
$(error Neither ifx nor ifort is found. Run "source /opt/intel/oneapi/setvars.sh" first)
endif
#OPT=-i-dynamic
VPATH = ./
#OPT0= -O3 -xAVX -mkl=parallel
#OPT0= -O3 -xCORE-AVX2 -mkl=parallel -lmpi
# -xAVX is used only with ifort / mpiifort 2026.11.25 CC
ifneq ($(findstring ifort,$(FC)),)
OPT0= -xAVX -fpe0 -check bounds -traceback -g -qmkl=parallel -qopenmp
else
OPT0= -fpe0 -check bounds -traceback -g -qmkl=parallel -qopenmp
endif

target1=test.exe
objs1=m_caltime.o m_matrix.o test.o

$(target1) : $(objs1)
	$(FC) $(OPT0) $(objs1) -o $(target1)


# rules
.SUFFIXES : .f90 .o
.f90.o:
	$(FC) -c $(OPT0) $< -o $@
.f.o:
	$(FC) -c $< -o $@

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.exe

