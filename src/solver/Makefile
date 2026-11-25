# Coded on August 15, 2015 by T. MINAMI

#FC=ifx  # commented out 2026.11.25 CC
# use ifx if ifx is available, otherwise ifort (ifort) 2026.11.25 CC
ifneq ($(shell command -v ifx 2>/dev/null),)
FC = ifx
else ifneq ($(shell command -v ifort 2>/dev/null),)
FC = ifort
else
$(error Neither ifx nor ifort is found. Run "source /opt/intel/oneapi/setvars.sh" first)
endif
#OPT=-i-dynamic
VPATH = ../common
# -xAVX is used only with ifort / mpiifort 2026.11.25 CC
ifneq ($(findstring ifort,$(FC)),)
OPT0= -O3 -xAVX -qmkl=parallel
else
OPT0= -O3 -qmkl=parallel
endif
#OPT0= -O3 -xAVX -qmkl=parallel # -xAVX is for Intel Xeon 
#OPT0= -O -fpe0 -check bounds -traceback -g  -mkl=parallel

target1=ebfem_bxyz.exe

all: $(target1)

objs1= m_caltime.o m_constants.o  m_param_ana.o dbesj1.o dbesj0.o intde1.o anahalfspace.o m_outerinnerproduct.o m_matrix.o m_fem_util.o m_fem_edge_util.o m_line_type.o  m_param.o m_obs_type.o m_outresp.o m_mesh_type.o  m_iccg_var_takuto.o solvePARDISO.o COCG13_single.o sort_index.o forward_bxyz.o m_triangle.o n_ebfem_bxyz.o
# bprimary2.o eprimary2.o
$(target1) : $(objs1)
	$(FC) $(OPT0) $(objs1) -o $(target1)


# rules
.SUFFIXES : .f90 .o
.f90.o:
	$(FC) -c $(OPT0) $< -o $@
.f.o:
	$(FC) -c $(OPT0) $< -o $@

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.exe
