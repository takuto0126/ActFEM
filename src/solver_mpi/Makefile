# Coded on 2017.05.30 by T. MINAMI
# compile in EIC with MPI and PARDISO

#FC=mpiifort  # commented out 2026.11.25 CC
# use mpiifx if ifx is available, otherwise mpiifort (ifort) 2026.11.25 CC
ifneq ($(shell command -v ifx 2>/dev/null),)
FC = mpiifx
else ifneq ($(shell command -v ifort 2>/dev/null),)
FC = mpiifort
else
$(error Neither ifx nor ifort is found. Run "source /opt/intel/oneapi/setvars.sh" first)
endif
MPIFC=mpif90
#OPT=-i-dynamic
VPATH = ../common ../solver ../common_mpi
#OPT0= -O3 -xAVX -mkl=parallel -lmpi
# -xAVX is used only with ifort / mpiifort 2026.11.25 CC
ifneq ($(findstring ifort,$(FC)),)
OPT0= -O -xAVX -fpe0 -check bounds -traceback -g -qmkl=parallel -lmpi
else
OPT0= -O -fpe0 -check bounds -traceback -g -qmkl=parallel -lmpi
endif

#2016.03.03 solver with many modules

target2=ebfem_bxyz_mpi.exe

all: $(target2)

objs2= m_constants.o  m_outerinnerproduct.o m_matrix.o m_fem_util.o m_fem_edge_util.o m_line_type.o  m_face_type.o m_param.o m_obs_type.o m_outresp.o m_mesh_type.o m_iccg_var_takuto.o solvePARDISO.o COCG13_single.o sort_index.o m_shareformpi_fwd.o forward_bxyz_mpi.o m_triangle.o n_ebfem_bxyz_mpi.o


$(target2) : $(objs2)
	$(FC) $(OPT0) $(objs2) -o $(target2)


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
