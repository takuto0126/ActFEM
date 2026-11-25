# Coded on Oct. 12, 2021 Takuto Minami

#FC=mpiifort  # commented out 2026.11.25 CC
# use mpiifx if ifx is available, otherwise mpiifort (ifort) 2026.11.25 CC
ifneq ($(shell command -v ifx 2>/dev/null),)
FC = mpiifx
else ifneq ($(shell command -v ifort 2>/dev/null),)
FC = mpiifort
else
$(error Neither ifx nor ifort is found. Run "source /opt/intel/oneapi/setvars.sh" first)
endif
#FC=mpiifx
#OPT=-i-dynamic
VPATH = ../common ../common_mpi/ ../solver  ../src_2D/ ../src_3DMT/
#OPT0= -O3 -xAVX -mkl=parallel  # commented out 2026.11.25 CC
# -xAVX is used only with ifort / mpiifort 2026.11.25 CC
ifneq ($(findstring ifort,$(FC)),)
OPT0= -O -xAVX -check bounds -traceback -g -qmkl=parallel
else
OPT0= -O -check bounds -traceback -g -qmkl=parallel
endif

#OPT0= -O3 -xCORE-AVX2 -mkl=parallel -lmpi

# solver with many modules
target1=n_ebfem_3DMT.exe
objs1= m_caltime.o m_constants.o  m_param_ana.o m_outerinnerproduct.o m_matrix.o m_fem_util.o m_fem_edge_util.o m_line_type.o  m_param.o m_obs_type.o m_outresp.o m_mesh_type.o  m_triangle.o m_face_type.o m_iccg_var_takuto.o sort_index.o m_surface_type.o solvePARDISO.o m_param_mt.o  forward_2DTM.o forward_3DMT.o m_shareformpi_mt.o n_ebfem_3DMT_mpi.o

$(target1) : $(objs1)
	$(FC) $(OPT0) $(objs1) -o $(target1)


# rules
.SUFFIXES : .f90 .o
.f90.o:
	$(FC) -c $(OPT0) $< -o $@
.f.o:
	$(FC) -c $< -o $@

clean:
	rm -f *.o    ../solver/*.o    ../src_2D/*.o   ../src_3DMT/*.o
	rm -f *.mod  ../solver/*.mod  ../src_2D/*.mod ../src_3DMT/*.mod
	rm -f *.exe  ../solver/*.exe  ../src_2D/*.exe ../src_3DMT/*.exe

