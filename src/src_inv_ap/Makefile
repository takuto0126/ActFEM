# Latest modification is on Dec 25, 2021
# Coded on August 15, 2015 by T. MINAMI

FC=mpiifort    # for poseidon
VPATH= ../common ../solver  ../src_inv/ ../common_mpi 

# for poseidon
OPT0= -fc=ifx -O3 -fpe0 -qmkl=parallel -traceback -g -lmpi

#2016.03.03 solver with many modules
target1=ebfem_inv_ap.exe
objs1= m_caltime.o m_constants.o  m_param_ana.o m_outerinnerproduct.o m_matrix.o m_fem_util.o m_fem_edge_util.o m_line_type.o  m_param.o m_obs_type.o m_outresp.o m_mesh_type.o  m_triangle.o m_face_type.o m_modelpart.o m_param_apinv.o m_iccg_var_takuto.o solvePARDISOinv.o m_solveCM_ap.o solvebeta.o sort_index.o m_jacobian_ap.o forward_apinv.o m_freq_mpi.o m_shareformpi_ap.o m_spectral.o m_modelroughness.o n_inv_ap.o

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
	rm -f ../src_mesh/*.o
	rm -f ../src_mesh/*.mod
	rm -f ../src_mesh/*.exe
	rm -f ../solver/*.o
	rm -f ../solver/*.mod
	rm -f ../solver/*.exe
