# slice.exe is added on 2017.12.20
# Coded on August 15, 2015 by T. MINAMI

FC=ifort
MPIFC=mpif90
#OPT=-i-dynamic
VPATH = ../common ../solver ../src_mesh ../src_inv
#OPT0= -O3 -xAVX -mkl=parallel
OPT0= -O -fpe0 -check bounds -traceback -g -qmkl=parallel

#2016.03.03 solver with many modules

target1=calprofile.exe
target2=slice.exe
target3=modmodel.exe
target4=n_model2cond.exe
target5=change_model.exe
target6=change_model2cond.exe
target7=model2vtk.exe

all: $(target1) $(target2) $(target3) $(target4) $(target5) $(target6) $(target7)

objs1= m_constants.o  m_outerinnerproduct.o m_matrix.o m_fem_util.o m_mesh_type.o  sort_index.o m_param.o m_modelpart.o m_post_param.o calprofile.o

objs2= m_constants.o  m_outerinnerproduct.o m_matrix.o m_fem_util.o m_mesh_type.o  sort_index.o m_param.o m_post_param.o m_slicesub.o slice.o

objs3= m_param.o m_fem_util.o m_outerinnerproduct.o sort_index.o m_cond_type.o m_mesh_type.o m_modelpart.o m_param_modmodel.o modmodel.o

objs4= m_param.o m_fem_util.o m_outerinnerproduct.o sort_index.o m_cond_type.o m_mesh_type.o m_modelpart.o n_model2cond.o

objs5=m_param.o m_fem_util.o m_outerinnerproduct.o sort_index.o m_cond_type.o m_mesh_type.o m_modelpart.o change_model.o

objs6=m_param.o m_fem_util.o m_outerinnerproduct.o sort_index.o m_cond_type.o m_mesh_type.o m_modelpart.o m_triangle.o m_cond_model_type.o change_model2cond.o

objs7= m_param.o m_fem_util.o m_outerinnerproduct.o sort_index.o m_cond_type.o m_mesh_type.o m_modelpart.o model2vtk.o

$(target1) : $(objs1)
	$(FC) $(OPT0) $(objs1) -o $(target1)

$(target2) : $(objs2)
	$(FC) $(OPT0) $(objs2) -o $(target2)

$(target3) : $(objs3)
	$(FC) $(OPT0) $(objs3) -o $(target3)

$(target4) : $(objs4)
	$(FC) $(OPT0) $(objs4) -o $(target4)

$(target5) : $(objs5)
	$(FC) $(OPT0) $(objs5) -o $(target5)

$(target6) : $(objs6)
	$(FC) $(OPT0) $(objs6) -o $(target6)

$(target7) : $(objs7)
	$(FC) $(OPT0) $(objs7) -o $(target7)


# rules
.SUFFIXES : .f90 .o
.f90.o:
	$(FC) -c $(OPT0) $< -o $@
.f.o:
	$(FC) -c $(OPT0) $< -o $@

clean:
	rm -f *.o    ../src_mesh/*.o   ../solver/*.o ../src_inv/*.o
	rm -f *.mod  ../src_mesh/*.mod  ../solver/*.mod ../src_inv/*.mod
	rm -f *.exe  ../src_mesh/*.exe ../solver/*.exe ../src_inv/*.exe
