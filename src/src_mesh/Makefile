# Coded on March 9, 2016 by T.MINAMI
# This makefile is for compiling the
# m_param.f90 3dgeo.f90
#

#FC=mpif90
#OPT=-i-dynamic
#FC=gfortran
#OPT= -O -fpe0 -check bounds -traceback -g  -mkl=parallel
FC=ifx
OPT= -O
OPT= -fpe0 -O -check bounds -traceback -g -qmkl=parallel
VPATH = ../common
all: 3dgeo mkline meshgen1 meshgen2  outmodel mkface compmodel compcond
#gridmdl2msh

# 3dgeo
objs= m_constants.o m_param.o m_horizontalresolution.o outbgmesh3d.o m_gmsh_geo.o 3dgeo.o

3dgeo: $(objs)
	$(FC) $(OPT) $(objs) -o 3dgeo.exe

# mkline
objs2= m_horizontalresolution.o m_constants.o m_param.o m_matrix.o m_outerinnerproduct.o m_fem_util.o m_fem_edge_util.o m_line_type.o sort_index.o m_mesh_type.o mkline.o
mkline: $(objs2)
	$(FC) $(OPT) $(objs2) -o mkline.exe

#-----------------------------------------------------
# meshgen1
objs3= m_constants.o m_param.o m_water_level.o m_topo_tool.o m_horizontalresolution.o outbgmesh2d.o outbgfield2d.o meshgen1.o
meshgen1: $(objs3)
	$(FC) $(OPT) $(objs3) -o meshgen1.exe

# meshgen2
objs4= m_constants.o m_param.o m_water_level.o m_topo_tool.o m_horizontalresolution.o outbgmesh3d.o m_outerinnerproduct.o sort_index.o m_fem_util.o m_mesh_type.o m_triangle.o calzobsr.o m_bell.o meshgen2.o
meshgen2: $(objs4)
	$(FC) $(OPT) $(objs4) -o meshgen2.exe


# gridmdl2msh
#objs5= m_constants.o m_param.o m_cond_type.o m_outerinnerproduct.o sort_index.o m_fem_util.o m_mesh_type.o gridmdl2msh.o
#gridmdl2msh: $(objs5)
#	$(FC) $(OPT) $(objs5) -o gridmdl2msh.exe

# outmodel
objs6= m_constants.o m_param.o m_outerinnerproduct.o sort_index.o m_fem_util.o m_mesh_type.o outmodel.o
outmodel: $(objs6)
	$(FC) $(OPT) $(objs6) -o outmodel.exe

#-----------------------------------------------------
# mkface
objs7= m_horizontalresolution.o m_constants.o m_param.o m_matrix.o m_outerinnerproduct.o m_fem_util.o m_fem_edge_util.o m_face_type.o sort_index.o m_mesh_type.o mkface.o
mkface: $(objs7)
	$(FC) $(OPT) $(objs7) -o mkface.exe

# compmodel
objs8= sort_index.o m_outerinnerproduct.o m_fem_util.o m_mesh_type.o compmodel.o
compmodel: $(objs8)
	$(FC) $(OPT) $(objs8) -o compmodel.exe

# compcond
objs9= compcond.o
compcond: $(objs9)
	$(FC) $(OPT) $(objs9) -o compcond.exe

# rules
.SUFFIXES : .f90 .o
.f90.o:
	$(FC) $(OPT) -c $< -o $@

clean:
	rm -f *.o
	rm -f *.mod
	rm -f *.exe
