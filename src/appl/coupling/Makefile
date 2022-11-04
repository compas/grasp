EXE=${GRASP}/bin/coupling
LIBS=-L ${GRASP}/lib/
FC_MODULES=

OBJS= \
	sqlsf.o \
	Coupling_constants.o \
	Coupling_structures.o \
	jj2lsj_data_3_C.o \
	jj2lsj_data_2_C.o \
	jj2lsj_data_1_C.o \
	Coupling_inside_shell.o \
	Coupling_data.o \
	Coupling_evaluation.o \
	Coupling_transform_cfg_LSjj2.o \
	Coupling_transform_LSjj2.o \
	Coupling_transform_cfg_LSjj1.o \
	Coupling_transform_LSjj1.o \
	Coupling_transform_cfg_LSLScjj.o \
	Coupling_transform_LSLScjj.o \
	Coupling_transform_cfg_LSjj3.o \
	Coupling_transform_LSjj3.o \
	Coupling_transform_cfg_LScLSJ3.o \
	Coupling_transform_LScLSJ3.o \
	Coupling_transform_cfg_LSJK3.o \
	Coupling_transform_LSJK3.o \
	Coupling_transform_cfg_LSLK3.o \
	Coupling_transform_LSLK3.o \
	Coupling_transform_cfg_LSLSJ3.o \
	Coupling_transform_LSLSJ3.o \
	Coupling_transform_cfg_LSLS3.o \
	Coupling_transform_LSLS3.o \
	Coupling_transform_cfg_JJJK.o \
	Coupling_transform_JJJK.o \
	Coupling_transform_cfg_LSLK.o \
	Coupling_transform_LSLK.o \
	Coupling_transform_cfg_LSJJ.o \
	Coupling_transform_LSJJ.o \
	Coupling_getdata_mchf.o \
	Coupling_main.o \
	Coupling.o

$(EXE): $(OBJS)
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f
	$(FC) -c $(FC_FLAGS) -o $@ $<

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
