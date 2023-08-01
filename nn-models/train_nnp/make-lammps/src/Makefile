TARGET = a.out

OBJECTS = allocarray.o \
          read_input_nnp3.o read_layer_node3.o read_layer_node_frex3.o \
          param_read_g2.o param_read_g5.o \
          main.o

MOD_FILES = allocarray.mod #parameters.mod

FC = ifort -O1 -shared-intel -mcmodel=large -heap-arrays -traceback
#FC = ifort -O1 -qopenmp -shared-intel -mcmodel=large
#Debag
#FC = ifort -O0 -qopenmp -check all -warn all -std -gen_interfaces -fpe0 -ftrapuv -traceback -shared-intel -mcmodel=large

#FFLAGS = -I/home/shimizu/intel/compilers_and_libraries_2018.1.163/linux/compiler/include
FFLAGS = #-I/home/shimizu/intel/compilers_and_libraries_2018.1.163/linux/compiler/lib/intel64_lin
LDFLAGS = #-openmp

#MKLROOT  = /opt/intel/composer_xe_2011/mkl
#FFLAGS  += -I${MKLROOT}/include/ia32 -I${MKLROOT}/include
#LDFLAGS += -L${MKLROOT}/lib/ia32 ${MKLROOT}/lib/ia32/libmklblas95.a

.SUFFIXES: .o .f90
.f90.o:
	${FC} -c $<

${TARGET}: ${OBJECTS}
	${FC} -o $@ ${OBJECTS} ${LDFLAGS} ${FFLAGS}

	cp ${TARGET} ../

.PHONY: clean
clean:
	${RM} ${TARGET} ${OBJECTS} ${MOD_FILES}
