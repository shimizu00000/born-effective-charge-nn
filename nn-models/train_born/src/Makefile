TARGET = a.out

OBJECTS = allocarray.o pre_random.o \
          read_input_nnp3.o read_io3.o read_layer_node3.o read_layer_node_frex3.o \
          energy_force_nnp3abc.o f_ef_deriv3abc.o \
          train3abc.o valid3abc.o \
          read_data_1.o read_data_2.o read_data_3.o \
          timer.o blas.o linpack.o lbfgsb.o \
          main_training.o

MOD_FILES = allocarray.mod

FC = mpiifort -O3 -march=core-avx2 -shared-intel -mcmodel=large -heap-arrays -traceback -fp-model strict

FFLAGS = 
LDFLAGS = 

MKLROOT  = /home/local/intel/compilers_and_libraries_2020.2.254/linux/mkl
MKLPATH  = ${MKLROOT}/lib/intel64

LDFLAGS = -L/$(MKLPATH) -Wl,--start-group -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -I/$(MKLROOT)/include/intel64 -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -Wl,--end-group

.SUFFIXES: .o .f90
.f90.o:
	${FC} -c $<

${TARGET}: ${OBJECTS}
	${FC} -o $@ ${OBJECTS} ${LDFLAGS} ${FFLAGS}

	cp ${TARGET} ../

.PHONY: clean
clean:
	${RM} ${TARGET} ${OBJECTS} ${MOD_FILES}
