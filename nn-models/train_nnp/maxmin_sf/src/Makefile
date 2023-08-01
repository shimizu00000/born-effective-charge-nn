TARGET = a.out

OBJECTS = read_input_nnp3.o \
          read_io_nnp3.o \
          read_io3.o \
          main_mean_disp.o 

MOD_FILES = #parameters.mod

FC = ifort -O1 -shared-intel -mcmodel=large

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
