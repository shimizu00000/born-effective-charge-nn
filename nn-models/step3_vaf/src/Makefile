TARGET = a.out

OBJECTS = allocarray.o \
          read_input_nnp4.o \
          param_read_g2.o \
          param_read_g5.o \
          read_io4.o \
          sf_1.o \
          sf_2_g2.o sf_2_g5.o sf_2.o \
          sf_3_g2.o sf_3_g5.o sf_3.o \
          sf_4_g2.o sf_4_g5.o sf_4.o \
          main_sf.o

MOD_FILES = allocarray.mod 

FC = mpiifort -O1 -shared-intel -mcmodel=large -heap-arrays -traceback -fp-model strict

FFLAGS = 
LDFLAGS = 

.SUFFIXES: .o .f90
.f90.o:
	${FC} -c $<

${TARGET}: ${OBJECTS}
	${FC} -o $@ ${OBJECTS} ${LDFLAGS} ${FFLAGS}

	cp ${TARGET} ../

.PHONY: clean
clean:
	${RM} ${TARGET} ${OBJECTS} ${MOD_FILES}
