TARGET = a.out

OBJECTS = allocarray.o \
          read_input_nnp3.o read_layer_node3.o read_layer_node_frex3.o \
          param_read_g2.o param_read_g5.o \
          main.o

MOD_FILES = allocarray.mod 

FC = ifort -O1 -shared-intel -mcmodel=large -heap-arrays -traceback

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
