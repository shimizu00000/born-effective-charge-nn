TARGET = a.out

OBJECTS = main_data.o 

MOD_FILES = 

FC = ifort -O1 -shared-intel -mcmodel=large -heap-arrays

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
