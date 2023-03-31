FC=gfortran
FFLAGS=-O3 -Wall -Wextra -fdefault-real-16
SRC=creatmodel.f90 modelbuildermod.f90 generalmodel.f90
OBJ=${SRC:.f90=.o}

%.o: %.f90
	${FC} ${FFLAGS} -o $@ -c $<

hssolatm: ${OBJ}
	${FC} $(FFLAGS) -o $@ $(OBJ)

clean:
	rm *.o *.mod solmodbuild

