FC=gfortran
FFLAGS=-O3 -Wall -Wextra -fdefault-real-16
SRC= modbuildermod.f90 generalmod.f90 createmodel.f90
OBJ=${SRC:.f90=.o}

%.o: %.f90
	${FC} ${FFLAGS} -o $@ -c $<

solmodbuild: ${OBJ}
	${FC} $(FFLAGS) -o $@ $(OBJ)

clean:
	rm *.o *.mod solmodbuild

