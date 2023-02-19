FC=gfortran
FFLAGS=-O3 -Wall -Wextra -fdefault-real-16
SRC=statbalancemod.f90 hydrostatsolatmos.f90 io.f90
OBJ=${SRC:.f90=.o}

%.o: %.f90
	${FC} ${FFLAGS} -o $@ -c $<

hssolatm: ${OBJ}
	${FC} $(FFLAGS) -o $@ $(OBJ)

clean:
	rm *.o *.mod hssolatm

