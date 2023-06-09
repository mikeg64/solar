FC=gfortran
FFLAGS=-O3 -Wall -Wextra -fdefault-real-16 -mcmodel=large -fPIC -rdynamic
SRC= typesmod.f90 generalmod.f90 modbuildermod.f90 fieldbuildermod.f90 createmodel.f90
SRCPARAM= typesmod.f90 generalmod.f90 modbuildermod.f90 fieldbuildermod.f90 createparametricmodel.f90
OBJ=${SRC:.f90=.o}
OBJPARAM=${SRCPARAM:.f90=.o}

%.o: %.f90
	${FC} ${FFLAGS} -o $@ -c $<

solmodbuild: ${OBJ}
	${FC} $(FFLAGS) -o $@ $(OBJ)

solparammodbuild: ${OBJPARAM}
	${FC} $(FFLAGS) -o $@ $(OBJPARAM)

clean:
	rm *.o *.mod solmodbuild solparammodbuild

