FLAGS=-fopenmp
LAPACK=/usr/lib/liblapack.so

windows: LAPACK=-llapack
windows: all

all: vars.o ranlux.o fmap.o
	gfortran -O3 -o fmap.x $^ $(LAPACK) $(FLAGS)

%.o: %.f95
	gfortran -O3 -c $^ $(FLAGS)

clean:
	@rm -fv *.x
	@rm -fv *.o
	@rm -fv *.out
	@rm -fv *.dat
	@rm -fv *.mod
	@rm -fv *.log
