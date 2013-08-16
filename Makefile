all: bin/fourier.lax bin/fourier.long bin/cn

.PHONY: clean

bin/fourier.lax: src/fourier.lax.c
	c99 -fopenmp -lm -lgsl -lgslcblas src/fourier.lax.c -o bin/fourier.lax

bin/fourier.long: src/fourier.long.double.c
	c99 -fopenmp -lm -lgsl -lgslcblas src/fourier.long.double.c -o bin/fourier.long

bin/cn: src/cn.c
	gcc -std=gnu99 -fopenmp -lm -lgsl -lgslcblas src/cn.c -o bin/cn

bin/fourier.bd: src/fourier.bd.c
	c99 -lm -lgsl -lgslcblas src/fourier.bd.c -o bin/fourier.bd

bin/boltzmann_c_solver: src/boltzmann_c_solver.c src/boltzmann_cli.c src/boltzmann_cli.h src/boltzmann.h src/boltzmann_solver.h
	gcc -std=gnu99 -lm -lgsl -lgslcblas $+ -o $@

clean: 
	rm -rf bin/fourier.lax bin/fourier.long bin/fourier.bd
