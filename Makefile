all: bch 

#CC = gcc
CC = clang

#-ftrapv: This option generates traps for signed overflow on addition, subtraction, multiplication operations. 
#COPTS = -O3 -fopenmp -ftrapv   -Wall -DUSE_INT128_T 

COPTS = -O3 -fopenmp -Wall -DUSE_INT128_T 

#COPTS = -g -Wall -DUSE_INT128_T 

bch: bch.h bch.c phi.c lie_series.c rightnormed.c goldberg.c
	$(CC) $(COPTS) phi.c goldberg.c rightnormed.c lie_series.c bch.c -o bch

bch_goldberg_30.txt: bch
	./bch goldberg_coefficients=1 N=30 > bch_goldberg_30.txt

bch_lyndon_20.txt: bch
	./bch N=20 M=14 verbosity_level=1 > bch_lyndon_20.txt

