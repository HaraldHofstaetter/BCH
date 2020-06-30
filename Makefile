all: bch 

CC = gcc
#CC = clang

#-ftrapv: This option generates traps for signed overflow on addition, subtraction, multiplication operations. 
#COPTS = -O3 -fopenmp -ftrapv   -Wall -DUSE_INT128_T 

COPTS = -O3 -fopenmp -Wall -DUSE_INT128_T 

#COPTS = -g -Wall -DUSE_INT128_T 

bch: bch.h bch.c phi.c lyndon.c goldberg.c
	$(CC) $(COPTS) phi.c goldberg.c lyndon.c bch.c -o bch
