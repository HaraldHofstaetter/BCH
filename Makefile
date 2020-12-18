all: bch 

CC = gcc 
#CC = clang 

#COPTS = -O3  -march=ivybridge  -fopenmp  -Wall -DUSE_INT128_T 
COPTS = -O3   -msse4.1   -fopenmp  -Wall -DUSE_INT128_T -DSIMD_VECTORIZED 
#COPTS = -O3  -msse4.1 -fopenmp -fsanitize=signed-integer-overflow -fsanitize=undefined -Wall -DUSE_INT128_T 
#COPTS = -g -Wall -DUSE_INT128_T 

bch: bch.h bch.c phi.c lie_series.c rightnormed.c goldberg.c convert_lyndon.c convert_rightnormed.c

	$(CC) $(COPTS) phi.c goldberg.c rightnormed.c lie_series.c bch.c -o bch

bch_goldberg_30.txt: bch
	./bch goldberg_coefficients=1 N=30 > bch_goldberg_30.txt

bch_lyndon_20.txt: bch
	./bch N=20 M=14 verbosity_level=1 > bch_lyndon_20.txt

bch_rightnormed_20.txt: bch
	./bch N=20 rightnormed_basis=1 verbosity_level=1 > bch_rightnormed_20.txt

