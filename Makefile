all: bch 

CC = gcc 
#CC = clang 

#CFLAGS = -O3 -fPIC -march=ivybridge  -fopenmp  -Wall 
CFLAGS = -O3 -fPIC  -msse4.1   -fopenmp  -Wall  
#CFLAGS = -O3 -fPIC -msse4.1 -fopenmp -fsanitize=signed-integer-overflow -fsanitize=undefined -Wall 
#CFLAGS = -g -fPIC  -Wall 

MAKE_SHARED_LIB = $(CC) -fopenmp -shared

SHARED_LIB = libbch.so

DEPS = bch.h khash.h
OBJS = phi.o lie_series.o lyndon.o rightnormed.o goldberg.o \
       convert_lyndon.o convert_rightnormed.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(SHARED_LIB): $(OBJS)
	$(MAKE_SHARED_LIB) -o $(SHARED_LIB) $(OBJS) 

bch: $(SHARED_LIB) bch.h 
	$(CC) $(CFLAGS) bch.c -o bch -L. -lbch

clean:
	rm -f *.o $(SHARED_LIB) bch

bch_goldberg_30.txt: bch
	./bch goldberg_coefficients=1 N=30 > bch_goldberg_30.txt

bch_lyndon_20.txt: bch
	./bch N=20 verbosity_level=1 > bch_lyndon_20.txt

bch_rightnormed_20.txt: bch
	./bch N=20 rightnormed_basis=1 verbosity_level=1 > bch_rightnormed_20.txt

