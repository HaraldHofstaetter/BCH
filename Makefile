all: bch 

CC = gcc 
#CC = clang 

CFLAGS = -O3 -fPIC -march=native  -fopenmp  -Wall 
#CFLAGS = -O3 -fPIC -march=native -fopenmp -fsanitize=signed-integer-overflow -fsanitize=undefined -Wall 

#CFLAGS = -g -fPIC -Wall 
#CFLAGS = -g -fPIC -fsanitize=address -fsanitize=signed-integer-overflow -fsanitize=undefined -Wall  

MAKE_SHARED_LIB = $(CC) -fopenmp -shared

SHARED_LIB = libbch.so

DEPS = bch.h khash.h
OBJS = phi.o lie_series.o lyndon.o rightnormed.o goldberg.o \
       convert_lyndon.o convert_rightnormed.o convert_hall.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(SHARED_LIB): $(OBJS)
	$(MAKE_SHARED_LIB) -o $(SHARED_LIB) $(OBJS) 

bch: $(SHARED_LIB) bch.h bch.c 
	$(CC) $(CFLAGS) bch.c -o bch -L. -lbch


clean:
	rm -f *.o $(SHARED_LIB) bch


tables:	bch_lyndon_20.txt bch_rightnormed_20.txt bch_hall_20.txt

bch_lyndon_20.txt: bch
	./bch N=20 verbosity_level=1 > bch_lyndon_20.txt

bch_rightnormed_20.txt: bch
	./bch N=20 basis=1 verbosity_level=1 > bch_rightnormed_20.txt

bch_hall_20.txt: bch
	./bch N=20 basis=2 verbosity_level=1 > bch_hall_20.txt



#Compile to WebAssembly (if the Emscripten SDK is properly installed):

wasm: bch.wasm bch.js bch.html

SRCS_WASM = bch.c phi.c lie_series.c lyndon.c rightnormed.c goldberg.c \
       convert_lyndon.c convert_rightnormed.c convert_hall.c

bch.wasm bch.js bch.html: $(DEPS) $(SRCS_WASM) shell_bch.html
	emcc $(SRCS_WASM) -O3 -s EXIT_RUNTIME=1 -s ALLOW_MEMORY_GROWTH=1 \
	    -o bch.html --shell-file ./shell_bch.html 
	

