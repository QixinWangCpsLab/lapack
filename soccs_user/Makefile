WORKSPACE := /home/parallels/workspace

CC := g++

CFLAGS := -g -Wall

HD_DIR := -I$(WORKSPACE)/lapack/LAPACKE/include

LDFLAGS := 	-llapacke -llapack -lrefblas \
						-lgfortran -lm

LDSOCCS :=	-lsoccs_sce_lapack

LIB_DIR := 	-L$(WORKSPACE)/lapack/libs

REF_DIR :=	-L$(WORKSPACE)/lapack/ref_libs

all: main ref_main

main: main.o
	$(CC) $^ $(CFLAGS) -o $@ $(LIB_DIR) $(LDFLAGS) $(LDSOCCS)

ref_main: main.o
	$(CC) $^ $(CFLAGS) -o $@ $(REF_DIR) $(LDFLAGS)

main.o: main.cpp
	$(CC) $^ -c $(CFLAGS) -o $@ $(HD_DIR)

.PHONY: clean
clean:
	rm -f main ref_main *.o

