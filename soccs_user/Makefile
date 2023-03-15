LAPACK_HOME := ..

CC := g++

RM := rm -rf

CFLAGS := -O3 -g -ggdb -std=c++11 -Wall \
																	-Wno-unused-result

HD_DIR := -I$(LAPACK_HOME)/LAPACKE/include

LDFLAGS := 	-llapacke -llapack -lrefblas \
						-lgfortran -lm

LDSOCCS :=	-lsoccs_sce_lapack -lrt

LIB_DIR := 	-L$(LAPACK_HOME)/libs

REF_DIR :=	-L$(LAPACK_HOME)/ref_libs

SOURCES := $(wildcard *.cpp)

SOCCS := $(SOURCES:%.cpp=%_soccs)

REFER := $(SOURCES:%_0.cpp=%_refer)

all: $(SOCCS) $(REFER) 

%_soccs: %.o $(LAPACK_HOME)/libs/libsoccs_sce_lapack.a
	$(CC) $^ $(CFLAGS) -o $@ $(LIB_DIR) $(LDFLAGS) $(LDSOCCS)

%_refer: %_0.o
	$(CC) $^ $(CFLAGS) -o $@ $(REF_DIR) $(LDFLAGS)

%.o: %.cpp
	$(CC) $^ -c $(CFLAGS) -o $@ $(HD_DIR)

.PRECIOUS: %.o

.PHONY: clean test_dppsv test_dgesvd loop

clean:
	$(RM) *_soccs *_refer *.o *.rs results/

test_dppsv:
	@if [ ! -d "results" ]; then mkdir results; fi
	$(MAKE) -B -s CFLAGS+=-DVERIFY
	@			./dppsv_refer data/dppsv_ap_200x200.data data/dppsv_bx_200x200.data \
											> results/refer_dppsv.rs \
		&& 	./dppsv_soccs data/dppsv_ap_200x200.data data/dppsv_bx_200x200.data \
											> results/soccs_dppsv.rs \
		&& 	diff results/refer_dppsv.rs results/soccs_dppsv.rs

test_dgesvd:
	# TODO

loop:
	$(MAKE) -B CFLAGS+=-DLOOP

