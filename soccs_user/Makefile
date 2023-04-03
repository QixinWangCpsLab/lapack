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

SOCCS := $(SOURCES:%_user.cpp=%_user_soccs)

REFER := $(SOURCES:%_user.cpp=%_user_refer)

STARTER := $(SOURCES:%_starter.cpp=%_starter)

all: $(SOCCS) $(REFER) $(STARTER)

%_soccs: %.o $(LAPACK_HOME)/libs/libsoccs_sce_lapack.a
	$(CC) $^ $(CFLAGS) -o $@ $(LIB_DIR) $(LDFLAGS) $(LDSOCCS)

%_refer: %.o
	$(CC) $^ $(CFLAGS) -o $@ $(REF_DIR) $(LDFLAGS)

%_starter: %_starter.o
	$(CC) $^ $(CFLAGS) -o $@

%.o: %.cpp
	$(CC) $^ -c $(CFLAGS) -o $@ $(HD_DIR)

.PRECIOUS: %.o %_starter.o

.PHONY: clean test_dppsv test_dgesvd loop

clean:
	$(RM) *_soccs *_refer *_starter *.o *.rs *.txt stat/ results/
	mkdir stat

test_dppsv:
	@if [ ! -d "results" ]; then mkdir results; fi
	$(MAKE) -B -s CFLAGS+=-DVERIFY
	@			./dppsv_user_refer 0 0 data/dppsv_ap_200x200.data data/dppsv_bx_200x200.data \
											> results/refer_dppsv.rs \
		&& 	./dppsv_user_soccs 0 0 data/dppsv_ap_200x200.data data/dppsv_bx_200x200.data \
											> results/soccs_dppsv.rs \
		&& 	diff results/refer_dppsv.rs results/soccs_dppsv.rs

test_dgesvd:
	# TODO

loop:
	$(MAKE) -B CFLAGS+=-DLOOP

