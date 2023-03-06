LAPACK_HOME := ..

CC := g++

RM := rm -f

CFLAGS := -O3 -g -ggdb -std=c++11 -Wall \
																	-Wno-unused-result

HD_DIR := -I$(LAPACK_HOME)/LAPACKE/include

LDFLAGS := 	-llapacke -llapack -lrefblas \
						-lgfortran -lm

LDSOCCS :=	-lsoccs_sce_lapack

LIB_DIR := 	-L$(LAPACK_HOME)/libs

REF_DIR :=	-L$(LAPACK_HOME)/ref_libs

SOURCES := $(wildcard *.cpp)

SOCCS := $(SOURCES:%.cpp=%_soccs)

REFER := $(SOURCES:%.cpp=%_refer)

all: $(SOCCS) $(REFER)

%_soccs: %.o
	$(CC) $^ $(CFLAGS) -o $@ $(LIB_DIR) $(LDFLAGS) $(LDSOCCS)

%_refer: %.o
	$(CC) $^ $(CFLAGS) -o $@ $(REF_DIR) $(LDFLAGS)

%.o: %.cpp
	$(CC) $^ -c $(CFLAGS) -o $@ $(HD_DIR)

.PRECIOUS: %.o

.PHONY: clean test_dppsv
clean:
	$(RM) *_soccs *_refer *.o *.rs

test_dppsv:
	make -B -s
	@			./dppsv_refer ap_200x200.data bx_200x200.data > dppsv_refer.rs \
		&& 	./dppsv_soccs ap_200x200.data bx_200x200.data > dppsv_soccs.rs \
		&& 	diff dppsv_refer.rs dppsv_soccs.rs

