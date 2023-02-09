WORKSPACE := /home/parallels/workspace

CFLAGS := -g -ggdb -std=c++11 -Wall
HD_DIR := -I$(WORKSPACE)/lapack/LAPACKE/include
LDFLAGS := -pthread -lrt
LIB_DIR := 

OBJS += $(patsubst %.cpp, %.o, $(wildcard *.cpp))

all: libsoccs_sce_lapack.a

libsoccs_sce_lapack.a: $(OBJS)
	$(AR) $(ARFLAGS) $@ $^

%.o: %.cpp
	$(CC) $^ -c $(CFLAGS) -o $@ $(HD_DIR) $(LDFLAGS)

.PHONY: clean
clean:
	rm -f *.a *.o
