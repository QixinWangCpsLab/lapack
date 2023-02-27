LAPACKHOME := ..

CC := g++

CFLAGS := -O3 -g -ggdb -std=c++11 $(DFLAG) -Wall \
															-Wno-address-of-packed-member \
															-Wno-unused-variable
HD_DIR := -I$(LAPACKHOME)/LAPACKE/include

LDFLAGS := -pthread -lrt

LIB_DIR := 

DFLAG :=

CLIENT_OBJS += 	lapack_dppsv_client.o \
								soccs_sce_client_server_circular_queue.o \
								soccs_sce_tools.o

SERVER_OBJS += 	lapack_dppsv_server.o \
								soccs_sce_client_server_circular_queue.o \
								soccs_sce_tools.o

# PACK CLIENT -> libsoccs_sce_lapack.a

all: server_soccs libsoccs_sce_lapack.a
	@rm -f /dev/shm/shm_SoCCS_SCE_LAPACK_DPPSV	

server_soccs: $(SERVER_OBJS)
	$(CC) $^ $(CFLAGS) -o $@ $(HD_DIR) $(LDFLAGS)

libsoccs_sce_lapack.a: $(CLIENT_OBJS)
	$(AR) $(ARFLAGS) $@ $^ && \
	cp $@ $(LAPACKHOME)/libs/

%.o: %.cpp
	$(CC) $^ -c $(CFLAGS) -o $@ $(HD_DIR) $(LDFLAGS)

.PHONY: clean debugging
clean:
	rm -f server_soccs *.a *.o 
	@rm -f /dev/shm/shm_SoCCS_SCE_LAPACK_DPPSV perf.data perf.data.old

debugging:
	make DFLAG=-DDEBUGGING all

