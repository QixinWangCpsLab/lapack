LAPACK_HOME := ..

CC := g++

RM := rm -f

CP := cp

CFLAGS 	:= 	-O3 -g -ggdb -std=c++11 $(DFLAG) -Wall \
															-Wno-address-of-packed-member \
															-Wno-unused-variable

HD_DIR 	:= 	-I$(LAPACK_HOME)/LAPACKE/include \
						-I./include

LDFLAGS := -pthread -lrt

DFLAG 	:=

COMMON_SRC	:=	$(wildcard ./commons/*.cpp)

CLIENT_SRC 	:= 	$(wildcard ./clients/*.cpp)

SERVER_SRC	:=	$(wildcard ./*.cpp)

CLIENT_OBJS += 	$(COMMON_SRC:%.cpp=%.o) \
								$(CLIENT_SRC:%.cpp=%.o)

SERVER_OBJS += 	$(SERVER_SRC:%.cpp=%.o) \
								$(COMMON_SRC:%.cpp=%.o)

# PACK CLIENTS -> libsoccs_sce_lapack.a

all: server_soccs libsoccs_sce_lapack.a
	@$(RM) 	/dev/shm/shm_SoCCS_SCE_LAPACK_DPPSV \
					/dev/shm/shm_SoCCS_SCE_LAPACK_DGESVD

server_soccs: $(SERVER_OBJS)
	$(CC) $^ $(CFLAGS) -o $@ $(HD_DIR) $(LDFLAGS)

libsoccs_sce_lapack.a: $(CLIENT_OBJS)
	$(AR) $(ARFLAGS) $@ $^ && \
	$(CP) $@ $(LAPACK_HOME)/libs/

%.o: %.cpp
	$(CC) $^ -c $(CFLAGS) -o $@ $(HD_DIR) $(LDFLAGS)

.PHONY: clean debugging
clean:
	$(RM) $(CLIENT_OBJS) $(SERVER_OBJS) server_soccs *.a
	@$(RM) perf.data perf.data.old \
				/dev/shm/shm_SoCCS_SCE_LAPACK_DPPSV \
				/dev/shm/shm_SoCCS_SCE_LAPACK_DGESVD

debugging:
	$(MAKE) -B DFLAG=-DDEBUGGING all

