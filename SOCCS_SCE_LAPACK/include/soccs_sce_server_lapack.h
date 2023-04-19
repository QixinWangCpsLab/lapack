#ifndef _SOCCS_SCE_SERVER_LAPACK_H_
#define _SOCCS_SCE_SERVER_LAPACK_H_

struct soccs_sce_server_pthread_argument {
  struct circular_queue *p_c2s_queue;
  struct circular_queue *p_s2c_queue;
} __attribute__((packed));

#endif
