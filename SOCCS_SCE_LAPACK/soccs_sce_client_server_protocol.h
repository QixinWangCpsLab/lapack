#ifndef _SOCCS_SCE_PROTOCOL_H_
#define _SOCCS_SCE_PROTOCOL_H_

#include <stdint.h>
#include "lapacke.h"
#include "soccs_sce_client_server_circular_queue.h"

#define MAX_REQEUST_PAYLOAD_LEN \
  (MAX_CIRCULAR_QUEUE_PACKET_SIZE - sizeof(struct request_packet_header))

#define MAX_REPLY_PAYLOAD_LEN \
  (MAX_CIRCULAR_QUEUE_PACKET_SIZE - sizeof(struct reply_packet_header))

enum function_id : uint16_t {
  LAPACK_DPPSV,
  // UNDEFINED  // undefined function id
};

struct request_packet_header {
  struct    circular_queue_packet_common_header common_header;
  uint32_t  client_pthread_id;        // client's pthread_t_to_uint32_t(pthread_self())
  uint64_t  transaction_id;           // 8 byte
  enum      function_id function_id;  // 2 byte, decides payload format
  uint32_t  payload_len;              // must be in [0, MAX_REQUEST_PAYLOAD_LEN]
} __attribute__((packed));

struct reply_packet_header {
  struct    circular_queue_packet_common_header common_header;
  uint32_t  client_pthread_id;        // client's pthread_t_to_uint32_t(pthread_self())
  uint64_t  transaction_id;           // 8 byte
  enum function_id function_id;       // 2 byte, decides payload format
  uint32_t payload_len;               // must be in [0, MAX_REPLY_PAYLOAD_LEN]
} __attribute__((packed));

struct request_lapack_dppsv_fixed {
  char uplo;
  lapack_int n;
  lapack_int nrhs;
  // double *ap;  flexible
  // double *b;   flexible
  lapack_int ldb;
} __attribute__((packed));

struct reply_lapack_dppsv_fixed {
  // double *ap;  flexible
  // double *b;   flexible
  lapack_int info;
} __attribute__((packed));

#define SHARED_MEM_NAME_LAPACK_DPPSV "/shm_SoCCS_SCE_LAPACK_DPPSV"

#endif

