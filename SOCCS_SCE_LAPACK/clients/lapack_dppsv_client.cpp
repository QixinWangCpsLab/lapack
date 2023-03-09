#include "soccs_sce_client_server_circular_queue.h"
#include "soccs_sce_client_server_protocol.h"
#include "soccs_sce_tools.h"

extern "C" {
  void 
    soccs_sce_lapack_dppsv (
        char *uplo, lapack_int *n, lapack_int *nrhs,
        double *ap, double *b, lapack_int *ldb,
        lapack_int *info);
}

#define NUM_CPUS 8
#define CLIENT_CORE_INDEX 0

static int fd_shm   = -1;
static void *p_shm  = nullptr;

static circular_queue *request_queue  = nullptr;
static circular_queue *reply_queue    = nullptr;

extern int ms_index;

void soccs_sce_lapack_dppsv (
    char *uplo, lapack_int *n, lapack_int *nrhs,
    double *ap, double *b, lapack_int *ldb,
    lapack_int *info) {

  /* client soccs_sce_dppsv */
  fprintf(stderr, "\033[31mcalling soccs_sce_lapack_dppsv\033[0m\n");

  uint32_t client_pthread_id = pthread_t_to_uint32_t(pthread_self());
  static uint64_t soccs_sce_dppsv_transaction_id = 0;
  uint64_t transaction_id = soccs_sce_dppsv_transaction_id++;

  if (fd_shm == -1) {
    char shm_name[32];
    if(ms_index == -1)
      strcpy(shm_name, SHARED_MEM_NAME_LAPACK);
    else
      sprintf(shm_name, "%s_%d", SHARED_MEM_NAME_LAPACK, ms_index);
    fprintf(stderr, "%s\n", shm_name);
    if((fd_shm = shm_open(shm_name, O_RDWR, 0666)) == -1)
      throw "shm_open failed.\n";

    struct stat tmp;
    if (fstat(fd_shm, &tmp) == -1)
      throw "fstat failed.\n";
    if (tmp.st_size != SHARED_MEM_SIZE)
      throw "shm size != SHARED_MEM_SIZE\n";
    if ((p_shm = mmap(NULL, SHARED_MEM_SIZE, PROT_READ | PROT_WRITE, 
            MAP_SHARED, fd_shm, 0)) == MAP_FAILED)
      throw "mmap failed.\n";
    request_queue = (struct circular_queue *) p_shm;
    reply_queue = request_queue + 1;
  }

#ifdef DEBUGGING
  fprintf(stderr, "Shared memory '%s'of %u bytes found and mmaped "
      "at virtual address %p.\n",
      SHARED_MEM_NAME_LAPACK, SHARED_MEM_SIZE, p_shm);
  fprintf(stderr, "c2s_queue of %lu bytes found at virtual address %p.\n",
      sizeof(struct circular_queue), request_queue);
  fprintf(stderr, "s2c_queue of %lu bytes found at virtual address %p.\n",
      sizeof(struct circular_queue), reply_queue);
#endif

  const lapack_int len_AP  = (*n) * ((*n) + 1) / 2;
  const lapack_int len_B   = (*ldb) * (*nrhs);

  int req_packet_payload_size = sizeof(struct request_lapack_dppsv_fixed) + 
    len_AP * sizeof(double) + len_B * sizeof(double);

  int req_packet_total_size = sizeof(struct request_packet_header) +
    req_packet_payload_size;

  /* acquire request queue lock */
  pthread_mutex_lock(&(request_queue->meta_info.mutex));
  while (!can_malloc_straight(request_queue, req_packet_total_size)) {
    /* wait until enough space in the request queue */
    pthread_cond_wait(&(request_queue->meta_info.can_produce),
        &(request_queue->meta_info.mutex));
  }

#ifdef DEBUGGING
  fprintf(stderr, "client: can enqueue now.\n");
#endif

  auto start0 = system_clock::now();

  uint8_t *raw_request_packet = 
    malloc_straight(request_queue, req_packet_total_size);
  if(raw_request_packet == nullptr)
    throw "raw_request_packet == nullptr\n";

  /* load packet header -> request_queue */
  struct request_packet_header *request_header = 
    (struct request_packet_header *) raw_request_packet;

  request_header->common_header.preamble    = (uint32_t) (PREAMBLE);
  request_header->common_header.packet_size = req_packet_total_size;
  request_header->client_pthread_id         = client_pthread_id;
  request_header->transaction_id            = transaction_id;
  request_header->function_id               = LAPACK_DPPSV;
  request_header->payload_len               = req_packet_payload_size;

  /* load fixed part -> request queue */
  struct request_lapack_dppsv_fixed *req_pkt_fixed = 
    reinterpret_cast<struct request_lapack_dppsv_fixed *> 
    (request_header + 1);

  req_pkt_fixed->uplo = *uplo;
  req_pkt_fixed->n    = *n;
  req_pkt_fixed->nrhs = *nrhs;
  req_pkt_fixed->ldb  = *ldb;

  /* load flexible part -> request_queue */
  double *req_pkt_flexible = reinterpret_cast<double *> (req_pkt_fixed + 1);
  double *req_pkt_flexible_pos = req_pkt_flexible;

  memcpy(req_pkt_flexible_pos, ap, len_AP * sizeof(double));
  req_pkt_flexible_pos += len_AP;
  memcpy(req_pkt_flexible_pos, b, len_B * sizeof(double));

  /* broadcast to servers and release lock */
  pthread_cond_broadcast(&(request_queue->meta_info.can_consume));

#ifdef DEBUGGING
  fprintf(stderr, "client: broadcast request queue -> can_consume.\n");
#endif

  pthread_mutex_unlock(&request_queue->meta_info.mutex);

#ifdef DEBUGGING
  fprintf(stderr, "PREAMBLE:%u\npreamble: %u\nsize: %d\n", 
      (uint32_t) PREAMBLE,
      request_header->common_header.preamble,
      request_header->common_header.packet_size);
  fprintf(stderr, "len_AP=%d, len_B=%d, n=%d, nrhs=%d, ldb=%d\n",
      len_AP, len_B, req_pkt_fixed->n, req_pkt_fixed->nrhs, req_pkt_fixed->ldb);
#endif

blocking_waiting_for_reply:

  /* acquire reply queue lock */
  pthread_mutex_lock(&(reply_queue->meta_info.mutex));
  auto start1 = system_clock::now();

  /* wait until reply_queue is not empty */
  while(!can_shallow_dequeue_straight(reply_queue)) {
    pthread_cond_wait(&(reply_queue->meta_info.can_consume),
        &(reply_queue->meta_info.mutex));
  }

  /* get packet header <- reply queue */
  int reply_packet_size = 0;
  struct reply_packet_header *reply_header = 
    (struct reply_packet_header *) 
    peek_shallow_dequeue_straight(reply_queue, &reply_packet_size);

  if(reply_packet_size <= 0 || reply_header == nullptr)
    throw "reply_packet_size <= 0 || reply_header == nullptr\n";

  if(reply_header->client_pthread_id != client_pthread_id) {
    /* reply packets for others */
    pthread_mutex_unlock(&(reply_queue->meta_info.mutex));
    sched_yield();
    goto blocking_waiting_for_reply;
  }

  /* reply packet for me */
  auto end1 = system_clock::now();
  if(reply_header->transaction_id != transaction_id)
    throw "wrong transaction_id\n";
  if(reply_header->function_id != LAPACK_DPPSV)
    throw "wrong function_id != LAPACK_DPPSV\n";

  /* get fixed part <- reply queue */
  struct reply_lapack_dppsv_fixed *rpl_pkt_fixed = 
    reinterpret_cast<struct reply_lapack_dppsv_fixed *> (reply_header + 1);
  *info = rpl_pkt_fixed->info;

  /* get flexible part <- reply queue */
  double *rpl_pkt_flexible =
    reinterpret_cast<double *> (rpl_pkt_fixed + 1);
  double *rpl_pkt_flexible_pos = rpl_pkt_flexible;

  memcpy(ap, rpl_pkt_flexible_pos, len_AP * sizeof(double));
  rpl_pkt_flexible_pos += len_AP;
  memcpy(b, rpl_pkt_flexible_pos, len_B * sizeof(double));

  auto end0 = system_clock::now();
  auto duration0 = duration_cast<nanoseconds> (end0 - start0);
  auto duration1 = duration_cast<nanoseconds> (end1 - start1);
  fprintf(stderr, "\033[032mC-S-C Cost\n%ld (s) : %-9ld (ns)\033[0m\n",
      duration0.count() / (long) 1e9, duration0.count() % (long) 1e9);
  fprintf(stderr, "\033[032mServer Cost\n%ld (s) : %-9ld (ns)\033[0m\n",
      duration1.count() / (long) 1e9, duration1.count() % (long) 1e9);

  /* free memory space in reply queue */
  reply_header = (struct reply_packet_header *)
    shallow_dequeue_straight(reply_queue, &reply_packet_size);

  pthread_cond_broadcast(&(reply_queue->meta_info.can_produce));
  pthread_mutex_unlock(&(reply_queue->meta_info.mutex));

#ifdef DEBUGGING
  fprintf(stderr, "client reached end.\n");
#endif

}

