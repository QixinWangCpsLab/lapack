#include "soccs_sce_client_server_circular_queue.h"
#include "soccs_sce_client_server_protocol.h"
#include "soccs_sce_tools.h"

extern "C" {
  void
    soccs_sce_lapack_dgesvd (
        char *jobu, char *jobvt,
        lapack_int *m, lapack_int *n,
        double *a, lapack_int *lda, double *s, double *u,
        lapack_int *ldu, double *vt, lapack_int *ldvt,
        double *work, double *lwork,
        lapack_int *info);
}

static int fd_shm   = -1;
static void *p_shm  = nullptr;

static circular_queue *request_queue  = nullptr;
static circular_queue *reply_queue    = nullptr;

extern int ms_index;

static bool terminated = false;
static void sigint_handler(int sig)
{
  if (sig != SIGINT)
    return;
  terminated = true;
  fprintf(stderr, "SIGINT received!\n");
}

void soccs_sce_lapack_dgesvd (
    char *jobu, char *jobvt,
    lapack_int *m, lapack_int *n,
    double *a, lapack_int *lda, double *s, double *u,
    lapack_int *ldu, double *vt, lapack_int *ldvt,
    double *work, double *lwork,
    lapack_int *info) {

  signal(SIGINT, sigint_handler);
  /* client soccs_sce_dgesvd */
#ifdef DEBUGGING
  fprintf(stderr, "\033[31mcalling soccs_sce_lapack_dgesvd\033[0m\n");
#endif

  uint32_t client_pthread_id = pthread_t_to_uint32_t(pthread_self());
  static uint64_t soccs_sce_dppsv_transaction_id = 0;
  uint64_t transaction_id = soccs_sce_dppsv_transaction_id++;

  if (fd_shm == -1) {
    char shm_name[32];
    if(ms_index == -1)
      strcpy(shm_name, SHARED_MEM_NAME_LAPACK);
    else
      sprintf(shm_name, "%s_%d", SHARED_MEM_NAME_LAPACK, ms_index);
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

  const lapack_int len_A  = (*lda) * (*n);
  const lapack_int len_S  = (*m) <= (*n) ? (*m) : (*n);
  const lapack_int len_U  = (*m) * (*ldu);
  const lapack_int len_VT = (*n) * (*ldvt);
  const lapack_int len_W  = (*lwork) >= 1 ? (*lwork) : 1;

  int req_packet_payload_size = sizeof(struct request_lapack_dgesvd_fixed) + 
    len_A * sizeof(double);

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
  request_header->function_id               = LAPACK_DGESVD;
  request_header->payload_len               = req_packet_payload_size;

  /* load fixed part -> request queue */
  struct request_lapack_dgesvd_fixed *req_pkt_fixed = 
    reinterpret_cast<struct request_lapack_dgesvd_fixed *> 
    (request_header + 1);

  req_pkt_fixed->jobu   = *jobu;
  req_pkt_fixed->jobvt  = *jobvt;
  req_pkt_fixed->m      = *m;
  req_pkt_fixed->n      = *n;
  req_pkt_fixed->lda    = *lda;
  req_pkt_fixed->ldu    = *ldu;
  req_pkt_fixed->ldvt   = *ldvt;
  req_pkt_fixed->lwork  = *lwork;

  /* load flexible part -> request_queue */
  double *req_pkt_flexible = reinterpret_cast<double *> (req_pkt_fixed + 1);

  memcpy(req_pkt_flexible, a, len_A * sizeof(double));

  /* broadcast to servers and release lock */
  pthread_cond_broadcast(&(request_queue->meta_info.can_consume));
  pthread_mutex_unlock(&request_queue->meta_info.mutex);

#ifdef DEBUGGING
  fprintf(stderr, "PREAMBLE: %u\npreamble: %u\nsize: %d\n", 
      (uint32_t) PREAMBLE,
      request_header->common_header.preamble,
      request_header->common_header.packet_size);
  fprintf(stderr, "len_A=%d, len_S=%d, len_U=%d, len_VT=%d, len_W=%d\n",
      len_A, len_S, len_U, len_VT, len_W);
  fprintf(stderr, "m=%d, n=%d, lda=%d, ldu=%d, ldvt=%d, lwork=%d\n",
      req_pkt_fixed->m, req_pkt_fixed->n, req_pkt_fixed->lda, 
      req_pkt_fixed->ldu, req_pkt_fixed->ldvt, req_pkt_fixed->lwork);
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
  if(reply_header->function_id != LAPACK_DGESVD)
    throw "wrong function_id != LAPACK_DGESVD\n";

  /* get fixed part <- reply queue */
  struct reply_lapack_dppsv_fixed *rpl_pkt_fixed = 
    reinterpret_cast<struct reply_lapack_dppsv_fixed *> (reply_header + 1);
  lapack_int r_info = rpl_pkt_fixed->info;

  /* get flexible part <- reply queue */
  double *rpl_pkt_flexible =
    reinterpret_cast<double *> (rpl_pkt_fixed + 1);

  double *rpl_pkt_flexible_pos = rpl_pkt_flexible;
  memcpy(a, rpl_pkt_flexible_pos, len_A * sizeof(double));
  rpl_pkt_flexible_pos += len_A;
  memcpy(s, rpl_pkt_flexible_pos, len_S * sizeof(double));
  rpl_pkt_flexible_pos += len_S;
  memcpy(u, rpl_pkt_flexible_pos, len_U * sizeof(double));
  rpl_pkt_flexible_pos += len_U;
  memcpy(vt, rpl_pkt_flexible_pos, len_VT * sizeof(double));
  rpl_pkt_flexible_pos += len_VT;
  memcpy(work, rpl_pkt_flexible_pos, len_W * sizeof(double));


  auto end0 = system_clock::now();
  auto duration0 = duration_cast<nanoseconds> (end0 - start0);
  auto duration1 = duration_cast<nanoseconds> (end1 - start1);

#ifdef TIMER
  fprintf(stderr, "\033[032mC-S-C Cost\n%ld (s) : %-9ld (ns)\033[0m\n",
      duration0.count() / (long) 1e9, duration0.count() % (long) 1e9);
  fprintf(stderr, "\033[032mServer Cost\n%ld (s) : %-9ld (ns)\033[0m\n",
      duration1.count() / (long) 1e9, duration1.count() % (long) 1e9);
#endif

  /* free memory space in reply queue */
  reply_header = (struct reply_packet_header *)
    shallow_dequeue_straight(reply_queue, &reply_packet_size);

  pthread_cond_broadcast(&(reply_queue->meta_info.can_produce));
  pthread_mutex_unlock(&(reply_queue->meta_info.mutex));

#ifdef DEBUGGING
  fprintf(stderr, "client reached end.\n");
#endif

  if(terminated == true) {
    fprintf(stderr, "EXIT\n");
    exit(EXIT_SUCCESS);
  }

}

