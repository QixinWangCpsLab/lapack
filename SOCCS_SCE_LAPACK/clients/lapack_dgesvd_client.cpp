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

#define NUM_CPUS 8
#define CLIENT_CORE_INDEX 0

static int fd_shm   = -1;
static void *p_shm  = nullptr;

static circular_queue *request_queue  = nullptr;
static circular_queue *reply_queue    = nullptr;

void soccs_sce_lapack_degsvd (
    char *jobu, char *jobvt,
    lapack_int *m, lapack_int *n,
    double *a, lapack_int *lda, double *s, double *u,
    lapack_int *ldu, double *vt, lapack_int *ldvt,
    double *work, double *lwork,
    lapack_int *info) {

  /* client soccs_sce_dppsv */
  fprintf(stderr, "\033[31mcalling soccs_sce_lapack_dppsv\033[0m\n");

  uint32_t client_pthread_id = pthread_t_to_uint32_t(pthread_self());
  static uint64_t soccs_sce_dppsv_transaction_id = 0;
  uint64_t transaction_id = soccs_sce_dppsv_transaction_id++;

  if (fd_shm == -1) {
    if((fd_shm = shm_open(SHARED_MEM_NAME_LAPACK_DPPSV, O_RDWR, 0666)) == -1)
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
      SHARED_MEM_NAME_LAPACK_DPPSV, SHARED_MEM_SIZE, p_shm);
  fprintf(stderr, "c2s_queue of %lu bytes found at virtual address %p.\n",
      sizeof(struct circular_queue), request_queue);
  fprintf(stderr, "s2c_queue of %lu bytes found at virtual address %p.\n",
      sizeof(struct circular_queue), reply_queue);
#endif



}
