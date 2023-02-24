#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <signal.h>
#include <errno.h>
#include <execinfo.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>

#include <time.h>
#include <thread>

#include "lapacke.h"
#include "/home/parallels/workspace/lapack/CBLAS/include/cblas.h"

#include "soccs_sce_client_server_circular_queue.h"
#include "soccs_sce_client_server_protocol.h"
#include "soccs_sce_tools.h"
#include "soccs_sce_server_lapack_dppsv.h"

static int fd_shm;
static void *p_shm = nullptr;
static struct circular_queue *request_queue = nullptr;
static struct circular_queue *reply_queue = nullptr;
static pthread_t server_pthread;
static bool terminate = false;

static void sigint_handler(int sig) {
  if(sig != SIGINT) return;
  terminate = true;
  fprintf(stderr, "SIGINT received!\n");
}

static void *server_soccs(void *p_arg) {

  /* server soccs_sce_lapack_dppsv */
  struct sched_param server_pthread_sched_param = {
    .sched_priority = sched_get_priority_max(SCHED_FIFO)
  };
  pthread_setschedparam(pthread_self(), SCHED_FIFO, 
      &server_pthread_sched_param);

  struct circular_queue *request_queue = 
    ((struct soccs_sce_server_pthread_argument *)p_arg)->p_c2s_queue;
  struct circular_queue *reply_queue =
    ((struct soccs_sce_server_pthread_argument *)p_arg)->p_s2c_queue;

  try {
    while (terminate == false) {
      /* work loop */

      // TODO timer

      /* acquire request queue lock */
      pthread_mutex_lock(&request_queue->meta_info.mutex);
      while (!can_shallow_dequeue_straight(request_queue)) {
        /* wait until any request packets or timeout */
        struct timespec abs_timeout = {
          .tv_sec = 0, .tv_nsec = 0
        };
        if(clock_gettime(CLOCK_REALTIME, &abs_timeout) < 0)
          throw "clock_gettime failed.\n";
        calc_timespec_sum(&wait_timeout, &abs_timeout);
        pthread_cond_timedwait(&request_queue->meta_info.can_consume,
            &(request_queue->meta_info.mutex), &abs_timeout);
        if(terminate == true) break;
      }

      if(terminate == true) {
        pthread_mutex_unlock(&(request_queue->meta_info.mutex));
        break;
      }

#ifdef DEBUGGING
      fprintf(stderr, "server: can shallow dequeue straightly now.\n");
#endif

      /* receive request packet header */
      int request_packet_size = 0;
      struct request_packet_header *request_header = 
        (struct request_packet_header *) shallow_dequeue_straight
        (request_queue, &request_packet_size);

      if(request_packet_size <= 0 || request_header == nullptr)
        throw "request_packet_size <= 0 || request_header == nullptr";

      /* work */
      switch (request_header->function_id)
      {
        case LAPACK_DPPSV:
          {

#ifdef DEBUGGING
            fprintf(stderr, "received packet from client %u\n", 
                request_header->client_pthread_id);
#endif
            struct request_lapack_dppsv_fixed *req_pkt_fixed = 
              reinterpret_cast<struct request_lapack_dppsv_fixed *>
              (request_header + 1);
            const char uplo       = req_pkt_fixed->uplo;
            const lapack_int n    = req_pkt_fixed->n;
            const lapack_int nrhs = req_pkt_fixed->nrhs;
            const lapack_int ldb  = req_pkt_fixed->ldb;
            const lapack_int len_AP  = n * (n + 1) / 2;
            const lapack_int len_B   = ldb * nrhs;
            lapack_int info          = 0;
            double *ap  = (double *) malloc (len_AP * sizeof(double));
            double *b   = (double *) malloc (len_B * sizeof(double));
/*
n = 4
ap = n*(n+1)/2 = 10
b = n = 4
ap[i][j] = ap[?]
*/
            double *req_pkt_flexible = 
              reinterpret_cast<double *> (req_pkt_fixed + 1);
            for(int i = 0; i < len_AP; i++)
              ap[i] = req_pkt_flexible[i];
            for(int i = 0; i < len_B; i++)
              b[i] = req_pkt_flexible[len_AP + i];

#ifdef DEBUGGING
            fprintf(stderr, "ap[%d]={", len_AP);
            for(int i = 0; i < len_AP; i++)
              fprintf(stderr, " %lf", ap[i]);
            fprintf(stderr, " }\n");
            fprintf(stderr, "b[%d]={", len_B);
            for(int i = 0; i < len_B; i++)
              fprintf(stderr, " %lf", b[i]);
            fprintf(stderr, " }\n");
#endif

            /* sanity checking */
            if (uplo != 'U' && uplo != 'u' && uplo != 'L' && uplo != 'l')
              info = -1;
            else if (n < 0)
              info = -2;
            else if (nrhs < 0)
              info = -3;
            else if (ldb < max(1, n))
              info = -6;
            if (info != 0)
              goto reply; // no need to calculate

            /* calculation starts from here */
            /* 190: dpptrf( uplo, n, ap) */
            int j, jc, jj;
            double ajj;
            if(uplo == 'U' || uplo == 'u') {
              /* 177: Compute the Choleskey factorization A = U**T**U */
              jj = 0;
              for(j = 1; j <= n; j++) {
                jc = jj + 1;
                jj = jj + j - 1;
                /* 184: Compute elements 0:j-2 of column j */ 
                if(j > 1) {
                  // dtpsv('U', 'T', 'N', j - 1, ap, ap + jc - 1, 1);
                  int kk_2 = 1, k_2, i_2, j_2;
                  double t_2;
                  for(j_2 = 1; j_2 <= j - 1; j_2++) {
                    t_2 = ap[jc + j_2 - 1];
                    k_2 = kk_2;
                    for(i_2 = 1; i_2 <= j_2 - 1; i_2++) {
                      t_2 = t_2 - ap[k_2 - 1] * ap[jc - 1 + i_2 - 1];
                      k_2 = k_2 + 1;
                    }
                    t_2 = t_2 / ap[kk_2 + j_2 - 1];
                    ap[jc - 1 + j_2 - 1] = t_2;
                    kk_2 = kk_2 + j_2;
                  }
                  // end of dtpsv():286-298-----------------------
                }
                // ddot(j - 1, AP[jc], 1, AP[jc], 1)
                double ddot = 0.0, dtemp = 0.0;
                int m_3 = (j - 1) % 5, i_3, mp1;
                if(m_3 != 0) {
                  for(i_3 = 1; i_3 <= m_3; i_3++)
                    dtemp = dtemp + ap[jc - 1 + i_3 - 1] * ap[jc - 1 + i_3 - 1];
                  if(j - 1 < 5) {
                    ddot = dtemp;
                    goto break_ddot;
                  }
                }
                mp1 = m_3 + 1;
                for(i_3 = mp1; i_3 <= j - 1; i_3 += 5) {
                  dtemp = dtemp + 
                    ap[jc + i_3 - 2] * ap[jc + i_3 - 2] +
                    ap[jc + i_3 - 1] * ap[jc + i_3 - 1] +
                    ap[jc + i_3 + 0] * ap[jc + i_3 + 0] +
                    ap[jc + i_3 + 1] * ap[jc + i_3 + 1] +
                    ap[jc + i_3 + 2] * ap[jc + i_3 + 2] ;
                }
                ddot = dtemp;
break_ddot:
                // end of ddot():113-127---------------------------

                /* 190: Compute U(j, j) and test for non-positive-definitenss */
                ajj = ap[jj] - ddot;
                if(ajj <= 0) {
                  ap[jj] = ajj;
                  info = j;
                  goto reply;
                }
                ap[jj] = sqrt(ajj);
              }
            }
            else {
              /* 206: Compute the Cholesky factorization A = L*L**T */
              // TODO
            }

            // end of dpptrf()--------------------------------------

            /* 193: Solve A * X = B, overwriting B with X */
            /* 195: dpptrs(uplo, n, nrhs, ap, b, ldb) */

            // end of dpptrs()--------------------------------------

reply:
            /* prepare reply packets */
#ifdef DEBUGGING
            fprintf(stderr, "server: ready to reply\n");
#endif
            int rpl_packet_payload_size = sizeof(reply_lapack_dppsv_fixed) +
              len_AP * sizeof(double) + len_B *sizeof(double);
            int rpl_packet_total_size = sizeof(struct reply_packet_header) + 
              rpl_packet_payload_size;

            // TODO timer

            pthread_mutex_lock(&(reply_queue->meta_info.mutex));
            while(!can_malloc_straight(reply_queue, rpl_packet_total_size)) {
              /* wait until there is space in reply queue or timeout */
              struct timespec abs_timeout = {
                .tv_sec = 0, .tv_nsec = 0
              };
              if(clock_gettime(CLOCK_REALTIME, &abs_timeout) < 0)
                throw "clock_gettime failed.\n";
              calc_timespec_sum(&wait_timeout, &abs_timeout);
              pthread_cond_timedwait(&(reply_queue->meta_info.can_produce),
                  &(reply_queue->meta_info.mutex), &abs_timeout);
              if(terminate == true) break;
            }

            if(terminate == true) {
              pthread_mutex_unlock(&(reply_queue->meta_info.mutex));
              pthread_mutex_unlock(&(request_queue->meta_info.mutex));
              break; // break switch
            }

#ifdef DEBUGGING
            fprintf(stderr, "server: can enqueue now.\n");
#endif

            uint8_t *raw_reply_packet = 
              malloc_straight(reply_queue, rpl_packet_total_size);
            if(raw_reply_packet == nullptr)
              throw "raw_reply_packet == nullptr";

            /* load packet header -> reply queue */
            struct reply_packet_header *reply_header = 
              (struct reply_packet_header *) raw_reply_packet;
            reply_header->common_header.preamble    = (uint32_t) (PREAMBLE);
            reply_header->common_header.packet_size = rpl_packet_total_size;
            reply_header->client_pthread_id = request_header->client_pthread_id;
            reply_header->transaction_id    = request_header->transaction_id;
            reply_header->function_id       = request_header->function_id;
            reply_header->payload_len       = rpl_packet_payload_size;

            /* load fixed part -> reply queue */
            struct reply_lapack_dppsv_fixed *rpl_pkt_fixed = 
              reinterpret_cast<reply_lapack_dppsv_fixed *>
              (reply_header + 1);
            rpl_pkt_fixed->info   = info;

            /* load flexible part -> reply queue */
            double *rpl_pkt_flexible = 
              reinterpret_cast<double *> (rpl_pkt_fixed + 1);
            for(int i = 0; i < len_AP; i++)
              rpl_pkt_flexible[i] = ap[i];
            for(int i = 0; i < len_B; i++)
              rpl_pkt_flexible[len_AP + i] = b[i];
            // TODO timer

            /* broadcast to clients and release lock */
            pthread_cond_broadcast(&(reply_queue->meta_info.can_consume));
            pthread_mutex_unlock(&(reply_queue->meta_info.mutex));

            pthread_cond_broadcast(&(request_queue->meta_info.can_produce));
            pthread_mutex_unlock(&(request_queue->meta_info.mutex));

#ifdef DEBUGGING
            fprintf(stderr, "finish 1 DPPSV.\n");
#endif

            break;
          }

      } // end of switch

      if(terminate == true) {
        pthread_mutex_unlock(&(reply_queue->meta_info.mutex));
        pthread_mutex_unlock(&(reply_queue->meta_info.mutex));
        break; // break while
      }

    } // end of while
  } catch (const char *msg) {
    fprintf(stderr, "worker thread exception: %s", msg);
    return nullptr;
  }

  return nullptr;
}

int main(int argc, char *argv[], char *envp[]) {

  /* server */
  if(argc != 3) {
    fprintf(stderr, "Usage: ./server <num_cores> <server_core_index>\n");
    exit(EXIT_FAILURE);
  }
  signal(SIGINT, sigint_handler);

  try {
    /* bind cpu core to server thread */
    int num_cpus = atoi(argv[1]);
    int server_core_index = atoi(argv[2]);
    cpu_set_t *cpuset = CPU_ALLOC(num_cpus);
    size_t cpuset_size = CPU_ALLOC_SIZE(num_cpus);
    CPU_ZERO_S(cpuset_size, cpuset);
    CPU_SET_S(server_core_index, cpuset_size, cpuset);
    pthread_attr_t server_thread_attr;
    pthread_attr_init(&server_thread_attr);
    pthread_attr_setaffinity_np(&server_thread_attr, cpuset_size, cpuset);

    if ((fd_shm = shm_open(SHARED_MEM_NAME_LAPACK_DPPSV, 
            O_RDWR | O_CREAT /*| O_EXCL*/, 0666)) == -1)
      throw "shm_open failed.\n";
    if (ftruncate(fd_shm, SHARED_MEM_SIZE) == -1)
      throw "ftrancate failed.\n";
    if ((p_shm = mmap(NULL, SHARED_MEM_SIZE, PROT_READ | PROT_WRITE, 
            MAP_SHARED, fd_shm, 0)) == MAP_FAILED)
      throw "mmap failed.\n";

    request_queue = (struct circular_queue *) p_shm;
    init_circular_queue(request_queue);
    reply_queue = (struct circular_queue *) 
      (((uint8_t *) p_shm) + sizeof(struct circular_queue));
    init_circular_queue(reply_queue);

    struct soccs_sce_server_pthread_argument pthread_arg = {
      .p_c2s_queue = request_queue,
      .p_s2c_queue = reply_queue
    };

    fprintf(stderr, "Shared memory '%s'of %u bytes "
        "created and mmaped at virtual address %p.\n",
        SHARED_MEM_NAME_LAPACK_DPPSV, SHARED_MEM_SIZE, p_shm);
    fprintf(stderr, "c2s_queue of %lu bytes "
        "allocated and initialized at virtual address %p.\n",
        sizeof(struct circular_queue), request_queue);
    fprintf(stderr, "s2c_queue of %lu bytes "
        "allocated and initialized at virtual address %p.\n",
        sizeof(struct circular_queue), reply_queue);

    /* create server thread */
    pthread_create(&server_pthread, &server_thread_attr, server_soccs, 
        (void*) &pthread_arg);
    fprintf(stderr, "Main server thread %u created.\n",
        pthread_t_to_uint32_t(server_pthread));

    /* wait for the server thread */
    pthread_join(server_pthread, NULL);
    munmap(p_shm, SHARED_MEM_SIZE);
    close(fd_shm);
    shm_unlink(SHARED_MEM_NAME_LAPACK_DPPSV);
    fprintf(stderr, "Main server thread: exited.\n");
  } catch (const char *msg) {
    fprintf(stderr, "\n Exception: %s", msg);
    fprintf(stderr, "\n errno: %s\n", strerror(errno));
  }

  fprintf(stderr, "\nOK\n");

  return 0;
}

