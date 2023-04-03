#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <chrono>
#include <unistd.h>
#include <sched.h>

#include "lapacke.h"

#define NUM_CPUS 8

using namespace std;
using namespace std::chrono;

static FILE *ap_data, *bx_data, *timer_log;
double *ap, *bx, *ap_, *bx_;
int layout;
lapack_int n, nrhs, ldb;
char uplo;

int cl_index, ms_index = -1;
char t_file[32];

void user_call_dppsv() {
  /* generate data, call dppsv */ 
  memcpy(ap_, ap, n * (n + 1) / 2 * sizeof(double));
  memcpy(bx_, bx, nrhs * ldb * sizeof(double));

  try {
    lapack_int info;

    auto start = system_clock::now();
    info = LAPACKE_dppsv(layout, uplo, n, nrhs, ap_, bx_, ldb);
    auto end = system_clock::now();
    auto duration = duration_cast<nanoseconds>(end - start);
    fprintf(stderr, "\033[032mTotal Cost\n%ld (s) : %-9ld (ns)\033[0m\n",
        duration.count() / (long) 1e9, duration.count() % (long) 1e9);
    timer_log = fopen(t_file, "a");
    fprintf(timer_log, "%d %ld\n", getpid(), duration.count());
    fclose(timer_log);

    fprintf(stderr, "info=%d\n", info);

#ifdef VERIFY
    for(int i = 0; i < n * (n + 1) / 2; i++) {
      fprintf(stdout, "%lf ", ap_[i]);
    }
    fprintf(stdout, "\n");

    for(int i = 0; i < nrhs * ldb; i++) {
      fprintf(stdout, "%lf ", bx_[i]);
    }
    fprintf(stdout, "\n");
#endif

  } catch (const char *msg) {
    fprintf(stderr, "client throws: %s\n", msg);
    fprintf(stderr, "errno: %s\n", strerror(errno));
    exit(EXIT_FAILURE);
  }
}

int main(int argc, char *argv[]) {

  if(argc != 4) {
    fprintf(stderr, "Args: [cl_index] [ms_index] [ap_data] [bx_data]\n");
    exit(EXIT_FAILURE);
  }

  cl_index = atoi(argv[1]);
  ms_index = atoi(argv[2]);
  fprintf(stderr, "Using core %d.\n", ms_index);

  cpu_set_t *cpuset = CPU_ALLOC(NUM_CPUS);
  size_t cpuset_size = CPU_ALLOC_SIZE(NUM_CPUS);
  CPU_ZERO_S(cpuset_size, cpuset);
  CPU_SET_S(cl_index, cpuset_size, cpuset);
  sched_setaffinity(getpid(), cpuset_size, cpuset);

  layout  =   LAPACK_COL_MAJOR;
  uplo    =   'L'; // 'U'
  n       =   200;
  nrhs    =   1;
  ldb     =   n;

  ap  = (double *) malloc (n * (n + 1) / 2 * sizeof(double));
  bx  = (double *) malloc (nrhs * ldb * sizeof(double));
  ap_ = (double *) malloc (n * (n + 1) / 2 * sizeof(double));
  bx_ = (double *) malloc (nrhs * ldb * sizeof(double));

  ap_data = fopen(argv[3], "r");
  bx_data = fopen(argv[4], "r");
  sprintf(t_file, "stat/log_dppsv_%d.txt", ms_index);
  timer_log = fopen(t_file, "w");
  fclose(timer_log);

  char sd[32];
  int k = 0;
  for(int i = 0; i < n; i++) 
  {
    for(int j = 0; j < n; j++) 
    {
      fscanf(ap_data, "%s", sd);
      if(j >= i) ap[k++] = strtod(sd, nullptr);
    }
  }

  for(int i = 0; i < nrhs * ldb; i++) {
    fscanf(bx_data, "%s", sd);
    bx[i] = strtod(sd, nullptr);
  }

#ifdef LOOP
  while(true)
#endif
    user_call_dppsv();

  return 0;
}

