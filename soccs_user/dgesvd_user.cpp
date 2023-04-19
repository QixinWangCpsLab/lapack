#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <chrono>
#include <unistd.h>
#include <sched.h>
#include <random>

#include "lapacke.h"

#define NUM_CPUS 8

static FILE *timer_log;
double *a, *s, *u, *vt, *work;
char jobu, jobvt;

int ms_index = -1;
int cl_pid;
int cl_index, cl_prior = -1;
char t_file[32];

void user_call_dgesvd();
void print_matrix(const char*, lapack_int, lapack_int, double*, lapack_int);

void user_call_dgesvd() { 

  auto start0 = std::chrono::steady_clock::now();

  const int  layout  = LAPACK_ROW_MAJOR;
  const char jobu    = 'A';
  const char jobvt   = 'A';
  const int  m       = 100;
  const int  n       = 100;
  const int  lda     = n;
  const int  ldu     = m;
  const int  ldvt    = n;
  const int  len_A   = lda * m;
  const int  len_S   = n;
  const int  len_U   = ldu * m;
  const int  len_VT  = ldvt * n;
  const int  len_W   = (m <= n ? m : n) - 1;

  a     = (double *) malloc (len_A * sizeof(double));
  s     = (double *) malloc (len_S * sizeof(double));
  u     = (double *) malloc (len_U * sizeof(double));
  vt    = (double *) malloc (len_VT * sizeof(double));
  work  = (double *) malloc (len_W * sizeof(double));

  const double data_min = -100.0;
  const double data_max =  100.0;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> uniform_double(data_min, data_max);

  for(int i = 0; i < len_A; i++)
    a[i] = uniform_double(gen);

  lapack_int info;

  auto start = std::chrono::steady_clock::now();
  /* CALL LAPACKE */
  info = LAPACKE_dgesvd(layout, jobu, jobvt, m, n,
      a, lda, s, u, ldu, vt, ldvt, work);

  auto end = std::chrono::steady_clock::now();
  auto duration = std::chrono::duration_cast
    <std::chrono::nanoseconds>(end - start);
  timer_log = fopen(t_file, "a");
  fprintf(timer_log, "%d %ld\n", getpid(), duration.count());
  fclose(timer_log);

  //  fprintf(stderr, "\033[032mTotal Cost\n%ld (s) : %-9ld (ns)\033[0m\n",
  //      duration.count() / (long) 1e9, duration.count() % (long) 1e9); 
  //  fprintf(stderr, "info=%d\n", info);


#ifdef VERIFY
  /* Print singular values */
  print_matrix( "Singular values", 1, n, s, 1);
  /* Print left singular vectors */
  print_matrix( "Left singular vectors", m, n, u, ldu);
  /* Print right singular vectors */
  print_matrix( "Right singular vectors", n, n, vt, ldvt);
#endif

  auto end0 = std::chrono::steady_clock::now();
  auto duration0 = std::chrono::duration_cast
    <std::chrono::nanoseconds>(end0 - start0);

  struct timespec to_sleep = {
    .tv_sec = 0, .tv_nsec = 10000000 - duration0.count()
  };
  nanosleep(&to_sleep, NULL);

  free(a);
  free(s);
  free(u);
  free(vt);
  free(work);
}

int main(int argc, char *argv[]) {  

  if(argc != 4) {
    fprintf(stderr, "Args: [cl_index] [ms_index] [cl_priority]\n");
    exit(EXIT_FAILURE);
  }

  cl_pid = getpid();
  cl_index = atoi(argv[1]);
  ms_index = atoi(argv[2]);
  cl_prior = atoi(argv[3]);

  /* set client CPU affinity */
  cpu_set_t *cpuset = CPU_ALLOC(NUM_CPUS);
  size_t cpuset_size = CPU_ALLOC_SIZE(NUM_CPUS);
  CPU_ZERO_S(cpuset_size, cpuset);
  CPU_SET_S(cl_index, cpuset_size, cpuset);
  sched_setaffinity(getpid(), cpuset_size, cpuset);
  /* set RR priority */
  struct sched_param rr_param = { 
    .sched_priority = cl_prior
  };
  sched_setscheduler(cl_pid, SCHED_FIFO, &rr_param);

  sprintf(t_file, "stat/log_dgesvd_%d.txt", ms_index);
  timer_log = fopen(t_file, "w");
  fclose(timer_log);

  sleep(1);
  while(true) {
//    try {
      user_call_dgesvd();
//    } catch (const char *msg) {
//      fprintf(stderr, "client throws: %s\n", msg);
//      fprintf(stderr, "errno: %s\n", strerror(errno));
//      exit(EXIT_FAILURE);
//    }
#ifndef LOOP
    break;
#endif
  }

  return 0;
}

void print_matrix(const char* desc, lapack_int m, lapack_int n,
    double* mat, lapack_int lda) {
  int i, j;
  fprintf(stdout, "\n %s\n", desc );
  for( i = 0; i < m; i++) {
    for( j = 0; j < n; j++) 
      fprintf(stdout, " %6.2f", mat[i * lda + j]);
    fprintf(stdout, "\n" );
  }
}

