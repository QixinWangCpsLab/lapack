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

static FILE *a_data, *timer_log;
double *a, *s, *u, *vt, *work;
double *a_, *s_, *u_, *vt_, *work_;
int layout;
char jobu, jobvt;
lapack_int m, n, lda, ldu, ldvt;
lapack_int len_A, len_S, len_U, len_VT, len_W;

int cl_index, ms_index = -1;
char t_file[32];

void print_matrix(const char*, lapack_int, lapack_int, double*, lapack_int);

void user_call_dgesvd() { 
  memcpy(a_, a, len_A * sizeof(double));
  memset(s_, 0x0, len_S * sizeof(double));
  memset(u_, 0x0, len_U * sizeof(double));
  memset(vt_, 0x0, len_VT * sizeof(double));
  memset(work_, 0x0, len_W * sizeof(double));

  try {
    lapack_int info;

    auto start = system_clock::now();
    /* CALL LAPACKE */
    info = LAPACKE_dgesvd(layout, jobu, jobvt, m, n,
        a_, lda, s_, u_, ldu, vt_, ldvt, work_);

    auto end = system_clock::now();
    auto duration = duration_cast<nanoseconds>(end - start);
    fprintf(stderr, "\033[032mTotal Cost\n%ld (s) : %-9ld (ns)\033[0m\n",
        duration.count() / (long) 1e9, duration.count() % (long) 1e9); 
    timer_log = fopen(t_file, "a");
    fprintf(timer_log, "%d %ld\n", getpid(), duration.count());
    fclose(timer_log);

    fprintf(stderr, "info=%d\n", info);

#ifdef VERIFY
    /* Print singular values */
    print_matrix( "Singular values", 1, n, s_, 1);
    /* Print left singular vectors */
    print_matrix( "Left singular vectors", m, n, u_, ldu);
    /* Print right singular vectors */
    print_matrix( "Right singular vectors", n, n, vt_, ldvt);
#endif
    struct timespec to_sleep = {
      .tv_sec = 0, .tv_nsec = 7500000 - duration.count()
    };
    nanosleep(&to_sleep, NULL);

  } catch (const char *msg) {
    fprintf(stderr, "client throws: %s\n", msg);
    fprintf(stderr, "errno: %s\n", strerror(errno));
    exit(EXIT_FAILURE);
  }
}

int main(int argc, char *argv[]) {  

  if(argc != 4) {
    fprintf(stderr, "Args: [cl_index] [ms_index] [a_data]\n");
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

  layout  = LAPACK_ROW_MAJOR;
  jobu    = 'A';
  jobvt   = 'A';
  m       = 100;
  n       = 100;
  lda     = n;
  ldu     = m;
  ldvt    = n;
  len_A   = lda * m;
  len_S   = n;
  len_U   = ldu * m;
  len_VT  = ldvt * n;
  len_W   = (m <= n ? m : n) - 1;

  a     = (double *) malloc (len_A * sizeof(double));
  s     = (double *) malloc (len_S * sizeof(double));
  u     = (double *) malloc (len_U * sizeof(double));
  vt    = (double *) malloc (len_VT * sizeof(double));
  work  = (double *) malloc (len_W * sizeof(double));
  a_    = (double *) malloc (len_A * sizeof(double));
  s_    = (double *) malloc (len_S * sizeof(double));
  u_    = (double *) malloc (len_U * sizeof(double));
  vt_   = (double *) malloc (len_VT * sizeof(double));
  work_ = (double *) malloc (len_W * sizeof(double));

  a_data = fopen(argv[3], "r");
  sprintf(t_file, "stat/log_dgesvd_%d.txt", ms_index);
  timer_log = fopen(t_file, "w");
  fclose(timer_log);

  char sd[32];
  for(int i = 0; i < len_A; i++) {
    fscanf(a_data, "%s", sd);
    a[i] = strtod(sd, nullptr);
  }

#ifdef LOOP
  while(true)
#endif
    user_call_dgesvd();

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

