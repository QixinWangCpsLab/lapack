#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <chrono>

#include "lapacke.h"

using namespace std;
using namespace std::chrono;

double *a, *s, *u, *vt, *work;
double *a_, *s_, *u_, *vt_, *work_;
int layout;
char jobu, jobvt;
lapack_int m, n, lda, ldu, ldvt;
lapack_int len_A, len_S, len_U, len_VT, len_W;

int ms_index = -1;
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
    /* CALL LAPACKE*/
    info = LAPACKE_dgesvd(layout, jobu, jobvt, m, n,
        a_, lda, s_, u_, ldu, vt_, ldvt, work_);

    auto end = system_clock::now();
    auto duration = duration_cast<nanoseconds>(end - start);
    fprintf(stderr, "\033[032mTotal Cost\n%ld (s) : %-9ld (ns)\033[0m\n",
        duration.count() / (long) 1e9, duration.count() % (long) 1e9);

    fprintf(stderr, "info=%d\n", info);
#define VERIFY
#ifdef VERIFY
    /* Print singular values */
    print_matrix( "Singular values", 1, n, s_, 1);
    /* Print left singular vectors */
    print_matrix( "Left singular vectors", m, n, u_, ldu);
    /* Print right singular vectors */
    print_matrix( "Right singular vectors", n, n, vt_, ldvt);
#endif

  } catch (const char *msg) {
    fprintf(stderr, "client throws: %s\n", msg);
    fprintf(stderr, "errno: %s\n", strerror(errno));
    exit(EXIT_FAILURE);
  }
}

int main(int argc, char *argv[]) {  

  if(argc != 2) {
    fprintf(stderr, "Args: [ms_index]\n");
    exit(EXIT_FAILURE);
  }

  ms_index = atoi(argv[1]);
  fprintf(stderr, "Using core %d.\n", ms_index);

  layout  = LAPACK_ROW_MAJOR;
  jobu    = 'A';
  jobvt   = 'A';
  m       = 6;
  n       = 5;
  lda     = n;
  ldu     = m;
  ldvt    = n;
  len_A   = lda * m;
  len_S   = n;
  len_U   = ldu * m;
  len_VT  = ldvt * n;
  len_W   = (m <= n ? m : n) - 1;

  a = (double *) malloc (len_A * sizeof(double));
  s = (double *) malloc (len_S * sizeof(double));
  u = (double *) malloc (len_U * sizeof(double));
  vt = (double *) malloc (len_VT * sizeof(double));
  work = (double *) malloc (len_W * sizeof(double));
  a_ = (double *) malloc (len_A * sizeof(double));
  s_ = (double *) malloc (len_S * sizeof(double));
  u_ = (double *) malloc (len_U * sizeof(double));
  vt_ = (double *) malloc (len_VT * sizeof(double));
  work_ = (double *) malloc (len_W * sizeof(double));

  double sa[len_A] = {
    8.79,   9.93,   9.83,   5.45,   3.16,
    6.11,   6.91,   5.04,  -0.27,   7.98,
   -9.15,  -7.93,   4.86,   4.85,   3.01,
    9.57,   1.64,   8.83,   0.74,   5.80,
   -3.49,   4.02,   9.80,   10.00,  4.27,
    9.84,   0.15,  -8.99,  -6.02,  -5.31
  };

  memcpy(a, sa, len_A * sizeof(double));

#ifdef LOOP
  while(true)
#endif
    user_call_dgesvd();

  return 0;
}

void print_matrix(const char* desc, lapack_int m, lapack_int n,
    double* mat, lapack_int lda) {
  int i, j;
  fprintf(stderr, "\n %s\n", desc );
  for( i = 0; i < m; i++) {
    for( j = 0; j < n; j++) 
      fprintf(stderr, " %6.2f", mat[i * lda + j]);
    fprintf(stderr, "\n" );
  }
}

