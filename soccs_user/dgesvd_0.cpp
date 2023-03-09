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

void print_matrix(const char*, lapack_int, lapack_int, double*, lapack_int);

void user_call_dgesvd() { 
  memcpy(a_, a, len_A * sizeof(double));
  memcpy(s_, s, len_S * sizeof(double));
  memcpy(u_, u, len_U * sizeof(double));
  memcpy(vt_, vt, len_VT * sizeof(double));
  memcpy(work_, work, len_W * sizeof(double));

  try {
    lapack_int info;

    auto start = system_clock::now();
    info = LAPACKE_dgesvd(layout, jobu, jobvt, m, n,
        a_, lda, s_, u_, ldu, vt_, ldvt, work_);
    auto end = system_clock::now();
    auto duration = duration_cast<nanoseconds>(end - start);
    fprintf(stderr, "\033[032mTotal Cost\n%ld (s) : %-9ld (ns)\033[0m\n",
        duration.count() / (long) 1e9, duration.count() % (long) 1e9);

    fprintf(stderr, "info=%d\n", info);

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
  /* 
     if(argc != 3) {
     fprintf(stderr, "Args: [ap_data] [bx_data]\n");
     exit(EXIT_FAILURE);
     }
   */

  layout  = LAPACK_COL_MAJOR;
  jobu    = 'A';
  jobvt   = 'A';
  m       = 6;
  n       = 5;
  lda     = m;
  ldu     = m;
  ldvt    = n;
  len_A   = lda * n;
  len_S   = m <= n ? m : n;
  len_U   = m * m;
  len_VT  = n * n;
  len_W   = m <= n ? m : n;

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
    8.79,  6.11, -9.15,  9.57, -3.49,  9.84,
    9.93,  6.91, -7.93,  1.64,  4.02,  0.15,
    9.83,  5.04,  4.86,  8.83,  9.80, -8.99,
    5.45, -0.27,  4.85,  0.74, 10.00, -6.02,
    3.16,  7.98,  3.01,  5.80,  4.27, -5.31
  };
  memcpy(a, sa, len_A * sizeof(double));

  /*
     ap_data = fopen(argv[1], "r");
     bx_data = fopen(argv[2], "r");

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
   */

#ifdef LOOP
  while(true)
#endif
    user_call_dgesvd();

  return 0;
}

void print_matrix(const char* desc, lapack_int m, lapack_int n,
    double* mat, lapack_int ld) {
  int i, j;
  fprintf(stderr, "\n %s\n", desc );
  for( i = 0; i < m; i++) {
    for( j = 0; j < n; j++) 
      fprintf(stderr, " %6.2f", mat[i + j * lda]);
    fprintf(stderr, "\n" );
  }
}

