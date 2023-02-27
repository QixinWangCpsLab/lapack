#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <chrono>

#include "lapacke.h"

using namespace std;
using namespace std::chrono;

static FILE *ap_data, *bx_data;

void user_call_dppsv() {
  int					layout	= 	LAPACK_COL_MAJOR;
  char 				uplo 		= 	'L'; // 'U'
  lapack_int 	n 			=		200;
  lapack_int 	nrhs		=		1;
  lapack_int	ldb 		=		n;

  double *ap = (double *) malloc (n * (n + 1) / 2 * sizeof(double));
  double *bx = (double *) malloc (nrhs * ldb * sizeof(double));

  char sd[32];
  int k = 0;
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      fscanf(ap_data, "%s", sd);
      if(j >= i) {
        ap[k++] = strtod(sd, nullptr);
      }
    }
  }

  for(int i = 0; i < nrhs * ldb; i++) {
    fscanf(bx_data, "%s", sd);
    bx[i] = strtod(sd, nullptr);
  }

  try {
    lapack_int info;

    auto start = system_clock::now();
    info = LAPACKE_dppsv(layout, uplo, n, nrhs, ap, bx, ldb);
    auto end = system_clock::now();
    auto duration = duration_cast<nanoseconds>(end - start);
    fprintf(stderr, "\033[032mTotal Cost\n%ld (s) : %-9ld (ns)\033[0m\n",
        duration.count() / (long) 1e9, duration.count() % (long) 1e9);

    fprintf(stderr, "info=%d\n", info);
    for(int i = 0; i < n * (n + 1) / 2; i++) {
      fprintf(stdout, "%lf ", ap[i]);
    }
    fprintf(stdout, "\n");

    for(int i = 0; i < nrhs * ldb; i++) {
      fprintf(stdout, "%lf ", bx[i]);
    }
    fprintf(stdout, "\n");

  } catch (const char *msg) {
    fprintf(stderr, "client throws: %s\n", msg);
    fprintf(stderr, "errno: %s\n", strerror(errno));
    exit(EXIT_FAILURE);
  }
}

int main(int argc, char *argv[]) {
  if(argc != 3) {
    fprintf(stderr, "make test\n");
    exit(EXIT_FAILURE);
  }
  ap_data = fopen(argv[1], "r");
  bx_data = fopen(argv[2], "r");

  user_call_dppsv();
  return 0;
}
