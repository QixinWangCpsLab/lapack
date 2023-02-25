#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "lapacke.h"

int main(int argc, char *argv[]) {

  int					layout	= 	LAPACK_COL_MAJOR;
  char 				uplo 		= 	'L'; // 'U'
  lapack_int 	n 			=		4;
  lapack_int 	nrhs		=		1;
  lapack_int	ldb 		=		4;

  double ap[10] = {
    4.16, -3.12, 0.56, -0.1, 5.03, -0.83, 1.18, 0.76, 0.34, 1.18
//    4, 2, 6, 17, -1, 14
  };

  double bx[4] = {
    8.70, -13.35, 1.89, -4.14
//    30, 19, 56
  }; 

  try {
    lapack_int info = LAPACKE_dppsv(layout, uplo, n, nrhs, ap, bx, ldb);
    
    printf("info=%d\n", info);
    for(int i = 0; i < (n*(n+1)/2); i++) {
      printf("%lf ", ap[i]);
    }
    printf("\n");

    for(int i = 0; i < n; i++) {
      printf("%lf ", bx[i]);
    }
    printf("\n");

  } catch (const char *msg) {
    fprintf(stderr, "client throws: %s", msg);
    fprintf(stderr, "errno: %s\n", strerror(errno));
  }

  return 0;
}
