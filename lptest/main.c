#include <stdio.h>
#include "lapacke.h"

int main(int argc, char *argv[]) {

	int					layout	= 	LAPACK_ROW_MAJOR;
	char 				uplo 		= 	'L'; // 'U'
	lapack_int 	n 			=		2;
	lapack_int 	nrhs		=		1;
	lapack_int	ldb 		=		2;

	double ap[10] = {
		4.16, -3.12, 5.03, 0.56, -0.83, 0.76, -0.10, 1.18, 0.34, 1.18
	};

	double bx[4] = {
		8.7, -13.35, 1.89, -4.14
	}; 

	lapack_int info = LAPACKE_dppsv(layout, uplo, n, nrhs, ap, bx, ldb);

	printf("info=%d\n", info);

	for(int i = 0; i < 10; i++) {
		printf("%lf ", ap[i]);
	}
	printf("\n");

	for(int i = 0; i < 4; i++) {
			printf("%lf ", bx[i]);
	}
	printf("\n");

	return 0;
}
