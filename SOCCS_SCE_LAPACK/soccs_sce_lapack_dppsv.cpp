#include <stdio.h>
#include "lapacke.h"

extern "C" {
	void 
		soccs_sce_lapack_dppsv (char *uplo, lapack_int *n, lapack_int *nrhs,
				double *ap, double *b, lapack_int *ld,
				lapack_int *info) {
			printf("calling soccs_sce_lapack_dppsv\n");

		}
}
