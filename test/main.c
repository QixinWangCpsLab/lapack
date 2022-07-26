#include "cblas.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(int argc, char *argv[]) {
	CBLAS_LAYOUT 	layout 	= CblasRowMajor;
	CBLAS_TRANSPOSE	TransA	= CblasNoTrans;
	CBLAS_INT		M		= 10;
	CBLAS_INT		N		= 20;
	float			alpha	= 1.0;
	float 			*A 		= (float*) malloc (sizeof(float) * M * N);
	CBLAS_INT 		lda		= M;
	float			*X		= (float*) malloc (sizeof(float) * N);
	CBLAS_INT		incX	= 1;
	float			beta	= 1.0;
	float			*Y		= (float*) malloc (sizeof(float) * M);
	CBLAS_INT		incY	= 1;
	/*
	float Ac[] = {	1.0, 1.0, 1.0,
					1.0, 1.0, 1.0,
					1.0, 1.0, 1.0,
					1.0, 1.0, 1.0,
					1.0, 1.0, 1.0};
	memcpy(A, Ac, sizeof(Ac));
	float Xc[] = {1.0, 1.0, 1.0};
	memcpy(X, Xc, sizeof(Xc));
	float Yc[] = {1.0, 1.0, 1.0, 1.0, 1.0};
	memcpy(Y, Yc, sizeof(Yc));
	*/
	int vi = 1;
	for(int i = 0; i < N * M; i++)
		A[i] = atof(argv[vi++]);
	for(int i = 0; i < N; i++)
		X[i] = atof(argv[vi++]);
	for(int i = 0; i < M; i++)
		Y[i] = atof(argv[vi++]);

	cblas_sgemv(layout, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);

	for(int i = 0; i < M; i++)
		printf("%f\n", Y[i]); 
	return 0;
}
