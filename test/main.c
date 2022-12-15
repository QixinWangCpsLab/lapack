#include "cblas.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(int argc, char *argv[]) {
	CBLAS_LAYOUT 	layout 	= CblasRowMajor;
	CBLAS_TRANSPOSE	TransA	= CblasNoTrans;
	FILE *input = fopen("test.in", "r");
	CBLAS_INT	M, N;
	fscanf(input, "%d", &M);
	fscanf(input, "%d", &N);

	float	alpha	= 1.0;
	float *A 	= (float*) malloc (sizeof(float) * M * N);
	CBLAS_INT lda	= N;
	int len_X = 0, len_Y = 0;

	if (TransA == CblasNoTrans) {
		len_X = N;
		len_Y = M;
	} else {
		len_X = M;
		len_Y = N;
	}

	float			*X		= (float*) malloc (sizeof(float) * len_X);
	CBLAS_INT		incX	= 1;
	float			beta	= 1.0;
	float			*Y		= (float*) malloc (sizeof(float) * len_Y);
	float			*TY		= (float*) malloc (sizeof(float) * len_Y);
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
	int vi = 3;
	for(int i = 0; i < N * M; i++) fscanf(input, "%f", &A[i]);
	for(int i = 0; i < len_X; i++) fscanf(input, "%f", &X[i]);
	for(int i = 0; i < len_Y; i++) fscanf(input, "%f", &Y[i]);
	while(1) {
		memcpy(TY, Y, sizeof(float) * len_Y);
		cblas_sgemv(layout, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
		printf("Finish 1 round!\n");
	}
	return 0;
}
