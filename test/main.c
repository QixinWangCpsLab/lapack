#include "cblas.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(int argc, char *argv[]) {
	CBLAS_LAYOUT 	layout 	= CblasRowMajor;
	CBLAS_TRANSPOSE	TransA	= CblasNoTrans;
	CBLAS_INT		M		= atoi(argv[1]);
	CBLAS_INT		N		= atoi(argv[2]);
	float			alpha	= 1.0;
	float 			*A 		= (float*) malloc (sizeof(float) * M * N);
	CBLAS_INT 		lda		= N;

	int len_X = 0, len_Y = 0;

	if (TransA == CblasNoTrans){
		len_X = N;
		len_Y = M;
	}else{
		len_X = M;
		len_Y = N;
	}

	float			*X		= (float*) malloc (sizeof(float) * len_X);
	CBLAS_INT		incX	= 1;
	float			beta	= 1.0;
	float			*Y		= (float*) malloc (sizeof(float) * len_Y);
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
	for(int i = 0; i < N * M; i++)
		A[i] = atof(argv[vi++]);
	for(int i = 0; i < len_X; i++)
		X[i] = atof(argv[vi++]);
	for(int i = 0; i < len_Y; i++)
		Y[i] = atof(argv[vi++]);

	cblas_sgemv(layout, TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);

	for(int i = 0; i < len_Y; i++)
		printf("%f\n", Y[i]); 
	return 0;
}
