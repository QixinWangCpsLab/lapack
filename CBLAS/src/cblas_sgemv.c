/*
 *
 * cblas_sgemv.c
 * This program is a C interface to sgemv.
 * Originally written by Keita Teranishi
 * 4/6/1998
 * Revised by Qixin Wang since
 * 12/7/2022
 *
 */
#include "cblas.h"
#include "cblas_f77.h"
#include <stdio.h>
#include <stdlib.h>

float max(const float a, const float b){
   return ((a) > (b))?(a):(b);
}

char lsame(const char c, const char l){
   //Start of concern >>>
   //Here we assume ASCII encoding. The reference BLAS implementation
   //assumes other encoding also.
   //<<< End of concern

   //Start of to do >>>
   //Please insert parameter sanity checks here.
   //e.g. assert('A' <= (*l) <= 'Z')
   //<<< End of to do

   if ((c) == (l)) return 1;
   else if ('a' <= (c) && (c) <= 'z' && ((c) - ('a' - 'A')) == (l)) return 1;
   else if ('A' <= (c) && (c) <= 'Z' && ((c) + ('a' - 'A')) == (l)) return 1;
   else return 0;
} 

void soccs_sce_sgemv(const CBLAS_TRANSPOSE trans, const CBLAS_INT m, const CBLAS_INT n, 
      const float alpha, const float *A /* assume row major storage ??? */, const CBLAS_INT lda, 
      const float *X, const CBLAS_INT inc_X, const float beta, 
      float *Y, const CBLAS_INT inc_Y){
   
   //floating point number 1 is not 1.0e+0, 0 is not 0.0e+0
   const float ONE = 1.0e+0, ZERO = 0.0e+0;

   //local scalars
   int len_X = 0, len_Y = 0, k_X = 0, k_Y = 0;

   //Start of to do >>>
   //Please insert parameter sanity checks here.
   //See reference implementation in BLAS/SRC/sgemv.f line 193-211.
   //<<< End of to do

   //Quick return if possible.
   if (m == 0 || n == 0 || ((alpha == ZERO) && (beta == ONE))) return;

   //Deprected old usage: if (lsame(trans, 'N')){
   if (trans == CblasNoTrans){
      len_X = n;
      len_Y = m;
   }else{
      len_X = m;
      len_Y = n;
   }

   //Take for example trans == 'n', if X is a column of a row-major order stored matrix,
   //inc_X > 0 is the number of elements in a row, 
   //-inc_X is the number of elements in a row when enumerating X from back to front.
   //
   //One case why inc_X is needed is that sgemv is called by a routine that calculates
   //
   //  $\mathblack{Y} = \alpha A \mathblack{X} + \beta \mathblack{Y}$,
   //
   //where $\mathblack{Y}$ is a $m \times inc_Y$ matrix, 
   //and $\mathblack{X}$ is a $n \times inc_X$ matrix.
   //The routine will then call sgemv('N', $m$, $n$, $\alpha$, $A$, $m$, &(\mathblack{X}[0][k]),
   //$inc_X$, $\beta$, &(\mathblack{Y}[0][k]),  $inc_Y$) to calculate the $k$th column of
   //$\mathblack{Y}$.
   if (inc_X > 0){
      k_X = 0;
   }else{
      k_X = 0 - (len_X - 1) * inc_X;
   }

   if (inc_Y > 0){
      k_Y = 0;
   }else{
      k_Y = 0 - (len_Y - 1) * inc_Y;
   }

   //First assign $Y := \beta Y$.
   if (beta != ONE){
      if (inc_Y == 1){
         if (beta == ZERO)
            for (int i = 0; i < len_Y; i++) Y[i] = ZERO;
         else
            for (int i = 0; i < len_Y; i++) Y[i] *= beta;
      }else{
         int i_Y = k_Y;
         if (beta == ZERO){
            for (int i = 0; i < len_Y; i++){
               Y[i_Y] = ZERO;
               i_Y += inc_Y;
            }
         }else{
            for (int i = 0; i < len_Y; i++){
               Y[i_Y] *= beta;
               i_Y += inc_Y;
            }
         }
      }
   }
   //Now $Y := \beta Y$.

   //Next assign $Y := \alpha A X + Y$ or $Y := \alpha A^T X + Y$.
   if (alpha == ZERO) return;
   //Deprecated old usage: if (lsame(trans, 'N')){
   if (trans == CblasNoTrans){
      //$Y := \alpha A X + Y$.
      int j_X = k_X;
      if (inc_Y == 1){
         for (int j = 0; j < n; j++){
            float temp = alpha * X[j_X];
            for (int i = 0; i < m; i++){
               Y[i] += temp * A[i * n + j]; //accessing A[i][j], assuming row major storage of A.
            }
            j_X += inc_X;
         }
      }else{
         for (int j = 0; j < n; j++){
            float temp = alpha * X[j_X];
            int i_Y = k_Y;
            for (int i = 0; i < m; i++){
               Y[i_Y] += temp * A[i * n + j]; //accessing A[i][j], assuming row major storage of A.
               i_Y += inc_Y;
            }
            j_X += inc_X;
         }
      }
   }else{
      //$Y := \alpha A^T X + Y$.
      int j_Y = k_Y;
      if (inc_X == 1){
         for (int j = 0; j < n; j++){
            float temp = ZERO;
            for (int i = 0; i < m; i++){
               temp += A[i * n + j] * X[i]; //accessing A^T[j][i], i.e. A[i][j], assuming row major storage of A.
            }
            Y[j_Y] += alpha * temp;
            j_Y += inc_Y;
         }
      }else{
         for (int j = 0; j < n; j++){
            float temp = ZERO;
            int i_X = k_X;
            for (int i = 0; i < m; i++){
               temp += A[i * n + j] * X[i_X]; //accessing A^T[j][i], i.e. A[i][j], assuming row major storage of A.
               i_X += inc_X;
            }
            Y[j_Y] += alpha * temp;
            j_Y += inc_Y;
         }
      }
   }

   return;
}

void cblas_sgemv(const CBLAS_LAYOUT layout,
				 const CBLAS_TRANSPOSE TransA, const CBLAS_INT M, const CBLAS_INT N, 
      			 const float alpha, const float *A , const CBLAS_INT lda, 
      			 const float *X, const CBLAS_INT incX, const float beta, 
      			 float *Y, const CBLAS_INT incY){
	if(layout == CblasRowMajor) 
		soccs_sce_sgemv(TransA, M, N, alpha, A, lda, X, incX, beta, Y, incY);
	else {
		//soccs_sce_sgemv(TransA, N, M, alpha, A, lda, X, incX, beta, Y, incY);
		perror("calling cblas_sgemv with layout == CblaslColMajor.");
		exit(-1);
	}

}
