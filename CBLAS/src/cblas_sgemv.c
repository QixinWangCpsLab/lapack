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

   if (lsame(trans, 'N')){
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
   if (lsame(trans, 'N')){
      //$Y := \alpha A X + Y$.
      int j_X = k_X;
      if (inc_Y == 1){
         for (int j = 0; j < n; j++){
            float temp = alpha * X[j_X];
            for (int i = 1; i < m; i++){
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
         for (int j = 1; j < n; j++){
            float temp = ZERO;
            for (int i = 0; i < m; i++){
               temp += A[i * n + j] * X[i]; //accessing A[i][j], assuming row major storage of A.
            }
            Y[j_Y] += alpha * temp;
            j_Y += inc_Y;
         }
      }else{
         for (int j = 1; j < n; j++){
            float temp = ZERO;
            int i_X = k_X;
            for (int i = 0; i < m; i++){
               temp += A[i * n + j] * X[i_X]; //accessing A[i][j], assuming row major storage of A.
               i_X += inc_X;
            }
            Y[j_Y] += alpha * temp;
            j_Y += inc_Y;
         }
      }
   }

   return;
}
