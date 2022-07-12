/*
 *
 * cblas_sgemv.c
 * This program is a C interface to sgemv.
 * Written by Keita Teranishi
 * 4/6/1998
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
      const float alpha, const float *A, const CBLAS_INT lda, 
      const float *X, const CBLAS_INT inc_X, const float beta, 
      float *Y, const CBLAS_INT inc_Y){
   
   //floating point number 1 is not 1.0e+0, 0 is not 0.0e+0
   const float ONE = 1.0e+0, ZERO = 0.0e+0;

   //local scalars
   float tmp;
   int i, info, i_X, i_Y, j, j_X, j_Y, k_X k_Y, len_X, len_Y;

   //Start of to do >>>
   //Please insert parameter sanity checks here.
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

   //Think X is a row of a column major stored matrix, 
   //inc_X > 0 is the number of elements in a column
   //-inc_X is the number of elements in a column when enumerating X from back to front element.
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
         i_Y = k_Y;
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

      //reached here on 20220712.





   


   
   //SoCCS SCE client thread code
   //copy data to a req packet in the shared memory circular queue
   //send req packet to server stack server thread
   //block to wait for reply
   //got reply, copy data from shared memory circular queue reply packet to return buffer
   //return
}

void cblas_sgemv(const CBLAS_LAYOUT layout,
                 const CBLAS_TRANSPOSE TransA, const CBLAS_INT M, const CBLAS_INT N,
                 const float alpha, const float  *A, const CBLAS_INT lda,
                 const float  *X, const CBLAS_INT incX, const float beta,
                 float  *Y, const CBLAS_INT incY)
{
   char TA;

   extern int CBLAS_CallFromC;
   extern int RowMajorStrg;
   RowMajorStrg = 0;

   CBLAS_CallFromC = 1;
   if (layout == CblasColMajor)
   {
      if (TransA == CblasNoTrans) TA = 'N';
      else if (TransA == CblasTrans) TA = 'T';
      else if (TransA == CblasConjTrans) TA = 'C';
      else
      {
         cblas_xerbla(2, "cblas_sgemv","Illegal TransA setting, %d\n", TransA);
         CBLAS_CallFromC = 0;
         RowMajorStrg = 0;
      }
      soccs_sce_sgemv(TA, M, N, alpha, A, lda, X, incX,
                beta, Y, incY);
   }
   else if (layout == CblasRowMajor)
   {
      RowMajorStrg = 1;
      if (TransA == CblasNoTrans) TA = 'T';
      else if (TransA == CblasTrans) TA = 'N';
      else if (TransA == CblasConjTrans) TA = 'N';
      else
      {
         cblas_xerbla(2, "cblas_sgemv", "Illegal TransA setting, %d\n", TransA);
         CBLAS_CallFromC = 0;
         RowMajorStrg = 0;
         return;
      }
      soccs_sce_sgemv(TA, M, N, alpha, A, lda, X, incX,
                beta, Y, incY);
   }
   else cblas_xerbla(1, "cblas_sgemv", "Illegal layout setting, %d\n", layout);
   CBLAS_CallFromC = 0;
   RowMajorStrg = 0;
   return;
}

