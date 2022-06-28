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

void soccs_sce_sgemv(&TA, &M, &N, &alpha, A, &lda, X, &incX,
                &beta, Y, &incY){
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
      soccs_sce_sgemv(&TA, &M, &N, &alpha, A, &lda, X, &incX,
                &beta, Y, &incY);
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
      soccs_sce_sgemv(&TA, &M, &N, &alpha, A, &lda, X, &incX,
                &beta, Y, &incY);
   }
   else cblas_xerbla(1, "cblas_sgemv", "Illegal layout setting, %d\n", layout);
   CBLAS_CallFromC = 0;
   RowMajorStrg = 0;
   return;
}

