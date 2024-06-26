/* SVE utility Library */
/* written by Victor Soria, BSC; Viktor K. Decyk, UCLA and Ricardo Fonseca, ISCTE */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>

/*--------------------------------------------------------------------*/
void sve_fallocate(float **s_f, int nsize, int *irc) {
/* allocate aligned float memory on SVE return pointer to C */
/* size is padded to be a multiple of the alignment length */
/* local data */
/* NV = vector length for 32 bit data */
#define NV             4
   int ns;
   void *sptr = NULL;
   ns = NV*((nsize - 1)/NV + 1);
   posix_memalign(&sptr, 4*NV, ns*sizeof(float));
   if (sptr==NULL) {
      printf("_mm_malloc float Error,len=%d\n",ns);
      *irc = 1;
   }
   *s_f = (float *)sptr;
   return;
#undef NV
}

/*--------------------------------------------------------------------*/
void sve_callocate(float complex **s_c, int nsize, int *irc) {
/* allocate aligned float complex memory on SVE return pointer to C */
/* size is padded to be a multiple of the alignment length          */
/* local data */
/* NV = vector length for 64 bit data */
#define NV             2
   int ns;
   void *sptr = NULL;
   ns = NV*((nsize - 1)/NV + 1);
   posix_memalign(&sptr,8*NV,ns*sizeof(float complex));
   if (sptr==NULL) {
      printf("_mm_malloc float complex Error,len=%d\n",ns);
      *irc = 1;
   }
   *s_c = (float complex *)sptr;
   return;
#undef NV
}

/*--------------------------------------------------------------------*/
void sve_iallocate(int **s_i, int nsize, int *irc) {
/* allocate aligned int memory on SVE, return pointer to C */
/* size is padded to be a multiple of the alignment length */
/* local data */
/* NV = vector length for 32 bit data */
#define NV             4
   int ns;
   void *sptr = NULL;
   ns = NV*((nsize - 1)/NV + 1);
   posix_memalign(&sptr,4*NV,ns*sizeof(int));
   if (sptr==NULL) {
      printf("_mm_malloc int Error,len=%d\n",ns);
      *irc = 1;
   }
   *s_i = (int *)sptr;
   return;
#undef NV
}

/*--------------------------------------------------------------------*/
void sve_deallocate(void *s_d) {
/* deallocate aligned memory on SVE */
   free(s_d);
   return;
}

/*--------------------------------------------------------------------*/
void csveiscan2(int *isdata, int nths) {
/* performs local prefix reduction of integer data shared by threads */
/* using binary tree method. */
/* requires SVE, isdata needs to be 16 byte aligned */
/* local data */
/* int j, ns, isum, ist;
   __m128i v_m1, v_m2, v_it, v_is, v_ioff;
   ns = 4*(nths/4);
   v_m1 = _mm_set_epi32(0,-1,0,-1);
   v_m2 = _mm_set_epi32(0,-1,-1,0);
   isum = 0;
   v_ioff = _mm_set1_epi32(isum);*/
/* vector loop over elements in blocks of 4 */
//   for (j = 0; j < ns; j+=4) {
/* load data */
//      v_it = _mm_load_si128((__m128i *)&isdata[j]);
/* first pass */
//      v_is = _mm_slli_si128(_mm_and_si128(v_it,v_m1),4);
//      v_it = _mm_add_epi32(v_is,v_it);
/* second pass */
//      v_is = _mm_shuffle_epi32(v_it,212);
//      v_is = _mm_slli_si128(_mm_and_si128(v_is,v_m2),4);
//      v_it = _mm_add_epi32(v_is,v_it);
/* add offset */
//      v_it = _mm_add_epi32(v_it,v_ioff);
/* next offset */
//      v_ioff = _mm_shuffle_epi32(v_it,255);
/* write data */
/*      _mm_store_si128((__m128i *)&isdata[j],v_it);
   }
   if (ns > 0)
      isum = isdata[ns-1];*/
/* loop over remaining elements */
/* for (j = ns; j < nths; j++) {
      ist = isdata[j];
      isum += ist;
      isdata[j] = isum;
   }
   return;*/
}

/* Interfaces to Fortran */

/*--------------------------------------------------------------------*/
void sve_deallocate_(void *sp_d) {
/* pointer in Fortran should also be nullified */
   sve_deallocate(sp_d);
   return;
}

/*--------------------------------------------------------------------*/
void csveiscan2_(int *isdata, int *nths) {
   csveiscan2(isdata,*nths);
   return;
}

void fcopyin_(float *f, float *g, int *n) {
   int j;
   for (j = 0; j < *n; j++) {
      f[j] = g[j];
   }
   return;
}
