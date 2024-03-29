/* 
 * Copyright (C) 2019 Recherche en Prevision Numerique
 * 
 * This is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this software; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "randomfunctions.h"
#include <mpi.h>
#if defined(PROFILE)
#define INSTRUMENT(A) A
#else
#define INSTRUMENT(A)
#endif
INSTRUMENT(extern unsigned int funcalls;)
INSTRUMENT(extern unsigned int funquick;)
INSTRUMENT(extern unsigned int funtails;)
INSTRUMENT(extern unsigned int funwedge;)
INSTRUMENT(extern unsigned int funloops;)
INSTRUMENT(extern unsigned int funused;)

int main(int argc, char **argv){
//   unsigned int lr;
  int i, j;
  double t0, t1, rval;
  double MPI_Wtime() ;
//   unsigned int ranbuf[1200000];
//   double ranbuf2[1200000];
//   int pos, neg, mask, postot, negtot ;
  double dmax, dmin, avg;
  unsigned long long *idmax, *idmin ;
  unsigned int maxpos, maxneg;
//   int gaussdist[10];
  int biggaussdist[2001];
  int index;
  generic_state *stream;
  uint32_t myseed = 123456;
  double x, p1, p2, prob, ptot;
  int est[2001];
  int mpi_rank = 0;
  int ierr;
  union{
    long l;
    double d;
  }v;
// check constants
  MPI_Init(&argc,&argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
  if(ierr != MPI_SUCCESS) {
    MPI_Finalize();
    return 1;
  }
  v.d = INVM31 ; if(mpi_rank == 0)  printf("%5i %16.16lx %24.20g\n",31,v.l,v.d);
  v.d = INVM32 ; if(mpi_rank == 0)  printf("%5i %16.16lx %24.20g\n",32,v.l,v.d);
  v.d = INVM33 ; if(mpi_rank == 0)  printf("%5i %16.16lx %24.20g\n",33,v.l,v.d);
  v.d = INVM63 ; if(mpi_rank == 0)  printf("%5i %16.16lx %24.20g\n",63,v.l,v.d);
  v.d = INVM64 ; if(mpi_rank == 0)  printf("%5i %16.16lx %24.20g\n",64,v.l,v.d);
  v.d = INVM65 ; if(mpi_rank == 0)  printf("%5i %16.16lx %24.20g\n",65,v.l,v.d);
//   void  RanNormalZigSetSeed(void *stream, int *piSeed, int cSeed)  ;
//   for(i=0 ; i<1200000 ; i++) ranbuf[i] = 0;
//   for(i=0 ; i<1200000 ; i++) ranbuf2[i] = 0.0;
  maxpos = 0x7FFFFFFF ;
  maxneg = 0x80000000 ;
  idmax = (unsigned long long *)&dmax;
  idmin = (unsigned long long *)&dmin;
  dmax = CVTDBL_32(maxpos) ;
  dmin = CVTDBL_32(maxneg) ;
  if(mpi_rank == 0)  printf("maxpos, maxneg transformed with CVTDBL_32  : %22.18f %22.18f , %16.16Lx, %16.16Lx\n",dmax,dmin,*idmax,*idmin);
  dmax = CVTDBLS_32(maxpos) ;
  dmin = CVTDBLS_32(maxneg) ;
  if(mpi_rank == 0)  printf("maxpos, maxneg transformed with CVTDBLS_32 : %22.18f %22.18f , %16.16Lx, %16.16Lx\n",dmax,dmin,*idmax,*idmin);

  stream = Ran_R250_new_stream(NULL, &myseed , 1);
  RanSetSeed_gaussian_stream(stream, &myseed, 1);
//   RanNormalZigSetSeed128(stream, &myseed, 1);
//   RanNormalZigSetSeed256(stream, &myseed, 1);

  dmin = 0.0 ; dmax = 0.0;
//   for( i=0 ; i < 10 ; i++) gaussdist[i] = 0;
  for( i=0 ; i < 2001 ; i++) biggaussdist[i] = 0;
  for(j=0; j<10 ; j++) ;
  avg = 0.0;
  for( i=0 ; i < 1000000000 ; i++) {
#if defined(TEST64)
//     rval = D64Ran_gaussian_stream(stream);      // use C entry point
    rval = F_D64Ran_gaussian_stream((statep *) &stream);  // use Fortran entry point
#else
//     rval = DRan_gaussian_stream(stream);        // use C entry point
    rval = F_DRan_gaussian_stream((statep *) &stream);    // use Fortran entry point
#endif
    avg = avg + rval ;
    dmin = (dmin < rval) ? dmin : rval ;
    dmax = (dmax > rval) ? dmax : rval ;
    if(rval > 10.0) rval = 10.0;
    if(rval < -10.0) rval = -10.0;
    index = 1001 + rval * 100;
    biggaussdist[index] ++;
  }
  if(mpi_rank == 0)  printf("for %d samples, min = %6.3f, max = %6.3f, avg = %10.7f\n",i,dmin,dmax,avg/i);
  ptot = 0.0;
  for( i=0 ; i < 2000 ; i++){
    x = (i-1001) / 100.0 ;
    p1 = exp(-0.5 * x * x)  * .01 ;
    x += .01;
    p2 =  exp(-0.5 * x * x) * .01 ;
    prob = .5 * (p1 + p2);
    ptot += prob ;
    p1 = 1000000000 * prob / 2.5066283;  // sqrt(2 * pi)
    est[i] = p1 * 1;
  }
  if(mpi_rank == 0)  printf("ptot = %g, expecting 2.5066283\n",ptot);
//   for( i=0 ; i < 10 ; i++) if(mpi_rank == 0)  printf("%9d ",gaussdist[i]);
//   for( i=0 ; i < 2000 ; i++) if(mpi_rank == 0)  printf("%9d %9d\n",i,biggaussdist[i]);
  if(mpi_rank == 0)  printf("%9s %9s %10s\n","slot","population","deviation(ppm) (center at slot 1000)");
  for( i=991 ; i <= 1010 ; i++) if(mpi_rank == 0)  printf("%9d %9d %10.0f\n",i,biggaussdist[i],1000000.0*(biggaussdist[i]-est[i])*1.0/est[i]);
  if(mpi_rank == 0)  printf("\n");
  MPI_Barrier(MPI_COMM_WORLD);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000000 ; i++) 
#if defined(TEST64)
    rval = F_D64Ran_gaussian_stream((statep *) &stream);  // time Fortran entry point (costlier)
#else
    rval = F_DRan_gaussian_stream((statep *) &stream);    // time Fortran entry point (costlier)
#endif
  t1 = MPI_Wtime();
  if(mpi_rank == 0)  printf("time for 1E+9 x 1 random DRan_gaussian_stream/R250 double value = %6.3f sec\n",t1-t0);  // DRan_gaussian_stream256

  t1 = 0 ; t0 = 1 ; 
  INSTRUMENT(t1 = funquick ; t0 = funcalls+1 ; )
  INSTRUMENT(if(mpi_rank == 0)  printf("%% quick calls in gaussian generator = %7.3f\n",t1 / t0 * 100.0);)
  INSTRUMENT(t1 = funtails ;)
  INSTRUMENT(if(mpi_rank == 0)  printf("%% tail calls in gaussian generator  = %7.3f\n",t1 / t0 * 100.0);)
  INSTRUMENT(t1 = funwedge ;)
  INSTRUMENT(if(mpi_rank == 0)  printf("%% wedge calls in gaussian generator = %7.3f\n",t1 / t0 * 100.0);)
  INSTRUMENT(t1 = funloops ;)
  INSTRUMENT(if(mpi_rank == 0)  printf("%% extra loops in gaussian generator = %7.3f\n",t1 / t0 * 100.0);)
  INSTRUMENT(t1 = funused ;)
  INSTRUMENT(if(mpi_rank == 0)  printf("uniform random values used per call = %7.3f\n",t1 / t0);)
  MPI_Finalize();
  return(0);
}
