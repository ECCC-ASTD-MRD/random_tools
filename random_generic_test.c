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

#define NSAMPLES 1000000000

#include <math.h>
#include <randomfunctions.h>
#include <mpi.h>
#include <sys/time.h>
#include <sys/resource.h>

int my_main(int argc, char **argv){
  unsigned int lr;
  int i, j;
  double t0, t1, rval;
  double MPI_Wtime() ;
  unsigned int ranbuf[1200000];
  double ranbuf2[1200000];
  int pos, neg, mask ;
  long long postot, negtot, walk ;
#if defined(CYCLIC_TEST)
  int ran;
  long long count, counts ;
#endif
  double dmax, dmin;
  unsigned long long *idmax, *idmin ;
  unsigned int maxpos, maxneg;
  generic_state *gen = NULL;
  uint32_t mySeed = 0;
  int mpi_rank = 0;
  int ierr;

  MPI_Init(&argc,&argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
  if(ierr != MPI_SUCCESS) {
    MPI_Finalize();
    return 1;
  }
  for(i=0 ; i<1200000 ; i++) ranbuf[i] = 0;
  for(i=0 ; i<1200000 ; i++) ranbuf2[i] = 0.0;
  maxpos = 0x7FFFFFFF ;
  maxneg = 0x80000000 ;
  idmax = (unsigned long long *)&dmax;
  idmin = (unsigned long long *)&dmin;
  dmax = CVTDBL_32(maxpos) ;
  dmin = CVTDBL_32(maxneg) ;
  if(mpi_rank == 0) printf("maxpos, maxneg transformed with CVTDBL_32  : %22.18f %22.18f , %16.16Lx, %16.16Lx\n",dmax,dmin,*idmax,*idmin);
  dmax = CVTDBLS_32(maxpos) ;
  dmin = CVTDBLS_32(maxneg) ;
  if(mpi_rank == 0) printf("maxpos, maxneg transformed with CVTDBLS_32 : %22.18f %22.18f , %16.16Lx, %16.16Lx\n",dmax,dmin,*idmax,*idmin);
  if(mpi_rank == 0) printf("time per value (nanoseconds)\n");
#if defined(TEST_R250)
  mySeed = 123456;
  gen = (generic_state *)  Ran_R250_new_stream(NULL, &mySeed, 1);
#endif
#if defined(TEST_MT19937)
  gen = (generic_state *)  Ran_MT19937_new_stream(NULL, &mySeed, 1);
#endif
#if defined(TEST_XSR128)
  gen = (generic_state *)  Ran_XSR128_new_stream(NULL, &mySeed, 1);
#endif
#if defined(TEST_XSR128R)
  gen = (generic_state *)  Ran_XSR128R_new_stream(NULL, &mySeed, 1);
#endif
#if defined(TEST_SHR3)
  mySeed = 123456;
  gen = (generic_state *)  Ran_SHR3_new_stream(NULL, &mySeed, 1);
#endif
#if defined(CYCLIC_TEST)
  ran = IRan_generic_stream(gen);
  counts = 0;
  fif(mpi_rank == 0) printf(stdout,"ran target = %d\n",ran);
  for(j=0 ; j<1000 ; j++){
    count = 0;
    while(ran != IRan_generic_stream(gen)) count ++ ;
    counts += count;
    fif(mpi_rank == 0) printf(stdout,"%5d-repeat after %12Ld , running average = %12Ld\n",j+1,count,counts / (j+1)) ;
    fflush(stdout);
  }
exit(0);
#endif

  for( i=0 ; i < 1000000 ; i++) lr = IRan_generic_stream(gen);  // prime the pump

  if(mpi_rank == 0) printf("random stream will be reseeded with defaults before each timing test\n");
  MPI_Barrier(MPI_COMM_WORLD);
  RanSetSeed_generic_stream(gen,NULL,0);
  t0 = MPI_Wtime();
  for( i=0 ; i < NSAMPLES ; i++) lr = IRan_generic_stream(gen);
  t1 = MPI_Wtime();
  if(mpi_rank == 0) printf("time for 1E+9 x 1    random generic scalar        integer value = %6.3f , last = %12d\n",t1-t0, lr);

  for( i=0 ; i < 1000000 ; i++) rval = DRan_generic_stream(gen);
  MPI_Barrier(MPI_COMM_WORLD);
  RanSetSeed_generic_stream(gen,NULL,0);
  t0 = MPI_Wtime();
  for( i=0 ; i < NSAMPLES ; i++) rval = DRan_generic_stream(gen);
  t1 = MPI_Wtime();
  if(mpi_rank == 0) printf("time for 1E+9 x 1    random generic scalar         double value = %6.3f , last = %12g\n",t1-t0,rval);

  for( i=0 ; i < 1000000 ; i++) rval = DRanS_generic_stream(gen);
  MPI_Barrier(MPI_COMM_WORLD);
  RanSetSeed_generic_stream(gen,NULL,0);
  t0 = MPI_Wtime();
  for( i=0 ; i < NSAMPLES ; i++) rval = DRanS_generic_stream(gen);
  t1 = MPI_Wtime();
  if(mpi_rank == 0) printf("time for 1E+9 x 1    random generic scalar  signed double value = %6.3f , last = %12g\n",t1-t0,rval);

  for( i=0 ; i < 10 ; i++) VecIRan_generic_stream(gen,ranbuf, 1000000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  RanSetSeed_generic_stream(gen,NULL,0);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000 ; i++){
    VecIRan_generic_stream(gen, ranbuf, 1000000) ;
  }
  t1 = MPI_Wtime();
  if(mpi_rank == 0) printf("time for 1E+3 x 1E+6 random generic Vector       integer values = %6.3f , last = %12d ",t1-t0,ranbuf[1000000-1]);

  postot = 0 ; negtot = 0;
  RanSetSeed_generic_stream(gen,NULL,0);
  for (j=0 ; j<100 ; j++) {
    VecIRan_generic_stream(gen,ranbuf, 1000000) ;
    mask = 1 ;
    while (mask) {
      pos = 0 ; neg = 0 ; 
      for( i=0 ; i < 1000000 ; i++) if(ranbuf[i] & mask) pos++ ; else neg++  ; 
      postot += pos ; negtot += neg ;
      mask <<= 1 ;//  if(mpi_rank == 0) printf("%5d ",pos-neg) ;
    }
  }
  if(mpi_rank == 0) printf(", pos - neg = %lld (random walk = %lld)\n", postot-negtot, walk = sqrtf(1.0f*(postot+negtot)));

  for( i=0 ; i < 10 ; i++) VecDRan_generic_stream(gen, ranbuf2, 1000000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  RanSetSeed_generic_stream(gen,NULL,0);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000 ; i++){
    VecDRan_generic_stream(gen, ranbuf2, 1000000) ;
  }
  t1 = MPI_Wtime();
  if(mpi_rank == 0) printf("time for 1E+3 x 1E+6 random generic Vector        double values = %6.3f,  last = %12g\n",t1-t0,ranbuf2[1000000-1]);

  for( i=0 ; i < 10 ; i++) VecDRanS_generic_stream(gen, ranbuf2, 1000000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  RanSetSeed_generic_stream(gen,NULL,0);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000 ; i++){
    VecDRanS_generic_stream(gen, ranbuf2, 1000000) ;
  }
  t1 = MPI_Wtime();
  if(mpi_rank == 0) printf("time for 1E+3 x 1E+6 random generic Vector signed double values = %6.3f,  last = %12g\n",t1-t0,ranbuf2[1000000-1]);

  for( i=0 ; i < 1000 ; i++) VecIRan_generic_stream(gen, &ranbuf[i], 1000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  RanSetSeed_generic_stream(gen,NULL,0);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000 ; i++){
    VecIRan_generic_stream(gen, &ranbuf[i], 1000) ;
  }
  t1 = MPI_Wtime();
  if(mpi_rank == 0) printf("time for 1E+6 x 1E+3 random generic Vector       integer values = %6.3f , last = %12d\n",t1-t0,ranbuf[i+1000-2]);

  for( i=0 ; i < 1000 ; i++) VecDRan_generic_stream(gen, &ranbuf2[i], 1000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  RanSetSeed_generic_stream(gen,NULL,0);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000 ; i++){
    VecDRan_generic_stream(gen, &ranbuf2[i], 1000) ;
  }
  t1 = MPI_Wtime();
  if(mpi_rank == 0) printf("time for 1E+6 x 1E+3 random generic Vector        double values = %6.3f,  last = %12g\n",t1-t0,ranbuf2[i+1000-2]);

  for( i=0 ; i < 1000 ; i++) VecDRanS_generic_stream(gen, &ranbuf2[i], 1000) ;
  MPI_Barrier(MPI_COMM_WORLD);
  RanSetSeed_generic_stream(gen,NULL,0);
  t0 = MPI_Wtime();
  for( i=0 ; i < 1000000 ; i++){
    VecDRanS_generic_stream(gen, &ranbuf2[i], 1000) ;
  }
  t1 = MPI_Wtime();
  if(mpi_rank == 0) printf("time for 1E+6 x 1E+3 random generic Vector signed double values = %6.3f,  last = %12g\n",t1-t0,ranbuf2[i+1000-2]);

  MPI_Finalize();
  return(0);
}

int main(int argc, char **argv){
  struct rlimit rlim ;

  getrlimit(RLIMIT_STACK, &rlim) ;
  printf("Stack limit size soft = %ld, hard = %ld\n", rlim.rlim_cur, rlim.rlim_max);
  if(rlim.rlim_cur < 128000000 && rlim.rlim_cur != -1) 
    rlim.rlim_cur = (rlim.rlim_max < 128000000) ? rlim.rlim_max : 128000000 ;
  printf("Stack limit size now soft = %ld, hard = %ld\n", rlim.rlim_cur, rlim.rlim_max);
  setrlimit(RLIMIT_STACK, &rlim) ;
  my_main(argc, argv) ;
}
