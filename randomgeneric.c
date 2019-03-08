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

// generic interface to random functions using streams  // !InTc!

#include <randomgeneric.h>

// Fortran interfaces for automated extraction
#if defined(NEVER_TO_BE_TRUE)

  type, bind(C) :: RANDOM_STREAM                                                          !InTf!
    type(C_PTR) :: p                                                                      !InTf!
  end type                                                                                !InTf!

! void F_RanSetSeed_generic_stream(statep *s   , int *piSeed, int cSeed)                  !InTf!
 interface                                                                                !InTf!
   subroutine RanSetSeed_generic_stream(stream, piSeed, cSeed) bind(C,name='F_RanSetSeed_generic_stream')  !InTf!
   import :: RANDOM_STREAM,C_INT                                                          !InTf!
   type(RANDOM_STREAM), intent(IN) :: stream                                              !InTf!
   integer(C_INT), intent(IN), value :: cSeed                                             !InTf!
   integer(c_INT), dimension(cSeed), intent(IN) :: piSeed                                 !InTf!
   end subroutine RanSetSeed_generic_stream                                               !InTf!
 end interface                                                                            !InTf!

! void F_RanSetSeed_gaussian_stream(statep *s   , int *piSeed, int cSeed)                  !InTf!
 interface                                                                                !InTf!
   subroutine RanSetSeed_gaussian_stream(stream, piSeed, cSeed) bind(C,name='F_RanSetSeed_gaussian_stream')  !InTf!
   import :: RANDOM_STREAM,C_INT                                                          !InTf!
   type(RANDOM_STREAM), intent(IN) :: stream                                              !InTf!
   integer(C_INT), intent(IN), value :: cSeed                                             !InTf!
   integer(c_INT), dimension(cSeed), intent(IN) :: piSeed                                 !InTf!
   end subroutine RanSetSeed_gaussian_stream                                              !InTf!
 end interface                                                                            !InTf!

! unsigned int F_IRan_generic_stream(statep *s   )                                        !InTf!
 interface                                                                                !InTf!
   function IRan_generic_stream(stream) result(ran) bind(C,name='F_IRan_generic_stream')  !InTf!
   import :: C_INT,RANDOM_STREAM                                                          !InTf!
   type(RANDOM_STREAM), intent(IN) :: stream                                              !InTf!
   integer(C_INT) :: ran                                                                  !InTf!
   end function IRan_generic_stream                                                       !InTf!
 end interface                                                                            !InTf!

! double F_DRan_generic_stream(statep *s   )                                              !InTf!
 interface                                                                                !InTf!
   function DRan_generic_stream(stream) result(ran) bind(C,name='F_DRan_generic_stream')  !InTf!
   import :: C_DOUBLE,RANDOM_STREAM                                                       !InTf!
   type(RANDOM_STREAM), intent(IN) :: stream                                              !InTf!
   real(C_DOUBLE) :: ran                                                                  !InTf!
   end function DRan_generic_stream                                                       !InTf!
 end interface                                                                            !InTf!

! double F_DRan_gaussian_stream(statep *s   )                                             !InTf!
 interface                                                                                !InTf!
   function DRan_gaussian_stream(stream) result(ran) bind(C,name='F_DRan_gaussian_stream')  !InTf!
   import :: C_DOUBLE,RANDOM_STREAM                                                       !InTf!
   type(RANDOM_STREAM), intent(IN) :: stream                                              !InTf!
   real(C_DOUBLE) :: ran                                                                  !InTf!
   end function DRan_gaussian_stream                                                      !InTf!
 end interface                                                                            !InTf!

! double F_D64Ran_gaussian_stream(statep *s   )                                              !InTf!
 interface                                                                                !InTf!
   function D64Ran_gaussian_stream(stream) result(ran) bind(C,name='F_D64Ran_gaussian_stream')  !InTf!
   import :: C_DOUBLE,RANDOM_STREAM                                                       !InTf!
   type(RANDOM_STREAM), intent(IN) :: stream                                              !InTf!
   real(C_DOUBLE) :: ran                                                                  !InTf!
   end function D64Ran_gaussian_stream                                                    !InTf!
 end interface                                                                            !InTf!

! double F_DRanS_generic_stream(statep *s   )                                             !InTf!
 interface                                                                                !InTf!
   function DRanS_generic_stream(stream) result(ran) bind(C,name='F_DRanS_generic_stream')  !InTf!
   import :: C_DOUBLE,RANDOM_STREAM                                                       !InTf!
   type(RANDOM_STREAM), intent(IN) :: stream                                              !InTf!
   real(C_DOUBLE) :: ran                                                                  !InTf!
   end function DRanS_generic_stream                                                      !InTf!
 end interface                                                                            !InTf!

! void F_VecIRan_generic_stream(statep *s   , unsigned int *ranbuf, int n)                !InTf!
 interface                                                                                !InTf!
   subroutine VecIRan_generic_stream(stream, ranbuf, n) bind(C,name='F_VecIRan_generic_stream') !InTf!
   import :: RANDOM_STREAM,C_INT                                                          !InTf!
   type(RANDOM_STREAM), intent(IN) :: stream                                              !InTf!
   integer(C_INT), intent(IN), value :: n                                                 !InTf!
   integer(c_INT), dimension(n), intent(OUT) :: ranbuf                                    !InTf!
   end subroutine VecIRan_generic_stream                                                  !InTf!
 end interface                                                                            !InTf!

! void F_VecDRanS_generic_stream(statep *s   , double *ranbuf, int n)                     !InTf!
 interface                                                                                !InTf!
   subroutine VecDRanS_generic_stream(stream, ranbuf, n) bind(C,name='F_VecDRanS_generic_stream') !InTf!
   import :: RANDOM_STREAM,C_INT,C_DOUBLE                                                 !InTf!
   type(RANDOM_STREAM), intent(IN) :: stream                                              !InTf!
   integer(C_INT), intent(IN), value :: n                                                 !InTf!
   real(C_DOUBLE), dimension(n), intent(OUT) :: ranbuf                                    !InTf!
   end subroutine VecDRanS_generic_stream                                                 !InTf!
 end interface                                                                            !InTf!

! void F_VecDRan_generic_stream(statep *s   , double *ranbuf, int n)                      !InTf!
 interface                                                                                !InTf!
   subroutine VecDRan_generic_stream(stream, ranbuf, n) bind(C,name='F_VecDRan_generic_stream') !InTf!
   import :: RANDOM_STREAM,C_INT,C_DOUBLE                                                 !InTf!
   type(RANDOM_STREAM), intent(IN) :: stream                                              !InTf!
   integer(C_INT), intent(IN), value :: n                                                 !InTf!
   real(C_DOUBLE), dimension(n), intent(OUT) :: ranbuf                                    !InTf!
   end subroutine VecDRan_generic_stream                                                  !InTf!
 end interface                                                                            !InTf!
#endif
#if defined(NEVER_TO_BE_TRUE)
//****f* librandom/Ran_MT19937_new_stream
// Synopsis
//    create a random generator "stream" of type MT19937
//
// ARGUMENTS
generic_state *Ran_MT19937_new_stream(void *clone, unsigned int *piSeed, int cSeed)
// subroutine Ran_MT19937_new_stream(stream, clone, piSeed, cSeed)
// INPUTS
//    clone      normally NULL (or pointer to same type stream to duplicate for testing)
//    piSeed     seed value array
//    cSeed      number of values in seed array (0 means default seeding)
// OUTPUTS
//    stream     pointer to linear stream of type MT19937
//****
//****f* librandom/Ran_R250_new_stream
// Synopsis
//    create a random generator "stream" of type R250
//
// ARGUMENTS
generic_state *Ran_R250_new_stream(void *clone, unsigned int *piSeed, int cSeed)
// subroutine Ran_R250_new_stream(stream, clone, piSeed, cSeed)
// INPUTS
//    clone      normally NULL (or pointer to same type stream to duplicate for testing)
//    piSeed     seed value array
//    cSeed      number of values in seed array (0 means default seeding)
// OUTPUTS
//    stream     pointer to linear stream of type R250
//****
//****f* librandom/Ran_SHR3_new_stream
// Synopsis
//    create a random generator "stream" of type SHR3
//
// ARGUMENTS
generic_state *Ran_SHR3_new_stream(void *clone, uint32_t *piSeed, int cSeed)
// subroutine Ran_SHR3_new_stream(stream, clone, piSeed, cSeed)
// INPUTS
//    clone      normally NULL (or pointer to same type stream to duplicate for testing)
//    piSeed     seed value array
//    cSeed      number of values in seed array (0 means default seeding)
// OUTPUTS
//    stream     pointer to linear stream of type SHR3
//****
//****f* librandom/Ran_XSR128_new_stream
// Synopsis
//    create a random generator "stream" of type XSR128
//
// ARGUMENTS
generic_state *Ran_XSR128_new_stream(void *clone, unsigned int *piSeed, int cSeed)
// subroutine Ran_XSR128_new_stream(stream, clone, piSeed, cSeed)
// INPUTS
//    clone      normally NULL (or pointer to same type stream to duplicate for testing)
//    piSeed     seed value array
//    cSeed      number of values in seed array (0 means default seeding)
// OUTPUTS
//    stream     pointer to linear stream of type XSR128
//****
//****f* librandom/Ran_XSR128R_new_stream
// Synopsis
//    create a random generator "stream" of type XSR128R
//
// ARGUMENTS
generic_state *Ran_XSR128R_new_stream(void *clone, unsigned int *piSeed, int cSeed)
// subroutine Ran_XSR128R_new_stream(stream, clone, piSeed, cSeed)
// INPUTS
//    clone      normally NULL (or pointer to same type stream to duplicate for testing)
//    piSeed     seed value array
//    cSeed      number of values in seed array (0 means default seeding)
// OUTPUTS
//    stream     pointer to linear stream of type XSR128R
//****
#endif
//****P* librandom/random value generators
// Synopsis
// generic access to multiple families of random number generators (NOT for cryptographic usage)
//
// the user will create a random number "stream" of a given type (see available generators)
// and call various functions/routines to produce random number sequences (scalar or vector).
// the stream generators are believed to be "thread safe" once they are initialized
//
// source code available at:
//   https://gitlab.com/mfvalin/random_tools
// 
// scalar values :
//   unsigned integer value      0 <=  value  < 2**32 - 1
//   unsigned double value     0.0 <   value  < 1.0
//   signed double value      -1.0 <   value  < 1.0
// vector of values : same 3 kinds as scalar values
//
// available generators  
// generic (linear) generators :
//    R250     (shift register sequence)     https://www.taygeta.com/rwalks/node2.html
//    MT19937  (Mersenne twister)            https://en.wikipedia.org/wiki/Mersenne_Twister
//    SHR3     (3-shift-register)
//    XSR128   (xor-shift)                   https://en.wikipedia.org/wiki/Xorshift
//    XSR128R  (xor-shift-rotate)            https://en.wikipedia.org/wiki/Xoroshiro128%2B
// gaussian (normal) generators :
//    NormalZig (ziggurat with 256 bins)     https://en.wikipedia.org/wiki/Ziggurat_algorithm
//
// performance:
//    the R250 linear generator is the fastest of all (fully vectorizable)
//
// this package contains C and Fortran functions
//     C function                                 equivalent Fortran function 
// RanSetSeed_generic_stream        subroutine RanSetSeed_generic_stream(stream, piSeed, cSeed)
// IRan_generic_stream              function   IRan_generic_stream(stream)  result(Iran)
// DRan_generic_stream              function   DRan_generic_stream(stream)  result(Dran)
// DRanS_generic_stream             function   DRanS_generic_stream(stream) result(Dran)
// VecIRan_generic_stream           subroutine VecIRan_generic_stream(stream , Ibuf, n)
// VecDRan_generic_stream           subroutine VecDRan_generic_stream(stream , Dbuf, n)
// VecDRanS_generic_stream          subroutine VecDRanS_generic_stream(stream, Dbuf, n)
// RanSetSeed_gaussian_stream       subroutine RanSetSeed_gaussian_stream(stream, piSeed, cSeed)
// DRan_gaussian_stream             function   DRan_gaussian_stream(stream) result(Dran)
// D64Ran_gaussian_stream           function   D64Ran_gaussian_stream(stream) result(Dran)
// 
// Stream creation functions
//     C function                                 equivalent Fortran function 
// Ran_R250_new_stream              subroutine Ran_R250_new_stream(newstream, clone, piSeed, cSeed)
// Ran_SHR3_new_stream              subroutine Ran_SHR3_new_stream(newstream, clone, piSeed, cSeed)
// Ran_MT19937_new_stream           subroutine Ran_MT19937_new_stream(newstream, clone, piSeed, cSeed)
// Ran_XSR128_new_stream            subroutine Ran_XSR128_new_stream(newstream, clone, piSeed, cSeed)
// Ran_XSR128R_new_stream           subroutine Ran_XSR128R_new_stream(newstream, clone, piSeed, cSeed)
//
// C functions and arguments description
//    generic_state *Ran_MT19937_new_stream(void *clone_in, unsigned int *piSeed, int cSeed)
//    generic_state *Ran_R250_new_stream(void *clone_in, unsigned int *piSeed, int cSeed)
//    generic_state *Ran_SHR3_new_stream(void *clone_in, uint32_t *piSeed, int cSeed)
//    generic_state *Ran_XSR128_new_stream(void *clone_in, unsigned int *piSeed, int cSeed)
//    generic_state *Ran_XSR128R_new_stream(void *clone_in, unsigned int *piSeed, int cSeed)
//    generic_state RanSetSeed_generic_stream(generic_state *stream, unsigned int *piSeed, int cSeed)
//    uint32_t IRan_generic_stream(generic_state *stream)
//    double DRan_generic_stream(generic_state *stream)
//    double DRanS_generic_stream(generic_state *stream)
//    void VecIRan_generic_stream(generic_state *stream, unsigned int *ibuf, int n)
//    void VecDRan_generic_stream(generic_state *stream, double *dbuf, int n)
//    void VecDRanS_generic_stream(generic_state *stream, double *dbuf, int n)
//    double DRan_gaussian_stream(generic_state *stream)
//    double D64Ran_gaussian_stream(generic_state *stream)
//
//    clone_in  : normally NULL (used only for tests)
//    cSeed     : number of seed values (see below for specifics)
//    piSeed    : array of seed values
//    stream    : pointer to a stream (obtained from the stream creeation routines)
//    n         : number of values
//    ibuf      : pointer to array of integers
//    dbuf      : pointer to array of doubles
//
// Fortran arguments description
//    type(RANDOM_STREAM), intent(IN)              :: stream
//    type(RANDOM_STREAM), intent(IN)              :: clone    ! normally RANDOM_STREAM(C_NULL_PTR)
//    type(RANDOM_STREAM), intent(OUT)             :: newstream
//    integer(c_INT), dimension(cSeed), intent(IN) :: piSeed   ! array of seed values
//    integer(C_INT), intent(IN), value            :: cSeed    ! number of seed values
//    integer(C_INT), intent(IN), value            :: n
//    integer(C_INT)                               :: Iran
//    real(C_DOUBLE)                               :: Dran
//    integer(C_INT), dimension(n), intent(OUT)    :: Ibuf
//    real(C_DOUBLE), dimension(n), intent(OUT)    :: Dbuf
//
// general notes:
//    before a stream can be used, it must be created by one of the stream creation functions
//    there is a specific creator for each type of stream
//
//    for the stream creation functions, cSeed = 0 forces a stream variety specific default seeding
//    cSeed < 0 is an ERROR and may produce UNEXPECTED RESULTS depending upon the stream variety
//
//    R250    accepts cSeed = 1 or cSeed = 250 (piSeed must contain cSeed positive values)
//    SHR3    accepts cSeed > 0  (one value will be used from piSeed)
//    MT19937 accepts cSeed > 0  (one value will be used from piSeed)
//    XSR128  accepts cSeed >= 4 (4 values will be used from piSeed)
//    XSR128R accepts cSeed >= 4 (4 values will be used from piSeed)
//
// C specific note:
//    to compile and load successfully, one must get the interface definitions from
//
//    #include <randomfunctions.h>
//
// Fortran specific note:
//    to compile and load successfully, one must get the interface definitions from
//
//    use ISO_C_BINDING
//    include 'randomfunctions.inc'
//
// ==================== a small Fortran demo program ====================
//
// program demo
//   use ISO_C_BINDING
//   implicit none
//   include 'randomfunctions.inc'
//   integer, parameter :: NI = 10
//   type(RANDOM_STREAM) :: s, clone
//   integer, dimension(NI) :: ibuf
//   integer :: iran, cSeed
//   integer, dimension(1) :: piSeed
//   real*8,  dimension(NI) :: dbuf, dsbuf
//   real *8 :: dran, dsran
//   cSeed = 0                         ! default initialization
//   clone = RANDOM_STREAM(C_NULL_PTR) ! null clone, not duplicating existing stream
//   call Ran_R250_new_stream(s, clone, piSeed, cSeed)   ! create R250 stream
//   iran  = IRan_generic_stream(s)                      ! get 1 integer value
//   dran  = DRan_generic_stream(s)                      ! get 1 real*8 value
//   dsran = DRanS_generic_stream(s)                     ! get 1 real*8 value
//   print *,'scalar ',iran, dran, dsran
//   call VecIRan_generic_stream(s, ibuf, NI)            ! get NI integer values
//   call VecDRan_generic_stream(s, dbuf, NI)            ! get NI real*8 values
//   call VecDRanS_generic_stream(s, dsbuf, NI)          ! get NI real*8 values
//   print *,'vector1',ibuf(1),dbuf(1),dsbuf(1)
//   print *,'vector2',ibuf(NI),dbuf(NI),dsbuf(NI)
//   call RanSetSeed_gaussian_stream(s,piSeed, cSeed)
//   dran = DRan_gaussian_stream(s)
//   print *,'gaussian ',dran
//   stop
// end
//
// ==================== a small C demo program ==================== 
//
// #include <stdio.h>
// #include <stdint.h>
// #include <randomfunctions.h>
// #define NI 10
// int main(){
//   uint32_t mySeed = 123456;
//   int cSeed = 0;
//   generic_state *s;
//   uint32_t iran, ibuf[NI];
//   double dran, dsran, dbuf[NI], dsbuf[NI];
//   s = (generic_state *)  Ran_R250_new_stream(NULL, &mySeed, cSeed);  // default seeding
//   iran  = IRan_generic_stream(s) ;                     // get 1 integer value
//   dran  = DRan_generic_stream(s) ;                     // get 1 double value
//   dsran = DRanS_generic_stream(s);                     // get 1 double value
//   printf("scalar  %10d %f %f\n",iran,dran,dsran);
//   VecIRan_generic_stream(s, ibuf, NI) ;                // get NI integer values
//   VecDRan_generic_stream(s, dbuf, NI) ;                // get NI double values
//   VecDRanS_generic_stream(s, dsbuf, NI) ;              // get NI double values
//   printf("vector1 %10d %f %f\n",ibuf[0],dbuf[0],dsbuf[0]);
//   printf("vector2 %10d %f %f\n",ibuf[NI-1],dbuf[NI-1],dsbuf[NI-1]);
//   RanSetSeed_gaussian_stream(s, &mySeed, cSeed);
//   dran = DRan_gaussian_stream(s);
//   printf("gaussian %f\n",dran);
//   return 0;
// }
//****

void RanNormalZigSetSeed(generic_state *stream, uint32_t *piSeed, int cSeed);
//****f* librandom/RanSetSeed_gaussian_stream
// Synopsis
//    reseed a random generator "stream" previously created by a call to
//    Ran_XXXXXX_new_stream (where XXXXXX is one of above mentioned linear generators)
//
// ARGUMENTS
void RanSetSeed_gaussian_stream(generic_state *stream, uint32_t *piSeed, int cSeed)  // !InTc!
// subroutine RanSetSeed_gaussian_stream(stream, piSeed, cSeed)
// INPUTS
//    stream     linear stream created by a Ran_XXXXXX_new_stream function
//    piSeed     not used  (for consistency with other reseeding functions)
//    cSeed      not used  (for consistency with other reseeding functions)
// OUTPUTS
//    none       the gaussian stream is initialized and its internal buffer is marked as empty
//****
{
  RanNormalZigSetSeed(stream, piSeed, cSeed);
}
void F_RanSetSeed_gaussian_stream(statep *s, uint32_t *piSeed, int cSeed)            // !InTc!
{
  RanNormalZigSetSeed(s->p, piSeed, cSeed);
}

double DRan_NormalZig_stream(generic_state *stream);
//****f* librandom/DRan_gaussian_stream
// Synopsis
//    generate a single 64bit positive floating point value, 
//    gaussian distribution with 32 significant bits in mantissa
// ARGUMENTS
double DRan_gaussian_stream(generic_state *stream)                              // !InTc!
// function DRan_gaussian_stream(stream) result(ran)
// INPUTS
//    stream     linear stream created by a Ran_XXXXXX_new_stream function
// OUTPUTS
//    64bit floating point number, gaussian distribution, 32 significant bits in mantissa
//****
{
  return DRan_NormalZig_stream(stream);
}
double F_DRan_gaussian_stream(statep *s)       // !InTc!
{
  return DRan_NormalZig_stream(s->p);
}

double D64Ran_NormalZig_stream(generic_state *stream);
//****f* librandom/D64Ran_gaussian_stream
// Synopsis
//    generate a single 64bit positive floating point value, 
//    gaussian distribution with 52 significant bits in mantissa
// ARGUMENTS
double D64Ran_gaussian_stream(generic_state *stream)                            // !InTc!
// function D64Ran_gaussian_stream(stream) result(ran)
// INPUTS
//    stream     linear stream created by a Ran_XXXXXX_new_stream function
// OUTPUTS
//    64bit floating point number, gaussian distribution, 52 significant bits in mantissa
//****
{
  return D64Ran_NormalZig_stream(stream);
}
double F_D64Ran_gaussian_stream(statep *s)       // !InTc!
{
  return D64Ran_NormalZig_stream(s->p);
}

//****f* librandom/RanSetSeed_generic_stream
// Synopsis
//    reseed a random generator "stream" previously created by a call to
//    Ran_XXXXXX_new_stream (where XXXXXX is one of above mentioned linear generators)
//
// ARGUMENTS
void RanSetSeed_generic_stream(generic_state *stream, unsigned int *piSeed, int cSeed)  // !InTc!
// subroutine RanSetSeed_generic_stream(stream, piSeed, cSeed)
// INPUTS
//    stream     linear stream created by a Ran_XXXXXX_new_stream function
//    piSeed     integer "seed" array  (NULL pointer means use default seed)
//    cSeed      size of "seed" array  (0 means use default seed)
// OUTPUTS
//    none       the stream is reseeded and its internal buffer is marked as empty
//****
{
  generic_state *state = stream ;
  state->seed(stream, piSeed, cSeed);
}
void F_RanSetSeed_generic_stream(statep *s, unsigned int *piSeed, int cSeed)  // Fortran interface using derived type
{
  RanSetSeed_generic_stream( s->p, piSeed, cSeed);
}

//****f* librandom/IRan_generic_stream
// Synopsis
//    generate a single unsigned 32bit integer value, according to the stream type
// ARGUMENTS
uint32_t IRan_generic_stream(generic_state *stream)       // !InTc!
// function IRan_generic_stream(stream) result(iran)
// INPUTS
//    stream     linear stream created by a Ran_XXXXXX_new_stream function
// OUTPUTS
//    unsigned 32bit integer
//****
{
//   generic_state *state = stream ;
//   return state->iran(stream);
  uint32_t value;
// fprintf(stderr,"IRan_generic_stream\n");
  if(stream->cur > stream->top) stream->refill(stream);
  value = stream->buf[stream->cur] ;
  stream->cur = stream->cur + 1;
  return value;
}
uint32_t F_IRan_generic_stream(statep *s)  // Fortran interface using derived type
{
  return(IRan_generic_stream(s->p));
}

//****f* librandom/DRan_generic_stream
// Synopsis
//    generate a single 64bit positive floating point value, according to the stream type
// ARGUMENTS
double DRan_generic_stream(generic_state *stream)       // !InTc!
// function DRan_generic_stream(stream) result(dran)
// INPUTS
//    stream     linear stream created by a Ran_XXXXXX_new_stream function
// OUTPUTS
//    64bit floating point number    0.0 < value < 1.0
//****
{
//   generic_state *state = stream ;
//   return state->dran(stream);
  uint32_t value;
  if(stream->cur > stream->top) stream->refill(stream);
  value = stream->buf[stream->cur] ;
  stream->cur = stream->cur + 1;
  return CVTDBL_32(value);
}
double F_DRan_generic_stream(statep *s)  // Fortran interface using derived type
{
  return(DRan_generic_stream(s->p));
}

//****f* librandom/DRanS_generic_stream
// Synopsis
//    generate a single 64bit signed floating point value, according to the stream type
// ARGUMENTS
double DRanS_generic_stream(generic_state *stream)       // !InTc!
// function DRanS_generic_stream(stream) result(dran)
// INPUTS
//    stream     linear stream created by a Ran_XXXXXX_new_stream function
// OUTPUTS
//    64bit floating point number    -1.0 < value < 1.0
//****
{
//   generic_state *state = stream ;
//   return state->drans(stream);
  uint32_t value;
  if(stream->cur > stream->top) stream->refill(stream);
  value = stream->buf[stream->cur] ;
  stream->cur = stream->cur + 1;
  return CVTDBLS_32(value);
}
double F_DRanS_generic_stream(statep *s)  // Fortran interface using derived type
{
  return(DRanS_generic_stream(s->p));
}

//****f* librandom/VecIRan_generic_stream
// Synopsis
//    generate multiple unsigned 32bit integer values, according to the stream type
// ARGUMENTS
void VecIRan_generic_stream(generic_state *stream, unsigned int *ranbuf, int n)       // !InTc!
// subroutine VecIRan_generic_stream(stream, ranbuf, n)
// INPUTS
//    stream     linear stream created by a Ran_XXXXXX_new_stream function
//    n          number of values to generate
// OUTPUTS
//    ranbuf     array of integer output values
//****
{
//   generic_state *state = stream ;
//   state->vec_iran(stream,ranbuf,n);
  int topp1 = stream->top + 1;
  int cur = stream->cur;
  uint32_t *buf = stream->buf;

  while(n > 0 && cur < topp1){
    *ranbuf = buf[cur] ; ranbuf++ ; cur++; n-- ;
  }
  if(cur < topp1) {
    stream->cur = cur ;
  }else{
    stream->refill(stream);   // this will set stream->cur
  }
  if(n <= 0) return ;         // done
  while(n > 0){
    if(n >= topp1){             // full buffers
      for(cur=0 ; cur<topp1 ; cur++) {ranbuf[cur] = buf[cur] ; }
      ranbuf += topp1;
      stream->refill(stream);
      n -= topp1;
    }else{                    // less than 1 buffer, no need for refill
      for(cur=0 ; cur<n ; cur++) {ranbuf[cur] = buf[cur] ; }
      stream->cur = cur ;     // store current pointer
      n = 0;
    }
  }
}
void F_VecIRan_generic_stream(statep *s, unsigned int *ranbuf, int n)  // Fortran interface using derived type
{
  VecIRan_generic_stream(s->p,ranbuf,n);
}

//****f* librandom/VecDRan_generic_stream
// Synopsis
//    generate multiple 64bit positive floating point values, according to the stream type
//    ( 0.0 < values < 1.0 )
// ARGUMENTS
void VecDRan_generic_stream(generic_state *stream, double *ranbuf, int n)       // !InTc!
// subroutine VecDRan_generic_stream(stream, ranbuf, n)
// INPUTS
//    stream     linear stream created by a Ran_XXXXXX_new_stream function
//    n          number of values to generate
// OUTPUTS
//    ranbuf     array of double(real*8) output values
//****
{
//   generic_state *state = stream ;
//   state->vec_dran(stream,ranbuf,n);
  int topp1 = stream->top + 1;
  int cur = stream->cur;
  uint32_t *buf = stream->buf;

  while(n > 0 && cur < topp1){
    *ranbuf = CVTDBL_32(buf[cur]) ; ranbuf++ ; cur++; n-- ;
  }
  if(cur < topp1) {
    stream->cur = cur ;
  }else{
    stream->refill(stream);   // this will set stream->cur
  }
  if(n <= 0) return ;         // done
  while(n > 0){
    if(n >= topp1){             // full buffers
      for(cur=0 ; cur<topp1 ; cur++) {ranbuf[cur] = CVTDBL_32(buf[cur]) ; }
      ranbuf += topp1;
      stream->refill(stream);
      n -= topp1;
    }else{                    // less than 1 buffer, no need for refill
      for(cur=0 ; cur<n ; cur++) {ranbuf[cur] = CVTDBL_32(buf[cur]) ; }
      stream->cur = cur ;     // store current pointer
      n = 0;
    }
  }
}
void F_VecDRan_generic_stream(statep *s, double *ranbuf, int n)  // Fortran interface using derived type
{
  VecDRan_generic_stream(s->p,ranbuf,n);
}

//****f* librandom/VecDRanS_generic_stream
// Synopsis
//    generate multiple 64bit positive floating point values, according to the stream type
//    ( -1.0 < values < 1.0 )
// ARGUMENTS
void VecDRanS_generic_stream(generic_state *stream, double *ranbuf, int n)       // !InTc!
// subroutine VecDRanS_generic_stream(stream, ranbuf, n)
// INPUTS
//    stream     linear stream created by a Ran_XXXXXX_new_stream function
//    n          number of values to generate
// OUTPUTS
//    ranbuf     array of double (real*8) output values
//****
{
//   generic_state *state = stream ;
//   state->vec_drans(stream,ranbuf,n);
  int topp1 = stream->top + 1;
  int cur = stream->cur;
  uint32_t *buf = stream->buf;

  while(n > 0 && cur < topp1){
    *ranbuf = CVTDBLS_32(buf[cur]) ; ranbuf++ ; cur++; n-- ;
  }
  if(cur < topp1) {
    stream->cur = cur ;
  }else{
    stream->refill(stream);   // this will set stream->cur
  }
  if(n <= 0) return ;         // done
  while(n > 0){
    if(n >= topp1){             // full buffers
      for(cur=0 ; cur<topp1 ; cur++) {ranbuf[cur] = CVTDBLS_32(buf[cur]) ; }
      ranbuf += topp1;
      stream->refill(stream);
      n -= topp1;
    }else{                    // less than 1 buffer, no need for refill
      for(cur=0 ; cur<n ; cur++) {ranbuf[cur] = CVTDBLS_32(buf[cur]) ; }
      stream->cur = cur ;     // store current pointer
      n = 0;
    }
  }
}
void F_VecDRanS_generic_stream(statep *s, double *ranbuf, int n)  // Fortran interface using derived type
{
  VecDRanS_generic_stream(s->p,ranbuf,n);
}


