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
//****P* librandom/random value generators
// Synopsis
// generic access to multiple families of random number generators (NOT for cryptographic usage)
//
// the user will create a random number "stream" of a given type (see available generators)
// and call various functions/routines to produce random number sequences (scalar or vector).
// the stream generators are believed to be "thread safe" once they are initialized
// 
// scalar values :
//   unsigned integer value      0 <=  value  < 2**32 - 1
//   unsigned double value     0.0 <   value  < 1.0
//   signed double value      -1.0 <   value  < 1.0
// vector of values : same 3 kinds as scalar values
//
// available generators  
// linear generators :
//    R250     (shift register sequence)     https://www.taygeta.com/rwalks/node2.html
//    MT19937  (Mersenne twister)            https://en.wikipedia.org/wiki/Mersenne_Twister
//    SHR3     (3-shift-register)
//    XSR128   (xor-shift)                   https://en.wikipedia.org/wiki/Xorshift
//    XSR128R  (xor-shift-rotate)            https://en.wikipedia.org/wiki/Xoroshiro128%2B
// non linear generators :
//    gaussian (ziggurat with 256 bins)      https://en.wikipedia.org/wiki/Ziggurat_algorithm
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
//
// Fortran arguments
//    type(RANDOM_STREAM), intent(IN)              :: stream
//    integer(c_INT), dimension(cSeed), intent(IN) :: piSeed
//    integer(C_INT), intent(IN), value            :: cSeed
//    integer(C_INT), intent(IN), value            :: n
//    integer(C_INT)                               :: Iran
//    real(C_DOUBLE)                               :: Dran
//    integer(C_INT), dimension(n), intent(OUT)    :: Ibuf
//    real(C_DOUBLE), dimension(n), intent(OUT)    :: Dbuf
//
// C note:
//    to compile and load successfully, one must get the interface definitions from
//
//    #include <randomfunctions.h>
//
// Fortran note:
//    to compile and load successfully, one must get the interface definitions from
//
//    use ISO_C_BINDING
//    include 'randomfunctions.inc'
//****
//****f* librandom/RanSetSeed_generic_stream
// Synopsis
//    reseed a random generator "stream" previously created by a call to
//    Ran_XXXXXX_new_stream (where XXXXXX is one of above mentioned linear generators)
//
// ARGUMENTS
void RanSetSeed_generic_stream(generic_state *stream, unsigned int *piSeed, int cSeed)  // !InTc!
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
// INPUTS
//    stream     linear stream created by a Ran_XXXXXX_new_stream function
// OUTPUTS
//    unsigned 32bit integer,   0 <= value < 1
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
//    ( 0 <= values < 1 )
// ARGUMENTS
void VecIRan_generic_stream(generic_state *stream, unsigned int *ranbuf, int n)       // !InTc!
// INPUTS
//    stream     linear stream created by a Ran_XXXXXX_new_stream function
//    n          number of values to generate
// OUTPUTS
//    ranbuf     array ot output values
//****
{
//   generic_state *state = stream ;
//   state->vec_iran(stream,ranbuf,n);
  int topp1 = stream->top + 1;
  int cur = stream->cur;
  int *buf = stream->buf;
  int i;
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
// INPUTS
//    stream     linear stream created by a Ran_XXXXXX_new_stream function
//    n          number of values to generate
// OUTPUTS
//    ranbuf     array ot output values
//****
{
//   generic_state *state = stream ;
//   state->vec_dran(stream,ranbuf,n);
  int topp1 = stream->top + 1;
  int cur = stream->cur;
  int *buf = stream->buf;
  int i;
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
// INPUTS
//    stream     linear stream created by a Ran_XXXXXX_new_stream function
//    n          number of values to generate
// OUTPUTS
//    ranbuf     array ot output values
//****
{
//   generic_state *state = stream ;
//   state->vec_drans(stream,ranbuf,n);
  int topp1 = stream->top + 1;
  int cur = stream->cur;
  int *buf = stream->buf;
  int i;
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


