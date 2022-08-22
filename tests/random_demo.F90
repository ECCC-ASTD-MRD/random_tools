program demo
  use ISO_C_BINDING
  implicit none
  include 'randomfunctions.inc'
  integer, parameter :: NI = 10
  type(RANDOM_STREAM) :: s, clone
  integer, dimension(NI) :: ibuf
  integer :: iran, cSeed
  integer, dimension(1) :: piSeed
  real*8,  dimension(NI) :: dbuf, dsbuf
  real *8 :: dran, dsran
  cSeed = 0                         ! default initialization
  clone = RANDOM_STREAM(C_NULL_PTR) ! null clone, not duplicating existing stream
  call Ran_R250_new_stream(s, clone, piSeed, cSeed)   ! create R250 stream
  iran  = IRan_generic_stream(s)                      ! get 1 integer value
  dran  = DRan_generic_stream(s)                      ! get 1 real*8 value
  dsran = DRanS_generic_stream(s)                     ! get 1 real*8 value
  print *,'scalar ',iran, dran, dsran
  call VecIRan_generic_stream(s, ibuf, NI)            ! get NI integer values
  call VecDRan_generic_stream(s, dbuf, NI)            ! get NI real*8 values
  call VecDRanS_generic_stream(s, dsbuf, NI)          ! get NI real*8 values
  print *,'vector1',ibuf(1),dbuf(1),dsbuf(1)
  print *,'vector2',ibuf(NI),dbuf(NI),dsbuf(NI)
  call RanSetSeed_gaussian_stream(s,piSeed, cSeed)
  dran = DRan_gaussian_stream(s)
  print *,'gaussian ',dran
  stop
end
