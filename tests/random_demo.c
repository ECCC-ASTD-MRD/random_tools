#include <stdio.h>
#include <stdint.h>
#include "randomfunctions.h"
#define NI 10
int main(){
  uint32_t mySeed = 123456;
  int cSeed = 0;
  generic_state *s;
  uint32_t iran, ibuf[NI];
  double dran, dsran, dbuf[NI], dsbuf[NI];
  s = (generic_state *)  Ran_R250_new_stream(NULL, &mySeed, cSeed);  // default seeding
  iran  = IRan_generic_stream(s) ;                     // get 1 integer value
  dran  = DRan_generic_stream(s) ;                     // get 1 double value
  dsran = DRanS_generic_stream(s);                     // get 1 double value
  printf("scalar  %10d %f %f\n",iran,dran,dsran);
  VecIRan_generic_stream(s, ibuf, NI) ;                // get NI integer values
  VecDRan_generic_stream(s, dbuf, NI) ;                // get NI double values
  VecDRanS_generic_stream(s, dsbuf, NI) ;              // get NI double values
  printf("vector1 %10d %f %f\n",ibuf[0],dbuf[0],dsbuf[0]);
  printf("vector2 %10d %f %f\n",ibuf[NI-1],dbuf[NI-1],dsbuf[NI-1]);
  RanSetSeed_gaussian_stream(s, &mySeed, cSeed);
  dran = DRan_gaussian_stream(s);
  printf("gaussian %f\n",dran);
  return 0;
}
