#ifndef USEHINTEST_H
#define USEHINTEST_H

#include <stddef.h>
#include <stdint.h>
#include <math.h>
#include "params.h"
#include "polyvec.h"
#include "poly.h"
#include "binary_int32.h"


void  masking_usinghint(uint32_t ring,uint32_t zoo);
void  num_masking_usinghint(uint32_t ring,uint32_t zoo);
void  noking_usinghint(uint32_t ring,uint32_t zoo);

void int32_masking_arithmetic_to_boolean_highbits(int32_t* pR1, int32_t * pR, int32_t base,pMtL*  pMR);
void int32_unmasking_arithmetic_to_boolean_highbits(int32_t* pR1, int32_t * pR, int32_t base,pMtL*  pMR);
int test_use_hint_both(int num);

#endif
