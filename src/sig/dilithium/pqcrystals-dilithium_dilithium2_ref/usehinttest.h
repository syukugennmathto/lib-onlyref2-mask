#ifndef USEHINTEST_H
#define USEHINTEST_H

#include <stddef.h>
#include <stdint.h>
#include "params.h"
#include "polyvec.h"
#include "poly.h"


void  masking_usinghint(uint32_t ring,uint32_t zoo);
void  num_masking_usinghint(uint32_t ring,uint32_t zoo);
void  noking_usinghint(uint32_t ring,uint32_t zoo);
int test_use_hint_both();

#endif
