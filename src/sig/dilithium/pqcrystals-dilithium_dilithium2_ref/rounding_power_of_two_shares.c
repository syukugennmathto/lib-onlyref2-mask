#include "binary.h"
#include "params.h"
#include "rounding_power_of_two_shares.h"

void masking_arithmetic_to_boolean_decompose(uint32_t* pR1, uint32_t* pR0, uint32_t * pR, uint32_t base)
{
    uint32_t mask = (1 << base) - 1;
	uint32_t b[shares],pR0p[shares],pR1p[shares],sd_1[shares] = {(mask>>1)+1,0};

    masking_arithmetic_to_boolean_convert(pR0p,pR);
	masking_boolean_lshift(pR0p,pR0p,32-base);
	masking_boolean_rshift(b,pR0p,31);
	masking_boolean_neg   (b,b);
	masking_boolean_lshift(b,b,base);
	masking_boolean_rshift(pR0p,pR0p,32-base);
	masking_boolean_xor   (pR0,pR0p,b);

	masking_arithmetic_add(b,pR,sd_1);
	masking_arithmetic_to_boolean_convert(pR1p,b);
	masking_boolean_rshift(pR1,pR1p,base);
}

void masking_arithmetic_to_boolean_highbits(uint32_t* pR1, uint32_t * pR, uint32_t base)
{
	uint32_t mask = (1 << base) - 1;
	uint32_t b[shares],pR1p[shares],sd_1[shares] = {(mask>>1)+1,0};
	masking_arithmetic_add(b,pR,sd_1);
	masking_arithmetic_to_boolean_convert(pR1p,b);
	masking_boolean_rshift(pR1,pR1p,base);
}

void masking_arithmetic_to_boolean_lowbits(uint32_t * pR0, uint32_t * pR, uint32_t base)
{
	uint32_t b[shares],pR0p[shares];
	masking_arithmetic_to_boolean_convert(pR0p,pR);
	masking_boolean_lshift(pR0p,pR0p,32-base);
	masking_boolean_rshift(b,pR0p,31);
	masking_boolean_neg   (b,b);
	masking_boolean_lshift(b,b,base);
	masking_boolean_rshift(pR0p,pR0p,32-base);
	masking_boolean_xor   (pR0,pR0p,b);
}


unsigned char masking_arithmetic_makeint(uint32_t* z, uint32_t* r, uint32_t base)
{
	uint32_t r_z[shares],r_z_1[shares],r_1[shares],t[shares],c;
	masking_arithmetic_add(r_z,r,z);
	masking_arithmetic_to_boolean_highbits(r_1,r,base);
	masking_arithmetic_to_boolean_highbits(r_z_1,r_z,base);
	masking_boolean_xor(t,r_1,r_z_1);
	masking_boolean_lshift(t,t,31);
	c = masking_boolean_fullxor(t);
	return (c >> 31);
}
