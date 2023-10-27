#ifndef _MASKING_ROUNDING_H_
#define _MASKING_ROUNDING_H_


// 関数リストを条件に従って書き換える


#define masking_arithmetic_to_boolean_decompose DILITHIUM_NAMESPACE(masking_arithmetic_to_boolean_decompose)
void masking_arithmetic_to_boolean_decompose(uint32_t* pR1, uint32_t* pR0, uint32_t * pR, uint32_t base);

#define masking_arithmetic_to_boolean_highbits DILITHIUM_NAMESPACE(masking_arithmetic_to_boolean_highbits)
void masking_arithmetic_to_boolean_highbits (uint32_t* pR1, uint32_t* pR , uint32_t base);

#define masking_arithmetic_to_boolean_lowbits DILITHIUM_NAMESPACE(masking_arithmetic_to_boolean_lowbits)
void masking_arithmetic_to_boolean_lowbits  (uint32_t* pR0, uint32_t* pR , uint32_t base);

#define masking_arithmetic_makeint DILITHIUM_NAMESPACE(masking_arithmetic_makeint)
unsigned char masking_arithmetic_makeint(uint32_t* z, uint32_t* r, uint32_t base);

#endif
