#ifndef ROUNDING_H
#define ROUNDING_H

#include <stdio.h>
#include <stdint.h>
#include "params.h"
#include "rounding.h"

#define power2round DILITHIUM_NAMESPACE(power2round)
int32_t power2round(int32_t *a0, int32_t a);

#define decompose_array DILITHIUM_NAMESPACE(decompose_array)
void decompose_array(int32_t *a0_array, int32_t *a1_array, const int32_t *a_array, size_t length);

#define decompose_array DILITHIUM_NAMESPACE(decompose_array)
int32_t decompose(int32_t *a0, int32_t a);

#define prime_decompose DILITHIUM_NAMESPACE(prime_decompose)
void prime_decompose(uint32_t *r1, uint32_t *r0, uint32_t r);

#define make_hint DILITHIUM_NAMESPACE(make_hint)
unsigned int make_hint(int32_t a0, int32_t a1);

#define use_hint DILITHIUM_NAMESPACE(use_hint)
int32_t use_hint(int32_t a, unsigned int hint);

#define highbits DILITHIUM_NAMESPACE(highbits)
void highbits(uint32_t *r1, uint32_t r, uint32_t base);

#define lowbits DILITHIUM_NAMESPACE(lowbits)
void lowbits(uint32_t *r0, uint32_t  r , uint32_t base);

#define masking_arithmetic_to_boolean_decompose DILITHIUM_NAMESPACE(masking_arithmetic_to_boolean_decompose)
void masking_arithmetic_to_boolean_decompose(uint32_t* pR1, uint32_t* pR0, uint32_t *pR, uint32_t base);

#define masking_arithmetic_to_boolean_highbits DILITHIUM_NAMESPACE(masking_arithmetic_to_boolean_highbits)
void masking_arithmetic_to_boolean_highbits (uint32_t* pR1, uint32_t* pR , uint32_t base);

#define masking_arithmetic_to_boolean_lowbits DILITHIUM_NAMESPACE(masking_arithmetic_to_boolean_lowbits)
void masking_arithmetic_to_boolean_lowbits  (uint32_t* pR0, uint32_t* pR , uint32_t base);

#define masking_arithmetic_makeint DILITHIUM_NAMESPACE(masking_arithmetic_makeint)
unsigned char masking_arithmetic_makeint(uint32_t* z, uint32_t* r, uint32_t base);

#define masking_boolean_refresh DILITHIUM_NAMESPACE(masking_boolean_refresh)
void masking_boolean_refresh(uint32_t* pC, uint32_t *pA);

#define masking_boolean_fullrefresh DILITHIUM_NAMESPACE(masking_boolean_fullrefresh)
void masking_boolean_fullrefresh(uint32_t* pC, uint32_t *pA);

#define masking_boolean_fullxor DILITHIUM_NAMESPACE(masking_boolean_fullxor)
uint32_t masking_boolean_fullxor(uint32_t* pA);

#define masking_boolean_generate DILITHIUM_NAMESPACE(masking_boolean_generate)
void masking_boolean_generate(uint32_t* pC, uint32_t a);

#define masking_boolean_xor DILITHIUM_NAMESPACE(masking_boolean_xor)
void masking_boolean_xor(uint32_t* pC, uint32_t *pA, uint32_t *pB);

#define masking_boolean_neg DILITHIUM_NAMESPACE(masking_boolean_neg)
void masking_boolean_neg(uint32_t* pC, uint32_t *pA);

#define masking_boolean_not DILITHIUM_NAMESPACE(masking_boolean_not)
void masking_boolean_not(uint32_t* pC, uint32_t *pA);

#define masking_boolean_lshift DILITHIUM_NAMESPACE(masking_boolean_lshift)
void masking_boolean_lshift(uint32_t* pC, uint32_t *pA, uint32_t shift);

#define masking_boolean_rshift DILITHIUM_NAMESPACE(masking_boolean_rshift)
void masking_boolean_rshift(uint32_t* pC, uint32_t *pA, uint32_t shift);

#define masking_boolean_and DILITHIUM_NAMESPACE(masking_boolean_and)
void masking_boolean_and(uint32_t* pC, uint32_t *pA, uint32_t *pB);

#define masking_boolean_add DILITHIUM_NAMESPACE(masking_boolean_add)
void masking_boolean_add(uint32_t* pC, uint32_t *pA, uint32_t *pB);

#define masking_boolean_unprotected_generate DILITHIUM_NAMESPACE(masking_boolean_unprotected_generate)
void masking_boolean_unprotected_generate(uint32_t* pC, uint32_t a);

#define masking_boolean_unprotected_recombine DILITHIUM_NAMESPACE(masking_boolean_unprotected_recombine)
uint32_t masking_boolean_unprotected_recombine(uint32_t* pA);

#define masking_arithmetic_add DILITHIUM_NAMESPACE(masking_arithmetic_add)
void masking_arithmetic_add(uint32_t* pC, uint32_t* pA, uint32_t* pB);

#define masking_arithmetic_generate DILITHIUM_NAMESPACE(masking_arithmetic_generate)
void masking_arithmetic_generate(uint32_t* pC, uint32_t a);

#define masking_arithmetic_unprotected_recombine DILITHIUM_NAMESPACE(masking_arithmetic_unprotected_recombine)
uint32_t masking_arithmetic_unprotected_recombine(uint32_t* pR);

#define masking_arithmetic_to_boolean_convert DILITHIUM_NAMESPACE(masking_arithmetic_to_boolean_convert)
void masking_arithmetic_to_boolean_convert(uint32_t* pC, uint32_t *pA);

#define masking_boolean_to_arithmetic_convert DILITHIUM_NAMESPACE(masking_boolean_to_arithmetic_convert)
void masking_boolean_to_arithmetic_convert(uint32_t* pC, uint32_t *pA);

#define random_uint32_t DILITHIUM_NAMESPACE(random_uint32_t)
uint32_t random_uint32_t();

#endif
