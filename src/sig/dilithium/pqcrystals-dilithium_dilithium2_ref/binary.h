#ifndef _MASKING_BASE_H_
#define _MASKING_BASE_H_
#include "random.h"
#include <stdint.h>
#include <stdio.h>




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

#endif // _MASKING_POLYNOMIALS_H_
