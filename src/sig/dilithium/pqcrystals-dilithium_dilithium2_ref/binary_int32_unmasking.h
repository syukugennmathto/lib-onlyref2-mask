#ifndef _MASKING_BASE_H_
#define _MASKING_BASE_H_
#include "random.h"
#include "rounding.h"
#include <stdint.h>
#include <stdio.h>




#define int_unmasking_boolean_refresh DILITHIUM_NAMESPACE(int_unmasking_boolean_refresh)
void int_unmasking_boolean_refresh(int32_t*pC, int32_t *pA,pMsm pMsx);

#define int_unmasking_boolean_generate DILITHIUM_NAMESPACE(int_unmasking_boolean_generate)
void int_unmasking_boolean_generate(int32_t* pC,int32_t a,svec *pVa);

#define int_unmasking_boolean_xor DILITHIUM_NAMESPACE(int_unmasking_boolean_xor)
void int_masking_boolean_xor(int32_t* pC, int32_t *pA, int32_t *pB);


#define int_unmasking_boolean_lshift DILITHIUM_NAMESPACE(int_unmasking_boolean_lshift)
void int_unmasking_boolean_lshift(int32_t*pC, int32_t *pA, int32_t shift);

#define int_unmasking_boolean_rshift DILITHIUM_NAMESPACE(int_unmasking_boolean_rshift)
void int_unmasking_boolean_rshift(int32_t*pC, int32_t *pA, int32_t shift);

#define int_unmasking_arithmetic_add DILITHIUM_NAMESPACE(int_unmasking_arithmetic_add)
void int_unmasking_arithmetic_add(int32_t* pC,int32_t* pA,int32_t* pB);

#define int_unmasking_boolean_and DILITHIUM_NAMESPACE(int_unmasking_boolean_and)
void int_unmasking_boolean_and(int32_t*pC, int32_t *pA, int32_t *pB,pMsm pMsx);

#define int_unmasking_boolean_add DILITHIUM_NAMESPACE(int_unmasking_boolean_add)
void int_unmasking_boolean_add(int32_t*pC, int32_t *pA, int32_t *pB,pMtM *pMx);


#define int_unmasking_arithmetic_to_boolean_convert DILITHIUM_NAMESPACE(int_unmasking_arithmetic_to_boolean_convert)
void int_unmasking_arithmetic_to_boolean_convert(int32_t*pC, int32_t *pA,pMtL*  pMR);


#endif // _MASKING_POLYNOMIALS_H_
