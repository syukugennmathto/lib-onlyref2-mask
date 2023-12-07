#ifndef ROUNDING_H
#define ROUNDING_H

#include <stdio.h>
#include <stdint.h>
#include "params.h"
//#include "rounding.h"


#define printBinary DILITHIUM_NAMESPACE(printBinary)
void printBinary(unsigned int num);

#define power2round DILITHIUM_NAMESPACE(power2round)
uint8_t power2round(int32_t *a0, int32_t a);

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

#define pt_decompose DILITHIUM_NAMESPACE(pt_decompose)
void pt_decompose(uint32_t *r1, uint32_t *r0, uint32_t  r , uint32_t base);

#define pt_highbits DILITHIUM_NAMESPACE(pt_highbits)
void pt_highbits(uint32_t *r1, uint32_t  r , uint32_t base);

#define pt_makehint DILITHIUM_NAMESPACE(pt_makehint)
unsigned char pt_makehint (uint32_t  z , uint32_t r, uint32_t base);

#define pt_usehint DILITHIUM_NAMESPACE(pt_usehint)
uint32_t pt_usehint  (uint32_t  y , uint32_t r, uint32_t base);




#define uuse_hint DILITHIUM_NAMESPACE(uuse_hint)
uint32_t uuse_hint(uint32_t r, unsigned int hint,uint32_t base);

#define arithmetic_to_boolean_use_hint DILITHIUM_NAMESPACE(arithmetic_to_boolean_use_hint)
void arithmetic_to_boolean_use_hint(uint32_t *r, unsigned int hint,uint32_t base);

#define hhighbbits DILITHIUM_NAMESPACE(hhighbbits)
uint32_t hhighbbits(uint32_t output, uint32_t input);


#define highbits DILITHIUM_NAMESPACE(highbits)
void highbits(uint32_t r1, uint32_t r, uint32_t base);

#define nos_highbits DILITHIUM_NAMESPACE(nos_highbits)
int nos_highbits(int32_t *input, int32_t base);

#define lowbits DILITHIUM_NAMESPACE(lowbits)
void lowbits(uint32_t r0, uint32_t  r , uint32_t base);

#define masking_arithmetic_to_boolean_decompose DILITHIUM_NAMESPACE(masking_arithmetic_to_boolean_decompose)
void masking_arithmetic_to_boolean_decompose(uint32_t* pR1, uint32_t* pR0, uint32_t *pR, uint32_t base);

#define masking_arithmetic_to_boolean_highbits DILITHIUM_NAMESPACE(masking_arithmetic_to_boolean_highbits)
void masking_arithmetic_to_boolean_highbits (uint32_t* pR1, uint32_t* pR , uint32_t base);

#define masking_arithmetic_to_boolean_lowbits DILITHIUM_NAMESPACE(masking_arithmetic_to_boolean_lowbits)
void masking_arithmetic_to_boolean_lowbits  (uint32_t* pR0, uint32_t* pR , uint32_t base);

#define masking_arithmetic_makeint DILITHIUM_NAMESPACE(masking_arithmetic_makeint)
unsigned char masking_arithmetic_makeint(uint32_t* z, uint32_t* r, uint32_t base);

#define masking_arithmetic_two_makeint DILITHIUM_NAMESPACE(masking_arithmetic_two_makeint)
unsigned char masking_arithmetic_two_makeint(uint32_t* z, uint32_t* r, uint32_t base);


#define new_masking_arithmetic_to_boolean_makehint DILITHIUM_NAMESPACE(new_masking_arithmetic_to_boolean_makehint)
void new_masking_arithmetic_to_boolean_makehint(int32_t* c,const int32_t* z, const int32_t* r,mask_point pd);



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


#define num_masking_arithmetic_to_boolean_decompose DILITHIUM_NAMESPACE(num_masking_arithmetic_to_boolean_decompose)
void num_masking_arithmetic_to_boolean_decompose(uint32_t pR1, uint32_t pR0, uint32_t pR, uint32_t base);

#define num_masking_arithmetic_to_boolean_highbits DILITHIUM_NAMESPACE(num_masking_arithmetic_to_boolean_highbits)
void num_masking_arithmetic_to_boolean_highbits(uint32_t pR1, uint32_t pR, uint32_t base);

#define num_masking_arithmetic_to_boolean_lowbits DILITHIUM_NAMESPACE(num_masking_arithmetic_to_boolean_lowbits)
void num_masking_arithmetic_to_boolean_lowbits(uint32_t pR0, uint32_t pR, uint32_t base);

#define num_masking_arithmetic_makehint DILITHIUM_NAMESPACE(num_masking_arithmetic_makehint)
unsigned char num_masking_arithmetic_makehint(uint32_t z, uint32_t r, uint32_t base);

#define num_masking_arithmetic_usehint DILITHIUM_NAMESPACE(num_masking_arithmetic_usehint)
unsigned char num_masking_arithmetic_usehint(unsigned char h, uint32_t r, uint32_t base);

#define num_masking_boolean_refresh DILITHIUM_NAMESPACE(num_masking_boolean_refresh)
void num_masking_boolean_refresh(uint32_t pC, uint32_t pA);

#define num_masking_boolean_fullrefresh DILITHIUM_NAMESPACE(num_masking_boolean_fullrefresh)
void num_masking_boolean_fullrefresh(uint32_t pC, uint32_t pA);

#define num_masking_boolean_fullxor DILITHIUM_NAMESPACE(num_masking_boolean_fullxor)
uint32_t num_masking_boolean_fullxor(uint32_t pA);

#define num_masking_boolean_generate DILITHIUM_NAMESPACE(num_masking_boolean_generate)
void num_masking_boolean_generate(uint32_t pC, uint32_t a);

#define num_masking_boolean_xor DILITHIUM_NAMESPACE(num_masking_boolean_xor)
void num_masking_boolean_xor(uint32_t pC, uint32_t pA, uint32_t pB);

#define num_masking_boolean_neg DILITHIUM_NAMESPACE(num_masking_boolean_neg)
void num_masking_boolean_neg(uint32_t pC, uint32_t pA);

#define num_masking_boolean_not DILITHIUM_NAMESPACE(num_masking_boolean_not)
void num_masking_boolean_not(uint32_t pC, uint32_t pA);

#define num_masking_boolean_lshift DILITHIUM_NAMESPACE(num_masking_boolean_lshift)
void num_masking_boolean_lshift(uint32_t pC, uint32_t pA, uint32_t shift);

#define num_masking_boolean_rshift DILITHIUM_NAMESPACE(num_masking_boolean_rshift)
void num_masking_boolean_rshift(uint32_t pC, uint32_t pA, uint32_t shift);

#define num_masking_boolean_and DILITHIUM_NAMESPACE(num_masking_boolean_and)
void num_masking_boolean_and(uint32_t pC, uint32_t pA, uint32_t pB);

#define num_masking_boolean_add DILITHIUM_NAMESPACE(num_masking_boolean_add)
void num_masking_boolean_add(uint32_t pC, uint32_t pA, uint32_t pB);

#define num_masking_boolean_unprotected_generate DILITHIUM_NAMESPACE(num_masking_boolean_unprotected_generate)
void num_masking_boolean_unprotected_generate(uint32_t pC, uint32_t a);

#define num_masking_boolean_unprotected_recombine DILITHIUM_NAMESPACE(num_masking_boolean_unprotected_recombine)
uint32_t num_masking_boolean_unprotected_recombine(uint32_t pA);

#define num_masking_arithmetic_add DILITHIUM_NAMESPACE(num_masking_arithmetic_add)
void num_masking_arithmetic_add(uint32_t pC, uint32_t pA, uint32_t pB);

#define num_masking_arithmetic_generate DILITHIUM_NAMESPACE(num_masking_arithmetic_generate)
void num_masking_arithmetic_generate(uint32_t pC, uint32_t a);

#define num_masking_arithmetic_unprotected_recombine DILITHIUM_NAMESPACE(num_masking_arithmetic_unprotected_recombine)
uint32_t num_masking_arithmetic_unprotected_recombine(uint32_t pR);

#define num_masking_arithmetic_to_boolean_convert DILITHIUM_NAMESPACE(num_masking_arithmetic_to_boolean_convert)
void num_masking_arithmetic_to_boolean_convert(uint32_t pC, uint32_t pA);

#define num_masking_boolean_to_arithmetic_convert DILITHIUM_NAMESPACE(num_masking_boolean_to_arithmetic_convert)
void num_masking_boolean_to_arithmetic_convert(uint32_t pC, uint32_t pA);


//#define masking_use_hint DILITHIUM_NAMESPACE(masking_use_hint)
//uint32_t masking_use_hint(uint32_t a, unsigned char hint,uint32_t base);


//#define prime_decompose DILITHIUM_NAMESPACE(prime_decompose)
//void prime_decompose(uint32_t *r1, uint32_t *r0, uint32_t r);

#define prime_lowbits DILITHIUM_NAMESPACE(prime_lowbits)
void prime_lowbits  (uint32_t *r0, uint32_t r);

#define prime_highbits DILITHIUM_NAMESPACE(prime_highbits)
void prime_highbits (uint32_t *r1, uint32_t r);

#define masking_arithmetic_to_boolean_prime_decompose DILITHIUM_NAMESPACE(masking_arithmetic_to_boolean_prime_decompose)
void masking_arithmetic_to_boolean_prime_decompose(uint32_t* pR1,uint32_t* pR0, uint32_t * pR,uint32_t aalpha);

#define masking_arithmetic_to_boolean_prime_lowbits DILITHIUM_NAMESPACE(masking_arithmetic_to_boolean_prime_lowbits)
void masking_arithmetic_to_boolean_prime_lowbits(uint32_t* pR0, uint32_t * pR);

#define masking_arithmetic_to_boolean_prime_highbits DILITHIUM_NAMESPACE(masking_arithmetic_to_boolean_prime_highbits)
void masking_arithmetic_to_boolean_prime_highbits(uint32_t* pR1, uint32_t * pR,uint32_t aalpha);

#define random_uint32_t_mono DILITHIUM_NAMESPACE(random_uint32_t_mono)
uint32_t random_uint32_t_mono(uint32_t r);

#define random_int32_t DILITHIUM_NAMESPACE(random_int32_t)
int32_t random_int32_t();

#define new_masking_decompose DILITHIUM_NAMESPACE(new_masking_decompose)
int32_t  new_masking_decompose(int32_t pR0,const int32_t pR);

#define new_masking_arithmetic_to_boolean_highbits DILITHIUM_NAMESPACE(new_masking_arithmetic_to_boolean_highbits)
void new_masking_arithmetic_to_boolean_highbits(int32_t* pR1,const int32_t* pR,mask_point* pD);

#define new_unsking_arithmetic_to_boolean_highbits DILITHIUM_NAMESPACE(new_unsking_arithmetic_to_boolean_highbits)
void new_unsking_arithmetic_to_boolean_highbits(int32_t* pR,int32_t* pR1,mask_point pD);


#define new_masking_arithmetic_to_boolean_lowbits DILITHIUM_NAMESPACE(new_masking_arithmetic_to_boolean_lowbits)
void new_masking_arithmetic_to_boolean_lowbits(int32_t* pR0,const int32_t* pR,mask_point pD);

#define new_nosking_arithmetic_to_boolean_highbits DILITHIUM_NAMESPACE(new_nosking_arithmetic_to_boolean_highbits)
void new_nosking_arithmetic_to_boolean_highbits(int32_t* pR1,const int32_t* pR,mask_point pD);


#define new_masking_use_hint DILITHIUM_NAMESPACE(new_masking_use_hint)
void new_masking_use_hint(int32_t* pR2,const int32_t* hint,mask_point pD);

#endif
