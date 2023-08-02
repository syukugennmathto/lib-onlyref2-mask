#include <stdio.h>
#include <stdint.h>
#include "params.h"
#include "rounding.h"
#include "random.h"
#include "randombytes.h"



/*************************************************
* Name:        power2round
*
* Description: For finite field element a, compute a0, a1 such that
*              a mod^+ Q = a1*2^D + a0 with -2^{D-1} < a0 <= 2^{D-1}.
*              Assumes a to be standard representative.
*
* Arguments:   - int32_t a: input element
*              - int32_t *a0: pointer to output element a0
*
* Returns a1.
**************************************************/
uint8_t power2round(int32_t *a0, int32_t a)  {
  int32_t a1;

  a1 = (a + (1 << (D-1)) - 1) >> D;
  *a0 = a - (a1 << D);
  return a1;
}

/*************************************************
* Name:        decompose
*
* Description: For finite field element a, compute high and low bits a0, a1 such
*              that a mod^+ Q = a1*ALPHA + a0 with -ALPHA/2 < a0 <= ALPHA/2 except
*              if a1 = (Q-1)/ALPHA where we set a1 = 0 and
*              -ALPHA/2 <= a0 = a mod^+ Q - Q < 0. Assumes a to be standard
*              representative.
*
* Arguments:   - int32_t a: input element
*              - int32_t *a0: pointer to output element a0
*
* Returns a1.
**************************************************/
void decompose_array(int32_t *a0_array, int32_t *a1_array, const int32_t *a_array, size_t length) {
    uint32_t mask = (1 << D) - 1;
    uint32_t d_1 = (mask >> 1) + 1;

    for (size_t i = 0; i < length; i++) {
        // Convert to unsigned 32-bit integers for intermediate calculations
        uint32_t a_unsigned = (uint32_t)a_array[i];

        // Calculate a1 and a0 in unsigned format
        uint32_t a1_unsigned = (a_unsigned + 127) >> 7;
    #if GAMMA2 == (Q-1)/32
        a1_unsigned = (a1_unsigned * 1025 + (1 << 21)) >> 22;
        a1_unsigned &= 15;
    #elif GAMMA2 == (Q-1)/88
        a1_unsigned = (a1_unsigned * 11275 + (1 << 23)) >> 24;
        a1_unsigned ^= ((43 - a1_unsigned) >> 31) & a1_unsigned;
    #endif

        uint32_t a0_unsigned = a_unsigned - a1_unsigned * 2 * GAMMA2;
        a0_unsigned -= (((Q - 1) / 2 - a0_unsigned) >> 31) & Q;

        // Convert back to signed 32-bit integer format and store in the arrays
        a0_array[i] = (int32_t)a0_unsigned;
        a1_array[i] = (int32_t)a1_unsigned;
    }
}


void prime_decompose(uint32_t *r1, uint32_t *r0, uint32_t r)
{
    uint32_t q     = 8380417;
    
    uint32_t r1p,r0p,m;

    r0p  = (r & 0x7FFFF) + ((r>>19)<<9);
    r0p -= alpha/2 + 1;
    r0p += (-(r0p >> 31)) & alpha;
    r0p -= alpha/2 - 1;

    r1p  = r - r0p;
    r1p  = (r1p >> 19) + 1 - ((r1p-1)>>31);

    m    = ~(-(r1p>>4));
    r1p  = m&r1p;
    r0p  = r0p + (m&1) - 1 + q;
    *r0  = r0p;
    *r1  = r1p;

}


/*************************************************
* Name:        decompose
*
* Description: For finite field element a, compute high and low bits a0, a1 such
*              that a mod^+ Q = a1*ALPHA + a0 with -ALPHA/2 < a0 <= ALPHA/2 except
*              if a1 = (Q-1)/ALPHA where we set a1 = 0 and
*              -ALPHA/2 <= a0 = a mod^+ Q - Q < 0. Assumes a to be standard
*              representative.
*
* Arguments:   - int32_t a: input element
*              - int32_t *a0: pointer to output element a0
*
* Returns a1.
**************************************************/
int32_t decompose(int32_t *a0, int32_t a) {
    uint32_t mask = (1 << D) - 1;
    uint32_t d_1 = (mask >> 1) + 1;
    uint32_t a0_unsigned, a1_unsigned;

    // Convert to unsigned 32-bit integers for intermediate calculations
    uint32_t a_unsigned = (uint32_t)a;

    // Calculate a1 and a0 in unsigned format
    a1_unsigned = (a_unsigned + 127) >> 7;
#if GAMMA2 == (Q-1)/32
    a1_unsigned = (a1_unsigned * 1025 + (1 << 21)) >> 22;
    a1_unsigned &= 15;
#elif GAMMA2 == (Q-1)/88
    a1_unsigned = (a1_unsigned * 11275 + (1 << 23)) >> 24;
    a1_unsigned ^= ((43 - a1_unsigned) >> 31) & a1_unsigned;
#endif

    a0_unsigned = a_unsigned - a1_unsigned * 2 * GAMMA2;
    a0_unsigned -= (((Q - 1) / 2 - a0_unsigned) >> 31) & Q;

    // Convert back to signed 32-bit integer format
    *a0 = (int32_t)a0_unsigned;

    // Convert a1 back to signed 32-bit integer format and return it
    return (int32_t)a1_unsigned;
}

/*************************************************
* Name:        make_hint
*
* Description: Compute hint bit indicating whether the low bits of the
*              input element overflow into the high bits.
*
* Arguments:   - int32_t a0: low bits of input element
*              - int32_t a1: high bits of input element
*
* Returns 1 if overflow.
**************************************************/
unsigned int make_hint(int32_t a0, int32_t a1) {
    
    uint32_t r1, v1;
    uint32_t mask = (1 << D) - 1;
    uint32_t d_1 = (mask >> 1) + 1;
    uint32_t a0_unsigned = (uint32_t)a0;
    uint32_t a1_unsigned = (uint32_t)a1;

    highbits(r1, a0_unsigned, D);
    highbits(v1, a0_unsigned + a1_unsigned, D);

    if (r1 == v1) {
        return 0;
    } else {
        return 1;
    }
}


/*************************************************
* Name:        use_hint
*
* Description: Correct high bits according to hint.
*
* Arguments:   - int32_t a: input element
*              - unsigned int hint: hint bit
*
* Returns corrected high bits.
**************************************************/
int32_t use_hint(int32_t a, unsigned int hint) {
    int32_t a0, a1;

    // コード2のdecompose関数を呼び出してa0とa1を計算
    a1 = decompose(&a0, a);

    if (hint == 1) {
        // hintが1の場合、条件によってa1を調整
#if GAMMA2 == (Q-1)/32
        if (a0 > 0) {
            a1 = (a1 + 1) & 15;
        } else {
            a1 = (a1 - 1) & 15;
        }
#elif GAMMA2 == (Q-1)/88
        if (a0 > 0) {
            a1 = (a1 == 43) ?  0 : a1 + 1;
        } else {
            a1 = (a1 ==  0) ? 43 : a1 - 1;
        }
#endif
    }

    return a1;
}



/*************************************************
* Name:        highbits
*
* Description: Compute the high bits of a finite field element `r` given a base value.
*              The base value determines the number of bits to shift the element `r` to the right.
*
* Arguments:   - uint32_t *r1: pointer to the output high bits of `r`
*              - uint32_t r: input element
*              - uint32_t base: base value used for shifting `r`
**************************************************/
void highbits(uint32_t output, uint32_t input, uint32_t base) {
    uint32_t mask = (1 << D) - 1;
    uint32_t d_1 = (mask >> 1) + 1;
    output = (input + d_1) >> base;
}


void lowbits(uint32_t output, uint32_t  input, uint32_t base)
{
    uint32_t r0p,b;
    //activate_trigger_aux();
    r0p = (input<<(32-base));
    //desactivate_trigger_aux();
    b   = (-(r0p>>(31)));
    b   = b<<base;
    r0p = (r0p>>(32-base));
    r0p = r0p^b;
    output = r0p;
}


void masking_arithmetic_to_boolean_decompose(uint32_t* pR1, uint32_t* pR0, uint32_t *pR, uint32_t base)
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

void masking_arithmetic_to_boolean_highbits(uint32_t* pR1, uint32_t* pR, uint32_t base)
{
    uint32_t mask = (1 << base) - 1;
    uint32_t b[shares],pR1p[shares],sd_1[shares] = {(mask>>1)+1,0};
    masking_arithmetic_add(b,pR,sd_1);
    masking_arithmetic_to_boolean_convert(pR1p,b);
    masking_boolean_rshift(pR1,pR1p,base);
}

void masking_arithmetic_to_boolean_lowbits(uint32_t* pR0, uint32_t* pR, uint32_t base)
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

void masking_boolean_refresh(uint32_t*pC, uint32_t *pA)
{
    for(uint32_t j = 0; j < shares; ++j)
    {
        pC[j] = pA[j];
    }

    for(uint32_t i = 0; i < shares;++i)
    {
        for(uint32_t j = i+1; j < shares;++j)
        {
            uint32_t r;
            r = random_uint32_t();
            pC[i] ^=r;
            pC[j] ^=r;
        }
    }
}

void masking_boolean_fullrefresh(uint32_t *pC, uint32_t *pA)
{
    for(uint32_t j = 0; j < shares; ++j)
    {
        pC[j] = pA[j];
    }

    for(uint32_t i = 0; i < shares;++i)
    {
        for(uint32_t j = 1; j < shares;++j)
        {
            uint32_t r;
            r = random_uint32_t();
            pC[0] ^=r;
            pC[j] ^=r;
        }
    }
}

uint32_t masking_boolean_fullxor(uint32_t *pA)
{
    uint32_t r ;
    uint32_t x [shares];
    masking_boolean_fullrefresh(x,pA);
    r = x[0];
    for(uint32_t j = 1; j < shares; ++j)
    {
        r ^=x[j];
    }
    return r;
}

void masking_boolean_generate(uint32_t* pC,uint32_t a)
{
    pC[0] = a;
    for(uint32_t i = 1; i < shares;++i)
    {
        pC[i] = random_uint32_t();
        pC[0] ^= pC[i];
    }
}

void masking_boolean_xor(uint32_t*pC, uint32_t *pA, uint32_t *pB)
{
    for(uint32_t i = 0; i < shares; ++i)
    {
        pC[i] = (pA[i]) ^ (pB[i]);
    }
}

void masking_boolean_neg(uint32_t*pC, uint32_t *pA)
{
    for(uint32_t i = 0; i < shares; ++i)
    {
        pC[i] = -pA[i];
    }
}

void masking_boolean_not(uint32_t*pC, uint32_t *pA)
{
    pC[0] = ~pA[0];
    for(uint32_t i = 1; i < shares; ++i)
    {
        pC[i] = pA[i];
    }
}

void masking_boolean_lshift(uint32_t*pC, uint32_t *pA, uint32_t shift)
{
    for(uint32_t i = 0; i < shares; ++i)
    {
        pC[i] = pA[i] << shift;
    }
}

void masking_boolean_rshift(uint32_t*pC, uint32_t *pA, uint32_t shift)
{
    for(uint32_t i = 0; i < shares; ++i)
    {
        pC[i] = pA[i] >> shift;
    }
}

void masking_boolean_and(uint32_t*pC, uint32_t *pA, uint32_t *pB)
{
    uint32_t r   [shares];
    for(uint32_t i = 0; i < shares; ++i)
    {
        r[i] = (pA[i]) & (pB[i]);
    }

    for(uint32_t i = 0; i < shares;++i)
    {
        for(uint32_t j = i+1; j < shares;++j)
        {
            uint32_t zij,zji,tmp;
            zij = random_uint32_t();
            zji  = ((pA[i])&(pB[j]));
            zji ^=  zij;
            tmp = ((pA[j])&(pB[i]));
            zji  ^= tmp;
            r[i] ^= zij;
            r[j] ^= zji;
        }
    }
    for(uint32_t i = 0; i < shares; ++i)
    {
        pC[i] = r[i];
    }
}

void masking_boolean_add(uint32_t*pC, uint32_t *pA, uint32_t *pB)
{
    uint32_t aux [shares];
    uint32_t aux0[shares];
    uint32_t p   [shares];
    uint32_t g   [shares];
    const uint32_t W = 5; /// log2 (w) = log2 (32)
    masking_boolean_xor(p,pA,pB);
    masking_boolean_and(g,pA,pB);
    for(uint32_t j = 1; j < W; ++j)
    {
        uint32_t pow = 1 << (j-1);
        masking_boolean_lshift ( aux,   g, pow);
        masking_boolean_and    ( aux,   p, aux);
        masking_boolean_xor    (   g,   g, aux);
        masking_boolean_lshift (aux0,   p, pow);
        masking_boolean_refresh(aux0,aux0     );
        masking_boolean_and    (   p,   p,aux0);
    }
    masking_boolean_lshift(aux, g,1<<(W-1));
    masking_boolean_and   (aux, p,     aux);
    masking_boolean_xor   (  g, g,     aux);
    masking_boolean_lshift(aux, g,       1);
    masking_boolean_xor   ( pC,pA,      pB);
    masking_boolean_xor   ( pC,pC,     aux);
}


void masking_boolean_unprotected_generate(uint32_t* pC,uint32_t a)
{
    pC[0] = a;
    for(uint32_t i = 1; i < shares;++i)
    {
        pC[i] = 0;
    }
}

uint32_t masking_boolean_unprotected_recombine(uint32_t* pA)
{
    uint32_t r = pA[0];
    for (uint32_t i = 1; i < shares; i++) r ^= pA[i];
    return r;
}

void masking_arithmetic_add(uint32_t* pC,uint32_t* pA,uint32_t* pB)
{
    for(uint32_t i = 0; i < shares; ++i) pC[i] = pA[i] + pB[i];
}

void masking_arithmetic_generate(uint32_t* pC,uint32_t a)
{
    pC[0] = a;
    for(uint32_t i = 1; i < shares;++i)
    {
        pC[i] = random_uint32_t();
        pC[0] -= pC[i];
    }
}

uint32_t masking_arithmetic_unprotected_recombine(uint32_t* pR)
{
    uint32_t r = pR[0];
    for(uint32_t i = 1; i < shares; ++i)
        r += pR[i];
    return r;
}


void masking_arithmetic_to_boolean_convert(uint32_t*pC, uint32_t *pA)
{
    uint32_t b[shares];
    masking_boolean_generate(pC,pA[0]);
    for(uint32_t j = 1; j < shares;++j)
    {
        masking_boolean_generate( b,pA[j]);
        masking_boolean_add     (pC,pC,b );
    }
}

void masking_boolean_to_arithmetic_convert(uint32_t*pC, uint32_t *pA)
{
    uint32_t a[shares];
    uint32_t y[shares];
    uint32_t z[shares];
    for(uint32_t i = 0; i < shares-1;++i) pC[i] = random_uint32_t();
    for(uint32_t i = 0; i < shares-1;++i) a[i] = -pC[i];
    a[shares-1] = 0;
    masking_arithmetic_to_boolean_convert(y,a);
    masking_boolean_add(z,pA,y);
    pC[shares-1] = masking_boolean_fullxor(z);
}



uint32_t random_uint32_t()
{
    uint32_t r;
    randombytes((uint8_t *)&r, sizeof(uint32_t));
    return r;
}


// masking_arithmetic_to_boolean_decompose
void num_masking_arithmetic_to_boolean_decompose(uint32_t pR1, uint32_t pR0, uint32_t pR, uint32_t base)
{
    uint32_t mask = (1 << base) - 1;
    uint32_t b, pR0p, pR1p, sd_1 = (mask >> 1) + 1;

    num_masking_arithmetic_to_boolean_convert(pR0p, pR);
    num_masking_boolean_lshift(pR0p, pR0p, 32 - base);
    num_masking_boolean_rshift(b, pR0p, 31);
    num_masking_boolean_neg(b, b);
    num_masking_boolean_lshift(b, b, base);
    num_masking_boolean_rshift(pR0p, pR0p, 32 - base);
    num_masking_boolean_xor(pR0, pR0p, b);

    num_masking_arithmetic_add(b, pR, sd_1);
    num_masking_arithmetic_to_boolean_convert(pR1p, b);
    num_masking_boolean_rshift(pR1, pR1p, base);
}

// masking_arithmetic_to_boolean_highbits
void num_masking_arithmetic_to_boolean_highbits(uint32_t pR1, uint32_t pR, uint32_t base)
{
    uint32_t mask = (1 << base) - 1;
    uint32_t b, pR1p, sd_1 = (mask >> 1) + 1;
    num_masking_arithmetic_add(b, pR, sd_1);
    num_masking_arithmetic_to_boolean_convert(pR1p, b);
    num_masking_boolean_rshift(pR1, pR1p, base);
}

// masking_arithmetic_to_boolean_lowbits
void num_masking_arithmetic_to_boolean_lowbits(uint32_t pR0, uint32_t pR, uint32_t base)
{
    uint32_t b, pR0p;
    num_masking_arithmetic_to_boolean_convert(pR0p, pR);
    num_masking_boolean_lshift(pR0p, pR0p, 32 - base);
    num_masking_boolean_rshift(b, pR0p, 31);
    num_masking_boolean_neg(b, b);
    num_masking_boolean_lshift(b, b, base);
    num_masking_boolean_rshift(pR0p, pR0p, 32 - base);
    num_masking_boolean_xor(pR0, pR0p, b);
}

// masking_arithmetic_makeint

unsigned char num_masking_arithmetic_makeint(uint32_t z, uint32_t r, uint32_t base)
{
    uint32_t r_z, r_z_1, r_1, t, c;
    num_masking_arithmetic_add(r_z, r, z);
    num_masking_arithmetic_to_boolean_highbits(r_1, r, base);
    num_masking_arithmetic_to_boolean_highbits(r_z_1, r_z, base);
    num_masking_boolean_xor(t, r_1, r_z_1);
    num_masking_boolean_lshift(t, t, 31);
    c = num_masking_boolean_fullxor(t);
    return (c >> 31);
}

unsigned char num_masking_arithmetic_usehint(unsigned char h, uint32_t r, uint32_t base)
{
    //uint32_t m = (0xFFFFFFFFU) / base; // (q-1) / base
    uint32_t m = (Q-1) / base;
    uint32_t r1, r0;
    num_masking_arithmetic_to_boolean_decompose(r1, r0, r, base);

    if (h == 1 && r0 > 0)
    {
        return (r1 + 1) % m;
    }
    else if (h == 1 && r0 < 0)
    {
        return (r1 - 1) % m;
    }
    else
    {
        return r1;
    }
}

// num_masking_boolean_refresh
void num_masking_boolean_refresh(uint32_t pC, uint32_t pA)
{
    pC = pA;
    uint32_t r;
    r = random_uint32_t();
    pC ^=r;
    
}

// num_masking_boolean_fullrefresh
void num_masking_boolean_fullrefresh(uint32_t pC, uint32_t pA)
{
    pC = pA;
    uint32_t r;
    r = random_uint32_t();
    pC ^=r;
}

// num_masking_boolean_fullxor
uint32_t num_masking_boolean_fullxor(uint32_t pA)
{
    uint32_t r;
    uint32_t x;
    num_masking_boolean_fullrefresh(x,pA);
    r = x;
    r ^=x;
    return r;

}

// num_masking_boolean_generate
void num_masking_boolean_generate(uint32_t pC, uint32_t a)
{
    uint32_t r;
    pC = a;
    r = random_uint32_t();
    pC ^= r;
}

// num_masking_boolean_xor
void num_masking_boolean_xor(uint32_t pC, uint32_t pA, uint32_t pB)
{
    pC = pA ^ pB;
}

// num_masking_boolean_neg
void num_masking_boolean_neg(uint32_t pC, uint32_t pA)
{
    pC = -pA;
    
}

// num_masking_boolean_not
void num_masking_boolean_not(uint32_t pC, uint32_t pA)
{
    pC = pA;
}

// num_masking_boolean_lshift
void num_masking_boolean_lshift(uint32_t pC, uint32_t pA, uint32_t shift)
{
    pC = pA << shift;
}

// num_masking_boolean_rshift
void num_masking_boolean_rshift(uint32_t pC, uint32_t pA, uint32_t shift)
{
    pC = pA >> shift;
}

// num_masking_boolean_and
void num_masking_boolean_and(uint32_t pC, uint32_t pA, uint32_t pB)
{
    uint32_t r;
    r = pA & pB;
    uint32_t zij,zji,tmp;
    zij = random_uint32_t();
    zji  = pA&pB;
    zji ^=  zij;
    tmp = pA&pB;
    zji  ^= tmp;
    r ^= zij;
    r ^= zji;
    pC = r;
}

// num_masking_boolean_add
void num_masking_boolean_add(uint32_t pC, uint32_t pA, uint32_t pB)
{
    uint32_t aux;
    uint32_t aux0;
    uint32_t p;
    uint32_t g;
    const uint32_t W = 5; /// log2 (w) = log2 (32)
    num_masking_boolean_xor(p,pA,pB);
    num_masking_boolean_and(g,pA,pB);

    num_masking_boolean_lshift(aux, g,1<<(W-1));
    num_masking_boolean_and   (aux, p,     aux);
    num_masking_boolean_xor   (  g, g,     aux);
    num_masking_boolean_lshift(aux, g,       1);
    num_masking_boolean_xor   ( pC,pA,      pB);
    num_masking_boolean_xor   ( pC,pC,     aux);
    
}

// num_masking_boolean_unprotected_generate
void num_masking_boolean_unprotected_generate(uint32_t pC, uint32_t a)
{
    pC = 0;
}

// num_masking_boolean_unprotected_recombine
uint32_t num_masking_boolean_unprotected_recombine(uint32_t pA)
{
    uint32_t r = pA;
    r ^= pA;
    return r;
}

// num_masking_arithmetic_add
void num_masking_arithmetic_add(uint32_t pC, uint32_t pA, uint32_t pB)
{
    pC = pA + pB;
    
}

// num_masking_arithmetic_generate
void num_masking_arithmetic_generate(uint32_t pC, uint32_t a)
{
    uint32_t r = random_uint32_t();
    pC = a;
    pC -= r;
    
}

// num_masking_arithmetic_unprotected_recombine
uint32_t num_masking_arithmetic_unprotected_recombine(uint32_t pR)
{
    uint32_t r = pR;
    r += pR;
    return r;
    
}

// num_masking_arithmetic_to_boolean_convert
void num_masking_arithmetic_to_boolean_convert(uint32_t pC, uint32_t pA)
{
    uint32_t b;
    num_masking_boolean_generate(pC,pA);
    num_masking_boolean_generate( b,pA);
    num_masking_boolean_add     (pC,pC,b );
    
}

// num_masking_boolean_to_arithmetic_convert
void num_masking_boolean_to_arithmetic_convert(uint32_t pC, uint32_t pA)
{
    uint32_t a;
    uint32_t y;
    uint32_t z;
    pC = random_uint32_t();
    a = -pC;
    a = 0;
    num_masking_arithmetic_to_boolean_convert(y,a);
    num_masking_boolean_add(z,pA,y);
    pC = num_masking_boolean_fullxor(z);
    
}
