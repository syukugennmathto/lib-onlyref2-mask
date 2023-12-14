#include <stdio.h>
#include <stdint.h>
#include "params.h"
#include "rounding.h"
#include "random.h"
#include "randombytes.h"

void printBinary(unsigned int num) {
    if (num > 1) {
        printBinary(num >> 1);
    }
    printf("%d", num & 1);
}


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
/*
#if GAMMA2 == (Q-1)/32
    a1_unsigned = (a1_unsigned * 1025 + (1 << 21)) >> 22;
    a1_unsigned &= 15;
#elif GAMMA2 == (Q-1)/88
    */
    a1_unsigned = (a1_unsigned * 11275 + (1 << 23)) >> 24;
    a1_unsigned ^= ((43 - a1_unsigned) >> 31) & a1_unsigned;
    a1_unsigned = a1_unsigned%44;
//#endif

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

    //alpha change D
    highbits(r1, a0_unsigned,alpha);
    highbits(v1, a0_unsigned + a1_unsigned, alpha);

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
        // m = (q-1)/alpha = (q-1)/2*gamma2 = 88*(q-1)/2*(q-1)= 44
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

uint32_t uuse_hint(uint32_t r, unsigned int hint,uint32_t base){
    uint32_t r0, r1;
    int32_t a0, a1,a;
    a = (int32_t)r;
    uint32_t m = (Q-1)/base;
    //r1 = decompose(&r0, r);
    printf("mask is %u\n",m);
    printf("\n");
    prime_decompose(&r1,&r0,r);
    printf("high is %u\n",r1);
    printf("lows is %u\n",r0);
    a1 = decompose(&a0,a);
    printf("\n");
    printf("high is %u\n",a1);
    
    if (hint == 1)
    {
        if (r0 > 0)
        {
            if (r1 > 0){
                return (r1 + 1) % m;
            }
            else{
                return (-1 * r1 + 1) % m;
            }
        }
        else
        {
            if (r1 > 0){
            return (r1 - 1) % m;
        }
        else{
            return (-1 * r1 - 1) % m;
        }
        }
    }
    return r1;
 
}
void pt_decompose(uint32_t *r1, uint32_t *r0, uint32_t  r , uint32_t base)
{
    uint32_t mask = (1 << base) - 1;
    uint32_t d_1 = (mask>>1)+1;
    uint32_t r0p,b;
    r0p = (r<<(32-base));
    b   = (-(r0p>>(31)));
    b   = b<<base;
    r0p = (r0p>>(32-base));

    *r0 = r0p^b;
    *r1  = ((r+d_1)>>base);
}
void pt_highbits(uint32_t *r1, uint32_t  r , uint32_t base)
{
    uint32_t mask = (1 << base) - 1;
    uint32_t d_1 = (mask>>1)+1;
    uint32_t r1p;
    r1p  = ((r+d_1)>>base);
    *r1 = r1p;
}

unsigned char pt_makehint(uint32_t  z , uint32_t r, uint32_t base)
{
    uint32_t r1,v1;
    pt_highbits(&r1, r  , base);
    pt_highbits(&v1, r+z, base);
    return (r1==v1)?0:1;
}

uint32_t pt_usehint(uint32_t  y , uint32_t r, uint32_t base)
{
    uint32_t r0,r1;
    pt_decompose(&r1,&r0,r,base);
    if( y == 1 )
    {
        if( r0&(1<<(base-1)) )
        {
            --r1;
        }else
        {
            ++r1;
        }
    }
    return r1;
}


void arithmetic_to_boolean_use_hint(uint32_t *r, unsigned int hint,uint32_t base){
    uint32_t r0[shares], r1[shares];
    uint32_t mask = 44;
    //uint32_t mask = (1 << base) - 1;
    //r1 = decompose(&r0, r);
    //printf("mask is %u\n",m);
    //printf("\n");
    masking_arithmetic_to_boolean_decompose(r1,r0,r,base);
        for(int i = 0;i<shares;i++)
        {
            if (hint == 1)
            {
                if (r0[i]&(1<<(base-1)))
                {
                    r1[i] = (r1[i] == 43) ?  0 : r1[i] + 1;
                }
                else
                {
                    r1[i] = (r1[i] ==  0) ? 43 : r1[i] - 1;
                }
            }
         else
         {
             r[i] = r1[i];
         }
        }
    

    for(int i = 0;i<shares;i++){
        printf("use_high_deco%d is %d \n",i,r1[i]);
        printf("use_lows_deco%d is %d \n",i,r0[i]);
        printf("use_ring_deco%d is %d \n",i,r[i]);
        printf("\n");
    }
}


uint32_t hhighbbits(uint32_t output, uint32_t input){
    uint32_t r0, r1;
    prime_decompose(&r1,&r0,input);
    printf("high is %u\n",r1);
    printf("lows is %u\n",r0);
    output = r1;
    return output;
}



/*************************************************
* Name:        highbitster
*
* Description: Compute the high bits of a finite field element `r` given a base value.
*              The base value determines the number of bits to shift the element `r` to the right.
*
* Arguments:   - uint32_t *r1: pointer to the output high bits of `r`
*              - uint32_t r: input element
*              - uint32_t base: base value used for shifting `r`
**************************************************/
void highbits(uint32_t output, uint32_t input, uint32_t base){
    uint32_t mask = (1 << D) - 1;
    uint32_t d_1 = (mask >> 1) + 1;
    output = ((input + d_1) >> base);
    //printf("output is %u\n",output);
}
/*************************************************
* Name:        highbitster
*
* Description: Compute the high bits of a finite field element `r` given a base value.
*              The base value determines the number of bits to shift the element `r` to the right.
*
* Arguments:   - uint32_t *r1: pointer to the output high bits of `r`
*              - uint32_t r: input element
*              - uint32_t base: base value used for shifting `r`
**************************************************/
int nos_highbits(int32_t *output0, int32_t input){
    int32_t output1;
    int32_t mask = (1 << 22) - 1;
    int32_t d_1 = (mask >> 1) + 1;
    //printf("intput is %u\n",(input + d_1));
    output1 = ((input + d_1) >> badominton);
    //printf("output is %d\n",output);

    return output1;
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
    printf("deco_mask is %u\n",mask);
    uint32_t b[shares],pR0p[shares],pR1p[shares],sd_1[shares] = {(mask>>1)+1,0};
    //lowbits
    masking_arithmetic_to_boolean_convert(pR0p,pR);
    masking_boolean_lshift(pR0p,pR0p,32-base);
    masking_boolean_rshift(b,pR0p,31);
    masking_boolean_neg(b,b);
    masking_boolean_lshift(b,b,base);
    masking_boolean_rshift(pR0p,pR0p,32-base);
    masking_boolean_xor(pR0,pR0p,b);

    //highbits
    masking_arithmetic_add(b,pR,sd_1);
    for(int i = 0;i<shares;i++){
        printf("b[%d] is %u\n",i,b[i]);
    }
    printf("\n");
    masking_arithmetic_to_boolean_convert(pR1p,b);
    for(int i = 0;i<shares;i++){
        printf("pR1p[%d] is %u\n",i,pR1p[i]);
    }
    printf("\n");
    masking_boolean_rshift(pR1,pR1p,base);
    for(int i = 0;i<shares;i++){
        printf("pR1p[%d] is %u\n",i,pR1p[i]);
    }
    printf("\n");

}

void masking_arithmetic_to_boolean_highbits(uint32_t* pR1, uint32_t* pR, uint32_t base)
{
    uint32_t mask = (1 << base) - 1;
    printf("high_mask is %u\n",mask);
    uint32_t b[shares],pR1p[shares],sd_1[shares] = {(mask>>1)+1,0};
    
    masking_arithmetic_add(b,pR,sd_1);
    for(int i = 0;i<shares;i++){
        printf("b[%d] is %d\n",i,b[i]);
    }
    printf("\n");
    masking_arithmetic_to_boolean_convert(pR1p,b);
    for(int i = 0;i<shares;i++){
        printf("pR1p[%d] is %d\n",i,pR1p[i]);
    }
    printf("\n");
    
    masking_boolean_rshift(pR1,pR1p,base);
    for(int i = 0;i<shares;i++){
        printf("pR1p[%d] is %d\n",i,pR1p[i]);
    }
    printf("\n");
}

void convers_arithmetic_to_boolean_highbits(uint32_t* pR1, uint32_t* pR, uint32_t base)
{
    uint32_t mask = (1 << base) - 1;
    printf("high_mask is %u\n",mask);
    uint32_t b[shares],pR1p[shares],sd_1[shares] = {(mask>>1)+1,0};
    
    masking_boolean_lshift(pR1,pR1p,base);
    for(int i = 0;i<shares;i++){
        printf("pR1p[%d] is %d\n",i,pR1p[i]);
    }
    printf("\n");
    
    masking_arithmetic_to_boolean_convert(pR1p,b);
    for(int i = 0;i<shares;i++){
        printf("pR1p[%d] is %d\n",i,pR1p[i]);
    }
    printf("\n");
    
    masking_arithmetic_add(b,pR,sd_1);
    for(int i = 0;i<shares;i++){
        printf("b[%d] is %d\n",i,b[i]);
    }
    printf("\n");


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

/*algorithm 9maskde algorithm of makehintQ with a prime mudulus Q. */
unsigned char masking_arithmetic_makeint(uint32_t* z, uint32_t* r, uint32_t base)
{
    uint32_t r_z[shares],r_z_1[shares],r_1[shares],t[shares],c;

    masking_arithmetic_to_boolean_highbits(r_1,r,base);
    masking_arithmetic_add(r_z,r,z);
    masking_arithmetic_to_boolean_highbits(r_z_1,r_z,base);

    masking_boolean_xor(t,r_1,r_z_1);

    masking_boolean_lshift(t,t,base);
 
    c = masking_boolean_fullxor(t);
    
    return (c >> 31);
}

unsigned char masking_arithmetic_two_makeint(uint32_t* z, uint32_t* r, uint32_t base)
{
    uint32_t r_z[shares],r_z_1[shares],r_1[shares],t[shares],c;
    masking_arithmetic_to_boolean_highbits(r_1,r,base); //1
    masking_arithmetic_add(r_z,r_1,z); //2
    masking_arithmetic_to_boolean_highbits(r_z_1,r_z,base); //3
    masking_boolean_xor(t,r_1,r_z_1); //4
    for(int i = 0;i<shares;i++){
        printf("t[%d] is %d\n",i,t[i]);
    }
    printf("\n");
    masking_boolean_lshift(t,t,base-1); //5
    for(int i = 0;i<shares;i++){
        printf("t[%d] is %d\n",i,t[i]);
    }
    printf("\n");
    c = masking_boolean_fullxor(t); //6

    return c; //7
}

void new_masking_arithmetic_to_boolean_makehint(int32_t* c,const int32_t* z,const int32_t* r,mask_point pd)
{
    int32_t r_z[N],r_z_1[N],r_1[N],t[N];
    //mask_point pd;
    
    //masking_arithmetic_to_boolean_highbits(r_1,r,base); //1
    new_nosking_arithmetic_to_boolean_highbits(r_1,r,pd);
    
    //masking_arithmetic_add(r_z,r_1,z); //2
    for(uint32_t i = 0; i < N; ++i)
    {
        r_z[i] = (r_1[i]+z[i])%Q;
    }
    //masking_arithmetic_to_boolean_highbits(r_z_1,r_z,base); //3
    new_nosking_arithmetic_to_boolean_highbits(r_z_1,r_z,pd);
    
    //masking_boolean_xor(t,r_1,r_z_1); //4
    for(uint32_t i = 0; i < N; ++i)
    {
        t[i] = (r_1[i]) ^ (r_z_1[i]);
    }

    
    //masking_boolean_lshift(t,t,base-1); //5
    for(uint32_t i = 0; i < N; ++i)
    {
        t[i] = t[i] << (D-1);

    }
    
    //c = masking_boolean_fullxor(t); //6

    c[0] = t[0];
    for(uint32_t j = 0; j < N; ++j)
    {
        c[j] =c[j] ^ t[j];
    }

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

unsigned char num_masking_arithmetic_makehint(uint32_t z, uint32_t r, uint32_t base)
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
    uint32_t m = (Q-1) / base; // (q-1) / base
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


/*************************************************
* Name:        use_hint
*
* Description: Correct high bits according to hint.
*
* Arguments:   - int32_t a: input element
*              - unsigned int hint: hint bit
*
* Returns corrected high bits.
*************************************************
uint32_t masking_use_hint(uint32_t a, unsigned char hint,uint32_t base) {
    uint32_t a0, a1;
    masking_arithmetic_to_boolean_decompose(&a1,&a0,&a,base);

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
}*/


/*void prime_decompose(uint32_t *r1, uint32_t *r0, uint32_t r)
{
    uint32_t q     = 8380417;
    uint32_t alpha = 523776;
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

}*/

void prime_lowbits(uint32_t *r0, uint32_t r)
{
    //uint32_t q     = 8380417;
    //uint32_t alpha = 523776;
    uint32_t r1p,r0p,m;

    r0p  = (r & 0x7FFFF) + ((r>>19)<<9);
    r0p -= alpha/2 + 1;
    r0p += (-(r0p >> 31)) & alpha;
    r0p -= alpha/2 - 1;

    r1p  = r - r0p;
    r1p  = (r1p >> 19) + 1 - ((r1p-1)>>31);

    m    = ~(-(r1p>>4));
    r0p  = r0p + (m&1) - 1;
    *r0  = r0p;
}

void prime_highbits(uint32_t *r1, uint32_t r)
{
    //uint32_t q     = 8380417;
    //uint32_t alpha = 523776;
    uint32_t r1p,r0p,m;

    r0p  = (r & 0x7FFFF) + ((r>>19)<<9);
    r0p -= alpha/2 + 1;
    r0p += (-(r0p >> 31)) & alpha;
    r0p -= alpha/2 - 1;

    r1p  = r - r0p;
    r1p  = (r1p >> 19) + 1 - ((r1p-1)>>31);

    m    = ~(-(r1p>>4));
    r1p  = m&r1p;
    *r1  = r1p;
}

void masking_arithmetic_to_boolean_prime_decompose(uint32_t* pR1,uint32_t* pR0, uint32_t * pR,uint32_t aalpha)
{
    //const uint32_t q     = 8380417;
    //const uint32_t aalpha = 523776;

    uint32_t pRp[shares],pR0p[shares],pR1p[shares];
    uint32_t m             [shares] = {0x7FFFF       ,0}; //二進数で11111111111111111111111 1
    uint32_t cALPHA        [shares] = {  aalpha       ,0}; // 2
    uint32_t cALPHA_BOUND_1[shares] = {-(aalpha/2 + 1),0}; // 3
    uint32_t cALPHA_BOUND_2[shares] = {-(aalpha/2 - 1),0}; // 4
    uint32_t cQ_MINUS_ONE  [shares] = {Q - 1         ,0}; // 5
    uint32_t cONE          [shares] = {1             ,0}; // 6


    /// Modular reduction by alpha ///
    masking_arithmetic_to_boolean_convert(pRp,pR); //7
    masking_boolean_and   (pR0p,pRp,m); //8
    masking_boolean_rshift(m,pRp,19); //9
    masking_boolean_lshift(m,m,9); //10
    masking_boolean_add   (pR0p,m,pR0p); //11
    masking_boolean_add   (pR0p,pR0p,cALPHA_BOUND_1); //12
    masking_boolean_rshift(m,pR0p,31); //13
    masking_boolean_neg   (m,m); //14
    masking_boolean_and   (m,m,cALPHA); //15
    masking_boolean_add   (pR0p,pR0p,m); //16
    masking_boolean_add   (pR0p,pR0p,cALPHA_BOUND_2); //17

    /// Computation of r1 ///

    masking_boolean_not   (pR1p,pR0p);
    masking_boolean_add   (pR1p,pR1p,pRp);
    masking_boolean_rshift(m,pR1p,31);
    masking_boolean_neg   (m,m);
    masking_boolean_add   (pR1p,pR1p,cONE);
    masking_boolean_rshift(pR1p,pR1p,19);
    masking_boolean_not   (m,m);
    masking_boolean_lshift(m,m,31);
    masking_boolean_rshift(m,m,31);
    masking_boolean_add   (pR1p,pR1p,m);

    /// distance test and correction if required ///
    masking_boolean_lshift(m,pR1p,27);
    masking_boolean_rshift(m,m,31);
    masking_boolean_neg   (m,m);
    masking_boolean_not   (m,m);
    masking_boolean_and   (pR1,pR1p,m);

    masking_boolean_lshift(m,m,31);
    masking_boolean_rshift(m,m,31);
    masking_boolean_add   (pR0p,pR0p,m);
    masking_boolean_add   (pR0,pR0p,cQ_MINUS_ONE);
}

void masking_arithmetic_to_boolean_prime_lowbits(uint32_t* pR0, uint32_t * pR)
{
    //const uint32_t q     = 8380417;
    //const uint32_t alpha = 523776;

    uint32_t pRp[shares],pR0p[shares],pR1p[shares];
    uint32_t m             [shares] = {0x7FFFF       ,0};
    uint32_t cONE          [shares] = {1             ,0};
    uint32_t cALPHA        [shares] = {  alpha       ,0};
    uint32_t cALPHA_BOUND_1[shares] = {-(alpha/2 + 1),0};
    uint32_t cALPHA_BOUND_2[shares] = {-(alpha/2 - 1),0};
    uint32_t cQ_MINUS_ONE  [shares] = {Q - 1         ,0};


    /// Modular reduction by alpha ///
    masking_arithmetic_to_boolean_convert(pRp,pR);
    masking_boolean_and   (pR0p,pRp,m);
    masking_boolean_rshift(m,pRp,19);
    masking_boolean_lshift(m,m,9);
    masking_boolean_add   (pR0p,m,pR0p);
    masking_boolean_add   (pR0p,pR0p,cALPHA_BOUND_1);
    masking_boolean_rshift(m,pR0p,31);
    masking_boolean_neg   (m,m);
    masking_boolean_and   (m,m,cALPHA);
    masking_boolean_add   (pR0p,pR0p,m);
    masking_boolean_add   (pR0p,pR0p,cALPHA_BOUND_2);

    /// Computation of r1 ///

    masking_boolean_not   (pR1p,pR0p);
    masking_boolean_add   (pR1p,pR1p,pRp);
    masking_boolean_rshift(m,pR1p,31);
    masking_boolean_neg   (m,m);
    masking_boolean_add   (pR1p,pR1p,cONE);
    masking_boolean_rshift(pR1p,pR1p,19);
    masking_boolean_not   (m,m);
    masking_boolean_lshift(m,m,31);
    masking_boolean_rshift(m,m,31);
    masking_boolean_add   (pR1p,pR1p,m);

    /// distance test and correction if required ///
    masking_boolean_lshift(m,pR1p,27);
    masking_boolean_rshift(m,m,31);
    masking_boolean_neg   (m,m);
    masking_boolean_not   (m,m);

    masking_boolean_lshift(m,m,31);
    masking_boolean_rshift(m,m,31);
    masking_boolean_add   (pR0p,pR0p,m);
    masking_boolean_add   (pR0,pR0p,cQ_MINUS_ONE);
}

void masking_arithmetic_to_boolean_prime_highbits(uint32_t* pR1, uint32_t * pR,uint32_t aalpha)
{
    //const uint32_t q     = 8380417;
    //const uint32_t aalpha = 523776;

    uint32_t pRp[shares],pR0p[shares],pR1p[shares];
    uint32_t m             [shares] = {0x7FFFF       ,0};
    uint32_t cONE          [shares] = {1             ,0};
    uint32_t cALPHA        [shares] = {  aalpha       ,0};
    uint32_t cALPHA_BOUND_1[shares] = {-(aalpha/2 + 1),0};
    uint32_t cALPHA_BOUND_2[shares] = {-(aalpha/2 - 1),0};
    uint32_t cQ_MINUS_ONE  [shares] = {Q - 1         ,0};


    /// Modular reduction by alpha ///
    printf("masking_arithmetic_to_boolean_convert\n");
    masking_arithmetic_to_boolean_convert(pRp,pR);
    for(int i = 0;i<shares;i++){
        printf("pRp%d is %d \n",i,pRp[i]);
    }
    printf("\n");
    printf("masking_boolean_and\n");
    masking_boolean_and   (pR0p,pRp,m);
    for(int i = 0;i<shares;i++){
        printf("pR0p%d is %d \n",i,pR0p[i]);
    }
    printf("\n");
    
    printf("masking_boolean_rshift\n");
    masking_boolean_rshift(m,pRp,19);
    for(int i = 0;i<shares;i++){
        printf("m%d is %d \n",i,m[i]);
    }
    printf("\n");
    
    printf("masking_boolean_lshift\n");
    masking_boolean_lshift(m,m,9);
    for(int i = 0;i<shares;i++){
        printf("m%d is %d \n",i,m[i]);
    }
    printf("\n");
    
    printf("masking_boolean_add\n");
    masking_boolean_add   (pR0p,m,pR0p);
    for(int i = 0;i<shares;i++){
        printf("pR0p%d is %d \n",i,pR0p[i]);
    }
    printf("\n");
    
    printf("masking_boolean_add\n");
    masking_boolean_add   (pR0p,pR0p,cALPHA_BOUND_1);
    for(int i = 0;i<shares;i++){
        printf("pR0p%d is %d \n",i,pR0p[i]);
    }
    printf("\n");
    
    printf("masking_boolean_rshift\n");
    masking_boolean_rshift(m,pR0p,31);
    for(int i = 0;i<shares;i++){
        printf("m%d is %d \n",i,m[i]);
    }
    printf("\n");
    
    printf("masking_boolean_neg\n");
    masking_boolean_neg   (m,m);
    for(int i = 0;i<shares;i++){
        printf("m%d is %d \n",i,m[i]);
    }
    printf("\n");
    
    printf("masking_boolean_add\n");
    masking_boolean_and   (m,m,cALPHA);
    for(int i = 0;i<shares;i++){
        printf("m%d is %d \n",i,m[i]);
    }
    printf("\n");
    
    printf("masking_boolean_add\n");
    masking_boolean_add   (pR0p,pR0p,m);
    for(int i = 0;i<shares;i++){
        printf("pR0p%d is %d \n",i,pR0p[i]);
    }
    
    printf("masking_boolean_add\n");
    
    printf("alpha is %d \n",aalpha);
    masking_boolean_add   (pR0p,pR0p,cALPHA_BOUND_2);
    for(int i = 0;i<shares;i++){
        printf("pR0p%d is %d \n",i,pR0p[i]);
    }
    printf("\n");
    /// Computation of r1 ///
    printf("alpha is %d \n",aalpha);
    printf("Computation of r1 \n");
    masking_boolean_not   (pR1p,pR0p);
    masking_boolean_add   (pR1p,pR1p,pRp);
    masking_boolean_rshift(m,pR1p,31);
    masking_boolean_neg   (m,m);
    masking_boolean_add   (pR1p,pR1p,cONE);
    masking_boolean_rshift(pR1p,pR1p,19);
    masking_boolean_not   (m,m);
    masking_boolean_lshift(m,m,31);
    masking_boolean_rshift(m,m,31);
    masking_boolean_add   (pR1p,pR1p,m);
    for(int i = 0;i<shares;i++){
        printf("pR1p%d is %d \n",i,pR1p[i]);
    }
    printf("\n");
    /// distance test and correction if required ///
    
    printf("distance test and correction if required \n");
    masking_boolean_lshift(m,pR1p,27);
    masking_boolean_rshift(m,m,31);
    masking_boolean_neg   (m,m);
    masking_boolean_not   (m,m);
    masking_boolean_and   (pR1,pR1p,m);
    for(int i = 0;i<shares;i++){
        printf("pR1p%d is %d \n",i,pR1p[i]);
    }
    printf("\n");
}

uint32_t random_uint32_t_mono(uint32_t r)
{
    randombytes((uint8_t *)&r, sizeof(uint32_t));
    //if(r <0){r = (-1)*r;}
    r = r%(1<<12);
    return r;
}

int32_t random_int32_t()
{
    int32_t r;
    randombytes((uint8_t *)&r, sizeof(int32_t));
    //if(r <0){r = (-1)*r;}
    r = r%(1<<16);
    return r;
}


int32_t  new_masking_decompose(int32_t pR0,const int32_t pR) {
    
    int32_t b=0,pR1p=0;
    int32_t pR1=0;
    int32_t pR1pt=0,pR2=0,pR3=0;
    int32_t mask = ((1 << 22) - 1);
    int32_t d = (mask>>1)+1;
    pR2 = (pR + d);
    pR3 ^=((43-pR2)>>25)&43;
    return pR3;
    
/*
 
 pR1 = (pR+127)>>7;
 //printf("pR1= %d\n",pR1);
 pR2 = (pR1*11275 +(1<<23))>>24;
 //printf("pR2= %d\n",pR2);
 pR3 ^=((43-pR2)>>31)&43;
// printf("pR3= %d\n",pR3);
 return pR3;
 
 
#if GAMMA2 == (Q-1)/32
    a1_unsigned = (a1_unsigned * 1025 + (1 << 21)) >> 22;
    a1_unsigned &= 15;
#elif GAMMA2 == (Q-1)/88
    
    a1_unsigned = (a1_unsigned * 11275 + (1 << 23)) >> 24;
    a1_unsigned ^= ((43 - a1_unsigned) >> 31) & a1_unsigned;
    //a1_unsigned = a1_unsigned%44;
//#endif

    a0_unsigned = a_unsigned - a1_unsigned * 2 * GAMMA2;
    a0_unsigned -= (((Q - 1) / 2 - a0_unsigned) >> 31) & Q;

    // Convert back to signed 32-bit integer format
    *a0 = (int32_t)a0_unsigned;

    // Convert a1 back to signed 32-bit integer format and return it
 */

}


void new_masking_arithmetic_to_boolean_highbits(int32_t* pR2,const int32_t* pR,mask_point* pD)
{
    int32_t based =alpha;
    int32_t mask = ((1 << 22) - 1);
    int32_t b[N]={0,0},pR1p[N]={0,0},pR1pt[N]={0,0},pR1[N]={0,0},pR11[N]={0,0},sd_1[N]={0,0};
    sd_1[0] = (mask>>1)+1;

    //printf("b[2] = %d\n",b[2]);
    
    for(uint32_t i = 0; i < N; ++i){b[i] = (pR[i] + sd_1[i]);}
    //start Masking function
    int32_t d[N]={0,0};
    
    //convert
    pR1p[0] = b[0];
    for(uint32_t i = 1; i < N;++i)
    {
        pD->outside[i] = random_int32_t();
        pR1p[i] = pD->outside[i]; // 同様にpDを逆参照
        pR1p[0] ^= pR1p[i];
    }
    //convert
    for(uint32_t j = 1; j < N;++j)
    {
        d[0] = b[j];
        for(uint32_t i = 1; i < N;++i)
        {
            pD->inside[j][i] = random_int32_t();
            d[i] = pD->inside[j][i]; // 同様にpDを逆参照
            d[0] ^= d[i];
        }
        pR1p[j]=pR1p[j]+d[j];
    }
    //printf("d[%d] = %d\n",3,d[3]);

    for(uint32_t i = 0; i < N; ++i){pR2[i]=pR1p[i]>>25;}
    
   //for(uint32_t i = 0; i < N; ++i){pR2[i]=new_masking_decompose(pR1[i],pR11[i]);}
    //for(uint32_t i = 0; i < N; ++i){printf("出力値____pR[%d] = %d\n",i,pR2[i]);}
}
void new_unmasking_arithmetic_to_boolean_highbits(int32_t* pR,const int32_t* pR1,mask_point pD)
{
    int32_t based =alpha;
    int32_t mask = ((1 << 22) - 1);

    int32_t b[N]={0,0},pR1p[N],sd_1[N]={0,0},d[N]={0,0},pR1pt[N]={0,0},pR2[N]={0,0},pR3[N]={0,0},pR11[N]={0,0},pR4[N]={0,0};
    sd_1[0] = (mask>>1)+1;

    for(uint32_t i = 1; i < N; ++i){b[i] = (pR1[i] + sd_1[i]);}
    //convert
    pR1p[0]=b[0];
    for(uint32_t i = 1; i < N;++i)
    {
        pR1p[i] = pD.outside[i];
        pR1p[0] ^= pR1p[i];
    }
    for(int32_t j = 1; j < N;++j)
    {
        d[0] = b[j];
        for(uint32_t i = 1; i < N;++i)
        {
            d[i] = pD.inside[j][i];
            d[0] ^= d[i];
        }
        pR1p[j]=pR1p[j]+d[j];
    }
    for(uint32_t i = 0; i < N; ++i){pR[i] ^= ((43-pR1p[i])>>25)&43;}
    //for(uint32_t i = 0; i < N; ++i){pR2[i]=pR1p[i]>>base;}
    //for(uint32_t i = 0; i < N; ++i){pR[i]=new_masking_decompose(pR3[i],pR2[i]);}
    
    //for(uint32_t i = 0; i < N; ++i)pR11[i] = (pR1pt[i]/11275 -(1<<23))<<24;
    //for(uint32_t i = 0; i < N; ++i)pR2[i] = (pR1[i]-127)<<7;
    //for(uint32_t i = 0; i < N; ++i)pR[i] ^=((43-pR2[i])>>31)&43;
    
    
    //for(uint32_t i = 0; i < N; ++i){pR[i]=new_masking_decompose(pR2[i],pR1pt[i]);}
}

void new_nosking_arithmetic_to_boolean_highbits(int32_t* pR1,const int32_t* pR,mask_point pD)
{
    int32_t based =alpha;
    int32_t mask = ((1 << 23) - 1);
    int32_t b[N]={0,0},pR1p[N]={0,0},sd_1[N],d[N]={0,0},pR1pt[N]={0,0},pR11[N]={0,0};
    for(uint32_t i = 0; i < N-2; ++i) sd_1[i] = (mask>>1)+1;
    sd_1[N-1]=0;
    for(uint32_t i = 0; i < N; ++i) {b[i] = (pR[i] + sd_1[i]);}
   // printf("b[2] = %d\n",b[2]);
    pR1p[0] = b[0];
    
    for(uint32_t i = 0; i < N; ++i){pR1[i]=new_masking_decompose(pR1pt[i],pR[i]);}
    
    //Masking function
    for(uint32_t i = 1; i < N;++i)
    {
        //random_uint32_t()
        pR1p[i] = pD.outside[i];
        pR1p[0] -= pR1p[i];
    }
    
    for(int32_t j = 1; j < N;++j)
    {
        d[0] = b[j];
        for(uint32_t i = 1; i < N;++i)
        {
            
            //random_uint32_t()
            d[i] = pD.inside[j][i];
            d[0]= d[i];
        }
        
        
    }
    
    for(uint32_t i = 0; i < N; ++i){pR1[i]=(pR1p[i]^d[i]);}
    
    
    

    //Masking function


    //printf("6pR1[2] = %d\n",pR1[2]);
}




void new_masking_arithmetic_to_boolean_lowbits(int32_t* pR0,const int32_t* pR,mask_point pD)
{
    int32_t mask = (1 << 22) - 1;
    int32_t b[N],pR0p[N]={0,0},sd_1[N] = {(mask>>1)+1,0};
     
     //masking_arithmetic_to_boolean_convert(pR0p,pR);
    int32_t d[N]={0,0};
     pR0p[0] = pR[0];

     for(uint32_t j = 1; j < N;++j)
     {
         d[0] = pR[j];
         for(uint32_t i = 1; i < N;++i)
         {
             
             //pD.e[j][i] = 16; //random_uint32_t()
             d[i] = pD.inside[j][i];
             d[0] -= d[i];
         }
         
         for(uint32_t i = 0; i < N; ++i){pR0p[i]=pR0p[i]+d[i];
         }
     }
     
    //masking_boolean_lshift(pR0p,pR0p,32-base);
    for(uint32_t i = 0; i < N; ++i)
    {
        pR0p[i] = pR0p[i] << (32-D);
    }
    
    //masking_boolean_rshift(b,pR0p,31);
    for(uint32_t i = 0; i < N; ++i)
    {
        b[i] = pR0p[i] >> 31;
    }
    
    //masking_boolean_neg   (b,b);
    for(uint32_t i = 0; i < N; ++i)
    {
        b[i] = -b[i];
    }
    
    //masking_boolean_lshift(b,b,base);
    for(uint32_t i = 0; i < N; ++i)
    {
        b[i] = b[i] << D;
    }
    
    //masking_boolean_rshift(pR0p,pR0p,32-base);
    for(uint32_t i = 0; i < N; ++i)
    {
        pR0p[i] = pR0p[i] >> (32-D);
    }
    
    //masking_boolean_xor   (pR0,pR0p,b);
    for(uint32_t i = 0; i < N; ++i)
    {
        pR0[i] = (pR0p[i]) ^ (b[i]);
    }
    
    for(uint32_t i = 0; i < N; ++i)
    {
        pR0[i] = pR0[i]%Q;
    }
    
}


void new_masking_use_hint(int32_t* pR2,const int32_t* hint,mask_point pD) {
    
    int32_t pR0p[N],pR2p[N];
    
    new_nosking_arithmetic_to_boolean_highbits(pR2p,pR2,pD);
    new_masking_arithmetic_to_boolean_lowbits(pR0p,pR2,pD);
    int32_t m = 44;
    
    for(uint32_t i = 0; i < N;++i)
    {
        if(hint[i]==0)
        {
            pR2[i]= pR2p[i];

        }
        else
        {
            if(pR0p[i]>0){pR2[i] = (pR2p[i] + 1)%m;}
            else{pR2[i] = (pR2p[i] - 1)%m;}
        }
        
    }

}
