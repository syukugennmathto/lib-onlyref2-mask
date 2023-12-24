#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <oqs/oqs.h>
#include <stdint.h>
#include <math.h>
#include "usehinttest.h"
#include "params.h"
#include "rounding.h"
#include "sign.h"
#include "packing.h"
#include "polyvec.h"
#include "poly.h"
#include "randombytes.h"
#include "symmetric.h"
#include "fips202.h"
#include "noise.h"
#include "binary_int32.h"
#include "binary_int32_unmasking.h"
#define type 10
#define GAMMA2 == (Q-1)/88

int32_t decompose(int32_t *lowbit, int32_t input) {
    int32_t a1;
    int32_t highbit;

  a1  = (input + 127) >> 7;
    
    
#if GAMMA2 == (Q-1)/32
  a1  = (a1*1025 + (1 << 21)) >> 22;
  a1 &= 15;
    
#elif GAMMA2 == (Q-1)/88
  a1  = (a1*11275 + (1 << 23)) >> 24;
  a1 ^= ((43 - a1) >> 31) & a1;

#endif
  
    highbit = a1;
  
    *a0  = input - highbit*2*GAMMA2;
    *a0 -= (((Q-1)/2 - *a0) >> 31) & Q;
  return highbit;
}

void u_unmasking_highbits(int32_t* output,int32_t* input,int32_t betaa,int32_t* randomness)
{
    int32_t d[type];
    int32_t r1[type];
    int32_t r[type];
    
    
    int32_t mask = betaa -1;
    
    // random xor mask number . here is arith::mask
    for(int i = 0;i<type;i++){d[i] =  ((mask >>1)+1) ^ randomness[i];}
    
    //  here is arith::lshift
    for(int i = 0;i<type;i++){r1[i] =  input[i] << log2(betaa);}
    
    //  here is arith::add
    for(int i = 0;i<type;i++){r[i] =  r1[i] - d[i];}
    
    //here is HighBit-function
    for(int i = 0;i<type;i++){r1[i] =  (r[i]+127)>>7;}
    for(int i = 0;i<type;i++){r1[i] =  (r1[i]*11275 +(1<<23))>>24;}
    for(int i = 0;i<type;i++){r1[i] ^= ((43-r1[i])>>31)&r1[i];}
    for(int i = 0;i<type;i++){output[i] = r1[i];}
    
}

void u_masking_highbits(int32_t* output,int32_t* input,int32_t betaa,int32_t* randomness)
{
    int32_t d[type];
    int32_t r1[type];
    
    
    int32_t mask = betaa -1;
    
    //making random numeber
    for(int i = 0;i<type;i++){randomness[i] = random_int32_t();}
    
    // random xor mask number . here is arith::mask
    for(int i = 0;i<type;i++){d[i] =  ((mask >>1)+1) ^ randomness;}
    
    //  here is arith::add
    for(int i = 0;i<type;i++){r1[i] =  input[i] + d[i];}

    //  here is arith::rshift
    for(int i = 0;i<type;i++){output[i] =  r1[i] >> log2(betaa);}
    
}

int main(){

    int32_t input1[type],input2[type],input3[type];
    int32_t output1[type],output2[type],output3[type],output4[type];
    
    int32_t input01[type],input02[type],input03[type];
    
    int32_t betaa = 16;

    for(int i = 0;i<type;i++){
        input01[i]  = (1<<(i+15))+(1<<(i+16))+(1<<(i+16))+(1<<(i+14));
        input02[i]  = input01[i] ;
        input03[i]  = input01[i] ;
    }
    
    for(int i = 0;i<shares;i++){
        
        printf("input[%d]...%d\n",i,input01[i]);
        
        //masking-HighBits
        u_masking_highbits(output1[i],input1[i],betaa,randomness[i]);
        printf("Masking number is ...%d\n",output1[i]);
        
        input2[i] = output1[i]

        //unmasking-HighBits
        u_masking_highbits(output2[i],input2[i],betaa,randomness[i]);
        printf("Masking Unmasking number is ...%d\n",output2[i]);
        
        //only HighBits
        output3[i] = decompose(output4[i],input02[i]);
        printf("only high bits.................%d\n",output3[i]);
        printf("\n");

    }
        return 0;
}

