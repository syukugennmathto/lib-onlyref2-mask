#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <oqs/oqs.h>
#include "pqcrystals-dilithium_dilithium2_ref/params.h"
#include "pqcrystals-dilithium_dilithium2_ref/sign.h"
//#include "pqcrystals-dilithium_dilithium2_ref/usehinttest.h"


int32_t  ew_masking_decompose(int32_t pR0,const int32_t pR) {
    
    int32_t mask = (1 << 23) - 1;
    int32_t b=0,pR1p=0,d_1 = (mask>>1)+1;
    int32_t pR1=0;
    int32_t pR1pt=0,pR2=0,pR3=0;
    
    pR1 = (pR+127)>>7;
    //printf("pR1= %d\n",pR1);
    pR2 = (pR1*11275 +(1<<23))>>24;
    //printf("pR2= %d\n",pR2);
    pR3 ^=((43-pR2)>>31)&43;
    // printf("pR3= %d\n",pR3);
    return pR3;
    
}

int main(){
    int32_t ringa,zooo,test;
    uint32_t crash;
    int32_t pR0,pR1[N],pR,pR2,pR3,pR4[N];
    ringa = 0;
    test =4*N;
    for(uint32_t i = 0;i<test;i++){
        pR1[i] = (1<<23)+(1<<30)+i;
    }
    
    for(uint32_t i = 0;i<test;i++){
        pR4[i]=ew_masking_decompose(pR0,pR1[i]);
        printf("HighBits (1<<%d)のとき%d\n",i,pR4[i]);
    }
    return 0;
}



