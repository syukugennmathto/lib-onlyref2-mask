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


int32_t  new_masking_lowbits(int32_t pR1,const int32_t pR) {

    

    int32_t mask = (1 << 23) - 1;

    int32_t b=0,pR1p=0,d_1 = (mask>>1)+1,basedd=4;

    int32_t pR0=0;

    int32_t pR1pt=0,pR2=0,pR3=0;

    

    //highbits

    pR1 = (pR+127)>>7;

    pR1 = (pR1*11275 +(1<<23))>>24;

    pR1 ^=((43-pR1)>>31)&pR1;

    

    //lowbits

    pR0  = pR - pR1*2*(Q-1)/88;

    //printf("pR0 = %d\n",pR0);

    pR0 -= (((Q-1)/2 - pR0) >> 31) & Q;

    //printf("pR0 = %d\n",pR0);

    

    return pR0;

}


void int32_masking_arithmetic_to_boolean_highbits(int32_t* pR1, int32_t * pR, int32_t base,pMtL*  pMR)
{
    int32_t mask = (1 << base) - 1;
    int32_t b[shares],pR1p[shares],sd_1[shares] = {(mask>>1)+1,0};
    int_masking_arithmetic_add(b,pR,sd_1);
    printf("pR1p...%d\n",b[4]);
    int_masking_arithmetic_to_boolean_convert(pR1p,b,pMR);

    int_masking_boolean_rshift(pR1,pR1p,base);
}

void int32_unmasking_arithmetic_to_boolean_highbits(int32_t* pR1, int32_t * pR, int32_t base,pMtL*  pMR)
{
    int32_t mask = (1 << base) - 1;
    int32_t b[shares],pR1p[shares],sd_1[shares] = {(mask>>1)+1,0};
    int_unmasking_arithmetic_add(b,pR,sd_1);
    printf("pR1p...%d\n",b[4]);
    int_unmasking_arithmetic_to_boolean_convert(pR1p,b,pMR);

    int_unmasking_boolean_rshift(pR1,pR1p,base);
}


int test_use_hint_both(int num){
    
    int32_t tips=28;
    int32_t time[tips],same[11][100];
    mask_point pD;
    pMtL pAM[shares];
    for (int i = 0;i<11;i++)
    {
        for (uint32_t round = 0;round<100;round++)
        {
            same[i][round]=0;;
        }
    }


        for(int i = 0;i<tips;i++){time[i]=0;}
    int32_t round=0;
                //printf("base値...%d\n",cycle);
                
                int32_t prime=0;
                int32_t ringa,zooo,bro1,bro2;
                
                uint32_t crash;
                int count=0;
                ringa = 0;
                crash = random_uint32_t_mono(1)%alpha;
                
                zooo = 546913;//2^32未満の値 546913 1071201
                int32_t ring1[shares];
                int32_t ring0[shares];
                int32_t ring2[shares];
                int32_t ring3[shares];
                int32_t ring4[shares];
                int32_t ring5[shares];
                int32_t ring6[shares];


                int32_t ring9[shares];
                int32_t ring10[shares];
                int32_t ring11[shares];
                int32_t ring12[shares];
                int32_t ring13[shares];
                int32_t ring14[shares];
                int32_t ring15[shares];
                int32_t ring16[shares];
                int32_t ring17[shares];
                int32_t ring18[shares];
                int32_t ring19[shares];
                int32_t ring20[shares];
                int32_t ring21[shares];
                int32_t ring22[shares];
                uint32_t ring7[shares];
                int32_t ring8[shares];
                int32_t zoo1[shares],zoo2[shares],zoo3[shares],zoo4[shares];
                //mask_point pD[shares],pE[shares];
                
                for(int i = 0;i<shares;i++){
                    ring0[i]=0;
                    ring1[i] = 0;
                    ring11[i] = 0;
                    ring2[i] = 0;
                    ring3[i] = 0;
                    ring4[i] = 0;
                    ring5[i] = 0;
                    ring6[i] = 0;
                    ring7[i] = 0;
                    ring8[i] = 0;
                    ring9[i] = 0;
                    ring10[i] = 0;
                    ring11[i] = 0;
                    ring12[i] = 0;
                    ring13[i] = 0;
                    ring14[i] = 0;
                    ring15[i] = 0;
                    ring16[i] = 0;
                    ring17[i] = 0;
                    ring18[i] = 0;
                    ring19[i] = 0;
                    ring20[i] = 0;
                    ring21[i] = 0;
                    ring22[i] = 0;
                    zoo1[i]  = (1<<(2*i+15))+(1<<(i+16))+(1<<(i+16))+(1<<(i+14));
                    //zoo1[i]  = zooo;
                    zoo2[i]  = zoo1[i];
                    zoo3[i]  = zoo1[i];
                    zoo4[i]  = zoo1[i];
                }
    
    /*
     ring3[i] = zoo1[i] ^ dec;
     //ring3[i] = new_masking_decompose(ring0[5],zoo1[i],10);
     ring2[i] = ring3[i] + d_1;
     ring5[i] = ring2[i] ;
     
     ring3[i] = new_masking_decompose(ring0[5],zoo1[i],10);
     ring2[i] = ring3[i]^ dec;
     ring5[i] = ring2[i]^d_1;
     */
    
    
    for(int i = 0;i<shares;i++){
        printf("input[%d]...%d\n",i,zoo1[i]);
        //printf("binary of input is...");
        //printBinary(zoo1[i]);printf("\n");
        //printf("\n");
        int32_t betaad = 16;
        int32_t mask=(1<<(1<<i))-(i*77),d_1=(mask>>1)+1;
        int32_t mask2 = betaad -1;
        int32_t dec = random_int32_t();
        printf("random[%d]...%d\n",i,dec);
        //masking-hb
        ring3[i] =  ((mask2 >>1)+1)^ dec;
        //ring3[i] = new_masking_decompose(ring0[5],zoo1[i],10);
        ring2[i] = zoo1[i] + ring3[i];
        ring5[i] = ring2[i] >> 4;
        printf("Masking number is ...%d\n",ring5[i]);

        
        //unmasking
        ring1[i] = ring5[i] << 4;
        ring4[i] = ring1[i] - ring3[i];
        ring7[i] = new_masking_decompose(ring0[5],ring4[i],10);
        printf("Masking Unmasking number is ...%d\n",ring7[i]);
        
        //hb
        ring1[i] = new_masking_decompose(ring0[5],zoo2[i],10);
        printf("only high bits.................%d\n",ring1[i]);
        
        new_masking_arithmetic_to_boolean_highbits(&ring20[i],&zoo3[i],&pD);
        new_unmasking_arithmetic_to_boolean_highbits(&ring21[i],&ring20[i],pD);
        ring0[5] = new_masking_lowbits(ring22[i],ring21[i]);
        printf("masking Highbits..............%d\n",ring22[i]);

        
        
        
        
        ring8[i] = new_masking_lowbits(ring11[i],zoo4[i]);
        printf("0olb...%d\n",ring8[i]);

        //masking-lb intput is pR

        ring9[i] = new_masking_lowbits(ring9[i],zoo3[i]);
        //convert
        ring10[i] = ring9[i] + dec;
        ring11[i] = ring10[i] ^ d_1;
        //ring11[i] = ring11[i] << 4;

        //ring7[i] =  ring9[i] + zoo3[i];

        //unmasking
        //ring11[i] = ring11[i] >> 4;
        ring12[i] = ring11[i] ^ d_1;
        ring13[i] = ring12[i] - dec;
        //ring7[i] = ring7[i]  -ring2[i];
        printf("3alb...%d\n",ring13[i]);
        
        printf("\n");
        
    }
    /*
    printf("first masking\n");
        for(int i = 0;i<shares -1;i++){
            int32_masking_arithmetic_to_boolean_highbits(&ring11[i],&zoo1[i],22,&pAM[i]);
            printf("Masking0[%d]...%d\n",i,ring11[i]);
            int32_unmasking_arithmetic_to_boolean_highbits(&ring3[i],&ring11[i],22,&pAM[i]);
            ring5[i] = new_masking_decompose(ringa,zoo2[i],22);
    }
    
    
    printf("i wanna know here\n");

    for(int i = 0;i<shares;i++){
        //ring6[i] = hard_masking_decompose(ringa,zoo3[i]);
        
    }
    for(int i = 0;i<shares;i++){
        //prime_highbits(&ring7[i],zoo3[i]);
    }
    for(int i = 0;i<shares -1;i++){
        int32_unmasking_arithmetic_to_boolean_highbits(&ring4[i],&ring11[i],22,&pAM[i]);
        int32_unmasking_arithmetic_to_boolean_highbits(&ring8[i],&zoo3[i],22,&pAM[i]);
        printf("\n");
    }
        //printf("HB,Masking\n");
              for(uint32_t s=0;s < shares;++s){
                    for(uint32_t i = 0; i < N;++i)
                    {
                        pE[s].outside[i] = pD[s].outside[i];
                    }
                    
                    for(uint32_t j = 0; j < N;++j)
                    {
                        for(uint32_t i = 0; i < N;++i)
                        {
                            pE[s].inside[j][i] = pD[s].inside[j][i];
                        }
                    }
                }
         
                //for(int i = 0;i<shares;i++){new_unmasking_arithmetic_to_boolean_highbits(&ring3[i],&ring1[i],pE[i]);}
                
    
        for(int i = 0;i<shares;i++){
            //printf("iの値...%d\n",i);
            printf("Masking01[%d]...%d\n",i,ring3[i]);
            printf("Masking02[%d]...%d\n",i,ring8[i]);
            printf("Masking11[%d]...%d\n",i,ring4[i]);
            printf("\n");
    }
                    
     printf("入力値...%d\n",zoo1[2]);
    printBinary(zoo1[2]);printf("\n");
    printf("Nesking1[%d]...%d\n",2,ring5[2]);
    printf("Hosking2[%d]...%d\n",2,ring6[2]);
    printf("Posking3[%d]...%d\n",2,ring7[2]);
    printf("\n");
    printBinary(ring5[2]);printf("\n");
    printBinary(ring6[2]);printf("\n");
    printf("\n");
    for(int i = 0;i<shares;i++){
        //printf("iの値...%d\n",i);
        printBinary(ring3[i]);printf("\n");
}
    
    
    for(int32_t j = 1; j < shares; ++j)
    {
        bro1 = bro1 ^ ring3[j];
        bro2 = bro2 ^ ring4[j];
    }
    
    printf("bro1...%d\n",bro1);
    printf("bro2...%d\n",bro2);
     */
                
                /*if(count>=0&&count<4){
                    printf("base is %dにおける不一致数が4以下...%d\n",25,count);
                    time[25] +=1;
                    
                }*/
                //printf("base is %dにおける不一致数...%d\n",cycle,count);}
                //if(count==0){if(prime>0){printf("2^%dにおける一致したときの43の個数...%d\n",21,prime);}}

            
            //printf("ここに至るまで %u 回かかった\n",round);
        //printf("\n");
            //printf("base is %dにおける不一致数...%d\n",25,time[25]);
            
        
        
    
        /*for (int i = 4;i<shares;i++)
        {
            for (uint32_t round = 0;round<10;round++)
            {
            printf("不一致だった箇所iの値%d,%d回目の個数...%d\n",i,round,same[i][round]);
            }
        }
    */
    
    return 0;
    
}
