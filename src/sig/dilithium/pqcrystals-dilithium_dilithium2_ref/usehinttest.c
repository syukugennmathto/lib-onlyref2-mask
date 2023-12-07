#include <stdint.h>
#include <math.h>
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






/*
 
 
 void  masking_usinghint(uint32_t ringa,uint32_t zooa)
{
    uint32_t hint,uusehint,uuusehint,hhighbits,ringzoo1,ring1,ring0;
    unsigned char hhint=0;
    uint32_t ringzoo[shares],ring[shares],zoo[shares],hhintt[shares];
    
    for(int i = 0;i<shares;i++){
        ring[i] = 10;
        zoo[i]  = 15;
        ringzoo[i] = ring[i] + zoo[i];
    }
    
    hhint = masking_arithmetic_makeint(zoo,ring,alpha);
    printf("ari_mask_hint is %u\n",hhint);
    
    uuusehint = uuse_hint(ring[1],hhint,alpha);
    printf("nosk_uusehint is %u \n",uuusehint);
    
    masking_arithmetic_to_boolean_decompose(&hhighbits,&ring0,ringzoo,alpha);
    printf("mask_decobits is %u \n",hhighbits);
    
    masking_arithmetic_to_boolean_highbits(&hhighbits,ringzoo, alpha);
    printf("mask_highbits is %u \n",hhighbits);

    uuusehint = masking_use_hint(ringzoo[1],hhint,alpha);

    printf("mask_uusehint is %u \n",uuusehint);
    printf("\n");
    printf("%u and %u are correct\n",uuusehint,hhighbits);
    
}

void  num_masking_usinghint(uint32_t ring,uint32_t zoo)
{
    uint32_t hint,uusehint,uuusehint,hhighbits,ringzoo,ring1,ring0;
    unsigned char hhint=0;


    ring = 10;
    zoo  = 15;
    hhint = 0;
    ringzoo = ring + zoo;
    printf("aaaaaaringzoo is %u\n",ringzoo);
    printf("ring is %u\n",ring);
    printf("zooo is %u\n",zoo);
    hhint = num_masking_arithmetic_makehint(zoo,ring,alpha);
    printf("num_mask_hint is %u\n",hhint);
    
    uusehint = num_masking_arithmetic_usehint(hhint,ring,alpha);
    printf("mask_uusehint is %u \n",uusehint);
    
    num_masking_arithmetic_to_boolean_decompose(hhighbits,ring0,ringzoo,alpha);
    printf("mask_decobits is %u \n",hhighbits);
    
    num_masking_arithmetic_to_boolean_highbits(hhighbits,ringzoo,alpha);
    printf("mask_highbits is %u \n",hhighbits);
    
    uuusehint = uuse_hint(ring,hhint,alpha);
    printf("nosk_uusehint is %u \n",uuusehint);
    printf("\n");
    printf("%u and %u are correct\n",uusehint,hhighbits);
    
    
}

void  noking_usinghint(uint32_t ring,uint32_t zoo)
{
    uint32_t uusehint,hhighbits,ringzoo,ring1,ring0,hhighbbitss,highb;
    unsigned int hhint;
    
    ring = alpha/2 + 3;
    zoo  = alpha/2 - 2;
    hhint = 0;
    ringzoo = ring + zoo;
    printf("aaaaaaringzoo is %u\n",ringzoo);
    printf("ring is %u\n",ring);
    printf("zooo is %u\n",zoo);
    
    hhint = make_hint(ring,zoo);
    printf("nosk_hintting is %u\n",hhint);
    printf("\n");
    uusehint = uuse_hint(ring,hhint,alpha);
    printf("nosk_uusehint is %u\n",uusehint);
    printf("\n");
    
    hhighbbits(hhighbits,ringzoo);
    printf("nosk_highbits is %u\n",hhighbits);
    highbits(highb,ringzoo,alpha);
    
    prime_decompose(&hhighbbitss,&ring0,ringzoo);
    printf("nosk_highdeco is %u\n",hhighbbitss);
    printf("comp_highbits is %u\n",highb);

    printf("%u and %u are correct\n",uusehint,hhighbits);
    
 
}
void  masking_prime_usinghint(uint32_t ringa,uint32_t zooa)
{
    uint32_t hint,uusehint,uuusehint,hhighbits,ringzoo1,ring1,ring0;
    unsigned char hhint;
    uint32_t ringzoo[shares];
    uint32_t hhintt[shares];
    uint32_t ring_a2[shares],ring_a1[shares],ring_a0[shares],zoo_a[shares],highbits_a[shares];
    
    
    
    //printf("nosk_delobits is %u \n",ring0);

    uint32_t base =alpha;
    
    //for(base = 12;base<32;base++){
        printf("\n");
        printf("base is %u\n",base);
        printf("\n");

    uint32_t ring[shares];
    uint32_t zoo[shares]; // r_1 の対応するビットと異なるビットが 1
    
    for(int i = 0;i<shares;i++){
        ring[i] = zooa;
        zoo[i]  = zooa;
        ringzoo[i] = ring[i] + zoo[i];
    }
        //printf("ringzoo is %u\n",ringzoo[1]);
        hhint = masking_arithmetic_two_makeint(zoo,ring,base);
        printf("ari_mask_hint is %d\n",hhint);
        printf("\n");

        arithmetic_to_boolean_use_hint(ring,hhint,base);
        
 for(int i = 0;i<shares;i++){
            ring_a2[i] = ring[i];
        }
        masking_arithmetic_to_boolean_decompose(ring_a1,ring_a0,ringzoo,base);
        masking_arithmetic_to_boolean_highbits(highbits_a,ringzoo,base);
        for(int i = 0;i<shares;i++){
            printf("mask_highbit%d is %x \n",i,highbits_a[i]);
            printf("mask_decombt%d is %x \n",i,ring_a1[i]);
            printf("arit_usehint%d is %x \n",i,ring[i]);
            printf("\n");
        }
    printf("\n");
    
    base =19;
    for(int i = 0;i<shares;i++){
        ring[i] = ringa;
        zoo[i]  = zooa;
        ringzoo[i] = ring[i] + zoo[i];
    }
    for(uint32_t i = 0; i < shares;++i){
        hhintt[i] = pt_makehint(zoo[i],ring[i],base);
        pt_highbits(&ring_a1[i],ringzoo[i],base);
        highbits_a[i] = pt_usehint(hhintt[i],ring[i],base);
        printf("nosk_makehin%d is %u \n",i,hhintt[i]);
        printf("nosk_highbit%d is %d \n",i,ring_a1[i]);
        printf("noit_usehint%d is %d \n",i,highbits_a[i]);
        printf("\n");
    }
    
    ring[1] = zoo[1];
    hhint = make_hint(ring[1],zoo[1]);
    printf("ari_mask_hint is %d\n",hhint);
    //}
    printf("alpha is %u\n",aalpha);
    masking_arithmetic_to_boolean_decompose(ring_a1,ring_a0,ringzoo,aalpha);
    masking_arithmetic_to_boolean_highbits(highbits_a,ringzoo,aalpha);
    for(int i = 0;i<shares;i++){
        printf("mask_decobit%d is %u \n",i,ring_a1[i]);
        printf("mask_highbit%d is %u \n",i,highbits_a[i]);
        printf("\n");
        
    }
    
    //printf("mask_delobits is %u \n",ring0);
    

    //prime_highbits(&hhighbits,ringzoo[1]);
    //printf("nosk_highbits is %u \n",hhighbits);
    


    

    

    
    //printf("nosk_hintting is %u\n",hhint);
    arithmetic_to_boolean_use_hint(ring,hint,base);
    uuusehint = uuse_hint(ring[1],hhint,alpha);
    printf("nosk_uusehint is %u \n",uuusehint);
    //uuusehint = masking_use_hint(ring[1],hhint,alpha);
    //printf("mask_uusehint is %u \n",uuusehint);
    printf("\n");

    printf("%u and %u are correct\n",uuusehint,hhighbits);
     
}


    
    crash = random_uint32_t_mono(1);
    printf("crash is %u\n",crash);
    crash = random_uint32_t_mono(1);
    printf("crash is %u\n",crash);
    
    
    //for(int t= 0; t<3;t++){

        //printf("\n");
        //printf("\n");
        //printf("numking_use hint\n");
        //num_masking_usinghint(10,15);
        //printf("nosking_use hint\n");
        //noking_usinghint(10,15);
       // printf("\n");
       // printf("\n");
       // printf("prime_masking_use_hint\n");
    //for(int i = 0;i<5;i++){
        //ring = (random_uint32_t() %(alpha/2)+alpha/2-1);

        //printf("ring is %u\n",ring);
       // masking_prime_usinghint(ring,zooo);
       // printf("\n");
       // printf("\n");
    //}

int test_use_hint_both(){
    
    int32_t ringa,zooo;
    uint32_t crash;
    ringa = 0;
    
    for(int s = 0; s<92;s++){
        printf("\n");
        printf("s = %d\n",s);
        zooo = -44 + s;
        int32_t ring1[shares];
        int32_t ring0[shares];
        int32_t ring11[shares];
        int32_t ring2[shares];
        int32_t ring3[shares];
        int32_t ring4[shares];
        int32_t ring5[shares];
        int32_t zoo1[shares],zoo2[shares];
        mask_point pD,pE;
        
        printf("1");
        for(int i = 0;i<shares;i++){
            ring0[i]=0;
            ring1[i] = 0;
            ring11[i] = 0;
            ring2[i] = 0;
            ring3[i] = 0;
            ring4[i] = 0;
            ring5[i] = 0;
            zoo1[i]  = zooo;
            zoo2[i]  = zooo;
        }
        
        
        printf("zoo = %d\n",zoo1[0]);
        
        for(int i = 0;i<shares;i++){
            new_masking_arithmetic_to_boolean_highbits(&ring1[i],&zoo1[i],&pD);
        }
        
        for(uint32_t i = 0; i < N;++i)
        {
            pE.outside[i] = pD.outside[i];
        }
        
        for(uint32_t j = 0; j < N;++j)
        {
            for(uint32_t i = 0; i < N;++i)
            {
                pE.inside[j][i] = pD.inside[j][i];
            }
        }
        
        for(int i = 0;i<shares;i++){
            printf("1ring[%d] = %d\n",i,ring1[i]);
        }
        
        printf("break\n");
        for(int i = 0;i<shares;i++){
            ring2[i] = decompose(&ring0[i],ring1[i]);
        }
        
        for(int i = 0;i<shares;i++){
            new_nosking_arithmetic_to_boolean_highbits(&ring3[i],&ring2[i],pE);
        }
        for(int i = 0;i<shares;i++){
            printf("3ring[%d] = %d\n",i,ring3[i]);
        }
        
        printf("break\n");
        for(int i = 0;i<shares;i++){
            ring4[i] = decompose(&ring0[i],ring3[i]);
            printf("4ring[%d] = %d\n",i,ring4[i]);
        }
        
        printf("break\n");
        for(int i = 0;i<shares;i++){
            ring5[i] = decompose(&ring0[i],zoo2[i]);
            printf("5ring[%d] = %d\n",i,ring5[i]);
        }
        
        
        printf("here is breaking point\n");
        for(int i = 0;i<shares;i++){
            ring1[i] = 0;
            ring11[i] = 0;
            ring2[i] = 0;
            ring3[i] = 0;
            ring4[i] = 0;
            ring5[i] = 0;
            zoo1[i]  = zooo;
            zoo2[i]  = zooo;
        }
        
        
        printf("zoo = %d\n",zoo1[0]);
        
        printf("break\n");
        for(int i = 0;i<shares;i++){
            ring1[i] = decompose(&ring0[i],zoo1[i]);
        }
        
        for(int i = 0;i<shares;i++){
            new_nosking_arithmetic_to_boolean_highbits(&ring2[i],&ring1[i],pE);
        }
        
        for(int i = 0;i<shares;i++){
            printf("1ring[%d] = %d\n",i,ring2[i]);
        }
        
        for(int i = 0;i<shares;i++){
            new_nosking_arithmetic_to_boolean_highbits(&ring3[i],&ring2[i],pE);
        }
        
        for(int i = 0;i<shares;i++){
            printf("3ring[%d] = %d\n",i,ring3[i]);
        }
        
        printf("break\n");
        for(int i = 0;i<shares;i++){
            ring4[i] = decompose(&ring0[i],ring3[i]);
            printf("4ring[%d] = %d\n",i,ring4[i]);
        }
        
        printf("break\n");
        for(int i = 0;i<shares;i++){
            ring5[i] = decompose(&ring0[i],zoo2[i]);
            printf("5ring[%d] = %d\n",i,ring5[i]);
        }
        
        
    }
        return 0;
        
    
}*/



int test_use_hint_both(){
    
    int32_t ringa,zooo;
    uint32_t crash;
    ringa = 0;
    crash = random_uint32_t_mono(1)%alpha;
        zooo = (1<<31)+Q-crash;
        int32_t ring1[shares];
        int32_t ring0[shares];
        int32_t ring11[shares];
        int32_t ring2[shares];
        int32_t ring3[shares];
        int32_t ring4[shares];
        int32_t ring5[shares];
        int32_t ring6[shares];
        int32_t zoo1[shares],zoo2[shares];
        mask_point pD,pE;
        
        printf("1");
        for(int i = 0;i<shares;i++){
            ring0[i]=0;
            ring1[i] = 0;
            ring11[i] = 0;
            ring2[i] = 0;
            ring3[i] = 0;
            ring4[i] = 0;
            ring5[i] = 0;
            zoo1[i]  = zooo;
            zoo2[i]  = zooo;
        }
        
        
        //printf("zoo = %d\n",zoo1[0]);
    

        

        
        for(int i = 0;i<shares;i++){
            //printf("1ring[%d] = %d\n",i,ring1[i]);
        }
        
        //printf("break\n");
        for(int i = 0;i<shares;i++){
        //    ring2[i] = new_masking_decompose(ring0[i],ring1[i]);
        }
        
        for(int i = 0;i<shares;i++){
       //     new_nosking_arithmetic_to_boolean_highbits(&ring3[i],&ring2[i],pE);
        }
        for(int i = 0;i<shares;i++){
            //printf("3ring[%d] = %d\n",i,ring3[i]);
        }
        
        //printf("break\n");
        for(int i = 0;i<shares;i++){
       //     ring4[i] = new_masking_decompose(ring0[i],ring3[i]);
            //printf("4ring[%d] = %d\n",i,ring4[i]);
        }
        
        //printf("break\n");
        for(int i = 0;i<shares;i++){
       //     ring5[i] = new_masking_decompose(ring0[i],zoo2[i]);
            //printf("5ring[%d] = %d\n",i,ring5[i]);
        }
        
        
       // printf("here is breaking point\n");
        for(int i = 0;i<shares;i++){
            ring1[i] = 0;
            ring11[i] = 0;
            ring2[i] = 0;
            ring3[i] = 0;
            ring4[i] = 0;
            ring5[i] = 0;
            zoo1[i]  = zooo;
            zoo2[i]  = zooo;
        }
        
        
       // printf("zoo = %d\n",zoo1[0]);
        
       // printf("break\n");
        for(int i = 0;i<shares;i++){
         //   ring1[i] = new_masking_decompose(ring0[i],zoo1[i]);
        }
        
        for(int i = 0;i<shares;i++){
        //    new_nosking_arithmetic_to_boolean_highbits(&ring2[i],&ring1[i],pE);
        }
        
        for(int i = 0;i<shares;i++){
      //      printf("1ring[%d] = %d\n",i,ring2[i]);
        }
        
        for(int i = 0;i<shares;i++){
         //   new_nosking_arithmetic_to_boolean_highbits(&ring3[i],&ring2[i],pE);
        }
        
        for(int i = 0;i<shares;i++){
      //      printf("3ring[%d] = %d\n",i,ring3[i]);
        }
        
        printf("break\n");
        for(int i = 0;i<shares;i++){
      //      ring4[i] = new_masking_decompose(ring0[i],ring3[i]);
       //     printf("4ring[%d] = %d\n",i,ring4[i]);
        }
        

        
    
   // printf("here is breaking point\n");
    for(int i = 0;i<shares;i++){
        ring0[i]=0;
        ring1[i] = 0;
        ring11[i] = 0;
        ring2[i] = 0;
        ring3[i] = 0;
        ring4[i] = 0;
        ring5[i] = 0;
        zoo1[i]  = zooo;
        zoo2[i]  = zooo;
    }
    
    
   // printf("zoo = %d\n",zoo1[0]);
    
   // printf("break\n");
    for(int i = 0;i<shares;i++){
     //   new_nosking_arithmetic_to_boolean_highbits(&ring2[i],&ring1[i],pE);
    }
    
    for(int i = 0;i<shares;i++){
  //      printf("1ring[%d] = %d\n",i,ring2[i]);
    }
  //  printf("break\n");
    for(int i = 0;i<shares;i++){
     //   new_nosking_arithmetic_to_boolean_highbits(&ring3[i],&ring2[i],pE);
    }
    
    for(int i = 0;i<shares;i++){
  //      printf("3ring[%d] = %d\n",i,ring3[i]);
    }
    
  //  printf("here is breaking point\n");
    for(int i = 0;i<shares;i++){
        ring0[i]=0;
        ring1[i] = 0;
        ring11[i] = 0;
        ring2[i] = 0;
        ring3[i] = 0;
        ring4[i] = 0;
        ring5[i] = 0;
        ring6[i] = 0;
        zoo1[i]  = zooo;
        zoo2[i]  = zooo;
    }
    
    
  //  printf("zoo = %d\n",zoo1[0]);
    for(int i = 0;i<shares;i++){
      new_masking_arithmetic_to_boolean_highbits(&ring1[i],&zoo1[i],&pD);
    //    printf("1ring[%d] = %d\n",i,ring1[i]);
    //    printBinary(ring1[i]);
    //    printf("\n");
    }



    for(uint32_t i = 0; i < N;++i)
    {
    pE.outside[i] = pD.outside[i];
    }
    
    for(uint32_t j = 0; j < N;++j)
    {
        for(uint32_t i = 0; i < N;++i)
        {
            pE.inside[j][i] = pD.inside[j][i];
        }
    }

    
  //  printf("break\n");
    for(int i = 0;i<shares;i++){
   //     new_nosking_arithmetic_to_boolean_highbits(&ring2[i],&zoo2[i],pE);
    }
 //   printf("ここからunmaskingした\n");
  //  printf("Masking\n");
    for(int i = 0;i<shares;i++){
        new_unsking_arithmetic_to_boolean_highbits(&ring3[i],&ring1[i],pE);
    }
  //  printf("HB,Masking,Unmasking,HBした値...\n");
    for(int i = 0;i<shares;i++){
  //      printf("3ring[%d] = %d\n",i,ring3[i]);
  //     printBinary(ring3[i]);
  //         printf("\n");
    }

    
   // printf("入力値は %d\n",zoo1[0]);
    

    
   // printf("実際に,HighBitsされた\n");
    for(int i = 0;i<shares;i++){
        ring5[i] = new_masking_decompose(ring0[i],zoo2[i]);
  //      printf("5ring[%d] = %d\n",i,ring5[i]);
  //     printBinary(ring5[i]);
  //         printf("\n");
    }
    printf("入力値...%d\n",zoo1[0]);
    printBinary(zoo1[0]);
    printf("\n");
    printf("HB,Masking,UnMaskingして最後にHBした値...%d\n",ring3[3]);
    printBinary(ring3[3]);
        printf("\n");
    printf("実際に,HighBitsされた値..................%d\n",ring5[3]);
    printBinary(ring5[3]);
        printf("\n");

        return 0;
        
    
}

