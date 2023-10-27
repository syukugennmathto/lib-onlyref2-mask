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
        /*for(int i = 0;i<shares;i++){
            ring_a2[i] = ring[i];
        }*/
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
    /*
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
     */
}

int test_use_hint_both(){
    
    uint32_t ring,zooo,crash;
    
    /*crash = random_uint32_t_mono(1);
    printf("crash is %u\n",crash);
    crash = random_uint32_t_mono(1);
    printf("crash is %u\n",crash);
    */
    
    //for(int t= 0; t<3;t++){

        //printf("\n");
        //printf("\n");
        //printf("numking_use hint\n");
        //num_masking_usinghint(10,15);
        //printf("nosking_use hint\n");
        //noking_usinghint(10,15);
        printf("\n");
        printf("\n");
        printf("prime_masking_use_hint\n");
    //for(int i = 0;i<5;i++){
        //ring = (random_uint32_t() %(alpha/2)+alpha/2-1);
        ring = Q-10000;
        zooo = alpha/2-20;
        //printf("ring is %u\n",ring);
        masking_prime_usinghint(ring,zooo);
        printf("\n");
        printf("\n");
    //}
    return 0;
    
}
