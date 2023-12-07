#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include "params.h"
#define QINV 58728449

/*
void masking_arithmetic_to_boolean_highbits(uint32_t* pR1, uint32_t* pR,mask_point* pD)
{
    uint32_t mask = ((1 << base) - 1);
    printf("mask is %u\n",mask);
    uint32_t b[shares]={0,0},pR1p[shares]={0,0},sd_1[shares] = {((mask>>1)+1),((mask>>1)+1),0};
    printf("sd_1 is %u\n",sd_1[0]);
    
    for(uint32_t i = 0; i < shares; ++i) b[i] = (pR[i] + sd_1[i])%Q;
    
    uint32_t d[shares]={0,0};
    pR1p[0] = b[0];
    for(uint32_t i = 1; i < shares;++i)
    {
        pD->d[i] = 64;  // pDをポインタ経由で逆参照してメンバーにアクセス
        
        pR1p[i] = pD->d[i]; // 同様にpDを逆参照
        
        pR1p[0] -= pR1p[i];
    }
    for(uint32_t j = 1; j < shares;++j)
    {
        d[0] = b[j];
        for(uint32_t i = 1; i < shares;++i)
        {
            pD->e[j][i] = 64; // pDを逆参照
            d[i] = pD->e[j][i]; // 同様にpDを逆参照
            d[0] -= d[i];
        }
        for(uint32_t i = 0; i < shares; ++i){pR1[i]=pR1p[i]+d[i];}
    }
    
    //for(uint32_t i = 0; i < shares; ++i){pR1[i] = pR1p[i] >> base;}
    for(uint32_t i = 0; i < shares; ++i){pR1[i] = pR1p[i] >> 4;}
    for(uint32_t i = 0; i < shares; ++i){pR1[i] = pR1[i]%Q;}
}

void nosking_arithmetic_to_boolean_highbits(uint32_t* pR1, uint32_t* pR,mask_point pD)
{
    uint32_t mask = ((1 << base) - 1);
    //printf("mask is %u\n",mask);
    uint32_t b[shares]={0,0},pR1p[shares]={0,0},sd_1[shares] = {((mask>>1)+1),((mask>>1)+1),0};
    //printf("sd_1 is %u\n",sd_1[0]);
    
    for(uint32_t i = 0; i < shares; ++i) b[i] = (pR[i] + sd_1[i])%Q;
    
    //masking_arithmetic_to_boolean_convert(pR0p,pR);
    uint32_t d[shares]={0,0};
    pR1p[0] = b[0];
    for(uint32_t i = 1; i < shares;++i)
    {
        //pD.d[i] = 10;  //random_uint32_t()
        pR1p[i] = pD.d[i];
        pR1p[0] -= pR1p[i];
    }
    
    for(uint32_t j = 1; j < shares;++j)
    {
        d[0] = b[j];
        for(uint32_t i = 1; i < shares;++i)
        {
            
            //pD.e[j][i] = 10; //random_uint32_t()
            d[i] = pD.e[j][i];
            d[0] -= d[i];
        }
        
        for(uint32_t i = 0; i < shares; ++i){pR1[i]=pR1p[i]+d[i];}
    }
    
    //for(uint32_t i = 0; i < shares; ++i){pR1[i] = pR1p[i] >> base;}
    for(uint32_t i = 0; i < shares; ++i){pR1[i] = pR1p[i] >> 4;}
    for(uint32_t i = 0; i < shares; ++i)
    {
        pR1[i] = pR1[i]%Q;
    }
}

void masking_arithmetic_to_boolean_lowbits(int32_t* pR0, uint32_t* pR,mask_point pD)
{
    uint32_t mask = (1 << base) - 1;
    uint32_t b[shares],pR0p[shares]={0,0},sd_1[shares] = {(mask>>1)+1,((mask>>1)+1),0};
     
     //masking_arithmetic_to_boolean_convert(pR0p,pR);
    uint32_t d[shares]={0,0};
     pR0p[0] = pR[0];
     for(uint32_t i = 1; i < shares;++i)
     {
     //pD.d[i] = 8;  //random_uint32_t()
     pR0p[i] = pD.d[i];
     pR0p[0] -= pR0p[i];
     }
     
     for(uint32_t j = 1; j < shares;++j)
     {
         d[0] = pR[j];
         for(uint32_t i = 1; i < shares;++i)
         {
             
             //pD.e[j][i] = 16; //random_uint32_t()
             d[i] = pD.e[j][i];
             d[0] -= d[i];
         }
         
         for(uint32_t i = 0; i < shares; ++i){pR0p[i]=pR0p[i]+d[i];
         }
     }
     
    //masking_boolean_lshift(pR0p,pR0p,32-base);
    for(uint32_t i = 0; i < shares; ++i)
    {
        pR0p[i] = pR0p[i] << (32-base);
    }
    
    //masking_boolean_rshift(b,pR0p,31);
    for(uint32_t i = 0; i < shares; ++i)
    {
        b[i] = pR0p[i] >> 31;
    }
    
    //masking_boolean_neg   (b,b);
    for(uint32_t i = 0; i < shares; ++i)
    {
        b[i] = -b[i];
    }
    
    //masking_boolean_lshift(b,b,base);
    for(uint32_t i = 0; i < shares; ++i)
    {
        b[i] = b[i] << base;
    }
    
    //masking_boolean_rshift(pR0p,pR0p,32-base);
    for(uint32_t i = 0; i < shares; ++i)
    {
        pR0p[i] = pR0p[i] >> (32-base);
    }
    
    //masking_boolean_xor   (pR0,pR0p,b);
    for(uint32_t i = 0; i < shares; ++i)
    {
        pR0[i] = (pR0p[i]) ^ (b[i]);
    }
    
    for(uint32_t i = 0; i < shares; ++i)
    {
        pR0[i] = pR0[i]%Q;
    }
    
}

void masking_arithmetic_two_makeint(uint32_t* c,uint32_t* z, uint32_t* r,mask_point pd)
{
    uint32_t r_z[shares],r_z_1[shares],r_1[shares],t[shares];
    //mask_point pd;
    
    //masking_arithmetic_to_boolean_highbits(r_1,r,base); //1
    nosking_arithmetic_to_boolean_highbits(r_1,r,pd);
    
    //masking_arithmetic_add(r_z,r_1,z); //2
    for(uint32_t i = 0; i < shares; ++i)
    {
        r_z[i] = (r_1[i]+z[i])%Q;
    }
    //masking_arithmetic_to_boolean_highbits(r_z_1,r_z,base); //3
    nosking_arithmetic_to_boolean_highbits(r_z_1,r_z,pd);
    
    //masking_boolean_xor(t,r_1,r_z_1); //4
    for(uint32_t i = 0; i < shares; ++i)
    {
        t[i] = (r_1[i]) ^ (r_z_1[i]);
    }

    
    //masking_boolean_lshift(t,t,base-1); //5
    for(uint32_t i = 0; i < shares; ++i)
    {
        t[i] = t[i] << (base-1);

    }
    
    //c = masking_boolean_fullxor(t); //6

    c[0] = t[0];
    for(uint32_t j = 0; j < shares; ++j)
    {
        c[j] =c[j] ^ t[j];
    }

}



void arithmetic_to_boolean_use_hint(uint32_t *r, unsigned int hint){
    uint32_t r0[shares], r1[shares];
    uint32_t mask = 44;
    //uint32_t mask = (1 << base) - 1;
    
    
    //masking_arithmetic_to_boolean_decompose(r1,r0,r,base);
    uint32_t mask = (1 << base) - 1;
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
    masking_arithmetic_to_boolean_highbits(uint32_t* pR1, uint32_t* pR,mask_point pD)

    
    
    
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
    
}*/




    /*
    int32_t pR0p[shares];
    uint32_t pR1p[shares],pR2p[shares];
    uint32_t pR1[shares],pR2[shares],pR3[shares],c[shares];
    mask_point pd,pe,pf;
    int32_t m = (Q-1)/base;
    
    for(uint32_t i = 0; i < shares+1;++i)
    {
        pR0p[i] =0;
        pR1p[i] =0;
        pR2p[i] =0;
        pR3[i] =0;
        c[i] = 0;
        pR1[i] = 2<<i*10;
        pR2[i] = 2>>i*3;
    }
    pR1[0] = base;
    for(uint32_t i = 0; i < shares+1;++i){
        pR1[i] = pR1[i] %Q;
        pR2[i] = pR2[i] %Q;
        pR3[i] = (pR1[i] + pR2[i])%Q;
    }
    
    for(uint32_t i = 0; i < shares+1;++i)
    {
        pd.d[i] =0;
        pe.d[i] =0;
        pf.d[i] =0;
        for(uint32_t j = 0; j < shares+1;++j)
        {
            pd.e[i][j] = 0;
            pe.e[i][j] = 0;
            pf.e[i][j] = 0;
        }
    }
    
    printf("masking_arithmetic_to_boolean_highbits\n");
    masking_arithmetic_to_boolean_highbits(pR1p,pR3,&pd);
    
    for(uint32_t i = 0; i < shares;++i)
    {
        pe.d[i] = pd.d[i];
        //printf("pd.d[%d] is %d \n",i,pd.d[i]);
        for(uint32_t j = 0; j < shares;++j)
        {
            pe.e[i][j] = pd.e[i][j];
            //printf("pd.e[%d][%d] is %d \n",i,j,pd.e[i][j]);
        }
    }
        for(uint32_t i = 0; i < shares;++i)
        {
            printf("pR1[%d] is %d \n",i,pR1p[i]);
        }
    printf("using_arithmetic_to_boolean_highbits\n");
    
    nosking_arithmetic_to_boolean_highbits(pR2p,pR2,pe);
    
    for(uint32_t i = 0; i < shares;++i)
    {
        pf.d[i] = pd.d[i];
        //printf("pd.d[%d] is %d \n",i,pe.d[i]);
        for(uint32_t j = 0; j < shares;++j)
        {
            pf.e[i][j] = pd.e[i][j];
            //printf("pd.e[%d][%d] is %d \n",i,j,pe.e[i][j]);
        }
    }
    
    masking_arithmetic_to_boolean_lowbits(pR0p,pR2,pf);
    
    for(uint32_t i = 0; i < shares;++i)
    {
        pf.d[i] = pd.d[i];
       //printf("pd.d[%d] is %d \n",i,pf.d[i]);
        for(uint32_t j = 0; j < shares;++j)
        {
            pf.e[i][j] = pd.e[i][j];
            //printf("pd.e[%d][%d] is %d \n",i,j,pf.e[i][j]);
        }
    }
    for(uint32_t i = 0; i < shares;++i){printf("pR0[%d] is %d \n",i,pR0p[i]);}
    printf("before computing\n");
    for(uint32_t i = 0; i < shares;++i){printf("pR2[%d] is %d \n",i,pR2p[i]);}
    
    
    masking_arithmetic_two_makeint(c,pR1,pR2,pd);
    
    for(uint32_t i = 0; i < shares;++i)
    {
        printf("c[%d] is %u \n",i,c[i]);
        if(c[i]==0)
        {
            pR2p[i]= pR2p[i];

        }
        else
        {
            if(pR0p[i]>0){pR2p[i] = (pR2p[i] + 1)%Q;}
            else{pR2p[i] = (pR2p[i] - 1)%Q;}
        }
        
    }
    
    printf("after computing\n");
    for(uint32_t i = 0; i < shares;++i){printf("pR2[%d] is %d \n",i,pR2p[i]);}
    
    printf("\n cheaking\n");
    for(uint32_t i = 0; i < shares;++i){printf("[%d]\n %d and %d\n",i,pR2p[i],pR1p[i]);}
     
     
     
     
     uint8_t sing[20],ct[10],zl[2][2],hat[2][2];
     
     unsigned int i,j;
     
     for(i=0; i < 20; ++i) sing[i]=0;
     for(i=0; i < 10; ++i) ct[i]=1;
     for(i=0; i < 2; ++i){for(j=0; j < 2; ++j)zl[i][j]=2;}
     for(i=0; i < 2; ++i){for(j=0; j < 2; ++j)hat[i][j]=2;}
     
     packtest_sig(sing,ct,zl,hat);
     
     for(uint32_t i=0; i < 20; ++i)printf("[%d] is %u\n",i,sing[i]);
    */

/*************************************************
* Name:        polyz_pack
*
* Description: Bit-pack polynomial with coefficients
*              in [-(GAMMA1 - 1), GAMMA1].
*
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            POLYZ_PACKEDBYTES bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
/*void polyz_pack(uint8_t *r, const uint8_t a[2]) {
  unsigned int i;
  uint32_t t[4];

  for(i = 0; i < 2; ++i) {
    t[0] = 3 - a[2*i+0];
    t[1] = 3 - a[2*i+1];

    r[5*i+0]  = t[0];
    r[5*i+1]  = t[0] >> 8;
    r[5*i+2]  = t[0] >> 16;
    r[5*i+2] |= t[1] << 4;
    r[5*i+3]  = t[1] >> 4;
    r[5*i+4]  = t[1] >> 12;
  }
}

void packtest_sig(uint8_t sig[20],
              const uint8_t c[10],
              const uint8_t z[2][2],
              const uint8_t h[2][2])
{
  unsigned int i, j, k;

  for(i=0; i < 9; ++i) sig[i] = c[i];
  sig += 9;

  for(i = 0; i < 3; ++i) polyz_pack(sig + i,z[i]);
  sig += 3;

  for(i = 0; i < 3 + 4; ++i) sig[i] = 0;

  k = 0;
  for(i = 0; i < 3; ++i) {
    for(j = 0; j < 3; ++j)
      if(h[i][j] != 0)
        sig[k++] = j;

    sig[3 + i] = k;
  
    
}
}*/


static const int32_t zetas[N] = {
         0,    25847, -2608894,  -518909,   237124,  -777960,  -876248,   466468,
   1826347,  2353451,  -359251, -2091905,  3119733, -2884855,  3111497,  2680103,
   2725464,  1024112, -1079900,  3585928,  -549488, -1119584,  2619752, -2108549,
  -2118186, -3859737, -1399561, -3277672,  1757237,   -19422,  4010497,   280005,
   2706023,    95776,  3077325,  3530437, -1661693, -3592148, -2537516,  3915439,
  -3861115, -3043716,  3574422, -2867647,  3539968,  -300467,  2348700,  -539299,
  -1699267, -1643818,  3505694, -3821735,  3507263, -2140649, -1600420,  3699596,
    811944,   531354,   954230,  3881043,  3900724, -2556880,  2071892, -2797779,
  -3930395, -1528703, -3677745, -3041255, -1452451,  3475950,  2176455, -1585221,
  -1257611,  1939314, -4083598, -1000202, -3190144, -3157330, -3632928,   126922,
   3412210,  -983419,  2147896,  2715295, -2967645, -3693493,  -411027, -2477047,
   -671102, -1228525,   -22981, -1308169,  -381987,  1349076,  1852771, -1430430,
  -3343383,   264944,   508951,  3097992,    44288, -1100098,   904516,  3958618,
  -3724342,    -8578,  1653064, -3249728,  2389356,  -210977,   759969, -1316856,
    189548, -3553272,  3159746, -1851402, -2409325,  -177440,  1315589,  1341330,
   1285669, -1584928,  -812732, -1439742, -3019102, -3881060, -3628969,  3839961,
   2091667,  3407706,  2316500,  3817976, -3342478,  2244091, -2446433, -3562462,
    266997,  2434439, -1235728,  3513181, -3520352, -3759364, -1197226, -3193378,
    900702,  1859098,   909542,   819034,   495491, -1613174,   -43260,  -522500,
   -655327, -3122442,  2031748,  3207046, -3556995,  -525098,  -768622, -3595838,
    342297,   286988, -2437823,  4108315,  3437287, -3342277,  1735879,   203044,
   2842341,  2691481, -2590150,  1265009,  4055324,  1247620,  2486353,  1595974,
  -3767016,  1250494,  2635921, -3548272, -2994039,  1869119,  1903435, -1050970,
  -1333058,  1237275, -3318210, -1430225,  -451100,  1312455,  3306115, -1962642,
  -1279661,  1917081, -2546312, -1374803,  1500165,   777191,  2235880,  3406031,
   -542412, -2831860, -1671176, -1846953, -2584293, -3724270,   594136, -3776993,
  -2013608,  2432395,  2454455,  -164721,  1957272,  3369112,   185531, -1207385,
  -3183426,   162844,  1616392,  3014001,   810149,  1652634, -3694233, -1799107,
  -3038916,  3523897,  3866901,   269760,  2213111,  -975884,  1717735,   472078,
   -426683,  1723600, -1803090,  1910376, -1667432, -1104333,  -260646, -3833893,
  -2939036, -2235985,  -420899, -2286327,   183443,  -976891,  1612842, -3545687,
   -554416,  3919660,   -48306, -1362209,  3937738,  1400424,  -846154,  1976782
};

int32_t montgomery_reduce(int64_t a) {
  int32_t t;

  t = (int64_t)(int32_t)a*QINV;
  t = (a - (int64_t)t*Q) >> 32;
  return t;
}

/*************************************************
* Name:        ntt
*
* Description: Forward NTT, in-place. No modular reduction is performed after
*              additions or subtractions. Output vector is in bitreversed order.
*
* Arguments:   - uint32_t p[N]: input/output coefficient array
**************************************************/
void ntt(int32_t *a) {
  unsigned int len, start, j, k;
  int32_t zeta, t;

  k = 0;
  for(len = 128; len > 0; len >>= 1) {
    for(start = 0; start < N; start = j + len) {
      zeta = zetas[++k];
      for(j = start; j < start + len; ++j) {
        t = montgomery_reduce((int64_t)zeta * a[j + len]);
        a[j + len] = a[j] - t;
        a[j] = a[j] + t;
      }
    }
  }
}

int main(){
    uint32_t x,t,sp;
    x = 1<<17;
    x = x -78;
    int a[N];
    sp=1<<3;
    for(int i=0;i<N;++i){
        a[i]=(1<<i)%(1<<30);
        
    }

    ntt(a);
    
    for(int i=0;i<N;++i){
        
        printf("%d is %d\n",a[i],i);
        printf("%d\n",x);
    }
    
}

