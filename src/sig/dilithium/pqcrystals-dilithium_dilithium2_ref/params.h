#ifndef PARAMS_H
#define PARAMS_H

#include <stdint.h>
#include "config.h"



#define N 256
#define Q 8380417
#define D 13 //D=13
#define SEEDBYTES 32
#define CRHBYTES 64
#define ROOT_OF_UNITY 1753
#define dilithium_n N
#define dilithium_gamma 19
#define large_noise_bound alpha
#define dilithium_nu 5
#define low_noise_bound dilithium_nu
#define shares 6
#define MAX_SIGN 1000
#define MAX_VALUE_2_13 ((1U << 13) - 1)
#define badominton 31




#if DILITHIUM_MODE == 2
#define K 4
#define L 4
#define ETA 2
#define TAU 39
#define BETA 78
#define GAMMA1 (1 << 17)
#define GAMMA2 ((Q-1)/88)
#define OMEGA 80
#define dilithium_k K
#define dilithium_l L
#define alpha 2*GAMMA2
#endif

#define POLYT1_PACKEDBYTES  320
#define POLYT0_PACKEDBYTES  416
#define POLYVECH_PACKEDBYTES (OMEGA + K)

#if GAMMA1 == (1 << 17)
#define POLYZ_PACKEDBYTES   576
#elif GAMMA1 == (1 << 19)
#define POLYZ_PACKEDBYTES   640
#endif

#if GAMMA2 == (Q-1)/88
#define POLYW1_PACKEDBYTES  192
#elif GAMMA2 == (Q-1)/32
#define POLYW1_PACKEDBYTES  128
#endif

#if ETA == 2
#define POLYETA_PACKEDBYTES  96
#elif ETA == 4
#define POLYETA_PACKEDBYTES 128
#endif

#define CRYPTO_PUBLICKEYBYTES (SEEDBYTES + K*POLYT1_PACKEDBYTES)
#define CRYPTO_SECRETKEYBYTES (3*SEEDBYTES \
                               + L*POLYETA_PACKEDBYTES \
                               + K*POLYETA_PACKEDBYTES \
                               + K*POLYT0_PACKEDBYTES)
#define CRYPTO_BYTES (SEEDBYTES + L*POLYZ_PACKEDBYTES + POLYVECH_PACKEDBYTES)
#define masking_CRYPTO_BYTES (SEEDBYTES + L*POLYZ_PACKEDBYTES + POLYVECH_PACKEDBYTES +N*K*(1+N))

typedef struct{
    uint8_t a[SEEDBYTES];
    uint8_t b[SEEDBYTES];
    uint8_t c[K*POLYW1_PACKEDBYTES];
    uint8_t d[K*POLYW1_PACKEDBYTES];
} cst_8t;

typedef struct{
    uint32_t outside[N];
    uint32_t inside[N][N];
} mask_point;

typedef struct {
  mask_point vec[K];
} maskpointveck;

typedef struct {
  int32_t coeffs[N];
} poly;

typedef struct {
  uint32_t coeffs[N];
} upoly;

#endif
