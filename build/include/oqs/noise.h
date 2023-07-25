#ifndef _NOISE_H_
#define _NOISE_H_

#include "../../dilithium/pqcrystals-dilithium_dilithium2_ref/random.h"
#include "../../dilithium/pqcrystals-dilithium_dilithium2_ref/poly.h"
#include "../../dilithium/pqcrystals-dilithium_dilithium2_ref/params.h"
#include "../../common/pqclean_shims/fips202.h"
#include <oqs/sha3.h>
#include <stdint.h>
#include <stdlib.h>
#include <stddef.h>

//#define SHAKE128_RATE 168
//#define SHAKE256_RATE 136

//typedef struct {OQS_SHA3_shake256_inc_ctx ctx;} shake256incctx;

//typedef struct {OQS_SHA3_shake128_inc_ctx ctx;} shake128incctx;

#define oqs_large_noise_rejection large_noise_rejection
#define oqs_large_noise_rejection_256 large_noise_rejection_256
#define oqs_shake128_absorb shake128_absorb
#define oqs_shake128_noise shake128_noise
#define oqs_shake128_squeezeblocks_noise shake128_squeezeblocks_noise
#define oqs_shake256_absorb shake256_absorb
#define oqs_shake256_noise shake256_noise
#define oqs_shake256_squeezeblocks_noise shake256_squeezeblocks_noise
#define oqs_small_noise_rejection small_noise_rejection
#define oqs_small_noise_rejection_256 small_noise_rejection_256



void sam(uint32_t* out_buffer, unsigned char *seed);

void small_bounded_noise_generation_256(uint32_t *out_buffer, uint8_t* seed, unsigned char nonce);

void small_bounded_noise_generation(uint32_t* out_buffer, unsigned char* seed, unsigned char nonce);

void large_bounded_noise_generation(poly *r, const uint8_t *seed, unsigned char nonce);

uint32_t large_noise_rejection_256(uint32_t *out_buffer, size_t out_len, unsigned char *buffer, size_t in_len);

void shake128_absorb(uint64_t *s,const unsigned char *input,unsigned long long inlen);


void shake128_squeezeblocks_noise(unsigned char *output,unsigned long nblocks,uint64_t *s);

void shake256_absorb(uint64_t *s,const unsigned char *input,unsigned long long inlen);

void shake256_squeezeblocks_noise(unsigned char *output,unsigned long nblocks,uint64_t *s);

void shake128_noise(unsigned char *output,unsigned long long outlen,const unsigned char *input,unsigned long long inlen);

void shake256_noise(unsigned char *output,unsigned long long outlen,const unsigned char *input,unsigned long long inlen);

static unsigned int rej_eta(int32_t *a, size_t len, const uint8_t *buf, size_t buflen);

#define shake128_inc_init OQS_SHA3_shake128_inc_init
#define shake128_inc_absorb OQS_SHA3_shake128_inc_absorb
#define shake128_inc_finalize OQS_SHA3_shake128_inc_finalize

#define shake256_inc_init OQS_SHA3_shake256_inc_init
#define shake256_inc_absorb OQS_SHA3_shake256_inc_absorb
#define shake256_inc_finalize OQS_SHA3_shake256_inc_finalize

//uint32_t small_noise_rejection(uint32_t* out_buffer, size_t out_len, unsigned char* buffer, size_t in_len);
//uint32_t large_noise_rejection(uint32_t* out_buffer, size_t out_len, unsigned char* buffer, size_t in_len);
//uint32_t large_noise_rejection_256(uint32_t* out_buffer, size_t out_len, unsigned char* buffer, size_t in_len);
//uint32_t small_noise_rejection_256(uint32_t* out_buffer, size_t out_len, unsigned char* buffer, size_t in_len);




#endif // _NOISE_H_
