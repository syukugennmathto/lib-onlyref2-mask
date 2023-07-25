#ifndef SYMMETRIC_H
#define SYMMETRIC_H
#include <oqs/common.h>

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "params.h"
#include "aes256ctr.h"
#include "noise.h"

#ifdef DILITHIUM_USE_AES


//void OQS_SHA3_shake256(uint8_t *output, size_t outlen, const uint8_t *input, size_t inplen);
//void OQS_SHA3_shake128(uint8_t *output, size_t outlen, const uint8_t *input, size_t inplen);

typedef aes256ctr_ctx stream128_state;
typedef aes256ctr_ctx stream256_state;

#define dilithium_aes256ctr_init DILITHIUM_NAMESPACE(dilithium_aes256ctr_init)
void dilithium_aes256ctr_init(aes256ctr_ctx *state,
                              const uint8_t key[32],
                              uint16_t nonce);

#define STREAM128_BLOCKBYTES AES256CTR_BLOCKBYTES
#define STREAM256_BLOCKBYTES AES256CTR_BLOCKBYTES

#define stream128_init(STATE, SEED, NONCE) \
        dilithium_aes256ctr_init(STATE, SEED, NONCE)
#define stream128_squeezeblocks(OUT, OUTBLOCKS, STATE) \
        aes256ctr_squeezeblocks(OUT, OUTBLOCKS, STATE)
#define stream128_release(STATE) \
        aes256_ctx_release(STATE)
#define stream256_init(STATE, SEED, NONCE) \
        dilithium_aes256ctr_init(STATE, SEED, NONCE)
#define stream256_squeezeblocks(OUT, OUTBLOCKS, STATE) \
        aes256ctr_squeezeblocks(OUT, OUTBLOCKS, STATE)
#define stream256_release(STATE) \
        aes256_ctx_release(STATE)

#else
typedef shake128incctx stream128_state;
typedef shake256incctx stream256_state;

/*
#define dilithium_shake128_stream_init DILITHIUM_NAMESPACE(dilithium_shake128_stream_init)
void dilithium_shake128_stream_init(shake128incctx *state,
                                    const uint8_t seed[SEEDBYTES],
                                    uint16_t nonce);

#define dilithium_shake256_stream_init DILITHIUM_NAMESPACE(dilithium_shake256_stream_init)
void dilithium_shake256_stream_init(shake256incctx *state,
                                    const uint8_t seed[CRHBYTES],
                                    uint16_t nonce);
 
 */

#define STREAM128_BLOCKBYTES 168
#define STREAM256_BLOCKBYTES 136

#define stream128_init(STATE, SEED, NONCE) \
        dilithium_shake128_stream_init(STATE, SEED, NONCE)
#define stream128_squeezeblocks(OUT, OUTBLOCKS, STATE) \
        shake128_squeezeblocks_noise(OUT, OUTBLOCKS, STATE)
#define stream128_release(STATE) shake128_inc_ctx_release(STATE)
#define stream256_init(STATE, SEED, NONCE) \
        dilithium_shake256_stream_init(STATE, SEED, NONCE)
#define stream256_squeezeblocks(OUT, OUTBLOCKS, STATE) \
        shake256_squeezeblocks_noise(OUT, OUTBLOCKS, STATE)
#define stream256_release(STATE) shake256_inc_ctx_release(STATE)

#endif

#endif
