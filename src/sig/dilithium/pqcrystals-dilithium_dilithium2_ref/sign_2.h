#ifndef SIGN_H
#define SIGN_H

#include <stddef.h>
#include <stdint.h>
#include "params.h"
#include "polyvec.h"
#include "poly.h"

#define masking_crypto_sign_keypair DILITHIUM_NAMESPACE(masking_crypto_sign_keypair)
int masking_crypto_sign_keypair(uint8_t *pk, uint8_t *sk);

#define masking_crypto_sign_signature DILITHIUM_NAMESPACE(masking_crypto_sign_signature)
int masking_crypto_sign_signature(uint8_t *sig,
                               size_t *siglen,
                               const uint8_t *m,
                               size_t mlen,
                               const uint8_t *sk);


#define masking_crypto_sign DILITHIUM_NAMESPACE(masking_crypto_sign)
int masking_crypto_sign(uint8_t *sm, size_t *smlen,const uint8_t *m, size_t mlen,const uint8_t *sk);

#define masking_crypto_sign_verify DILITHIUM_NAMESPACE(masking_crypto_sign_verify)
int masking_crypto_sign_verify(const uint8_t *sig, size_t siglen,const uint8_t *m, size_t mlen,const uint8_t *pk);


int crypto_sign_open(uint8_t *m, size_t *mlen,
                     const uint8_t *sm, size_t smlen,
                     const uint8_t *pk);

int number_chknorm(int32_t a, int32_t B);

#endif
