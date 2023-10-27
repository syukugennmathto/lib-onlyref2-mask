#ifndef PQCLEAN_MCELIECE348864_AVX2_AES256CTR_H
#define PQCLEAN_MCELIECE348864_AVX2_AES256CTR_H

#include <stddef.h>
#include <stdint.h>



#include "namespace.h"
#define AESCTR_NONCEBYTES 12
#define AES256_KEYBYTES 32

#define aes256ctr CRYPTO_NAMESPACE(aes256ctr)

void aes256ctr(
    uint8_t *out,
    size_t outlen,
    const uint8_t nonce[AESCTR_NONCEBYTES],
    const uint8_t key[AES256_KEYBYTES]
);

#endif
