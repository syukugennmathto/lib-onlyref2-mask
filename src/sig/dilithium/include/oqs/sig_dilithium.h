// SPDX-License-Identifier: MIT

#ifndef OQS_SIG_DILITHIUM_H
#define OQS_SIG_DILITHIUM_H

#include <oqs/oqs.h>



#ifdef OQS_ENABLE_SIG_dilithium_2
#define OQS_SIG_dilithium_2_length_public_key 1312
#define OQS_SIG_dilithium_2_length_secret_key 2528
#define OQS_SIG_dilithium_2_length_signature 2420

//OQS_SIG *OQS_SIG_dilithium_2_new(void);
//OQS_API OQS_STATUS OQS_SIG_keypair(uint8_t *public_key, uint8_t *secret_key);
//OQS_API OQS_STATUS OQS_SIG_sign(uint8_t *signature, size_t *signature_len, const uint8_t *message, size_t message_len, const uint8_t *secret_key);
//OQS_API OQS_STATUS OQS_SIG_verify(const uint8_t *message, size_t message_len, const uint8_t *signature, size_t signature_len, const uint8_t *public_key);
#endif

#endif
