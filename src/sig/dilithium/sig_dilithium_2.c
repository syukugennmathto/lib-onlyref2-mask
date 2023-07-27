// SPDX-License-Identifier: MIT

#include <stdlib.h>

#include <oqs/sig_dilithium.h>
#include "pqcrystals-dilithium_dilithium2_ref/sign.h"

#if defined(OQS_ENABLE_SIG_dilithium_2)

OQS_SIG *OQS_SIG_dilithium_2_new(void) {

	OQS_SIG *sig = malloc(sizeof(OQS_SIG));
	if (sig == NULL) {
		return NULL;
	}
	sig->method_name = OQS_SIG_alg_dilithium_2;
	sig->alg_version = "https://github.com/pq-crystals/dilithium/commit/d9c885d3f2e11c05529eeeb7d70d808c972b8409";

	sig->claimed_nist_level = 2;
	sig->euf_cma = true;

	sig->length_public_key = OQS_SIG_dilithium_2_length_public_key;
	sig->length_secret_key = OQS_SIG_dilithium_2_length_secret_key;
	sig->length_signature = OQS_SIG_dilithium_2_length_signature;

	sig->keypair = OQS_SIG_keypair;
	sig->sign = OQS_SIG_sign;
	sig->verify = OQS_SIG_verify;

	return sig;
}

extern int pqcrystals_dilithium2_ref_keypair(uint8_t *pk, uint8_t *sk);
extern int pqcrystals_dilithium2_ref_signature(uint8_t *sig, size_t *siglen, const uint8_t *m, size_t mlen, const uint8_t *sk);
extern int pqcrystals_dilithium2_ref_verify(const uint8_t *sig, size_t siglen, const uint8_t *m, size_t mlen, const uint8_t *pk);


OQS_API OQS_STATUS OQS_SIG_keypair(uint8_t *public_key, uint8_t *secret_key) {
return (OQS_STATUS) pqcrystals_dilithium2_ref_keypair(public_key, secret_key);
}

OQS_API OQS_STATUS OQS_SIG_sign(uint8_t *signature, size_t *signature_len, const uint8_t *message, size_t message_len, const uint8_t *secret_key) {
return (OQS_STATUS) pqcrystals_dilithium2_ref_signature(signature, signature_len, message, message_len, secret_key);
}

OQS_API OQS_STATUS OQS_SIG__verify(const uint8_t *message, size_t message_len, const uint8_t *signature, size_t signature_len, const uint8_t *public_key) {
return (OQS_STATUS) pqcrystals_dilithium2_ref_verify(signature, signature_len, message, message_len, public_key);
}

#endif


