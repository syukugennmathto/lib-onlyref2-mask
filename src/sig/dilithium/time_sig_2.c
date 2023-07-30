// SPDX-License-Identifier: MIT

#define OQS_DISABLE_SIG_dilithium_2
#define OQS_ENABLE_NO_DEPRECATED

// #include <openssl/opensslv.h>  // この行をコメントアウトします

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <oqs/oqs.h>
#include "params.h"

#if defined(OQS_USE_RASPBERRY_PI)
#define _OQS_RASPBERRY_PI
#endif
#if defined(OQS_SPEED_USE_ARM_PMU)
#define SPEED_USE_ARM_PMU
#endif
#include "ds_benchmark.h"
//#include "system_info.c"
#include "pqcrystals-dilithium_dilithium2_ref/sign.h"

#if OQS_USE_PTHREADS_IN_TESTS
#include <pthread.h>
#endif

#ifdef OQS_ENABLE_TEST_CONSTANT_TIME
#include <valgrind/memcheck.h>
#define OQS_TEST_CT_CLASSIFY(addr, len)  VALGRIND_MAKE_MEM_UNDEFINED(addr, len)
#define OQS_TEST_CT_DECLASSIFY(addr, len)  VALGRIND_MAKE_MEM_DEFINED(addr, len)
#else
#define OQS_TEST_CT_CLASSIFY(addr, len)
#define OQS_TEST_CT_DECLASSIFY(addr, len)
#endif


typedef struct magic_s {
    uint8_t val[31];
} magic_t;

static OQS_STATUS sig_test_correctness(const char *method_name) {

    OQS_SIG *sig = NULL;
    uint8_t *public_key = NULL;
    uint8_t *secret_key = NULL;
    uint8_t *message = NULL;
    size_t message_len = 100;
    uint8_t *signature = NULL;
    size_t signature_len;
    OQS_STATUS rc, ret = OQS_ERROR;
    
    sig = OQS_SIG_new(method_name);
    if (sig == NULL) {
        return OQS_SUCCESS;
    }
    size_t *sig_len = sig_len;

    public_key = malloc(sig->length_public_key);
    secret_key = malloc(sig->length_secret_key);
    message = malloc(message_len);
    signature = malloc(sig->length_signature);

    //The magic numbers are random values.
    //The length of the magic number was chosen to be 31 to break alignment
    magic_t magic;
    OQS_randombytes(magic.val, sizeof(magic_t));

    sig = OQS_SIG_new(method_name);
    if (sig == NULL) {
        fprintf(stderr, "ERROR: OQS_SIG_new failed\n");
        goto err;
    }

    printf("================================================================================\n");
    printf("Sample computation for signature %s\n", sig->method_name);
    printf("================================================================================\n");

    public_key = malloc(sig->length_public_key + 2 * sizeof(magic_t));
    secret_key = malloc(sig->length_secret_key + 2 * sizeof(magic_t));
    message = malloc(message_len + 2 * sizeof(magic_t));
    signature = malloc(sig->length_signature + 2 * sizeof(magic_t));

    if ((public_key == NULL) || (secret_key == NULL) || (message == NULL) || (signature == NULL)) {
        fprintf(stderr, "ERROR: malloc failed\n");
        goto err;
    }

    //Set the magic numbers before
    memcpy(public_key, magic.val, sizeof(magic_t));
    memcpy(secret_key, magic.val, sizeof(magic_t));
    memcpy(message, magic.val, sizeof(magic_t));
    memcpy(signature, magic.val, sizeof(magic_t));

    public_key += sizeof(magic_t);
    secret_key += sizeof(magic_t);
    message += sizeof(magic_t);
    signature += sizeof(magic_t);

    // and after
    memcpy(public_key + sig->length_public_key, magic.val, sizeof(magic_t));
    memcpy(secret_key + sig->length_secret_key, magic.val, sizeof(magic_t));
    memcpy(message + message_len, magic.val, sizeof(magic_t));
    memcpy(signature + sig->length_signature, magic.val, sizeof(magic_t));

    OQS_randombytes(message, message_len);
    OQS_TEST_CT_DECLASSIFY(message, message_len);

    rc = mask_crypto_sign_keypair(public_key, secret_key);
    OQS_TEST_CT_DECLASSIFY(&rc, sizeof rc);
    if (rc != OQS_SUCCESS) {
        fprintf(stderr, "ERROR: OQS_SIG_keypair failed\n");
        goto err;
    }

    rc = mask_crypto_sign(signature, &signature_len, message, message_len, secret_key);
    OQS_TEST_CT_DECLASSIFY(&rc, sizeof rc);
    if (rc != OQS_SUCCESS) {
        fprintf(stderr, "ERROR: OQS_SIG_sign failed\n");
        goto err;
    }

    OQS_TEST_CT_DECLASSIFY(public_key, sig->length_public_key);
    OQS_TEST_CT_DECLASSIFY(signature, signature_len);
    rc = mask_crypto_sign_verify((const uint8_t *)sig, *sig_len,message, message_len,public_key);
    OQS_TEST_CT_DECLASSIFY(&rc, sizeof rc);
    if (rc != OQS_SUCCESS) {
        fprintf(stderr, "ERROR: OQS_SIG_verify failed\n");
        goto err;
    }

    /* modify the signature to invalidate it */
    OQS_randombytes(signature, signature_len);
    OQS_TEST_CT_DECLASSIFY(signature, signature_len);
    rc = mask_crypto_sign_verify((const uint8_t *)sig, *sig_len,message, message_len,public_key);
    OQS_TEST_CT_DECLASSIFY(&rc, sizeof rc);
    if (rc != OQS_ERROR) {
        fprintf(stderr, "ERROR: OQS_SIG_verify should have failed!\n");
        goto err;
    }

#ifndef OQS_ENABLE_TEST_CONSTANT_TIME
    /* check magic values */
    int rv = memcmp(public_key + sig->length_public_key, magic.val, sizeof(magic_t));
    rv |= memcmp(secret_key + sig->length_secret_key, magic.val, sizeof(magic_t));
    rv |= memcmp(message + message_len, magic.val, sizeof(magic_t));
    rv |= memcmp(signature + sig->length_signature, magic.val, sizeof(magic_t));
    rv |= memcmp(public_key - sizeof(magic_t), magic.val, sizeof(magic_t));
    rv |= memcmp(secret_key - sizeof(magic_t), magic.val, sizeof(magic_t));
    rv |= memcmp(message - sizeof(magic_t), magic.val, sizeof(magic_t));
    rv |= memcmp(signature - sizeof(magic_t), magic.val, sizeof(magic_t));
    if (rv) {
        fprintf(stderr, "ERROR: Magic numbers do not mtach\n");
        goto err;
    }
#endif

    printf("verification passes as expected\n");
    ret = OQS_SUCCESS;
    goto cleanup;

err:
    ret = OQS_ERROR;

cleanup:
    if (secret_key) {
        OQS_MEM_secure_free(secret_key - sizeof(magic_t), sig->length_secret_key + 2 * sizeof(magic_t));
    }
    if (public_key) {
        OQS_MEM_insecure_free(public_key - sizeof(magic_t));
    }
    if (message) {
        OQS_MEM_insecure_free(message - sizeof(magic_t));
    }
    if (signature) {
        OQS_MEM_insecure_free(signature - sizeof(magic_t));
    }
    OQS_SIG_free(sig);

    return ret;
}

#ifdef OQS_ENABLE_TEST_CONSTANT_TIME
static void TEST_SIG_randombytes(uint8_t *random_array, size_t bytes_to_read) {
    // We can't make direct calls to the system randombytes on some platforms,
    // so we have to swap out the OQS_randombytes provider.
    OQS_randombytes_switch_algorithm("system");
    OQS_randombytes(random_array, bytes_to_read);
    OQS_randombytes_custom_algorithm(&TEST_SIG_randombytes);

    // OQS_TEST_CT_CLASSIFY tells Valgrind's memcheck tool to issue a warning if
    // the program branches on any byte that depends on random_array. This helps us
    // identify timing side-channels, as these bytes often contain secret data.
    OQS_TEST_CT_CLASSIFY(random_array, bytes_to_read);
}
#endif

#if OQS_USE_PTHREADS_IN_TESTS
struct thread_data {
    char *alg_name;
    OQS_STATUS rc;
};

void *test_wrapper(void *arg) {
    struct thread_data *td = arg;
    td->rc = sig_test_correctness(td->alg_name);
    return NULL;
}
#endif

static void fullcycle(uint8_t *sig, uint8_t *public_key, uint8_t *secret_key, uint8_t *signature, size_t signature_len, uint8_t *message, size_t message_len) {
    
    unsigned int n;
    size_t *sig_len = sig_len;
    
	if (mask_crypto_sign_keypair(public_key, secret_key) != OQS_SUCCESS) {
		printf("keygen error. Exiting.\n");
		exit(-1);
	}
	if (mask_crypto_sign(signature, &signature_len, message, message_len, secret_key) != OQS_SUCCESS) {
		printf("sign error. Exiting.\n");
		exit(-1);
	}
	if (mask_crypto_sign_verify(sig, *sig_len,message, message_len,public_key) != OQS_SUCCESS) {
		printf("verify error. Exiting.\n");
		exit(-1);
	}
}

static OQS_STATUS sig_speed_wrapper(const char *method_name, uint64_t duration, bool printInfo, bool doFullCycle,unsigned int n) {

	OQS_SIG *sig = NULL;
    size_t *sig_len = sig_len;
	uint8_t *public_key = NULL;
	uint8_t *secret_key = NULL;
	uint8_t *message = NULL;
	uint8_t *signature = NULL;
	size_t message_len = 50;
	size_t signature_len = 0;
	OQS_STATUS ret = OQS_ERROR;

	sig = OQS_SIG_new(method_name);
	if (sig == NULL) {
		return OQS_SUCCESS;
	}
    
    

	public_key = malloc(sig->length_public_key);
	secret_key = malloc(sig->length_secret_key);
	message = malloc(message_len);
	signature = malloc(sig->length_signature);

    // 各メモリ領域のバイト数を出力
    printf("Public Key Size: %zu bytes\n", sig->length_public_key);
    printf("Secret Key Size: %zu bytes\n", sig->length_secret_key);
    printf("Message Size: %zu bytes\n", message_len);
    printf("Signature Size: %zu bytes\n", sig->length_signature);
    
    if ((public_key == NULL) || (secret_key == NULL) || (message == NULL) || (signature == NULL)) {
		fprintf(stderr, "ERROR: malloc failed\n");
		goto err;
	}

	OQS_randombytes(message, message_len);

	printf("%-36s | %10s | %14s | %15s | %10s | %25s | %10s\n", sig->method_name, "", "", "", "", "", "");
	if (!doFullCycle) {
        printf("MAX is %d. \n",MAX_SIGN);
		TIME_OPERATION_SECONDS(mask_crypto_sign_keypair(public_key, secret_key), "keypair", duration)
        printf("keypair has done.\n");
		TIME_OPERATION_SECONDS(mask_crypto_sign(signature, &signature_len, message, message_len, secret_key), "sign", duration)
        printf("sign has done.\n");
		TIME_OPERATION_SECONDS(mask_crypto_sign_verify((const uint8_t *)sig, *sig_len,message, message_len,public_key), "verify", duration)
        printf("verify has done.\n");
	} else {
		TIME_OPERATION_SECONDS(fullcycle((uint8_t *)sig, public_key, secret_key, signature, signature_len, message, message_len), "fullcycle", duration)
	}


	if (printInfo) {
		printf("public key bytes: %zu, secret key bytes: %zu, signature bytes: %zu\n", sig->length_public_key, sig->length_secret_key, sig->length_signature);
		if (signature_len != sig->length_signature) {
			printf("   Actual signature length returned (%zu) less than declared maximum signature length (%zu)\n", signature_len, sig->length_signature);
		}
	}

	ret = OQS_SUCCESS;
	goto cleanup;

err:
	ret = OQS_ERROR;

cleanup:
	if (sig != NULL) {
		OQS_MEM_secure_free(secret_key, sig->length_secret_key);
	}
	OQS_MEM_insecure_free(public_key);
	OQS_MEM_insecure_free(signature);
	OQS_MEM_insecure_free(message);
	OQS_SIG_free(sig);

	return ret;
}

static OQS_STATUS printAlgs(void) {
	for (size_t i = 0; i < OQS_SIG_algs_length; i++) {
		OQS_SIG *sig = OQS_SIG_new(OQS_SIG_alg_identifier(i));
		if (sig == NULL) {
			printf("%s (disabled)\n", OQS_SIG_alg_identifier(i));
		} else {
			printf("%s\n", OQS_SIG_alg_identifier(i));
		}
		OQS_SIG_free(sig);
	}
	return OQS_SUCCESS;
}

int main(int argc, char **argv) {
    unsigned int global_n = 0;
	int ret = EXIT_SUCCESS;
	OQS_STATUS rc;

	bool printUsage = false;
	uint64_t duration = 3;
	bool printSigInfo = false;
	bool doFullCycle = false;
    unsigned int n;

	OQS_SIG *single_sig = NULL;

	OQS_init();
	OQS_randombytes_switch_algorithm(OQS_RAND_alg_openssl);
    


	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "--algs") == 0) {
			rc = printAlgs();
			if (rc == OQS_SUCCESS) {
				OQS_destroy();
				return EXIT_SUCCESS;
			} else {
				OQS_destroy();
				return EXIT_FAILURE;
			}
		} else if ((strcmp(argv[i], "--duration") == 0) || (strcmp(argv[i], "-d") == 0)) {
			if (i < argc - 1) {
				duration = (uint64_t)strtol(argv[i + 1], NULL, 10);
				if (duration > 0) {
					i += 1;
					continue;
				}
			}
		} else if ((strcmp(argv[i], "--help") == 0) || (strcmp(argv[i], "-h") == 0)) {
			printUsage = true;
			break;
		} else if ((strcmp(argv[i], "--info") == 0) || (strcmp(argv[i], "-i") == 0)) {
			printSigInfo = true;
			continue;
		} else if ((strcmp(argv[i], "--fullcycle") == 0) || (strcmp(argv[i], "-f") == 0)) {
			doFullCycle = true;
			continue;
		} else {
			single_sig = OQS_SIG_new(argv[i]);
			if (single_sig == NULL) {
				printUsage = true;
				break;
			}
		}
	}

	if (printUsage) {
		fprintf(stderr, "Usage: speed_sig <options> <alg>\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "<options>\n");
		fprintf(stderr, "--algs             Print supported algorithms and terminate\n");
		fprintf(stderr, "--duration n\n");
		fprintf(stderr, " -d n              Run each speed test for approximately n seconds, default n=3\n");
		fprintf(stderr, "--help\n");
		fprintf(stderr, " -h                Print usage\n");
		fprintf(stderr, "--info\n");
		fprintf(stderr, " -i                Print info (sizes, security level) about each SIG\n");
		fprintf(stderr, "--fullcycle\n");
		fprintf(stderr, " -f                Test full keygen-sign-verify cycle of each SIG\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "<alg>              Only run the specified SIG method; must be one of the algorithms output by --algs\n");
		OQS_destroy();
		return EXIT_FAILURE;
	}

	//print_system_info();

	printf("Speed test\n");
	printf("==========\n");

	PRINT_TIMER_HEADER
	if (single_sig != NULL) {
		rc = sig_speed_wrapper(single_sig->method_name, duration, printSigInfo, doFullCycle,n);
		if (rc != OQS_SUCCESS) {
			ret = EXIT_FAILURE;
		}
		OQS_SIG_free(single_sig);
        printf("n is %u \n", n);

	} else {
		for (size_t i = 0; i < OQS_SIG_algs_length; i++) {
			rc = sig_speed_wrapper(OQS_SIG_alg_identifier(i), duration, printSigInfo, doFullCycle,n);
			if (rc != OQS_SUCCESS) {
				ret = EXIT_FAILURE;
                printf("n is %u \n", n);
			}
		}
	}
	PRINT_TIMER_FOOTER
	OQS_destroy();

	return ret;
    
    
}
