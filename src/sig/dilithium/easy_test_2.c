#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <oqs/oqs.h>
#include "pqcrystals-dilithium_dilithium2_ref/params.h"
#include "pqcrystals-dilithium_dilithium2_ref/sign.h"


typedef struct magic_s {
    uint8_t val[31];
} magic_t;

static OQS_STATUS sig_test_correctness(const char *method_name) {
    OQS_SIG *sig = NULL;
    size_t *sig_len = sig_len;
    uint8_t *public_key = NULL;
    uint8_t *secret_key = NULL;
    uint8_t *message = NULL;
    size_t message_len = 5;
    uint8_t *signature = NULL;
    OQS_STATUS rc, ret = OQS_ERROR;
    size_t signature_len;

    cst_8t result;
    
    
    // マジックナンバーを初期化
    magic_t magic;
    OQS_randombytes(magic.val, sizeof(magic_t));

    sig = OQS_SIG_new("dilithium2");
    if (sig == NULL) {
        fprintf(stderr, "ERROR: OQS_SIG_new failed\n");
        goto err;
    }

    public_key = malloc(sig->length_public_key);
    secret_key = malloc(sig->length_secret_key);
    message = malloc(message_len);
    signature = malloc(sig->length_signature);
    
    printf("Public Key Size: %zu bytes\n", sig->length_public_key);
    printf("Secret Key Size: %zu bytes\n", sig->length_secret_key);
    printf("Message Size: %zu bytes\n", message_len);
    printf("Signature Size: %zu bytes\n", sig->length_signature);

    if ((public_key == NULL) || (secret_key == NULL) || (message == NULL) || (signature == NULL)) {
        fprintf(stderr, "ERROR: malloc failed\n");
        goto err;
    }

    OQS_randombytes(message, message_len);
    printf("message  is %s\n",message);
    
    rc = mask_crypto_sign_keypair(public_key, secret_key);
    if (rc != OQS_SUCCESS) {
        fprintf(stderr, "ERROR: OQS_SIG_keypair failed\n");
        goto err;
    }
    uint8_t *original_secret_key = secret_key;
    uint8_t *original_public_key = public_key;
    
    printf("sk is ");
    for (size_t i = 0; i < sig->length_secret_key; i++) {
        printf("%02x", original_secret_key[i]);
    }
    printf("\n");

    printf("pk is ");
    for (size_t i = 0; i < sig->length_public_key; i++) {
        printf("%02x", original_public_key[i]);
    }
    printf("\n");
    


    rc = mask_crypto_sign(signature, &signature_len, message, message_len, secret_key,public_key);
    printf("test num  is %u\n",rc);
    if (rc != OQS_SUCCESS) {
        fprintf(stderr, "ERROR: OQS_SIG_sign failed\n");
        goto err;
    }

    printf("sig is ");
    for (size_t i = 0; i < sig->length_signature; i++) {
        printf("%02x", signature[i]);
    }
    
    printf("\n");

    //sigは署名技術のステータスを示すものでsignatureに値が入っているのでそちらを入れるべきでは?//
    result = mask_crypto_sign_verify((uint8_t *)sig, *sig_len, message, message_len, public_key);
    
    for(int i = 0; i < SEEDBYTES; ++i){
        printf("c[%d]  is   %d /",i, result.a[i]);

        if(i%6 == 0) printf("\n");
    }
    printf("\n");
    for(int i = 0; i < SEEDBYTES; ++i){

        printf("c2[%d] is   %d /",i, result.b[i]);
        if(i%6 == 0) printf("\n");
    }
    printf("\n");
    printf("misstake is ");
    for(int i = 0; i < SEEDBYTES; ++i) {if(result.a[i] != result.b[i]) {printf("%d/",i);}}
    printf("\n");

    printf("\n");
    for(int i = 0; i < K*POLYW1_PACKEDBYTES; ++i){
        printf("sig[%d] is %d /",i, result.c[i]);

        if(i%6 == 0) printf("\n");

    }
    printf("\n");
    for(int i = 0; i < K*POLYW1_PACKEDBYTES; ++i){

        printf("buf[%d] is %d /",i, result.d[i]);
        if(i%6 == 0) printf("\n");

    }
    printf("\n");
    printf("misstake is ");
    for(int i = 0; i < K*POLYW1_PACKEDBYTES; ++i) {if(result.c[i] != result.d[i]) {printf("%d/",i);}}
    printf("\n");
    
    if(result.a[0] == result.b[0]) {
        rc = 0;
    }
    else rc = -1;
    
    
    
    if (rc != OQS_SUCCESS) {
        fprintf(stderr, "ERROR: OQS_SIG_verify failed\n");

        goto err;
    }

    ret = OQS_SUCCESS;
    goto cleanup;

err:
    ret = OQS_ERROR;

cleanup:
    free(secret_key);
    free(public_key);
    free(message);
    free(signature);
    OQS_SIG_free(sig);
    return ret;
}

int main(int argc, char **argv) {
    OQS_init();
    
    if (argc != 2) {
        fprintf(stderr, "Usage: test_sig algname\n");
        OQS_destroy();
        return EXIT_FAILURE;
    }
    
    char *alg_name = argv[1];
    if (!OQS_SIG_alg_is_enabled(alg_name)) {
        printf("Signature algorithm %s not enabled!\n", alg_name);
        OQS_destroy();
        return EXIT_FAILURE;
    }
    
    OQS_STATUS rc = sig_test_correctness(alg_name);
    if (rc != OQS_SUCCESS) {
        OQS_destroy();
        printf("sig hasn't done");
        return EXIT_FAILURE;
    }
    else{printf("sig has done");}
    
    OQS_destroy();
    
    return EXIT_SUCCESS;
    
    
}