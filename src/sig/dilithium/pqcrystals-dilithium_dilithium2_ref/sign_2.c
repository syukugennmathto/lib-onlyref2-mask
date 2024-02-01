#include <stdint.h>
#include <stdio.h>

#include "params.h"
#include "sign_2.h"
#include "packing.h"
#include "polyvec.h"
#include "poly.h"
#include "randombytes.h"
#include "symmetric.h"
#include "fips202.h"

/*************************************************
* Name:        crypto_sign_keypair
*
* Description: Generates public and private key.
*
* Arguments:   - uint8_t *pk: pointer to output public key (allocated
*                             array of CRYPTO_PUBLICKEYBYTES bytes)
*              - uint8_t *sk: pointer to output private key (allocated
*                             array of CRYPTO_SECRETKEYBYTES bytes)
*
* Returns 0 (success)
**************************************************/
int masking_crypto_sign_keypair(uint8_t *pk, uint8_t *sk) {
    uint8_t seedbuf[2*SEEDBYTES + CRHBYTES] = {0}; // すべての要素を0で初期化
    uint8_t rho2[CRYPTO_PUBLICKEYBYTES];
  uint8_t tr[SEEDBYTES],k;
    uint32_t i;
    uint16_t Ll= L;
  const uint8_t *rho, *rhoprime,*rhoprime2,*key;
  polyvecl mat[K];
  polyvecl s11, s12,s1hat;
  polyveck s21,s22, t1, t0;

  /* Get randomness for rho, rhoprime and key */
  randombytes(seedbuf, SEEDBYTES);
  shake256(seedbuf, 2*SEEDBYTES + CRHBYTES, seedbuf, SEEDBYTES);
  rho = seedbuf;
  rhoprime = rho + SEEDBYTES;
  rhoprime2 =rho + SEEDBYTES;
  key = rhoprime + CRHBYTES;

  /* Expand matrix */
  polyvec_matrix_expand(mat, rho);

  /* Sample short vectors s1 and s2 */
    //after masking
    polyvecl_uniform_eta(&s11, rhoprime, 0);
    //printf("gene1\n");
    polyvecl_small_bounded_noise_generation_256(&s12,(const uint8_t *) rhoprime, 0);
    //printf("gene2\n");
  polyveck_uniform_eta(&s21, rhoprime, L);
    
    for (size_t i = 0; i < L; ++i) {
        for (size_t j = 0; j < N; ++j) {
            s11.vec[i].coeffs[j] = s11.vec[i].coeffs[j]%3;
            //printf("[%zu][%zu] = %d and %d\n", i,j ,s11.vec[i].coeffs[j],s12.vec[i].coeffs[j]);
            
        }
    }
    
    polyveck_small_bounded_noise_generation_256(&s22,(const uint8_t *)rhoprime2, Ll);
    //printf("gene3\n");
  /* Matrix-vector multiplication */
  s1hat = s11;
    //printf("gene4\n");
  polyvecl_ntt(&s1hat);
  polyvec_matrix_pointwise_montgomery(&t1, mat, &s1hat);
  polyveck_reduce(&t1);
  polyveck_invntt_tomont(&t1);

  /* Add error vector s2 */
    //printf("gene5\n");
  polyveck_add(&t1, &t1, &s21);

  /* Extract t1 and write public key */
  polyveck_caddq(&t1);
  polyveck_power2round(&t1, &t0, &t1);
    //printf("gene6\n");
    for (size_t k = 0; k < 2*SEEDBYTES + CRHBYTES; ++k) {
      //  printf("seedbuf[%zu] = %u\n", k, seedbuf[k]);

    }
    
    for (size_t k = 0; k <2*SEEDBYTES + CRHBYTES; ++k)rho2[k] = seedbuf[k];


  pack_pk(pk,rho2, &t1);

  /* Compute H(rho, t1) and write secret key */
    //printf("gene7\n");
  shake256(tr, SEEDBYTES, pk, CRYPTO_PUBLICKEYBYTES);
    //printf("gene8\n");
  pack_sk(sk, rho, tr, key, &t0, &s11, &s21);

  return 0;
}

/*************************************************
* Name:        crypto_sign_signature
*
* Description: Computes signature.
*
* Arguments:   - uint8_t *sig:   pointer to output signature (of length CRYPTO_BYTES)
*              - size_t *siglen: pointer to output length of signature
*              - uint8_t *m:     pointer to message to be signed
*              - size_t mlen:    length of message
*              - uint8_t *sk:    pointer to bit-packed secret key
*
* Returns 0 (success)
**************************************************/
int masking_crypto_sign_signature(uint8_t *sig,
                          size_t *siglen,
                          const uint8_t *m,
                          size_t mlen,
                          const uint8_t *sk)
{
    int n=0;
  unsigned int i,j;
    int rc;
  uint8_t seedbuf[3*SEEDBYTES + 2*CRHBYTES];
  uint8_t *rho, *tr, *key, *mu, *rhoprime;
  uint16_t nonce = 0;
  polyvecl mat[K], s1, y, z,z1;
  polyveck t0, s2,w4, w3,w2,w1, w0, h,h1,h2,h3;
  poly cp;
  shake256incctx state;
  maskpointveck pD;
   // printf("here0\n");
  rho = seedbuf;
  tr = rho + SEEDBYTES;
  key = tr + SEEDBYTES;
  mu = key + SEEDBYTES;
  rhoprime = mu + CRHBYTES;
  unpack_sk(rho, tr, key, &t0, &s1, &s2, sk);

  /* Compute CRH(tr, msg) */
  shake256_inc_init(&state);
  shake256_inc_absorb(&state, tr, SEEDBYTES);
  shake256_inc_absorb(&state, m, mlen);
  shake256_inc_finalize(&state);
  shake256_inc_squeeze(mu, CRHBYTES, &state);

#ifdef DILITHIUM_RANDOMIZED_SIGNING
  randombytes(rhoprime, CRHBYTES);
#else
  shake256(rhoprime, CRHBYTES, key, SEEDBYTES + CRHBYTES);
#endif

  /* Expand matrix and transform vectors */
  polyvec_matrix_expand(mat, rho);
  polyvecl_ntt(&s1);
  polyveck_ntt(&s2);
  polyveck_ntt(&t0);
   // printf("here1\n");
rej:
   // printf("this is res\n");
  /* Sample intermediate vector y */
  //polyvecl_uniform_gamma1(&y, rhoprime, nonce++);
    //after masking
    polycecl_large_bounded_noise_generation(&y, rhoprime, nonce++);
    
  /* Matrix-vector multiplication */
  z = y;
  polyvecl_ntt(&z);
  polyvec_matrix_pointwise_montgomery(&w1, mat, &z);
  polyveck_reduce(&w1);
  polyveck_invntt_tomont(&w1);

  /* Decompose w and call the random oracle */
  polyveck_caddq(&w1);
  //  printf("here2\n");
  //polyveck_decompose(&w1, &w0, &w1);
  //highbits is w1,lowbits is w0
    w4=w1;
    polyveck_masking_decompose(&w0, &w1, &w1,&pD);
    polyveck_masking_arithmetic_to_boolean_highbits(&w1,&w4,&pD);
    //for(int i = 0;i<K;++i){for(int j = 0;j<N;++j){printf("w is %u\n",w1.vec[i].coeffs[j]);}}
    //new_nosking_arithmetic_to_boolean_highbits(int32_t* pR1,const int32_t* pR,mask_point pD)
  //  printf("this is decom\n");
    polyveck_unmasking_arithmetic_to_boolean_highbits(&w4,&w1,&pD);
    polyveck_pack_w1(sig, &w4);

  shake256_inc_ctx_reset(&state);
  shake256_inc_absorb(&state, mu, CRHBYTES);
  shake256_inc_absorb(&state, sig, K*POLYW1_PACKEDBYTES);
  shake256_inc_finalize(&state);
  shake256_inc_squeeze(sig, SEEDBYTES, &state);
  poly_challenge(&cp, sig);
  poly_ntt(&cp);
    /*uint32_t time =0;
    for(int k=0;k<N;++k){
        printf("cp[%d] is %d\n",k,cp.coeffs[k]);
        for(int kt=0;kt<L;++kt){
            printf("s1[%d][%d] is %d\n",k,kt,s1.vec[k].coeffs[kt]);
            printf("times is %d\n",cp.coeffs[k]*s1.vec[k].coeffs[kt]);
            if((cp.coeffs[k]*s1.vec[k].coeffs[kt] < GAMMA1 - BETA)&&(-1*cp.coeffs[k]*s1.vec[k].coeffs[kt] < GAMMA1 - BETA)){
                printf("this is out[%u]\n",time);
                time +=1;
            }
        }
        
    }
    return 1;
    */
  /* Compute z, reject if it reveals secret */
  polyvecl_pointwise_poly_montgomery(&z, &cp, &s1);
  polyvecl_invntt_tomont(&z);
  polyvecl_add(&z, &z, &y);
  polyvecl_reduce(&z);
    //ここの機構を関数化する
    for(int i = 0;i<L;++i){for(int j = 0;j<N;++j){z1.vec[i].coeffs[j]=z.vec[i].coeffs[j]%(GAMMA1 - BETA-1);}}
//    printf("this is z\n");
  //  printf("bound is %d\n",GAMMA1 - BETA);
    for(int i = 0;i<L;++i){for(int j = 0;j<N;++j){//printf("z is %u\n",z1.vec[i].coeffs[j]);
        
    }}
 //   printf("rc is %d\n",rc);
    
    if(polyvecl_chknorm(&z1, GAMMA1 - BETA)) goto rej;

    

  /* Check that subtracting cs2 does not change high bits of w and low bits
   * do not reveal secret information */
  polyveck_pointwise_poly_montgomery(&h, &cp, &s2);
 //   printf("this is 1w0\n");
    h1 = h;
  polyveck_invntt_tomont(&h1);
 //   printf("this is 2w0\n");
  polyveck_sub(&w0, &w0, &h1);
 //   printf("this is 3w0\n");
  polyveck_reduce(&w0);
 //   printf("this is 4w0\n");
    for(int i = 0;i<K;++i){for(int j = 0;j<N;++j){w1.vec[i].coeffs[j]=w0.vec[i].coeffs[j]%(GAMMA2 - BETA-1);}}
    for(int i = 0;i<K;++i){for(int j = 0;j<N;++j){//printf("w is %u\n",w1.vec[i].coeffs[j]);
        
    }}
  if(polyveck_chknorm(&w1, GAMMA2 - BETA)) goto rej;
   // printf("this is w0\n");
    
    
  /* Compute hints for w1 */
  polyveck_pointwise_poly_montgomery(&h1, &cp, &t0);
  polyveck_invntt_tomont(&h1);
  polyveck_reduce(&h1);
    
    
    for(int i = 0;i<K;++i){for(int j = 0;j<N;++j){h2.vec[i].coeffs[j]=h1.vec[i].coeffs[j]%GAMMA2;}}
    for(int i = 0;i<K;++i){for(int j = 0;j<N;++j){//printf("h is %u\n",h2.vec[i].coeffs[j]);
        
    }}
    
  if(polyveck_chknorm(&h2, GAMMA2))goto rej;
 //   printf("this is h\n");
  polyveck_add(&w3, &w1, &h2);
  //n = polyveck_make_hint(&h, &w0, &w1);
    polyveck_masking_make_hint(&h2, &w3, &w1,&pD);
    
    for(int i = 0;i<K;++i){for(int j = 0;j<N;++j){//printf("h is %u\n",h2.vec[i].coeffs[j]);
        
    }}
    for(i = 0; i < K; ++i){
        for(j=0;j<N;++j){
            if(h2.vec[i].coeffs[j]!=0){
                n +=1;
            }
        }
    }
   // printf("n and OMEGA are %d and %d\n",n,OMEGA);
    
    //if(n > OMEGA) goto rej;
  //  printf("this is n\n");
  shake256_inc_ctx_release(&state);

  /* Write signature */
  pack_sig(sig, sig, &z1, &h2);
  *siglen = masking_CRYPTO_BYTES;
  return 0;
}

/*************************************************
* Name:        crypto_sign
*
* Description: Compute signed message.
*
* Arguments:   - uint8_t *sm: pointer to output signed message (allocated
*                             array with CRYPTO_BYTES + mlen bytes),
*                             can be equal to m
*              - size_t *smlen: pointer to output length of signed
*                               message
*              - const uint8_t *m: pointer to message to be signed
*              - size_t mlen: length of message
*              - const uint8_t *sk: pointer to bit-packed secret key
*
* Returns 0 (success)
**************************************************/
int masking_crypto_sign(uint8_t *sm,
                size_t *smlen,
                const uint8_t *m,
                size_t mlen,
                const uint8_t *sk)
{
  size_t i;

  for(i = 0; i < mlen; ++i)
    sm[CRYPTO_BYTES + mlen - 1 - i] = m[mlen - 1 - i];
  masking_crypto_sign_signature(sm, smlen, sm + masking_CRYPTO_BYTES, mlen, sk);
  *smlen += mlen;
  return 0;
}

/*************************************************
* Name:        crypto_sign_verify
*
* Description: Verifies signature.
*
* Arguments:   - uint8_t *m: pointer to input signature
*              - size_t siglen: length of signature
*              - const uint8_t *m: pointer to message
*              - size_t mlen: length of message
*              - const uint8_t *pk: pointer to bit-packed public key
*
* Returns 0 if signature could be verified correctly and -1 otherwise
**************************************************/
int masking_crypto_sign_verify(const uint8_t *sig,
                       size_t siglen,
                       const uint8_t *m,
                       size_t mlen,
                       const uint8_t *pk)
{
  unsigned int i;
  uint8_t buf[K*POLYW1_PACKEDBYTES];
  uint8_t rho[SEEDBYTES];
  uint8_t mu[CRHBYTES];
  uint8_t c[SEEDBYTES];
  uint8_t c2[SEEDBYTES];
  poly cp;
  polyvecl mat[K], z;
  polyveck t1, w1, h;
  shake256incctx state;
  maskpointveck pD;
    printf("verify1\n");
  if(siglen != CRYPTO_BYTES)
    return -1;
    printf("verify2\n");
  unpack_pk(rho, &t1, pk);
  if(unpack_sig(c, &z, &h,sig))
    return -1;
    printf("verify3\n");
  if(polyvecl_chknorm(&z, GAMMA1 - BETA))
    return -1;

  /* Compute CRH(H(rho, t1), msg) */
  shake256(mu, SEEDBYTES, pk, CRYPTO_PUBLICKEYBYTES);
  shake256_inc_init(&state);
  shake256_inc_absorb(&state, mu, SEEDBYTES);
  shake256_inc_absorb(&state, m, mlen);
  shake256_inc_finalize(&state);
  shake256_inc_squeeze(mu, CRHBYTES, &state);

  /* Matrix-vector multiplication; compute Az - c2^dt1 */
  poly_challenge(&cp, c);
  polyvec_matrix_expand(mat, rho);

  polyvecl_ntt(&z);
  polyvec_matrix_pointwise_montgomery(&w1, mat, &z);

  poly_ntt(&cp);
  polyveck_shiftl(&t1);
  polyveck_ntt(&t1);
  polyveck_pointwise_poly_montgomery(&t1, &cp, &t1);

  polyveck_sub(&w1, &w1, &t1);
  polyveck_reduce(&w1);
  polyveck_invntt_tomont(&w1);

  /* Reconstruct w1 */
  polyveck_caddq(&w1);
  polyveck_use_hint(&w1, &w1, &h);
    for(int i = 0;i<K;++i){for(int j = 0;j<N;++j){printf("w is %u\n",w1.vec[i].coeffs[j]);}}
    //polyveck_masking_use_hint(&w1,&h,&pD);
    
    polyveck_pack_w1(buf, &w1);

  /* Call random oracle and verify challenge */
  shake256_inc_ctx_reset(&state);
  shake256_inc_absorb(&state, mu, CRHBYTES);
  shake256_inc_absorb(&state, buf, K*POLYW1_PACKEDBYTES);
  shake256_inc_finalize(&state);
  shake256_inc_squeeze(c2, SEEDBYTES, &state);
  shake256_inc_ctx_release(&state);
  for(i = 0; i < SEEDBYTES; ++i)
    if(c[i] != c2[i])
      return -1;

  return 0;
}

/*************************************************
* Name:        crypto_sign_open
*
* Description: Verify signed message.
*
* Arguments:   - uint8_t *m: pointer to output message (allocated
*                            array with smlen bytes), can be equal to sm
*              - size_t *mlen: pointer to output length of message
*              - const uint8_t *sm: pointer to signed message
*              - size_t smlen: length of signed message
*              - const uint8_t *pk: pointer to bit-packed public key
*
* Returns 0 if signed message could be verified correctly and -1 otherwise
**************************************************/
int crypto_sign_open(uint8_t *m,
                     size_t *mlen,
                     const uint8_t *sm,
                     size_t smlen,
                     const uint8_t *pk)
{
  size_t i;

  if(smlen < CRYPTO_BYTES)
    goto badsig;

  *mlen = smlen - CRYPTO_BYTES;
  if(masking_crypto_sign_verify(sm, CRYPTO_BYTES, sm + CRYPTO_BYTES, *mlen, pk))
    goto badsig;
  else {
    /* All good, copy msg, return 0 */
    for(i = 0; i < *mlen; ++i)
      m[i] = sm[CRYPTO_BYTES + i];
    return 0;
  }

badsig:
  /* Signature verification failed */
  *mlen = -1;
  for(i = 0; i < smlen; ++i)
    m[i] = 0;

  return -1;
}
