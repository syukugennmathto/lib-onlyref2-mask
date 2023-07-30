#include <stdint.h>
#include "params.h"
#include "rounding.h"
#include "sign.h"
#include "packing.h"
#include "polyvec.h"
#include "poly.h"
#include "randombytes.h"
#include "symmetric.h"
#include "fips202.h"
#include "noise.h"


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
int mask_crypto_sign_keypair(uint8_t *pk, uint8_t *sk) {
    uint8_t seedbuf[2*SEEDBYTES + CRHBYTES];
    uint8_t tr[SEEDBYTES];
    uint8_t *rho, *rhoprime, *key;
    polyvecl mat[K];
    polyvecl s1, s1hat;
    polyveck s2, t1, t0;
    

    /* Get randomness for rho, rhoprime and key */
    randombytes(seedbuf, SEEDBYTES);
    shake256_noise(seedbuf, 2*SEEDBYTES + CRHBYTES, seedbuf, SEEDBYTES);
    rho = seedbuf;
    rhoprime = rho + SEEDBYTES;
    key = rhoprime + CRHBYTES;

    /* Expand matrix */
    polyvec_matrix_expand(mat, rho);

    /* Sample short vectors s1 and s2 */
    polyvecl_uniform_eta(&s1, rhoprime, 0);
    polyveck_uniform_eta(&s2, rhoprime, L);

    // mask the s1 and s2 vectors with small_bounded_noise_generation_256
    unsigned char nonce = 0;
    for (uint32_t i = 0; i < L; ++i)
        small_bounded_noise_generation_256((uint32_t*)s1.vec[i].coeffs, rhoprime, nonce++);

    nonce = 0;
    for (uint32_t i = 0; i < K; ++i)
        small_bounded_noise_generation_256((uint32_t*)s2.vec[i].coeffs, rhoprime, nonce++);

    /* Matrix-vector multiplication */
    s1hat = s1;
    polyvecl_ntt(&s1hat);
    polyvec_matrix_pointwise_montgomery(&t1, mat, &s1hat);
    polyveck_reduce(&t1);
    polyveck_invntt_tomont(&t1);

    /* Add error vector s2 */
    polyveck_add(&t1, &t1, &s2);

    /* Extract t1 and write public key */
    polyveck_caddq(&t1);
    polyveck_power2round(&t1, &t0, &t1);
        pack_pk(pk, rho, &t1);

    /* Compute H(rho, t1) and write secret key */
    shake256(tr, SEEDBYTES, pk, CRYPTO_PUBLICKEYBYTES);
    pack_sk(sk, rho, tr, key, &t0, &s1, &s2);

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
*     r0ってなんですか？
* Returns 0 (success)
**************************************************/
int mask_crypto_sign_signature(uint8_t *sig,
                          size_t *siglen,
                          const uint8_t *m,
                          size_t mlen,
                          const uint8_t *sk)
{
 
    unsigned int n;
    unsigned int nxt = 0;

  uint8_t seedbuf[3*SEEDBYTES + 2*CRHBYTES];
  uint8_t *rho, *tr, *key, *mu, *rhoprime;
  uint16_t nonce = 0;
  polyvecl mat[K], s1, y, z;
  polyveck t0, s2, w1, w0, h;
  poly cp;
  shake256incctx state;
    long long int cycle[8] = {0,0,0,0,0,0,0,0},MAX = MAX_SIGN;
    int times =0;
    int ret = 0;
    
    uint32_t a0_w1[K * N];
    uint32_t a1_w1[K * N];
    uint32_t a0_w0[K * N];
    uint32_t a1_w0[K * N];

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



  /* Expand matrix and transform vectors */
  polyvec_matrix_expand(mat, rho);
  polyvecl_ntt(&s1);
  polyveck_ntt(&s2);
  polyveck_ntt(&t0);
  
#ifdef DILITHIUM_RANDOMIZED_SIGNING
  randombytes(rhoprime, CRHBYTES);
#else
    shake256_noise(rhoprime, CRHBYTES, key, SEEDBYTES + CRHBYTES);
#endif
    for (uint32_t i = 0; i < dilithium_l; ++i)large_bounded_noise_generation(&s1.vec[i], rhoprime, nonce++);
    for (uint32_t i = 0; i < dilithium_k; ++i)large_bounded_noise_generation(&s2.vec[i], rhoprime, nonce++);
 
    

rej:
    
   
    times = times + 1;
    cycle[6] = 0;
    for(int i =0;i<7;++i){cycle[6] = cycle[6] + cycle[i];}

    

    
  /*Sample intermediate vector y line5 */
  polyvecl_uniform_gamma1(&y, rhoprime, nonce++);

  /*Matrix-vector multiplication line6 */
  z = y;
  polyvecl_ntt(&z);
  polyvec_matrix_pointwise_montgomery(&w1, mat, &z);
  polyveck_reduce(&w1);
  polyveck_invntt_tomont(&w1);
    /* whta is ? */
    polyveck_caddq(&w1);
    
    
    /*ここにhigh bits関数が必要 w1,2/gamma line7 */
    for (uint32_t i = 0; i < K; ++i)     masking_arithmetic_to_boolean_highbits((uint32_t *)(w1.vec[i].coeffs), (uint32_t *)(w1.vec[i].coeffs), 2 * GAMMA2);
    
    
    /*masking-function befor polyveck_decompose(&w1, &w0, &w1);*/
    // forループを使用してw1の各要素に対してdecompose関数を呼び出す
    /* for (uint32_t i = 0; i < K; ++i) {
         cycle[0] = cycle[0] + 1;
            decompose_array(a0_w1, a1_w1, w1.vec[i].coeffs, K * N);
    }
    
    
    for (uint32_t i = 0; i < K; ++i) {
        cycle[0] = cycle[0] + 1;
        prime_decompose(&a0_w1[i * N], &a1_w1[i * N], w1.vec[i].coeffs[0]);
    }
     for (uint32_t i = 0; i < K; ++i) {
         cycle[0] = cycle[0] + 1;
         masking_arithmetic_to_boolean_decompose(&a0_w1[i * N], &a1_w1[i * N], (uint32_t*)w1.vec[i].coeffs, 32);
         masking_arithmetic_to_boolean_decompose(&a0_w0[i * N], &a1_w0[i * N], (uint32_t*)w0.vec[i].coeffs, 32);
 // ここでmasking_arithmetic_to_boolean_decomposeを呼び出す
     }
     */
    

    
    cycle[6] = 0;
    for(int i =0;i<7;++i){cycle[6] = cycle[6] + cycle[i];}
    if(cycle[6] > MAX){ret = 1; goto fish;}
    
    // forループを使用してw0の各要素に対してdecompose関数を呼び出す
    /*for (uint32_t i = 0; i < K; ++i) {
        
            decompose_array(a0_w0, a1_w0, w0.vec[i].coeffs, K * N);
    }
   
    
    for (uint32_t i = 0; i < K; ++i) {
        cycle[0] = cycle[0] + 1;
        prime_decompose(&a0_w0[i * N], &a1_w0[i * N], w0.vec[i].coeffs[0]);
    }
     */
    


    
    cycle[6] = 0;
    for(int i =0;i<7;++i){cycle[6] = cycle[6] + cycle[i];}
    if(cycle[6] > MAX){ret = 1; goto fish;}
    
    cycle[1] = cycle[1] + 1;
    
  polyveck_pack_w1(sig, &w1);

  /*hash function  Hfunction line8? */
  shake256_inc_ctx_reset(&state);
  shake256_inc_absorb(&state, mu, CRHBYTES);
  shake256_inc_absorb(&state, sig, K*POLYW1_PACKEDBYTES);
  shake256_inc_finalize(&state);
  shake256_inc_squeeze(sig, SEEDBYTES, &state);
  poly_challenge(&cp, sig);
  poly_ntt(&cp);
  /* Compute z, reject if it reveals secret line9 */
  polyvecl_pointwise_poly_montgomery(&z, &cp, &s1);
  polyvecl_invntt_tomont(&z);
  polyvecl_add(&z, &z, &y);
  polyvecl_reduce(&z);
    

    
    cycle[6] = 0;
    for(int i =0;i<7;++i){cycle[6] = cycle[6] + cycle[i];}
    if(cycle[6] > MAX){ret = 1; goto fish;}
    
    /* line11 */
    if(polyvecl_chknorm(&z, GAMMA1 - BETA)){ goto rej;}
    


  /* Check that subtracting cs2 does not change high bits of w and low bits
   * do not reveal secret information */
  polyveck_pointwise_poly_montgomery(&h, &cp, &s2);
  polyveck_invntt_tomont(&h);
  polyveck_sub(&w0, &w0, &h);
  polyveck_reduce(&w0);
    
    
    /*lowbits line10*/
    for (uint32_t i = 0; i < K; ++i) masking_arithmetic_to_boolean_lowbits((uint32_t *)(w0.vec[i].coeffs), (uint32_t *)(w0.vec[i].coeffs), 2 * GAMMA2);
    
    cycle[2] = cycle[2] + 1;
    
    cycle[6] = 0;
    for(int i =0;i<7;++i){cycle[6] = cycle[6] + cycle[i];}
    if(cycle[6] > MAX){ret = 1; goto fish;}
    
    /* line12 */
    if(polyveck_chknorm(&w0, GAMMA2 - BETA)){ goto rej;}
    
    /* line12 */
    if(polyveck_chknorm(&h, GAMMA2)){ goto rej;}
    

  /* Compute hints for w1 */
  polyveck_pointwise_poly_montgomery(&h, &cp, &t0);
  polyveck_invntt_tomont(&h);
  polyveck_reduce(&h);
    
    cycle[3] = cycle[3] + 1;
    
    cycle[6] = 0;
    for(int i =0;i<7;++i){cycle[6] = cycle[6] + cycle[i];}
    if(cycle[6] > MAX){ret = 1; goto fish;}
    


   polyveck_add(&w0, &w0, &h);
   nxt = polyveck_make_hint(&h, &w0, &w1);

    
    
    /* line13 */
    if(nxt > OMEGA){goto rej;}
    

    

    
    cycle[6] = 0;
    for(int i =0;i<7;++i){cycle[6] = cycle[6] + cycle[i];}
    if(cycle[6] > MAX){ret = 1; goto fish;}
    
    cycle[4] = cycle[4] + 1;
    

    
    cycle[5] = cycle[5] + 1;

  shake256_inc_ctx_release(&state);

  /* Write signature */
  pack_sig(sig, sig, &z, &h);
  *siglen = CRYPTO_BYTES;
    
fish:
    return ret;
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
int mask_crypto_sign(uint8_t *sm,size_t *smlen,const uint8_t *m,size_t mlen,const uint8_t *sk)
{
  size_t i;

    for(i = 0; i < mlen; ++i){
        sm[CRYPTO_BYTES + mlen - 1 - i] = m[mlen - 1 - i];
    }
    mask_crypto_sign_signature(sm, smlen, sm + CRYPTO_BYTES, mlen, sk);
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
int mask_crypto_sign_verify(const uint8_t *sig,
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

  if(siglen != CRYPTO_BYTES)
    return -1;

  unpack_pk(rho, &t1, pk);
  if(unpack_sig(c, &z, &h, sig))
    return -1;
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
  if(mask_crypto_sign_verify(sm, CRYPTO_BYTES, sm + CRYPTO_BYTES, *mlen, pk))
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
