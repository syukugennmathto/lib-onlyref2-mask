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
    uint8_t seedbuf1[SEEDBYTES];
    uint8_t tr[SEEDBYTES];
    uint8_t *rho, *rhoprime, *key;
    polyvecl mat[K];
    polyvecl s1, s1hat;
    polyveck s2, t1, t0;
    

    /* Get randomness for rho, rhoprime and key */
    randombytes(seedbuf, SEEDBYTES);
    randombytes(seedbuf1, SEEDBYTES);
    
    shake256_noise(seedbuf, 2*SEEDBYTES + CRHBYTES, seedbuf, SEEDBYTES);
    shake256_noise(seedbuf1, SEEDBYTES, seedbuf1, SEEDBYTES);
    rho = seedbuf1; //rho = seedbuf;
    rhoprime = rho; //rhoprime = rho + SEEDBYTES
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
    pack_pk(pk, rho, &t1); //pack_pk(pk, rho, &t1);

    /* Compute H(rho, t1)=tr and write secret key */
    shake256(tr, SEEDBYTES, pk, CRYPTO_PUBLICKEYBYTES);
    pack_sk(sk, rhoprime, tr, key, &t0, &s1, &s2); //pack_pk(pk, rho, &t1);

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
*
* Returns 0 (success)
**************************************************/
int mask_crypto_sign_signature(uint8_t *sig,
                          size_t *siglen,
                          const uint8_t *m,
                          size_t mlen,
                          const uint8_t *sk,
                          const uint8_t *pk)
{
 
    unsigned int n;
    uint32_t nxt = 0;
    uint8_t buf[K*POLYW1_PACKEDBYTES];
  uint8_t seedbuf[3*SEEDBYTES + 2*CRHBYTES];
  uint8_t *rho, *tr, *key, *mu, *rhoprime;
  uint32_t nonce = -1,ret = 0;
  polyvecl mat[K], s1, y, z;
  polyveck t3,t2,t0, s2, w5,w4,w3,w2,w1, w0, h;
  poly cp;
  shake256incctx state;
    long long int cycle[8] = {0,0,0,0,0,0,0,0},MAX = MAX_SIGN;
    int times =0;
    uint8_t c[SEEDBYTES];
    uint8_t c2[SEEDBYTES];

    
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
  unpack_pk(rho, &t2, pk);
 

  /* Compute CRH(tr, message) */
  shake256_inc_init(&state);
  shake256_inc_absorb(&state, tr, SEEDBYTES);
  shake256_inc_absorb(&state, m, mlen);
  shake256_inc_finalize(&state);
  shake256_inc_squeeze(mu, CRHBYTES, &state);

#ifdef DILITHIUM_RANDOMIZED_SIGNING
  randombytes(rhoprime, CRHBYTES);
#else
    shake256_noise(rhoprime, CRHBYTES, key, SEEDBYTES + CRHBYTES);
#endif


  /* Expand matrix and transform vectors */
  polyvec_matrix_expand(mat, rho);
  polyvecl_ntt(&s1);
  polyveck_ntt(&s2);
  polyveck_ntt(&t0);
    
    /* t0=t-t1*2^d line3*/
  

    //for (uint32_t i = 0; i < dilithium_l; ++i)large_bounded_noise_generation(&s1.vec[i], rhoprime, nonce++);
    //for (uint32_t i = 0; i < dilithium_k; ++i)large_bounded_noise_generation(&s2.vec[i], rhoprime, nonce++);
    /* line4 */
    
    

rej:
    nonce = nonce +1;
    if(nonce > MAX){
        ret = 1;
        goto fish;
        
    }
   


  //Sample intermediate vector y line5
  polyvecl_uniform_gamma1(&y, rhoprime, nonce);

  //Matrix-vector multiplication line6
  z = y;
  polyvecl_ntt(&z);
  polyvec_matrix_pointwise_montgomery(&w1, mat, &z);
  polyveck_reduce(&w1);
  polyveck_invntt_tomont(&w1);

    polyveck_caddq(&w1);

    
    //ここにhigh bits関数が必要 w1,2 gamma line7
    for (uint32_t i = 0; i < K; ++i) {
        for (size_t j = 0; j < N; j++){
            num_masking_arithmetic_to_boolean_highbits(w2.vec[i].coeffs[j],w1.vec[i].coeffs[j], alpha);
        }
    }
    

    
  polyveck_pack_w1(sig, &w2);

  /*hash function  Hfunction line8? */
  shake256_inc_init(&state);
  shake256_inc_ctx_reset(&state);
  shake256_inc_absorb(&state, mu, CRHBYTES);
  shake256_inc_absorb(&state, sig, K*POLYW1_PACKEDBYTES);
  shake256_inc_finalize(&state);
  shake256_inc_squeeze(sig, SEEDBYTES, &state);
  poly_challenge(&cp, sig);
  poly_ntt(&cp);
    
    /*hash function  Hfunction line8? */
    shake256_inc_init(&state);
    shake256_inc_ctx_reset(&state);
    shake256_inc_absorb(&state, mu, CRHBYTES);
    shake256_inc_absorb(&state, sig, K*POLYW1_PACKEDBYTES); // bufを追加
    shake256_inc_finalize(&state); // stateを最終化
    shake256_inc_squeeze(c, SEEDBYTES, &state);
    
  /* Compute z, reject if it reveals secret line9 */
  polyvecl_pointwise_poly_montgomery(&z, &cp, &s1);
  polyvecl_invntt_tomont(&z);
  polyvecl_add(&z, &z, &y);
  polyvecl_reduce(&z);
    

    
    /* line11 */
    if(polyvecl_chknorm(&z, GAMMA1 - BETA)){ goto rej;}
    
    /*for (uint32_t j = 0; j < N; ++j){
        for (uint32_t i = 0; i < L; ++i){
            if(number_chknorm(z.vec[i].coeffs[j], (GAMMA1 - BETA))){ goto rej;
            }
        }
    }*/
    


  /* Check that subtracting cs2 does not change high bits of w and low bits
   * do not reveal secret information */
  polyveck_pointwise_poly_montgomery(&h, &cp, &s2);
  polyveck_invntt_tomont(&h);
  polyveck_sub(&w0, &w0, &h);
  polyveck_reduce(&w0);
    
    
    /*lowbits line10*/
    for (uint32_t i = 0; i < K; ++i) {
        for (size_t j = 0; j < N; j++){
            num_masking_arithmetic_to_boolean_lowbits(w0.vec[i].coeffs[j],w0.vec[i].coeffs[j], alpha);
            if(number_chknorm(w0.vec[i].coeffs[j],GAMMA2)){ goto rej;}
        }
    }

    
    /* line12 */
    if(polyveck_chknorm(&w0, GAMMA2 - BETA)){ goto rej;}
    
    /* line12 */
    if(polyveck_chknorm(&h, GAMMA2)){ goto rej;}
    

  /* Compute hints for w0=c*t0*/
  polyveck_pointwise_poly_montgomery(&h, &cp, &t0);
  polyveck_invntt_tomont(&h);
  polyveck_reduce(&h);
    /* Compute hints for w0=w-cs2+ct0*/
  polyveck_add(&w0, &w0, &h);
    



    for (size_t i = 0; i < K; i++) {
        for (size_t j = 0; j < N; j++) {
            h.vec[i].coeffs[j] = num_masking_arithmetic_makeint(w0.vec[i].coeffs[j],w1.vec[i].coeffs[j], alpha);
        }
    }


    
    for (uint32_t j = 0; j < N; j++){
        for (uint32_t i = 0; i < K; i++){
            nxt = nxt + h.vec[i].coeffs[j];
        }
    }
    
    /* line13 */
    if(nxt > OMEGA){goto rej;}
    
    /*計算ミス防止用start*/
    /*Az =w4*/
    polyvec_matrix_pointwise_montgomery(&w4, mat, &z);
      
      
    
    /*times 2^D*/
    polyveck_shiftl(&t2);
    /*c*t1*/
    polyveck_ntt(&t2);
    polyveck_pointwise_poly_montgomery(&t2, &cp, &t2);
    /*Az-ct1*2^D*/
    polyveck_sub(&w4, &w1, &t2);
    polyveck_reduce(&w4);
    polyveck_invntt_tomont(&w4);
    
    for (size_t i = 0; i < K; i++) {
        for (size_t j = 0; j < N; j++) {
            w5.vec[i].coeffs[j] = num_masking_arithmetic_usehint(h.vec[i].coeffs[j],w4.vec[i].coeffs[j],alpha);
        }
    }
   polyveck_pack_w1(buf, &w5);
    
    shake256_inc_init(&state);
    shake256_inc_ctx_reset(&state); // stateをリセット
    shake256_inc_absorb(&state, mu, CRHBYTES); // muを追加
    shake256_inc_absorb(&state, buf, K*POLYW1_PACKEDBYTES); // bufを追加
    shake256_inc_finalize(&state); // stateを最終化
    shake256_inc_squeeze(c2, SEEDBYTES, &state); // c2を得る
    /*計算ミス防止用end*/
               //for(int i = 0; i < K*POLYW1_PACKEDBYTES; i++){if(buf[i] != sig[i]){goto rej;}}
               //for(int i = 0; i < K*POLYW1_PACKEDBYTES; i++){if(buf[i] != sig[i]){return -1;}}
    for(uint32_t i = 0; i < SEEDBYTES; ++i)if(c2[i] != c[i])goto rej;

    //for(int i = 0; i < K*POLYW1_PACKEDBYTES; i++)sig[i]=buf[i];


  shake256_inc_ctx_release(&state);

fish:
  /* Write signature */
  pack_sig(sig,sig, &z, &h);
  *siglen = CRYPTO_BYTES;
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
int mask_crypto_sign(uint8_t *sm,size_t *smlen,const uint8_t *m,size_t mlen,const uint8_t *sk,const uint8_t *pk)
{
  size_t i;
    uint32_t xt = 0;

    for(i = 0; i < mlen; ++i){
        sm[CRYPTO_BYTES + mlen - 1 - i] = m[mlen - 1 - i];
    }
   xt = mask_crypto_sign_signature(sm, smlen, sm + CRYPTO_BYTES, mlen, sk,pk);
  *smlen += mlen;
  return xt;
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
cst_8t mask_crypto_sign_verify(const uint8_t *sig,
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
  polyveck t1, w5,w4,w3,w2,w1, w0, h;
  shake256incctx state;
  cst_8t result;
    uint32_t nxt = 0;

  //if(siglen != CRYPTO_BYTES)return siglen;

  unpack_pk(rho, &t1, pk);
  if(unpack_sig(c, &z, &h, sig)) return result;
    
    /*line28 */
    for (uint32_t j = 0; j < N; ++j){
        for (uint32_t i = 0; i < L; ++i){
            if(number_chknorm(z.vec[i].coeffs[j], GAMMA1 - BETA)){
                return result;
            }
        }
    }
    for (uint32_t i = 0; i < SEEDBYTES; ++i){
        c2[i] = 0;
    }

  /* Compute CRH(H(rho, t1), msg) line28 */
  shake256_inc_init(&state);
  shake256(mu, SEEDBYTES, pk, CRYPTO_PUBLICKEYBYTES);

  shake256_inc_absorb(&state, mu, SEEDBYTES);
  shake256_inc_absorb(&state, m, mlen);
  shake256_inc_finalize(&state);
  shake256_inc_squeeze(mu, CRHBYTES, &state);

  /* Matrix-vector multiplication; compute Az - c2^dt1 in line27 */
  //poly_challenge(&cp, c);
  polyvec_matrix_expand(mat, rho);

    /*Az =w1*/
    polyvecl_ntt(&z);
    polyvec_matrix_pointwise_montgomery(&w4, mat, &z);
    poly_ntt(&cp);
    /*times 2^D*/
    polyveck_shiftl(&t1);
    polyveck_ntt(&t1);
    polyveck_pointwise_poly_montgomery(&t1, &cp, &t1);
    polyveck_sub(&w4, &w4, &t1);
    polyveck_reduce(&w4);
    polyveck_invntt_tomont(&w4);
    polyveck_caddq(&w4);
    for (size_t i = 0; i < L; i++) {
        for (size_t j = 0; j < N; j++) {
            w5.vec[i].coeffs[j] = num_masking_arithmetic_usehint(h.vec[i].coeffs[j],w4.vec[i].coeffs[j],alpha);
        }
    }
   polyveck_pack_w1(buf, &w5);
    
    
    


    
    for(i = 0; i < K*POLYW1_PACKEDBYTES; ++i){
        result.c[i] =sig[i];
        result.d[i] =buf[i];
    }
    

    
    /* Call random oracle and verify challenge line28*/
    shake256_inc_init(&state);
    shake256_inc_ctx_reset(&state); // stateをリセット
    shake256_inc_absorb(&state, mu, CRHBYTES); // muを追加
    shake256_inc_absorb(&state, buf, K*POLYW1_PACKEDBYTES); // bufを追加
    shake256_inc_finalize(&state); // stateを最終化
    shake256_inc_squeeze(c2, SEEDBYTES, &state); // c2を得る
     
    
    c[0]=0;
    c2[0]=0;
    for(i = 0; i < SEEDBYTES; ++i){
        result.a[i] = c[i];
        result.b[i] = c2[i];
    }


    
    for(i = 0; i < SEEDBYTES; ++i)
      if(c[i] != c2[i])
        return result;

    return result;
  /*for(i = 0; i < SEEDBYTES; ++i){
      
      result.b[i] = 0;
  }

  for(i = 0; i < K*POLYW1_PACKEDBYTES; ++i){
      result.c[i] = 0;
      result.d[i] = 0;
  }*/
    //for(i = 0; i < SEEDBYTES; ++i){if(c[i] != c2[i]) {result.a[0] = -4;result.b[0] = -5;}}
    

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
  /*if(mask_crypto_sign_verify(sm, CRYPTO_BYTES, sm + CRYPTO_BYTES, *mlen, pk))
    goto badsig;
  else {
    // All good, copy msg, return 0
    for(i = 0; i < *mlen; ++i)
      m[i] = sm[CRYPTO_BYTES + i];
    return 0;
  }
    */
badsig:
  /* Signature verification failed */
  *mlen = -1;
  for(i = 0; i < smlen; ++i)
    m[i] = 0;

  return -1;
}

int number_chknorm(int32_t a, int32_t B){
  unsigned int i;
  int32_t t;

  if(B > (Q-1)/8)
    return 1;

  for(i = 0; i < N; ++i) {
    /* Absolute value */
    t = a >> 31;
    t = a - (t & 2*a);

    if(t >= B) {
      return 1;
    }
  }
 return 0;
}
