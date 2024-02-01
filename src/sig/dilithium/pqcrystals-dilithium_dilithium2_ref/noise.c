#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include "../../common/pqclean_shims/fips202.h"
#include "noise.h"
#include "random.h"
#include "params.h"
#include "poly.h"

#define BLOCKS 3
#define NROUNDS 24
#define ROL(a, offset) ((a << offset) ^ (a >> (64-offset)))


 //#define DBENCH_START()
//#define DBENCH_STOP(t)



void sam(uint32_t* out_buffer, unsigned char *seed) /// len in 32 bits
{
    unsigned char* it = (unsigned char*)out_buffer;
    shake128_noise(it, dilithium_n * dilithium_k * dilithium_l * sizeof(uint32_t), seed, 32);
}

uint32_t small_noise_rejection_256(uint32_t* out_buffer, size_t out_len, uint8_t* buffer, size_t in_len)
{
    uint32_t ctr, pos;
    unsigned char t0, t1;

    ctr = pos = 0;
    while(ctr < out_len) {
        t0 = buffer[pos] & 0x0F;
        t1 = buffer[pos++] >> 4;

        if(t0 <= 2*low_noise_bound)
        out_buffer[ctr++] = low_noise_bound - t0;
        if((t1 <= (2*low_noise_bound)) && (ctr < out_len))
        out_buffer[ctr++] = low_noise_bound - t1;
        if(pos >= in_len)
        break;
    }
    return ctr;
}

void small_bounded_noise_generation_256(int32_t *out_buffer, const uint8_t *seed,uint16_t nonce)
{    unsigned int i;
    uint32_t ctr;
    unsigned char inbuf[SEEDBYTES + 1];
    unsigned char outbuf[2 * SHAKE256_RATE];
    uint64_t state[25];
    
    
    for (i = 0; i < SEEDBYTES; ++i) inbuf[i] = (unsigned char)seed[i];
    inbuf[SEEDBYTES] = nonce;
    shake256_absorb(state, inbuf, SEEDBYTES + 1);
    shake256_squeezeblocks_noise(outbuf, 2, state);
    ctr = small_noise_rejection_256((uint32_t *)out_buffer, N, outbuf, 2 * SHAKE256_RATE);

    if (ctr < N) {
        shake256_squeezeblocks_noise(outbuf, 1, state);
        small_noise_rejection_256((uint32_t *)out_buffer + ctr, N - ctr, outbuf, SHAKE256_RATE);
    }
}




// 修正：large_noise_rejection_256関数を修正
uint32_t large_noise_rejection_256(uint32_t *out_buffer, size_t out_len, unsigned char *buffer, size_t in_len){
    
    uint32_t ctr, pos;
    uint32_t t0, t1;

    ctr = pos = 0;
    while (ctr < out_len)
    {
        t0 = buffer[pos];
        t0 |= (uint32_t)buffer[pos + 1] << 8;
        t0 |= (uint32_t)buffer[pos + 2] << 16;
        t0 &= MAX_VALUE_2_13;

        t1 = buffer[pos + 2] >> 4;
        t1 |= (uint32_t)buffer[pos + 3] << 4;
        t1 |= (uint32_t)buffer[pos + 4] << 12;
        t1 &= MAX_VALUE_2_13;

        pos += 5;

        if (t0 <= (2 * large_noise_bound))
            out_buffer[ctr++] = large_noise_bound - t0;
        if ((t1 <= (2 * large_noise_bound)) && (ctr < out_len))
            out_buffer[ctr++] = large_noise_bound - t1;

        if (pos >= in_len)
            break;
    }
    return ctr;
}




// 修正：large_bounded_noise_generation関数を修正
void large_bounded_noise_generation(poly *r, const uint8_t *seed,uint16_t nonce)
{
    unsigned int i, ctr;
    unsigned char inbuf[SEEDBYTES + CRHBYTES + 2];
    unsigned char outting_buf[5 * 136];
    uint64_t state[25];

    for (i = 0; i < SEEDBYTES + CRHBYTES; ++i)
        inbuf[i] = seed[i];
    inbuf[SEEDBYTES + CRHBYTES + 0] = nonce & 0xFF;
    inbuf[SEEDBYTES + CRHBYTES + 1] = nonce >> 8;
    shake256_absorb(state, inbuf, SEEDBYTES + CRHBYTES + 2);
    shake256_squeezeblocks_noise(outting_buf, 5, state); // change : outbuf -> outting_buf

    ctr = large_noise_rejection_256((uint32_t *)r->coeffs, 256, outting_buf, 5 * SHAKE256_RATE); // 修正：large_noise_rejection -> large_noise_rejection_256
    if (ctr < 256)
    {
        shake256_squeezeblocks_noise(outting_buf, 1, state); // 修正：outbuf -> outting_buf
        ctr += large_noise_rejection_256((uint32_t *)(outting_buf + ctr), 256 - ctr, outting_buf, SHAKE256_RATE); // 修正：large_noise_rejection -> large_noise_rejection_256
    }
}









/// --------------------------------------------------------------------------- ///

uint32_t small_noise_rejection(uint32_t* out_buffer, size_t out_len, unsigned char* buffer, size_t in_len)
{
    uint32_t ctr, pos;
    unsigned char t0, t1;

    ctr = pos = 0;
    while(ctr < out_len) {
        t0 = buffer[pos] & 0x0F;
        t1 = buffer[pos++] >> 4;

        if(t0 <= 2*low_noise_bound)
        out_buffer[ctr++] = low_noise_bound - t0;
        if((t1 <= (2*low_noise_bound)) && (ctr < out_len))
        out_buffer[ctr++] = low_noise_bound - t1;

        if(pos >= in_len)
        break;
    }
    return ctr;
}




uint32_t large_noise_rejection(uint32_t* out_buffer, size_t out_len, unsigned char* buffer, size_t in_len)
{
    uint32_t ctr, pos;
    uint32_t t0,t1;

    ctr = pos = 0;
    while(ctr < out_len) {

        t0  = buffer[pos];
        t0 |= (uint32_t)buffer[pos + 1] << 8;
        t0 |= (uint32_t)buffer[pos + 2] << 16;
        t0 &= 0x7FFFF;

        t1  = buffer[pos + 2] >> 4;
        t1 |= (uint32_t)buffer[pos + 3] << 4;
        t1 |= (uint32_t)buffer[pos + 4] << 12;
        t1 &= 0x7FFFF;

        pos += 5;

        if(t0 <= (2*large_noise_bound))
        out_buffer[ctr++] = large_noise_bound - t0;
        if((t1 <= (2*large_noise_bound)) && (ctr < out_len))
        out_buffer[ctr++] = large_noise_bound - t1;

        if(pos >= in_len)
            break;
    }
    return ctr;
}

static uint64_t load64(const unsigned char *x) {
  unsigned int i;
  uint64_t r = 0;

  for (i = 0; i < 8; ++i)
    r |= (uint64_t)x[i] << 8*i;

  return r;
}

static void store64(unsigned char *x, uint64_t u) {
  unsigned int i;

  for(i = 0; i < 8; ++i)
    x[i] = u >> 8*i;
}

static const uint64_t KeccakF_RoundConstants[NROUNDS] = {
  (uint64_t)0x0000000000000001ULL,
  (uint64_t)0x0000000000008082ULL,
  (uint64_t)0x800000000000808aULL,
  (uint64_t)0x8000000080008000ULL,
  (uint64_t)0x000000000000808bULL,
  (uint64_t)0x0000000080000001ULL,
  (uint64_t)0x8000000080008081ULL,
  (uint64_t)0x8000000000008009ULL,
  (uint64_t)0x000000000000008aULL,
  (uint64_t)0x0000000000000088ULL,
  (uint64_t)0x0000000080008009ULL,
  (uint64_t)0x000000008000000aULL,
  (uint64_t)0x000000008000808bULL,
  (uint64_t)0x800000000000008bULL,
  (uint64_t)0x8000000000008089ULL,
  (uint64_t)0x8000000000008003ULL,
  (uint64_t)0x8000000000008002ULL,
  (uint64_t)0x8000000000000080ULL,
  (uint64_t)0x000000000000800aULL,
  (uint64_t)0x800000008000000aULL,
  (uint64_t)0x8000000080008081ULL,
  (uint64_t)0x8000000000008080ULL,
  (uint64_t)0x0000000080000001ULL,
  (uint64_t)0x8000000080008008ULL
};

void KeccakF1600_StatePermute(uint64_t *state)
{
        int round;

        uint64_t Aba, Abe, Abi, Abo, Abu;
        uint64_t Aga, Age, Agi, Ago, Agu;
        uint64_t Aka, Ake, Aki, Ako, Aku;
        uint64_t Ama, Ame, Ami, Amo, Amu;
        uint64_t Asa, Ase, Asi, Aso, Asu;
        uint64_t BCa, BCe, BCi, BCo, BCu;
        uint64_t Da, De, Di, Do, Du;
        uint64_t Eba, Ebe, Ebi, Ebo, Ebu;
        uint64_t Ega, Ege, Egi, Ego, Egu;
        uint64_t Eka, Eke, Eki, Eko, Eku;
        uint64_t Ema, Eme, Emi, Emo, Emu;
        uint64_t Esa, Ese, Esi, Eso, Esu;

        //copyFromState(A, state)
        Aba = state[ 0];
        Abe = state[ 1];
        Abi = state[ 2];
        Abo = state[ 3];
        Abu = state[ 4];
        Aga = state[ 5];
        Age = state[ 6];
        Agi = state[ 7];
        Ago = state[ 8];
        Agu = state[ 9];
        Aka = state[10];
        Ake = state[11];
        Aki = state[12];
        Ako = state[13];
        Aku = state[14];
        Ama = state[15];
        Ame = state[16];
        Ami = state[17];
        Amo = state[18];
        Amu = state[19];
        Asa = state[20];
        Ase = state[21];
        Asi = state[22];
        Aso = state[23];
        Asu = state[24];

        for( round = 0; round < NROUNDS; round += 2 )
        {
            //    prepareTheta
            BCa = Aba^Aga^Aka^Ama^Asa;
            BCe = Abe^Age^Ake^Ame^Ase;
            BCi = Abi^Agi^Aki^Ami^Asi;
            BCo = Abo^Ago^Ako^Amo^Aso;
            BCu = Abu^Agu^Aku^Amu^Asu;

            //thetaRhoPiChiIotaPrepareTheta(round  , A, E)
            Da = BCu^ROL(BCe, 1);
            De = BCa^ROL(BCi, 1);
            Di = BCe^ROL(BCo, 1);
            Do = BCi^ROL(BCu, 1);
            Du = BCo^ROL(BCa, 1);

            Aba ^= Da;
            BCa = Aba;
            Age ^= De;
            BCe = ROL(Age, 44);
            Aki ^= Di;
            BCi = ROL(Aki, 43);
            Amo ^= Do;
            BCo = ROL(Amo, 21);
            Asu ^= Du;
            BCu = ROL(Asu, 14);
            Eba =   BCa ^((~BCe)&  BCi );
            Eba ^= (uint64_t)KeccakF_RoundConstants[round];
            Ebe =   BCe ^((~BCi)&  BCo );
            Ebi =   BCi ^((~BCo)&  BCu );
            Ebo =   BCo ^((~BCu)&  BCa );
            Ebu =   BCu ^((~BCa)&  BCe );

            Abo ^= Do;
            BCa = ROL(Abo, 28);
            Agu ^= Du;
            BCe = ROL(Agu, 20);
            Aka ^= Da;
            BCi = ROL(Aka,  3);
            Ame ^= De;
            BCo = ROL(Ame, 45);
            Asi ^= Di;
            BCu = ROL(Asi, 61);
            Ega =   BCa ^((~BCe)&  BCi );
            Ege =   BCe ^((~BCi)&  BCo );
            Egi =   BCi ^((~BCo)&  BCu );
            Ego =   BCo ^((~BCu)&  BCa );
            Egu =   BCu ^((~BCa)&  BCe );

            Abe ^= De;
            BCa = ROL(Abe,  1);
            Agi ^= Di;
            BCe = ROL(Agi,  6);
            Ako ^= Do;
            BCi = ROL(Ako, 25);
            Amu ^= Du;
            BCo = ROL(Amu,  8);
            Asa ^= Da;
            BCu = ROL(Asa, 18);
            Eka =   BCa ^((~BCe)&  BCi );
            Eke =   BCe ^((~BCi)&  BCo );
            Eki =   BCi ^((~BCo)&  BCu );
            Eko =   BCo ^((~BCu)&  BCa );
            Eku =   BCu ^((~BCa)&  BCe );

            Abu ^= Du;
            BCa = ROL(Abu, 27);
            Aga ^= Da;
            BCe = ROL(Aga, 36);
            Ake ^= De;
            BCi = ROL(Ake, 10);
            Ami ^= Di;
            BCo = ROL(Ami, 15);
            Aso ^= Do;
            BCu = ROL(Aso, 56);
            Ema =   BCa ^((~BCe)&  BCi );
            Eme =   BCe ^((~BCi)&  BCo );
            Emi =   BCi ^((~BCo)&  BCu );
            Emo =   BCo ^((~BCu)&  BCa );
            Emu =   BCu ^((~BCa)&  BCe );

            Abi ^= Di;
            BCa = ROL(Abi, 62);
            Ago ^= Do;
            BCe = ROL(Ago, 55);
            Aku ^= Du;
            BCi = ROL(Aku, 39);
            Ama ^= Da;
            BCo = ROL(Ama, 41);
            Ase ^= De;
            BCu = ROL(Ase,  2);
            Esa =   BCa ^((~BCe)&  BCi );
            Ese =   BCe ^((~BCi)&  BCo );
            Esi =   BCi ^((~BCo)&  BCu );
            Eso =   BCo ^((~BCu)&  BCa );
            Esu =   BCu ^((~BCa)&  BCe );

            //    prepareTheta
            BCa = Eba^Ega^Eka^Ema^Esa;
            BCe = Ebe^Ege^Eke^Eme^Ese;
            BCi = Ebi^Egi^Eki^Emi^Esi;
            BCo = Ebo^Ego^Eko^Emo^Eso;
            BCu = Ebu^Egu^Eku^Emu^Esu;

            //thetaRhoPiChiIotaPrepareTheta(round+1, E, A)
            Da = BCu^ROL(BCe, 1);
            De = BCa^ROL(BCi, 1);
            Di = BCe^ROL(BCo, 1);
            Do = BCi^ROL(BCu, 1);
            Du = BCo^ROL(BCa, 1);

            Eba ^= Da;
            BCa = Eba;
            Ege ^= De;
            BCe = ROL(Ege, 44);
            Eki ^= Di;
            BCi = ROL(Eki, 43);
            Emo ^= Do;
            BCo = ROL(Emo, 21);
            Esu ^= Du;
            BCu = ROL(Esu, 14);
            Aba =   BCa ^((~BCe)&  BCi );
            Aba ^= (uint64_t)KeccakF_RoundConstants[round+1];
            Abe =   BCe ^((~BCi)&  BCo );
            Abi =   BCi ^((~BCo)&  BCu );
            Abo =   BCo ^((~BCu)&  BCa );
            Abu =   BCu ^((~BCa)&  BCe );

            Ebo ^= Do;
            BCa = ROL(Ebo, 28);
            Egu ^= Du;
            BCe = ROL(Egu, 20);
            Eka ^= Da;
            BCi = ROL(Eka, 3);
            Eme ^= De;
            BCo = ROL(Eme, 45);
            Esi ^= Di;
            BCu = ROL(Esi, 61);
            Aga =   BCa ^((~BCe)&  BCi );
            Age =   BCe ^((~BCi)&  BCo );
            Agi =   BCi ^((~BCo)&  BCu );
            Ago =   BCo ^((~BCu)&  BCa );
            Agu =   BCu ^((~BCa)&  BCe );

            Ebe ^= De;
            BCa = ROL(Ebe, 1);
            Egi ^= Di;
            BCe = ROL(Egi, 6);
            Eko ^= Do;
            BCi = ROL(Eko, 25);
            Emu ^= Du;
            BCo = ROL(Emu, 8);
            Esa ^= Da;
            BCu = ROL(Esa, 18);
            Aka =   BCa ^((~BCe)&  BCi );
            Ake =   BCe ^((~BCi)&  BCo );
            Aki =   BCi ^((~BCo)&  BCu );
            Ako =   BCo ^((~BCu)&  BCa );
            Aku =   BCu ^((~BCa)&  BCe );

            Ebu ^= Du;
            BCa = ROL(Ebu, 27);
            Ega ^= Da;
            BCe = ROL(Ega, 36);
            Eke ^= De;
            BCi = ROL(Eke, 10);
            Emi ^= Di;
            BCo = ROL(Emi, 15);
            Eso ^= Do;
            BCu = ROL(Eso, 56);
            Ama =   BCa ^((~BCe)&  BCi );
            Ame =   BCe ^((~BCi)&  BCo );
            Ami =   BCi ^((~BCo)&  BCu );
            Amo =   BCo ^((~BCu)&  BCa );
            Amu =   BCu ^((~BCa)&  BCe );

            Ebi ^= Di;
            BCa = ROL(Ebi, 62);
            Ego ^= Do;
            BCe = ROL(Ego, 55);
            Eku ^= Du;
            BCi = ROL(Eku, 39);
            Ema ^= Da;
            BCo = ROL(Ema, 41);
            Ese ^= De;
            BCu = ROL(Ese, 2);
            Asa =   BCa ^((~BCe)&  BCi );
            Ase =   BCe ^((~BCi)&  BCo );
            Asi =   BCi ^((~BCo)&  BCu );
            Aso =   BCo ^((~BCu)&  BCa );
            Asu =   BCu ^((~BCa)&  BCe );
        }

        //copyToState(state, A)
        state[ 0] = Aba;
        state[ 1] = Abe;
        state[ 2] = Abi;
        state[ 3] = Abo;
        state[ 4] = Abu;
        state[ 5] = Aga;
        state[ 6] = Age;
        state[ 7] = Agi;
        state[ 8] = Ago;
        state[ 9] = Agu;
        state[10] = Aka;
        state[11] = Ake;
        state[12] = Aki;
        state[13] = Ako;
        state[14] = Aku;
        state[15] = Ama;
        state[16] = Ame;
        state[17] = Ami;
        state[18] = Amo;
        state[19] = Amu;
        state[20] = Asa;
        state[21] = Ase;
        state[22] = Asi;
        state[23] = Aso;
        state[24] = Asu;
}

static void keccak_absorb(uint64_t *s,
                          unsigned int r,
                          const unsigned char *m,
                          unsigned long long mlen,
                          unsigned char p)
{
  unsigned int i;
  unsigned char t[200];

  for(i = 0; i < 25; ++i)
    s[i] = 0;

  while(mlen >= r) {
    for(i = 0; i < r/8; ++i)
      s[i] ^= load64(m + 8*i);

    KeccakF1600_StatePermute(s);
    mlen -= r;
    m += r;
  }

  for(i = 0; i < r; ++i)
    t[i] = 0;
  for(i = 0; i < mlen; ++i)
    t[i] = m[i];
  t[i] = p;
  t[r-1] |= 128;
  for(i = 0; i < r/8; ++i)
    s[i] ^= load64(t + 8*i);
}

static void keccak_squeezeblocks(unsigned char *h,
                                 unsigned long nblocks,
                                 uint64_t *s,
                                 unsigned int r)
{
  unsigned int i;

  while(nblocks > 0) {
    KeccakF1600_StatePermute(s);
    for(i=0; i < (r >> 3); i++) {
      store64(h + 8*i, s[i]);
    }
    h += r;
    --nblocks;
  }
}

void shake128_absorb(uint64_t *s,
                     const unsigned char *input,
                     unsigned long long inlen)
{
  keccak_absorb(s, 168, input, inlen, 0x1F);
}

void shake128_squeezeblocks_noise(unsigned char *output,
                            unsigned long nblocks,
                            uint64_t *s)
{
  keccak_squeezeblocks(output, nblocks, s, 168);
}

void shake256_absorb(uint64_t *s,
                     const unsigned char *input,
                     unsigned long long inlen)
{
  keccak_absorb(s, SHAKE256_RATE, input, inlen, 0x1F);
}

void shake256_squeezeblocks_noise(unsigned char *output,
                            unsigned long nblocks,
                            uint64_t *s)
{
  keccak_squeezeblocks(output, nblocks, s, SHAKE256_RATE);
}


void shake128_noise(unsigned char *output,
              unsigned long long outlen,
              const unsigned char *input,
              unsigned long long inlen)
{
  unsigned int i;
  unsigned long nblocks = outlen/168;
  unsigned char t[168];
  uint64_t s[25];

  shake128_absorb(s, input, inlen);
    shake128_squeezeblocks_noise(output, nblocks, s);

  output += nblocks*168;
  outlen -= nblocks*168;

  if(outlen) {
      shake128_squeezeblocks_noise(t, 1, s);
    for(i = 0; i < outlen; ++i)
      output[i] = t[i];
  }
}

void shake256_noise(unsigned char *output,
              unsigned long long outlen,
              const unsigned char *input,
              unsigned long long inlen)
{
  unsigned int i;
  unsigned long nblocks = outlen/SHAKE256_RATE;
  unsigned char t[SHAKE256_RATE];
  uint64_t s[25];

  shake256_absorb(s, input, inlen);
  shake256_squeezeblocks_noise(output, nblocks, s);

  output += nblocks*SHAKE256_RATE;
  outlen -= nblocks*SHAKE256_RATE;

  if(outlen) {
    shake256_squeezeblocks_noise(t, 1, s);
    for(i = 0; i < outlen; ++i)
      output[i] = t[i];
  }
}



/*************************************************
* Name:        rej_eta
*
* Description: Sample uniformly random coefficients in [-ETA, ETA] by
*              performing rejection sampling on array of random bytes.
*
* Arguments:   - int32_t *a: pointer to output array (allocated)
*              - unsigned int len: number of coefficients to be sampled
*              - const uint8_t *buf: array of random bytes
*              - unsigned int buflen: length of array of random bytes
*
* Returns number of sampled coefficients. Can be smaller than len if not enough
* random bytes were given.
**************************************************/
static unsigned int rej_eta(int32_t *a, size_t len, const uint8_t *buf, size_t buflen)
{
    unsigned int ctr, pos;
    
  //DBENCH_START();
    uint32_t t0, t1;


  ctr = pos = 0;
  while(ctr < len && pos < buflen) {
    t0 = buf[pos] & 0x0F;
    t1 = buf[pos++] >> 4;

#if ETA == 2
    if(t0 < 15) {
      t0 = t0 - (205*t0 >> 10)*5;
      a[ctr++] = 2 - t0;
    }
    if(t1 < 15 && ctr < len) {
      t1 = t1 - (205*t1 >> 10)*5;
      a[ctr++] = 2 - t1;
    }
#elif ETA == 4
    if(t0 < 9)
      a[ctr++] = 4 - t0;
    if(t1 < 9 && ctr < len)
      a[ctr++] = 4 - t1;
#endif
  }

  //DBENCH_STOP(*tsample);
  return ctr;
}
