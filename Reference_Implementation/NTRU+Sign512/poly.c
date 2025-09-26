#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "ntt.h"
#include "symmetric.h"
#include "sampler.h"
#include "gaussian.h"
#include <stdlib.h>
#include "fft.h"

#ifdef PARAM_1
#define QINV 62209 // q^(-1) mod 2^16
#endif

#ifdef PARAM_2
#define QINV 57857 // q^(-1) mod 2^16
#endif

static inline int16_t barrett_reduce(int16_t a)
{
  int16_t t;
  const int16_t v = ((1<<26) + NTRUPLUS_Q/2)/NTRUPLUS_Q;

  t  = ((int32_t)v*a + (1<<25)) >> 26;
  t *= NTRUPLUS_Q;
  return a - t;
}

static inline int16_t caddq(int16_t a) {
    a += (a >> 15) & NTRUPLUS_Q;
    return a;
}

static inline int16_t csubq(int16_t a) {
    int16_t t = ((a-NTRUPLUS_Q) >> 15);
    a -= ~t & NTRUPLUS_Q;
    return a;
}

static inline int16_t caddp(int16_t a) {
    a += (a >> 15) & NTRUPLUS_P;
    return a;
}

static inline int16_t csubp(int16_t a) {
    int16_t t = ((a-NTRUPLUS_P) >> 15);
    a -= ~t & NTRUPLUS_P;
    return a;
}

static inline int max(int16_t a, int16_t b) {
  int16_t t = (a - b) >> 15;
  return (a & ~t) | (b & t);
}

static inline int16_t freezep(int16_t a) {
    a += ((a + halfp) >> 15) & NTRUPLUS_P;
    a -= ((halfp - a) >> 15) & NTRUPLUS_P;
    return a;
}

void poly_cbd1(poly *r, const unsigned char seed[CRHBYTES], uint16_t nonce)
{
    
    uint8_t buf[NTRUPLUS_N/4];
    stream128_state state;

    stream128_init(&state, seed, nonce);
    shake128_squeeze(buf, NTRUPLUS_N/4, &state);

    for(int i = 0; i < NTRUPLUS_N/8; i++)
    {   
        r->coeffs[8*i]   =  (buf[2*i]       & 0x1) -  (buf[2*i+1]       & 0x1);
        r->coeffs[8*i+1] = ((buf[2*i] >> 1) & 0x1) - ((buf[2*i+1] >> 1) & 0x1);
        r->coeffs[8*i+2] = ((buf[2*i] >> 2) & 0x1) - ((buf[2*i+1] >> 2) & 0x1);
        r->coeffs[8*i+3] = ((buf[2*i] >> 3) & 0x1) - ((buf[2*i+1] >> 3) & 0x1);
        r->coeffs[8*i+4] = ((buf[2*i] >> 4) & 0x1) - ((buf[2*i+1] >> 4) & 0x1);
        r->coeffs[8*i+5] = ((buf[2*i] >> 5) & 0x1) - ((buf[2*i+1] >> 5) & 0x1);
        r->coeffs[8*i+6] = ((buf[2*i] >> 6) & 0x1) - ((buf[2*i+1] >> 6) & 0x1);
        r->coeffs[8*i+7] = ((buf[2*i] >> 7) & 0x1) - ((buf[2*i+1] >> 7) & 0x1);
    }
}

uint64_t innerproduct(poly *a1, poly *a2)
{
    int res = 0;
    
    for (int i = 0; i < NTRUPLUS_N ; ++i)
    {
        res += a1->coeffs[i] * a2->coeffs[i];
    }
           
    return res;
}

void poly_ntt(poly *r)
{
    ntt(r->coeffs);
}

void poly_invntt(poly *r)
{
    invntt(r->coeffs);
}

void poly_add(poly *r, const poly *a, const poly *b)
{
    for(int i = 0; i < NTRUPLUS_N; ++i)
        r->coeffs[i] = a->coeffs[i] + b->coeffs[i];
}

void poly_add_qhat(poly *r, const poly *a)
{
    for(int i = 0; i < NTRUPLUS_N/4; ++i)
    {
        r->coeffs[4*i  ] = a->coeffs[4*i] + qhat;
        r->coeffs[4*i+1] = a->coeffs[4*i+1];
        r->coeffs[4*i+2] = a->coeffs[4*i+2];
        r->coeffs[4*i+3] = a->coeffs[4*i+3];
    }
}

void poly_caddq(poly *v)
{
    for (int i = 0; i < NTRUPLUS_N; ++i)
        v->coeffs[i] = caddq(v->coeffs[i]);
}

void poly_csubq(poly *v)
{
    for (int i = 0; i < NTRUPLUS_N; ++i)
        v->coeffs[i] = csubq(v->coeffs[i]);
}

void poly_caddp(poly *v)
{
    for (int i = 0; i < NTRUPLUS_N; ++i)
        v->coeffs[i] = caddp(v->coeffs[i]);
}

void poly_csubp(poly *v)
{
    for (int i = 0; i < NTRUPLUS_N; ++i)
        v->coeffs[i] = csubp(v->coeffs[i]);
}

void poly_freezep(poly *v)
{
    for (int i = 0; i < NTRUPLUS_N; ++i)
        v->coeffs[i] = freezep(v->coeffs[i]);
}


uint16_t poly_eqcheck(const poly *a, const poly *b){
    
    uint16_t ret = 0;

    for (int i = 0; i < NTRUPLUS_N; ++i)
    {
        ret |= a->coeffs[i] - b->coeffs[i];
    }
    return ret;
}

uint16_t c_eqcheck(const uint8_t c[CRYPTO_C_BYTES], const uint8_t c2[CRYPTO_C_BYTES]){
    
    uint16_t ret = 0;

    for (int i = 0; i < CRYPTO_C_BYTES; ++i)
    {
        ret |= c[i] - c2[i];
    }

    return ret;
}


void poly_bc(poly *v, const poly *c, uint8_t b)
{
    b = 1-b;
    for (int i = 0; i < NTRUPLUS_N; ++i)
    {
        v->coeffs[i] = c->coeffs[i] & b;
    }
}


void poly_cqhat(poly *v)
{
    for (int i = 0; i < NTRUPLUS_N; ++i)
        v->coeffs[i] *= qhat;
}

void poly_barrett_reduce(poly *v)
{
    for (int i = 0; i < NTRUPLUS_N; ++i)
        v->coeffs[i] = barrett_reduce(v->coeffs[i]);

}

void poly_sub(poly *r, const poly *a, const poly *b)
{
    for(int i = 0; i < NTRUPLUS_N; ++i)
        r->coeffs[i] = a->coeffs[i] - b->coeffs[i];
}

void poly_cneg(poly *r, const uint8_t b)
{
    for (int i = 0; i < NTRUPLUS_N; i++)
    {
        r->coeffs[i] *= 1 - 2 * b;
    }
}

void poly_neg(poly *r)
{
    for (int i = 0; i < NTRUPLUS_N; i++)
    {
        r->coeffs[i] *= -1;
    }
}

void poly_shift(poly *r, const poly *a)
{
    for (int i = 0; i < NTRUPLUS_N; i++)
    {
        r->coeffs[i] = a->coeffs[i] << NTRUPLUS_d;
    }
}

uint64_t poly_sqnorm2(const poly *a, const poly *b)
{
    
    uint64_t ret = 0;

    for (int i = 0; i < NTRUPLUS_N; ++i)
    {
            ret += (uint64_t)a->coeffs[i] * a->coeffs[i] + b->coeffs[i] * b->coeffs[i];
    }
    return ret;
}

uint16_t poly_inftynorm(const poly *a, const poly *b)
{
    
    uint16_t ret = 0;

    for (int i = 0; i < NTRUPLUS_N; ++i) 
            ret = max(ret, abs(a->coeffs[i]));
    
    for (int i = 0; i < NTRUPLUS_N; ++i) 
            ret = max(ret, abs(b->coeffs[i]));

    return ret;
}

void poly_copy(poly *r, const poly *a)
{
    for (int i = 0; i < NTRUPLUS_N; ++i) 
            r->coeffs[i] = a->coeffs[i];   
}

void poly_basemul(poly *c, const poly *a, const poly *b)
{
    for(int i = 0; i < NTRUPLUS_N/8; i++)
    {
        basemul(c->coeffs + 8*i,     a->coeffs + 8*i,     b->coeffs + 8*i,      zetas[ntt_level + i]);
        basemul(c->coeffs + 8*i + 4, a->coeffs + 8*i + 4, b->coeffs + 8*i + 4, -zetas[ntt_level + i]);
    }
}

int poly_baseinv(poly *b, const poly *a)
{
    int result = 0;

    for(int i = 0; i < NTRUPLUS_N/8; ++i)
    {
        result += baseinv(b->coeffs + 8*i, a->coeffs + 8*i, zetas[ntt_level + i]);
    }

    return result;
}

#ifdef PARAM_1
void polyq_pack(uint8_t *r, const poly *a)
{
    int16_t t0, t1;

    for(int i = 0; i < NTRUPLUS_N/2; i++)
    {
        t0 = a->coeffs[2*i];
        t0 += (t0 >> 15) & NTRUPLUS_Q;
        t1 = a->coeffs[2*i+1];
        t1 += (t1 >> 15) & NTRUPLUS_Q;
        r[3*i] = t0;
        r[3*i+1] = (t0 >> 8) | (t1 << 4);
        r[3*i+2] = t1 >> 4;
  }
}

void polyq_unpack(poly *r, const uint8_t *a)
{
    int i;

    for(i=0;i<NTRUPLUS_N/2;i++){
        r->coeffs[2*i]   = a[3*i]        | ((uint16_t)a[3*i+1] & 0x0f) << 8;
        r->coeffs[2*i+1] = a[3*i+1] >> 4 | ((uint16_t)a[3*i+2] & 0xff) << 4;
    }
}
#endif

#ifdef PARAM_2
void polyq_pack(uint8_t *r, const poly *a)
{
    int i;

    for(i=0;i<NTRUPLUS_N/8;i++){
        int16_t t[8];

        for(int j = 0; j < 8 ; j++)\
        {
            t[j] = barrett_reduce(a->coeffs[8*i+j]);
            t[j] += (t[j] >> 15) & NTRUPLUS_Q;
        }

    
        r[13*i+ 0] =  t[0];
        r[13*i+ 1] = (t[0] >>  8) | ((t[1] & 0x07) << 5);
        r[13*i+ 2] = (t[1] >>  3);
        r[13*i+ 3] = (t[1] >> 11) | (t[2] << 2);
        r[13*i+ 4] = (t[2] >>  6) | (t[3] << 7);
        r[13*i+ 5] = (t[3] >>  1);
        r[13*i+ 6] = (t[3] >>  9) | (t[4] << 4);
        r[13*i+ 7] = (t[4] >>  4);
        r[13*i+ 8] = (t[4] >> 12) | (t[5] << 1);
        r[13*i+ 9] = (t[5] >>  7) | (t[6] << 6);
        r[13*i+10] = (t[6] >>  2);
        r[13*i+11] = (t[6] >> 10) | (t[7] << 3);
        r[13*i+12] = (t[7] >>  5);
    }

}

void polyq_unpack(poly *r, const uint8_t *a)
{
    int i;
    for(i=0;i<NTRUPLUS_N/8;i++)
    {
        r->coeffs[8*i+0] =  a[13*i+ 0]       | (((uint16_t)a[13*i+ 1] & 0x1f) << 8);
        r->coeffs[8*i+1] = (a[13*i+ 1] >> 5) | (((uint16_t)a[13*i+ 2]       ) << 3) | (((uint16_t)a[13*i+ 3] & 0x03) << 11);
        r->coeffs[8*i+2] = (a[13*i+ 3] >> 2) | (((uint16_t)a[13*i+ 4] & 0x7f) << 6);
        r->coeffs[8*i+3] = (a[13*i+ 4] >> 7) | (((uint16_t)a[13*i+ 5]       ) << 1) | (((uint16_t)a[13*i+ 6] & 0x0f) <<  9);
        r->coeffs[8*i+4] = (a[13*i+ 6] >> 4) | (((uint16_t)a[13*i+ 7]       ) << 4) | (((uint16_t)a[13*i+ 8] & 0x01) << 12);
        r->coeffs[8*i+5] = (a[13*i+ 8] >> 1) | (((uint16_t)a[13*i+ 9] & 0x3f) << 7);
        r->coeffs[8*i+6] = (a[13*i+ 9] >> 6) | (((uint16_t)a[13*i+10]       ) << 2) | (((uint16_t)a[13*i+11] & 0x07) << 10);
        r->coeffs[8*i+7] = (a[13*i+11] >> 3) | (((uint16_t)a[13*i+12]       ) << 5);
    }
}
#endif



void polyeta_pack(uint8_t *r, const poly *a)
{
    unsigned int i;
    uint8_t t[8];

    for (i = 0; i < NTRUPLUS_N / 4; ++i)
{
        t[0] = 1 - a->coeffs[4 * i + 0];
        t[1] = 1 - a->coeffs[4 * i + 1];
        t[2] = 1 - a->coeffs[4 * i + 2];
        t[3] = 1 - a->coeffs[4 * i + 3];
        r[i] = t[0] >> 0;
        r[i] |= t[1] << 2;
        r[i] |= t[2] << 4;
        r[i] |= t[3] << 6;
    }
}

void polyeta_unpack(poly *r, const uint8_t *a)
{
    unsigned int i;
    for (i = 0; i < NTRUPLUS_N / 4; ++i)
{
        r->coeffs[4 * i + 0] = a[i] >> 0;
        r->coeffs[4 * i + 0] &= 0x3;

        r->coeffs[4 * i + 1] = a[i] >> 2;
        r->coeffs[4 * i + 1] &= 0x3;

        r->coeffs[4 * i + 2] = a[i] >> 4;
        r->coeffs[4 * i + 2] &= 0x3;

        r->coeffs[4 * i + 3] = a[i] >> 6;
        r->coeffs[4 * i + 3] &= 0x3;

        r->coeffs[4 * i + 0] = 1 - r->coeffs[4 * i + 0];
        r->coeffs[4 * i + 1] = 1 - r->coeffs[4 * i + 1];
        r->coeffs[4 * i + 2] = 1 - r->coeffs[4 * i + 2];
        r->coeffs[4 * i + 3] = 1 - r->coeffs[4 * i + 3];
    }
}

void poly_pack_highbits(uint8_t *buf, const poly *a)
{
    unsigned int i;
    for (i = 0; i < NTRUPLUS_N / 4; i++)
{
        buf[3 * i + 0] =   a->coeffs[4 * i + 0];
        buf[3 * i + 0] |= (a->coeffs[4 * i + 1] << 6) & 0xff;

        buf[3 * i + 1] =  (a->coeffs[4 * i + 1] >> 2) & 0x0f;
        buf[3 * i + 1] |= (a->coeffs[4 * i + 2] << 4) & 0xff;

        buf[3 * i + 2] =  (a->coeffs[4 * i + 2] >> 4) & 0x03;
        buf[3 * i + 2] |= (a->coeffs[4 * i + 2] << 2) & 0xff;
    }
}

#ifdef PARAM_1
void poly_challenge(poly *c, const uint8_t highbits[POLY_HIGHBITS_PACKEDBYTES], const uint8_t mu[SEEDBYTES])
{
    unsigned int i, b, pos = 0;
    uint8_t buf[XOF256_BLOCKBYTES];
    xof256_state state;

    // H(HighBits(A * y mod 2q), LSB(round(y0) * j), M)
    xof256_absorbe_twice(&state, highbits,
                         POLY_HIGHBITS_PACKEDBYTES, mu,
                         SEEDBYTES);
    xof256_squeezeblocks(buf, 1, &state);

    for (i = 0; i < NTRUPLUS_N; ++i)
        c->coeffs[i] = 0;
    for (i = NTRUPLUS_N - NTRUPLUS_N; i < NTRUPLUS_N; ++i)
    {
        do {
            if (pos >= XOF256_BLOCKBYTES)
            {
                xof256_squeezeblocks(buf, 1, &state);
                pos = 0;
            }
            b = ((buf[pos] & 0x1)<<8)  | (buf[pos+1]); //need to check
            pos += 2;
        } while (b > i); //?????

        c->coeffs[i] = c->coeffs[b];
        c->coeffs[b] = 1;
    }

}
#endif

#ifdef PARAM_2
void poly_challenge(poly *c, const uint8_t highbits[POLY_HIGHBITS_PACKEDBYTES], const uint8_t mu[SEEDBYTES])
{
    unsigned int i, b, pos = 0;
    uint8_t buf[XOF256_BLOCKBYTES];
    xof256_state state;

    // H(HighBits(A * y mod 2q), LSB(round(y0) * j), M)
    xof256_absorbe_twice(&state, highbits,
                         POLY_HIGHBITS_PACKEDBYTES, mu,
                         SEEDBYTES);
    xof256_squeezeblocks(buf, 1, &state);

    for (i = 0; i < NTRUPLUS_N; ++i)
        c->coeffs[i] = 0;
    for (i = NTRUPLUS_N - NTRUPLUS_TAU; i < NTRUPLUS_N; ++i)
    {
        do {
            if (pos >= XOF256_BLOCKBYTES)
            {
                xof256_squeezeblocks(buf, 1, &state);
                pos = 0;
            }
            b = ((buf[pos] & 0x3)<<8)  | (buf[pos+1]); //need to check
            pos += 2;
        } while (b > i);

        c->coeffs[i] = c->coeffs[b];
        c->coeffs[b] = 1;
    }

}
#endif

#ifdef PARAM_1
void poly_decomposed_pack(uint8_t *buf, const poly *a)
{
    unsigned int i;
    for (i = 0; i < NTRUPLUS_N/8; i++)
{
        buf[7 * i + 0] = ((a->coeffs[8 * i + 1] & 0x01) << 7) | ((a->coeffs[8 * i + 0]) & 0x7f);
        buf[7 * i + 1] = ((a->coeffs[8 * i + 2] & 0x03) << 6) | ((a->coeffs[8 * i + 1] >> 1) & 0x3f);
        buf[7 * i + 2] = ((a->coeffs[8 * i + 3] & 0x07) << 5) | ((a->coeffs[8 * i + 2] >> 2) & 0x1f);
        buf[7 * i + 3] = ((a->coeffs[8 * i + 4] & 0x0f) << 4) | ((a->coeffs[8 * i + 3] >> 3) & 0x0f);
        buf[7 * i + 4] = ((a->coeffs[8 * i + 5] & 0x1f) << 3) | ((a->coeffs[8 * i + 4] >> 4) & 0x07);
        buf[7 * i + 5] = ((a->coeffs[8 * i + 6] & 0x3f) << 2) | ((a->coeffs[8 * i + 5] >> 5) & 0x03);
        buf[7 * i + 6] = ((a->coeffs[8 * i + 7] & 0x7f) << 1) | ((a->coeffs[8 * i + 6] >> 6) & 0x01);

    }
}

void poly_decomposed_unpack(poly *a, const uint8_t *buf)
{
    for(int i=0;i<NTRUPLUS_N/8;i++)
    {
        a->coeffs[8 * i + 0] =   buf[7 * i + 0] & 0x7f;  
        a->coeffs[8 * i + 1] = ((buf[7 * i + 1] & 0x3f) << 1) | ((buf[7 * i + 0] >> 7) & 0x01);
        a->coeffs[8 * i + 2] = ((buf[7 * i + 2] & 0x1f) << 2) | ((buf[7 * i + 1] >> 6) & 0x03); 
        a->coeffs[8 * i + 3] = ((buf[7 * i + 3] & 0x0f) << 3) | ((buf[7 * i + 2] >> 5) & 0x07);
        a->coeffs[8 * i + 4] = ((buf[7 * i + 4] & 0x07) << 4) | ((buf[7 * i + 3] >> 4) & 0x0f);
        a->coeffs[8 * i + 5] = ((buf[7 * i + 5] & 0x03) << 5) | ((buf[7 * i + 4] >> 3) & 0x1f);
        a->coeffs[8 * i + 6] = ((buf[7 * i + 6] & 0x01) << 6) | ((buf[7 * i + 5] >> 2) & 0x3f);
        a->coeffs[8 * i + 7] =                                   (buf[7 * i + 6] >> 1) & 0x7f;
        }
}
#endif

#ifdef PARAM_2
void poly_decomposed_pack(uint8_t *buf, const poly *a)
{
    for (int i = 0; i < NTRUPLUS_N; i++)
    {
        buf[i] = a->coeffs[i];
    }
}

void poly_decomposed_unpack(poly *a, const uint8_t *buf)
{
    for(int i = 0; i < NTRUPLUS_N; i++)
    {
        a->coeffs[i] = buf[i];
    }
}
#endif

void poly_compose(poly *a, const poly *ha, const poly *la)
{
    for (int i = 0; i < NTRUPLUS_N; ++i)
    {
        a->coeffs[i] = (ha->coeffs[i] << NTRUPLUS_d) + la->coeffs[i];
        a->coeffs[i] -= ((~(la->coeffs[i] - HALF_ALPHA_HINT)) & ALPHA_HINT);
    }
}

static inline void minmax(int64_t *x, int64_t *y) // taken from djbsort
{
    int64_t a = *x;
    int64_t b = *y;
    int64_t ab = b ^ a;
    int64_t c = b - a;
    c ^= ab & (c ^ b);
    c >>= 63;
    c &= ab;
    *x = a ^ c;
    *y = b ^ c;
}

int64_t poly_sqsing_value(const poly *s1, const poly *s2)
{
    __int128_t res = 0;
    complex_fp64_32 input[NTRUPLUS_N] = {0};

    int64_t sum[NTRUPLUS_N] = {0};
    int64_t bestm[NTRUPLUS_M + 1] = {0};

    fft_init(input, s1);
    fft(input);
	
    for(int i = 0; i < NTRUPLUS_N; i++)
    {
        sum[i] = complex_fp_sqabs(input[i]);
    }

    fft_init(input, s2);
    fft(input);

    for(int i = 0; i < NTRUPLUS_N; i++)
    {
        sum[i] += complex_fp_sqabs(input[i]);
    }

	// compute max m
    for(int i = 0; i < NTRUPLUS_M + 1; ++i)
    {
        bestm[i] = sum[i];
    }

    for(int i = NTRUPLUS_M + 1; i < NTRUPLUS_N; i++)
    {
        for(int j = 0; j < NTRUPLUS_M + 1; j++)
        {
            minmax(&sum[i], &bestm[j]);
        }
    }

    for(int i = NTRUPLUS_M; i > 0; i--)
    {
        minmax(&bestm[i-1], &bestm[i]);
    }

    for(int i = 1; i < NTRUPLUS_M + 1; i++)
    {
		bestm[i] = (bestm[i] + (1 << 15)) >> 16;	
        res += bestm[i];
    }

	res = NTRUPLUS_TAU * res + NTRUPLUS_R * ((bestm[0] + (1 << 15)) >> 16) + (1 << 15);
	res = res >> 16;

    return res;
}
//#endif