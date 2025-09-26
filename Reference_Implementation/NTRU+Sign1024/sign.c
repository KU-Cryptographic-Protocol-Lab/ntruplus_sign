#include "api.h"

#include "sign.h"
#include "packing.h"
#include "params.h"
#include "poly.h"
#include "randombytes.h"
#include "encoding.h"
#include "decompose.h"
#include "symmetric.h"
#include "sampler.h"
#include "gaussian.h"
#include "info.h"
#include "fft.h"
#include <inttypes.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>

unsigned long long ctr_sign;

int crypto_sign_keypair(unsigned char *pk, unsigned char *sk)
{
    int64_t r;
    poly f, tmp;
    poly g, ginv;
    poly a;
    int64_t squared_singular_value; 
    uint8_t seedbuf[2*CRHBYTES] = {0};
    uint16_t nonce = 0;
    const uint8_t *sigma, *key;
    
    do {
        r = 0;
        randombytes(seedbuf, 2*CRHBYTES);
    
        sigma = seedbuf;
        key = sigma + CRHBYTES;
    
        // Sample secret vectors f and g
        poly_cbd1(&g, seedbuf, nonce++);
        poly_cbd1(&f, seedbuf, nonce++);
     
        squared_singular_value = poly_sqsing_value(&g, &f);
        r = ~(squared_singular_value - ((int64_t) NTRUPLUS_GAMMA2N)) >> 63;
        //if(squared_singular_value > maxg) maxg = squared_singular_value;

        poly_ntt(&g);
        r |= poly_baseinv(&ginv, &g);
    } while(r);

    //compute a = (f+qhat)/g in NTT form
    poly_ntt(&f);
    poly_add_qhat(&tmp,&f);
    poly_basemul(&a,&ginv,&tmp);

    pack_pk(pk, &a);
    pack_sk(sk, pk, &g, &f, key);

    return 0;
}

int crypto_sign(unsigned char *sm, unsigned long long *smlen,
            const unsigned char *m, unsigned long long mlen,
            const unsigned char *sk)
{
    ctr_sign=0;
    uint8_t buf[POLY_HIGHBITS_PACKEDBYTES] = {0};
    uint8_t seedbuf[CRHBYTES+16] = {0}, key[SEEDBYTES] = {0};
    uint8_t mu[SEEDBYTES] = {0};

    uint8_t b, nonce =0;                  
    
    uint8_t reject;
    
    poly a, ay1, u;
    poly s1, cs1;
    poly s2, cs2;
    poly highbits, highbits2;
    poly y1, y2;
    poly z1, z2;
    poly h, hshift, htmp; 
    poly c, chat, ccheck;
    poly bc;
    xof256_state state;

    int64_t val;
    int64_t res;
    int b1, b2, b3;

    uint8_t randbytes[17] = {0};
    uint8_t tmp[CRHBYTES + 2];
    
    // Unpack secret key
    unpack_sk(&a, &s1, &s2, key, sk);

    xof256_absorbe_twice(&state, sk, CRYPTO_PUBLICKEYBYTES, m, mlen);
    xof256_squeeze(mu, SEEDBYTES, &state);
    xof256_absorbe_twice(&state, key, SEEDBYTES, mu, SEEDBYTES);
    xof256_squeeze(seedbuf, CRHBYTES+16, &state);
    
loop_reject:
           
    ctr_sign++;

    //poly y3, y4;
    sample_y(&y1, &y2, seedbuf, nonce++);

    for (int i = 0; i < CRHBYTES; i++)
        tmp[i] = seedbuf[i];

    tmp[CRHBYTES + 0] = nonce >> 0;     
    tmp[CRHBYTES + 1] = nonce >> 8;
    shake256(randbytes, 17, tmp, CRHBYTES+2);
    b = randbytes[0] & 0x1;

    //compute u = ay1+y2
    poly_ntt(&y1);  
    poly_basemul(&ay1, &a, &y1);
    poly_invntt(&y1);
    poly_invntt(&ay1);

    //compute [u]_d 
    poly_add(&u, &ay1, &y2);    
    poly_caddq(&u);
    poly_highbits_hint(&highbits, &u);

    //get c = H([u]_d, mu)
    poly_pack_highbits(buf, &highbits);

    uint8_t c_test[CRYPTO_C_BYTES] = {0};
    hash_H(c_test, buf, mu);
    encode_c(&c, c_test);

    for (int i = 0; i < NTRUPLUS_N; i++)
        chat.coeffs[i] = (-c.coeffs[i]) & Rmodq;
    
    ccheck = c;
    poly_ntt(&chat);
    poly_basemul(&cs1, &chat, &s1);
    poly_invntt(&cs1);
    poly_basemul(&cs2, &s2, &chat);
    poly_invntt(&cs2);

    //compute z = y + (-1)^b cs
    poly_cneg(&cs1, b);
    poly_cneg(&cs2, 1-b); //if s2 = -f, then the way to compute cs2 is same with the way getting cs1, but s2 is set as f for convenience.
    poly_add(&z1, &y1, &cs1);
    poly_add(&z2, &y2, &cs2);


    val = B_scsquare - poly_sqnorm2(&cs1, &cs2);

    b1 = SampleBernExp(val, &randbytes[1]);
    res = innerproduct(&z1, &cs1) + innerproduct(&z2, &cs2);
    res <<= 1;
    b2 = SampleBernCosh(res, &randbytes[9]);
    b3 = b1*b2;
    
    if (1-b3) { 
        goto loop_reject;
    }

    //compute htmp = u-z2+(1-b)c
    poly_sub(&htmp, &u, &z2);
    poly_bc(&bc, &c, b);
    poly_add(&htmp, &htmp, &bc);
    poly_caddq(&htmp);
    poly_csubq(&htmp);
    
    //compute h = [u]_d - [u-z2+(1-b)c]_d mod+- p
    poly_highbits_hint(&htmp, &htmp);
    poly_sub(&h, &highbits, &htmp); 
    poly_freezep(&h);

    //hcheck = h*2^d
    poly_shift(&hshift, &h);
    
    //Bound check B2, Binfty
    reject = (((uint64_t)B2bound - poly_sqnorm2(&z1, &hshift)) >> 63) & 1;
    reject |= ((Binftybound - poly_inftynorm(&z1, &hshift)) >> 15)& 1;

       
    if (reject) { //required to add rejection sampling 
        goto loop_reject;
    }
    //equality check: [u]_d= [[u+(-1)^b c]_d]
    poly_cneg(&ccheck, b);
    poly_add(&highbits2, &u, &ccheck); 
    poly_highbits_hint(&highbits2, &highbits2);
    reject |= poly_eqcheck(&highbits, &highbits2);
    
    if (reject) { //required to add rejection sampling 
        goto loop_reject;
    }

    poly z1_high, z1_low;
    poly_high_row_bits(&z1_high, &z1_low, &z1);
            
    if (pack_sig(sm, c_test, &z1_high, &z1_low, &h)){
        goto loop_reject;
    }

    *smlen = CRYPTO_BYTES;
    return 0;
}

int crypto_sign_verify(unsigned char *m, unsigned long long mlen,
                 const unsigned char *sm, unsigned long long smlen,
                 const unsigned char *pk)
{
  
    poly a, c, z1, z1_high, z1_low, h, az1, chat, hshift;
    uint8_t buf[POLY_HIGHBITS_PACKEDBYTES] = {0};
    uint8_t c_test[CRYPTO_C_BYTES];
    uint8_t mu[SEEDBYTES];
    xof256_state state;
    uint8_t cprime[CRYPTO_C_BYTES] = {0};

    // Check signature length
    if (smlen != CRYPTO_BYTES) {
        return -1;
    }

    //unpack_sig_test(&c, &z1_high, &z1_low, &h, sm);
    if (unpack_sig(c_test, &z1_high, &z1_low, &h, sm)) {
        return -1;
    }

    encode_c(&c, c_test);
    poly_compose(&z1, &z1_high, &z1_low);

    //unpack_sig(&c, &z1, &h, sm);
    unpack_pk(&a, pk);
    
    //Bound check B2, Binfty
    poly_shift(&hshift, &h);
    if (((((uint64_t)B2bound - poly_sqnorm2(&z1, &hshift)) >> 63) & 1) | (((Binftybound - poly_inftynorm(&z1, &hshift)) >> 15)& 1)) {
        return -2;
    }
    
    //compute az1 
    poly_ntt(&z1);
    poly_basemul(&az1, &a, &z1);
    poly_invntt(&az1);
 
    //compute [az1+cqhat]_d
    chat = c;
    poly_cqhat(&chat);
    poly_add(&az1, &az1, &chat);
    poly_caddq(&az1);
    //poly_csubq(&az1);
    poly_highbits_hint(&az1, &az1);

    //compute h+[az1+cqhat]_d 
    poly_add(&az1, &az1, &h);
    poly_caddp(&az1);
    poly_csubp(&az1);
 
    //check if c = H(h+[az1+cqhat]_d, mu) 
    poly_pack_highbits(buf, &az1);
    xof256_absorbe_twice(&state, pk, CRYPTO_PUBLICKEYBYTES, m, mlen);
    xof256_squeeze(mu, SEEDBYTES, &state);

    hash_H(cprime, buf, mu);    
    
    if(c_eqcheck(c_test, cprime)){
        return -3;
    }
    
    return 0;
    
}
