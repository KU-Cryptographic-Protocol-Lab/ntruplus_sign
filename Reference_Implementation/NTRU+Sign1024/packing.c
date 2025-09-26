#include "packing.h"
#include "params.h"
#include "poly.h"
#include <string.h>
#include "encoding.h"

#include <stdint.h>
#include <stdlib.h>
void pack_pk(uint8_t pk[CRYPTO_PUBLICKEYBYTES], const poly *a) {

    polyq_pack(pk, a);
    
}

void unpack_pk(poly *a, const uint8_t pk[CRYPTO_PUBLICKEYBYTES]) {

    polyq_unpack(a, pk);
}

void pack_sk(uint8_t sk[CRYPTO_SECRETKEYBYTES],
             const uint8_t pk[CRYPTO_PUBLICKEYBYTES], const poly *s0,
             const poly *s1, const uint8_t key[SEEDBYTES]) {

    memcpy(sk, pk, CRYPTO_PUBLICKEYBYTES);

    sk += CRYPTO_PUBLICKEYBYTES;
    polyq_pack(sk, s0);

    sk += CRYPTO_PUBLICKEYBYTES;

    polyq_pack(sk, s1);
    sk += CRYPTO_PUBLICKEYBYTES;
    memcpy(sk, key, SEEDBYTES);
}

void unpack_sk(poly *a, poly *s0, poly *s1, uint8_t *key,
               const uint8_t sk[CRYPTO_SECRETKEYBYTES]) {

    polyq_unpack(a, sk);
    sk += CRYPTO_PUBLICKEYBYTES;

    polyq_unpack(s0, sk);
    sk += CRYPTO_PUBLICKEYBYTES;

    polyq_unpack(s1, sk);
    sk += CRYPTO_PUBLICKEYBYTES;

    memcpy(key, sk, SEEDBYTES);
}

int pack_sig(uint8_t sig[CRYPTO_BYTES], const uint8_t c[CRYPTO_C_BYTES], const poly *z1_high, const poly *z1_low, const poly *h) {

    uint16_t size_enc_hb_z1, size_enc_h;
    uint8_t encoded_hb_z1[NTRUPLUS_N], encoded_h[NTRUPLUS_N];
    uint8_t offset_enc_h, offset_enc_hb_z1;
        
    // init/padding with zeros:
    memset(sig, 0, CRYPTO_BYTES);
    
    // encode challenge
    for (size_t i = 0; i < CRYPTO_C_BYTES; i++)
        sig[i] = c[i];
    
    sig += CRYPTO_C_BYTES;

    poly_decomposed_pack(sig, z1_low);
    sig += NTRUPLUS_d * NTRUPLUS_N / 8;
    
    size_enc_hb_z1 = encode_hb_z1(encoded_hb_z1, z1_high);
    size_enc_h = encode_h(encoded_h, h);
    
    //printf("sum = %d, hb_z1 = %d, h = %d \n", size_enc_h + size_enc_hb_z1, size_enc_hb_z1, size_enc_h);s

    if(size_enc_hb_z1 + size_enc_h > ENC_bound){
        return 1;
    }

    if (size_enc_h == 0 || size_enc_hb_z1 == 0) {
        return 1; // encoding failed
    }

    // The size of the encoded h and HB(z1) does not always fit in one byte,
    // thus we output a one byte offset to a fixed baseline
    if (size_enc_h < BASE_ENC_H || size_enc_hb_z1 < BASE_ENC_HB_Z1 ||
        size_enc_h > BASE_ENC_H + 255 ||
        size_enc_hb_z1 > BASE_ENC_HB_Z1 + 255) {
        return 1; // encoding size offset out of range
    }

    offset_enc_hb_z1 = size_enc_hb_z1 - BASE_ENC_HB_Z1;
    offset_enc_h = size_enc_h - BASE_ENC_H;

    if (NTRUPLUS_d * NTRUPLUS_N/8 + 2 + size_enc_hb_z1 + size_enc_h + CRYPTO_C_BYTES > CRYPTO_BYTES) {
        return 1; // signature too big
    }

    sig[0] = offset_enc_hb_z1;
    sig[1] = offset_enc_h;
    sig += 2;

    memcpy(sig, encoded_hb_z1, size_enc_hb_z1);
    sig += size_enc_hb_z1;

    memcpy(sig, encoded_h, size_enc_h);
    sig += size_enc_h;

    return 0;
}

int unpack_sig(uint8_t c[CRYPTO_C_BYTES], poly *z1_high, poly *z1_low, poly *h, const uint8_t sig[CRYPTO_BYTES]) {

    unsigned int i;
    uint16_t size_enc_hb_z1, size_enc_h;
    
    for (i = 0; i < CRYPTO_C_BYTES; i++)
    {
      c[i] = sig[i];
    }
    sig += CRYPTO_C_BYTES;

    poly_decomposed_unpack(z1_low, sig);
    sig += NTRUPLUS_d * NTRUPLUS_N / 8;
    
    size_enc_hb_z1 = (uint16_t)sig[0] + BASE_ENC_HB_Z1;
    size_enc_h = (uint16_t)sig[1] + BASE_ENC_H;
    sig += 2;

    if (CRYPTO_BYTES < (NTRUPLUS_d * NTRUPLUS_N/8 + 2 + CRYPTO_C_BYTES + size_enc_hb_z1 + size_enc_h))
        return 1; // invalid size_enc_h and/or size_enc_hb_z1


    if (decode_hb_z1(z1_high, sig, size_enc_hb_z1)) {
        return 1; // decoding failed
    }

    sig += size_enc_hb_z1;

    if (decode_h(h, sig, size_enc_h)) {
        return 1; // decoding failed
    }

    sig += size_enc_h;

    return 0;
}