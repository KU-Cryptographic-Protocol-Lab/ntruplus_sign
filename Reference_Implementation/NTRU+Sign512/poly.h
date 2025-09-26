#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include "params.h"

typedef struct{
  int16_t coeffs[NTRUPLUS_N];
} poly;

void poly_cbd1(poly *r, const unsigned char seed[CRHBYTES], uint16_t nonce);
void poly_sample_y2(poly *a1, poly *a2);

uint64_t innerproduct(poly *a1, poly *a2);

void poly_ntt(poly *r);
void poly_invntt(poly *r);

void poly_uniform_eta(poly *a, const uint8_t seed[CRHBYTES], uint16_t nonce);

void poly_add(poly *r, const poly *a, const poly *b);
void poly_add_qhat(poly *r, const poly *a);
void poly_caddq(poly *v);
void poly_csubq(poly *v);
void poly_caddp(poly *v);
void poly_csubp(poly *v);
void poly_freezep(poly *v);
uint16_t poly_eqcheck(const poly *a, const poly *b);
uint16_t c_eqcheck(const uint8_t c[CRYPTO_C_BYTES], const uint8_t c2[CRYPTO_C_BYTES]);
void poly_cqhat(poly *v);
void poly_copy(poly *r, const poly *a);
void poly_sub(poly *r, const poly *a, const poly *b);
void poly_cneg(poly *r, const uint8_t b);
void poly_neg(poly *r);
void poly_barrett_reduce(poly *v);
void poly_shift(poly *r, const poly *a);
void poly_bc(poly *v, const poly *c, uint8_t b);
uint64_t poly_sqnorm2(const poly *a, const poly *b);
uint16_t poly_inftynorm(const poly *a, const poly *b);

void poly_basemul(poly *r, const poly *a, const poly *b);
void poly_basemul2(poly *c, const poly *a, const poly *b);
int poly_baseinv(poly *b, const poly *a);

void poly_pack_highbits(uint8_t *buf, const poly *a);
void polyq_pack(uint8_t *r, const poly *a);
void polyq_unpack(poly *r, const uint8_t *a);
void polyeta_pack(uint8_t *r, const poly *a);
void polyeta_unpack(poly *r, const uint8_t *a);

void poly_challenge(poly *c, const uint8_t highbits[POLY_HIGHBITS_PACKEDBYTES], const uint8_t mu[SEEDBYTES]);

void poly_h_pack(uint8_t *buf, const poly *a);
void poly_h_unpack(poly *r, const uint8_t *a);
void poly_z1_pack(uint8_t *buf, const poly *a);
void poly_z1_unpack(poly *r, const uint8_t *a);

void poly_decomposed_pack(uint8_t *buf, const poly *a);
void poly_decomposed_unpack(poly *a, const uint8_t *buf);
void poly_compose(poly *a, const poly *ha, const poly *la);

int64_t poly_sqsing_value(const poly *s1, const poly *s2);
#endif