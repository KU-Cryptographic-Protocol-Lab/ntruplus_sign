// clang-format off

#include "params.h"
#include "api.h"
#include <stdint.h>
#include "poly.h"

void pack_pk(uint8_t pk[CRYPTO_PUBLICKEYBYTES], const poly *a);
void unpack_pk(poly *a, const uint8_t pk[CRYPTO_PUBLICKEYBYTES]);
void pack_sk(uint8_t sk[CRYPTO_SECRETKEYBYTES],
             const uint8_t pk[CRYPTO_PUBLICKEYBYTES], const poly *s0,
             const poly *s1, const uint8_t key[SEEDBYTES]);
void unpack_sk(poly *a, poly *s0, poly *s1, uint8_t *key,
               const uint8_t sk[CRYPTO_SECRETKEYBYTES]);
/*void pack_sig(uint8_t sig[CRYPTO_BYTES], const poly *c, const poly *z1, const poly *h);
void unpack_sig(poly *c, poly *z1, poly *h, const uint8_t sig[CRYPTO_BYTES]);*/

int pack_sig(uint8_t sig[CRYPTO_BYTES], const uint8_t c[CRYPTO_C_BYTES], const poly *z1_high, const poly *z1_low, const poly *h);
int unpack_sig(uint8_t c[CRYPTO_C_BYTES], poly *z1_high, poly *z1_low, poly *h, const uint8_t sig[CRYPTO_BYTES]);