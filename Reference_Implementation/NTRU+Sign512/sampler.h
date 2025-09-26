/*#include "fixpoint.h"*/
#include "params.h"
#include "poly.h"
#include <stdint.h>

void hash_H(uint8_t *c_bin, const uint8_t highbits[POLY_HIGHBITS_PACKEDBYTES], const uint8_t mu[SEEDBYTES]);
void encode_c(poly *c, const uint8_t *c_bin);