/*#include "fixpoint.h"*/
#include "params.h"
#include <stdint.h>
#include "poly.h"
#include <string.h>
#include "fips202.h"

void hash_H(uint8_t *c_bin, const uint8_t highbits[POLY_HIGHBITS_PACKEDBYTES], const uint8_t mu[SEEDBYTES]);
void encode_c(poly *c, const uint8_t *c_bin);