#ifndef HAETAE_ROUNDING_H
#define HAETAE_ROUNDING_H

#include "params.h"

#include "poly.h"
#include <stdint.h>

void poly_highbits_hint(poly *w, const poly *v);
void decompose_hint(int16_t *highbits, const int16_t r);

void poly_high_row_bits(poly *w, poly *w2, const poly *v);
void decompose_high_row(int16_t *highbits, int16_t *lowbits, const int16_t r);

void poly_highbits_hint_test2(poly *w, poly *w2, const poly *v);
void decompose_hint_test2(int16_t *highbits, int16_t *lowbits, const int16_t r);

void poly_highbits_hint_test3(poly *w, poly *w2, const poly *v);
void decompose_hint_test3(int16_t *highbits, int16_t *lowbits, const int16_t r);
#endif