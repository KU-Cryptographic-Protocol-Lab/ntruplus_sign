#include "decompose.h"
#include "params.h"
#include <stdint.h>

/*************************************************
 * Name:        decompose_hint
 *
 * Description: For finite field element r, compute highbits
 *              hb, lb such that r = hb * b + lb with -b/4 < lb <= b/4.
 *
 * Arguments:   - int32_t r: input element
 *              - int32_t *highbits: pointer to output element hb
 **************************************************/


void poly_highbits_hint(poly *w, const poly *v) {
    for (int i = 0; i < NTRUPLUS_N; i++) {
            decompose_hint(&w->coeffs[i], v->coeffs[i]);
    }
}
#if defined(PARAM_1) || defined(PARAM_2)
void decompose_hint(int16_t *highbits, const int16_t r) {
    int16_t hb, edgecase;

    hb = (r + HALF_ALPHA_HINT) >> NTRUPLUS_d;
    edgecase =
        ((NTRUPLUS_Q - 1) / ALPHA_HINT - (hb + 1)) >> 15; // if hb == (DQ-2)/ALPHA
    hb -= (NTRUPLUS_Q - 1) / ALPHA_HINT & edgecase;       // hb = 0

    *highbits = hb;
}
#endif

void poly_high_row_bits(poly *w, poly *w2, const poly *v) {
    for (int i = 0; i < NTRUPLUS_N; i++) {
        decompose_high_row(&w->coeffs[i], &w2->coeffs[i], v->coeffs[i]);
    }
}

#if defined(PARAM_1) || defined(PARAM_2)
void decompose_high_row(int16_t *highbits, int16_t *lowbits, const int16_t r) {
    int16_t hb, edgecase, lb;
    int16_t modALPHA_HINT = ALPHA_HINT - 1;

    hb = (r + HALF_ALPHA_HINT) >> NTRUPLUS_d;
    edgecase =
        ((NTRUPLUS_Q - 1) / ALPHA_HINT - (hb + 1)) >> 15; // if hb == (DQ-2)/ALPHA
    hb -= (NTRUPLUS_Q - 1) / ALPHA_HINT & edgecase;       // hb = 0
    lb = r & modALPHA_HINT;

    *highbits = hb;
    *lowbits = lb;
}
#endif
