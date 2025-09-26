#ifndef FFT__H
#define FFT__H

#include <stdint.h>
#include "poly.h"
#include "params.h"

typedef struct {
    int64_t real;
    int64_t imag;
} complex_fp64_32;

void fft(complex_fp64_32 data[NTRUPLUS_N]);
void fft_init(complex_fp64_32 r[NTRUPLUS_N], const poly *a);
int64_t complex_fp_sqabs(complex_fp64_32 x);

#endif
