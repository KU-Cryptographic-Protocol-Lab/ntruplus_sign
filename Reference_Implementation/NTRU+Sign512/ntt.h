#ifndef NTT_H
#define NTT_H

#include <stdint.h>
#include "params.h"
#ifdef PARAM_1
extern const int16_t zetas[128];
#endif

#ifdef PARAM_2
extern const int16_t zetas[256];
#endif
void ntt(int16_t r[NTRUPLUS_N]);
void invntt(int16_t r[NTRUPLUS_N]);

#ifdef PARAM_1
void basemul(int16_t r[4], const int16_t a[4], const int16_t b[4], int16_t zeta);
int  baseinv(int16_t r[8], const int16_t a[8], int16_t zeta);
#endif

#ifdef PARAM_2
void basemul(int16_t r[4], const int16_t a[4], const int16_t b[4], int16_t zeta);
int  baseinv(int16_t r[8], const int16_t a[8], int16_t zeta);
#endif

#endif