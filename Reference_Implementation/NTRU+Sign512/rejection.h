#ifndef REJ_H
#define REJ_H

#include "info.h"
#include "params.h"
//#include "poly.h"

#define LSBMASK(c)	(-((c)&1))
#define	CMUX(x,y,c)	(((x)&(LSBMASK(c)))^((y)&(~LSBMASK(c))))

int SampleBernExpSimple(uint64_t x, const unsigned char *buf);

int SampleBernExp(uint64_t x, const unsigned char *buf);
int SampleBernCosh(uint64_t x, const unsigned char *buf);

uint64_t PExpG(uint64_t xin);
uint64_t PExpR(uint64_t xin);
uint64_t exp_extended(int64_t x, int* exponent);

#endif
