#include "rejection.h"

static inline uint64_t myMul64Shift(uint64_t a, uint64_t b, int shift)
{
	__int128 a128 = (__int128)a;
	__int128 b128 = (__int128)b;
	__int128 c128;

	c128 = (a128 * b128 + (1 << (shift -1))) >> shift;

	return (uint64_t)c128;
}

static inline uint64_t myMul64Shift_NoRounding(uint64_t a, uint64_t b, int shift)
{
	__int128 a128 = (__int128)a;
	__int128 b128 = (__int128)b;
	__int128 c128;

	c128 = (a128 * b128) >> shift;

	return (uint64_t)c128;
}

uint64_t PExpG(uint64_t xin)
{
	uint64_t r;

	r = myMul64Shift(Coef_G[0], xin, expDiff_G[0]);
	for (int i = 1; i < DegG; i++) {
		r = myMul64Shift(Coef_G[i] - r, xin, expDiff_G[i]);
	}
	r = (Coef_G[DegG] - r);

	return r;
}

uint64_t PExpR(uint64_t xin)
{
	int64_t r;
	r = myMul64Shift(Coef_R[0], xin, expDiff_R[0]);
	for (int i = 1; i < DegR; i++) {
		r = myMul64Shift(Coef_R[i] - r, xin, expDiff_R[i]);
	}
	r = (Coef_R[DegR] - r);

	return r;
}

static void ConDiv(uint64_t x, uint64_t c, uint64_t* k, uint64_t* y, int diffE)
{
	uint64_t b;

	(*k) = 0;
	(*y) = x >> diffE;			

	for (int i = diffE - 1; i >= 0; i--) {
		(*y) <<= 1;
		(*y) |= (x >> i) & 1;
		b = LSBMASK((*y) >= c);
		(*k) |= (1ULL << i) & b;
		(*y) -= c & b;
	}
}


int SampleBernExpSimple(uint64_t x, const unsigned char *buf)
{
    uint64_t val, ran;

    val = PExpG(x);

    // 60비트 난수 생성 (PARAM_1 기준)

	#ifdef PARAM_1
	ran = ((uint64_t)(buf[0] & 0x07) << 56) | 
          ((uint64_t)buf[1] << 48) |
          ((uint64_t)buf[2] << 40) |
          ((uint64_t)buf[3] << 32) |
          ((uint64_t)buf[4] << 24) |
          ((uint64_t)buf[5] << 16) |
          ((uint64_t)buf[6] << 8)  | 
          (uint64_t)buf[7];
	#endif

    #ifdef PARAM_2
    ran = ((uint64_t)(buf[0] & 0x0f) << 56) | 
          ((uint64_t)buf[1] << 48) |
          ((uint64_t)buf[2] << 40) |
          ((uint64_t)buf[3] << 32) |
          ((uint64_t)buf[4] << 24) |
          ((uint64_t)buf[5] << 16) |
          ((uint64_t)buf[6] << 8)  | 
          (uint64_t)buf[7];  
    #endif

    return (ran <= val) ? 1 : 0;
}

uint64_t exp_extended(int64_t x, int* exponent)
{
	uint64_t absx, q, r, temp;
	int signx = (x >= 0);

	absx = CMUX(x, -x, (int64_t)signx) << ThetaChat;

	ConDiv(absx, ValChat, &q, &r, BitGamma);

	*exponent = CMUX(vartheta2 + q, vartheta2 - (q + 1), (int64_t)signx);
	r = CMUX(r, ValChat - r, (int64_t)signx);

	r = PExpR(r);
	temp = CMUX(q * (uint64_t)EpsVal, (q+1) * (uint64_t)EpsVal, (int64_t)signx);	
	temp = myMul64Shift(r, temp, ThetaEps);											
	r = CMUX(r - temp, r + temp, (int64_t)signx);

	return r;
}

int SampleBernExp(uint64_t x, const unsigned char *buf)
{
	uint64_t val, ran;
	int e;

	val = exp_extended(x, &e);
	val >>= (e - vartheta2);

#ifdef PARAM_1
	ran = ((uint64_t)(buf[0] & 0x0f) << 56) | 
          ((uint64_t)buf[1] << 48) |
          ((uint64_t)buf[2] << 40) |
          ((uint64_t)buf[3] << 32) |
          ((uint64_t)buf[4] << 24) |
          ((uint64_t)buf[5] << 16) |
          ((uint64_t)buf[6] << 8)  | 
          (uint64_t)buf[7]; 
#endif

#ifdef PARAM_2
	ran = ((uint64_t)(buf[0] & 0x07) << 56) | 
          ((uint64_t)buf[1] << 48) |
          ((uint64_t)buf[2] << 40) |
          ((uint64_t)buf[3] << 32) |
          ((uint64_t)buf[4] << 24) |
          ((uint64_t)buf[5] << 16) |
          ((uint64_t)buf[6] << 8)  | 
          (uint64_t)buf[7];
#endif

	return (ran <= val);
}

int SampleBernCosh(uint64_t x, const unsigned char *buf)
{
	uint64_t val1, val2;
	uint64_t val, ran;
	int e1, e2;

	uint64_t absx;
	int signx = (x >= 0);
	absx = CMUX(x, -1 * x, (int64_t)signx);

	val1 = exp_extended(-1 * absx, &e1);	// exp(+)
	val2 = exp_extended(     absx, &e2);	// exp(-)

	e2 -= e1;
	val2 = CMUX(val2 >> e2, 0, e2 <= vartheta2);
	val1 += val2;
	e1 += 1;

#ifdef PARAM_1
	ran = ((uint64_t)(buf[0] & 0x0f) << 56) | 
          ((uint64_t)buf[1] << 48) |
          ((uint64_t)buf[2] << 40) |
          ((uint64_t)buf[3] << 32) |
          ((uint64_t)buf[4] << 24) |
          ((uint64_t)buf[5] << 16) |
          ((uint64_t)buf[6] << 8)  | 
          (uint64_t)buf[7]; 
#endif

#ifdef PARAM_2
	ran = ((uint64_t)(buf[0] & 0x07) << 56) | 
          ((uint64_t)buf[1] << 48) |
          ((uint64_t)buf[2] << 40) |
          ((uint64_t)buf[3] << 32) |
          ((uint64_t)buf[4] << 24) |
          ((uint64_t)buf[5] << 16) |
          ((uint64_t)buf[6] << 8)  | 
          (uint64_t)buf[7];
#endif

	val = myMul64Shift_NoRounding(ran, val1, vartheta2);

	return (val < (1ULL << e1)) ? 1 : 0;
}