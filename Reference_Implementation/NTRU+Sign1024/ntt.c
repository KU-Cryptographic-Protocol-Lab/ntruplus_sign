#include "params.h"
#include "ntt.h"

#ifdef PARAM_1
#define QINV 62209 // q^(-1) mod 2^16
#endif

#ifdef PARAM_2
#define QINV 57857 // q^(-1) mod 2^16
#endif

#ifdef PARAM_1
const int16_t zetas[128] = {
	2285, 2571, 2970, 1812, 1493, 1422,  287,  202, 
	3158,  622, 1577,  182,  962, 2127, 1855, 1468, 
	 573, 2004,  264,  383, 2500, 1458, 1727, 3199, 
	2648, 1017,  732,  608, 1787,  411, 3124, 1758,
	1223,  652, 2777, 1015, 2036, 1491, 3047, 1785,
	 516, 3321, 3009, 2663, 1711, 2167,  126, 1469,
	2476, 3239, 3058,  830,  107, 1908, 3082, 2378,
	2931,  961, 1821, 2604,  448, 2264,  677, 2054,
	2226,  430,  555,  843, 2078,  871, 1550,  105,
	 422,  587,  177, 3094, 3038, 2869, 1574, 1653,
	3083,  778, 1159, 3182, 2552, 1483, 2727, 1119,
	1739,  644, 2457,  349,  418,  329, 3173, 3254,
	 817, 1097,  603,  610, 1322, 2044, 1864,  384,
	2114, 3193, 1218, 1994, 2455,  220, 2142, 1670,
	2144, 1799, 2051,  794, 1819, 2475, 2459,  478,
	3221, 3021,  996,  991,  958, 1869, 1522, 1628};
#endif

#ifdef PARAM_2
const int16_t zetas[256] = {
	-3593,  3777, -3182,  3625, -3696, -1100,  2456,  2194,
	  121, -2250,   834, -2495, -2319,  2876, -1701,  1414,
	 2816, -2088, -2237,  1986, -1599,  1993,  3706, -2006,
	-1525, -2557,  1296,  1483, -2830,  3364,   617,  1921,
	-3689, -1738,  3266, -3600,   810,  1887,  -638,    -7,
	 -438,  -679, -1305, -1760,   396, -3174, -3555, -1881,
	 3772, -2535, -2440, -2555,  1535,  -549,  3153,  2310,
	-1399,  1321,   514, -2956,  -103,  2804, -2043, -1431,
	-1054,  1698, -3456,  1166,  2426,  3831,   915,    -2,
	-3417,  -194,  2919,  2789,  3405,  2385, -2113, -2732,
	 2175,   373,  3692,  -730, -1756,  3135, -2391,   660,
	-1497,  2572, -3145,  1350, -2224, -3588, -1681,  2883,
	-1390,  1598,  3750,  2762,  2835,  2764, -2233,  3816,
	-1533,  1464,  -727,  1521,  1386, -3428,  -921, -2743,
	-2160,  2649,  -859,  2579,  1532,  1919,  -486,   404,
	-1056,   783,  1799, -2665,  3480,  2133, -3310, -1168,
	  -17,  3744,  2422,  2001,  1278,   929, -1348, -2230,
	 -179, -1242, -2059, -1070,  2161,  1649,  2072,  3177,
	-2071,  1121,  -436,   236,   715,   670,  -658, -1476,
	-2378,  2767,  3542,  -226,  1203,  1181,  -151, -3794,
	 1712,  -222,  2786,  -451, -3547,  1779, -1151,  -434,
	 3568, -3693,  3581, -1586,  1509,  2918,  2339, -1407,
	 3434, -3550,  2340,  2891,  2998, -3314,  3461, -2719,
	-2247, -2589,  1144,  1072,  1295, -2815, -3770,  3450,
	 3781, -2258,   796,  3163, -3208,  -589,  2963,  -124,
	 3214,  3334, -3366, -3745,  3723,  1931,  -429,  -402,
	-3408,    83, -1526,   826, -1338,  2345, -2303,  2515,
	 -642, -1837, -2965,  -791,   370,   293,  3312,  2083,
	-1689,  -777,  2070,  2262,  -893,  2386,  -188, -1519,
	-2874, -1404,  1012,  2130,  1441,  2532, -3335, -1084,
	-3343,  2937,   509, -1403,  2812,  3763,   592,  2005,
	 3657,  2460, -3677,  3752,   692,  1669,  2167, -3287
};
#endif


static inline int16_t montgomery_reduce(int32_t a)
{
  int16_t t;

  t = (int16_t)a*QINV;
  t = (a - (int32_t)t*NTRUPLUS_Q) >> 16;
  return t;
}

static inline int16_t barrett_reduce(int16_t a)
{
  int16_t t;
  const int16_t v = ((1<<26) + NTRUPLUS_Q/2)/NTRUPLUS_Q;

  t  = ((int32_t)v*a + (1<<25)) >> 26;
  t *= NTRUPLUS_Q;
  return a - t;
}

/*************************************************
* Name:        fqmul
*
* Description: Multiplication followed by Montgomery reduction
*
* Arguments:   - int16_t a: first factor
*              - int16_t b: second factor
*
* Returns 16-bit integer congruent to a*b*R^{-1} mod q
**************************************************/
static int16_t fqmul(int16_t a, int16_t b)
{
  return montgomery_reduce((int32_t)a*b);
}

/*************************************************
* Name:        fqinv
*
* Description: Inversion
*
* Arguments:   - int16_t a: first factor a = x * R mod q
*
* Returns 16-bit integer congruent to x^{-1} * R^3 mod q
**************************************************/

#ifdef PARAM_1
static inline int16_t fqinv(int16_t a)
{
	//q -2 =  3327 = 1100 11 11 1111
	int16_t t1,t2,t3;

	t1 = fqmul(a, a);    //10
	t2 = fqmul(t1, a);   //11
	t1 = fqmul(t2, t2);  //110
	t3 = fqmul(t1, t1);  //1100
	t1 = fqmul(t3, t3);  //11000
	t1 = fqmul(t1, t1);  //110000
	t1 = fqmul(t1, t2);  //110011
	t2 = fqmul(t1, t3);  //111111

	t1 = fqmul(t1, t1);  //1100110
	t1 = fqmul(t1, t1);  //11001100
	t1 = fqmul(t1, t1);  //110011000
	t1 = fqmul(t1, t1);  //1100110000
	t1 = fqmul(t1, t1);  //11001100000
	t1 = fqmul(t1, t1);  //110011000000
	t1 = fqmul(t1, t2);  //110011111111

	return t1;
}
#endif

#ifdef PARAM_2
static int16_t fqinv(int16_t a)
{
	//q -2 =  7679 = 1110111111111
	int16_t t1,t2,t3;

	t1 = fqmul(a, a);    //10
	t2 = fqmul(t1, t1);  //100
	
	t3 = fqmul(t1, a);   //11
	t1 = fqmul(t2, t3);  //111
	t1 = fqmul(t1, t1);  //1110
	t1 = fqmul(t1, t1);  //11100
	t1 = fqmul(t1, t1);  //111000
	t1 = fqmul(t1, t3);  //111011

	t2 = fqmul(t1, t2);  //111111

	t1 = fqmul(t1, t1);  //1110110
	t1 = fqmul(t1, t1);  //11101100
	t1 = fqmul(t1, t1);  //111011000
	t1 = fqmul(t1, t1);  //1110110000
	t1 = fqmul(t1, t1);  //11101100000
	t1 = fqmul(t1, t1);  //111011000000
	t1 = fqmul(t1, t2);  //111011111111
	t1 = fqmul(t1, t1);  //1110111111110
	t1 = fqmul(t1, a);   //1110111111111

	return t1;
}
#endif

/*************************************************
* Name:        ntt
*
* Description: Inplace number-theoretic transform (NTT) in Rq.
*              input is in standard order, output is in bitreversed order
*
* Arguments:   - int16_t r[512]: pointer to input/output vector of elements of Zq
**************************************************/

#ifdef PARAM_1
void ntt(int16_t r[NTRUPLUS_N])
{
	int16_t t;
	int16_t zeta;
	int k = 1;

	for(int step = NTRUPLUS_N/2; step >= 4; step >>= 1)
	{
		for(int start = 0; start < NTRUPLUS_N; start += (step << 1))
		{
			zeta = zetas[k++];

			for(int i = start; i < start + step; i++)
			{
				t = fqmul(zeta, r[i + step]);
				
				r[i + step] = barrett_reduce(r[i] - t);
				r[i       ] = barrett_reduce(r[i] + t);
			}
		}
	}
}
#endif

#ifdef PARAM_2
void ntt(int16_t r[NTRUPLUS_N])
{
	int16_t t;
	int16_t zeta;
	int k = 1;

	for(int step = 512; step >= 128; step >>= 1)
	{
		for(int start = 0; start < NTRUPLUS_N; start += (step << 1))
		{
			zeta = zetas[k++];

			for(int i = start; i < start + step; i++)
			{
				t = fqmul(zeta, r[i + step]);
				
				r[i + step] = (r[i] - t);
				r[i       ] = (r[i] + t);
			}
		}
	}

	for(int start = 0; start < NTRUPLUS_N; start += 128)
	{
		zeta = zetas[k++];

		for(int i = start; i < start + 64; i++)
		{
			t = fqmul(zeta, r[i + 64]);
			
			r[i + 64] = barrett_reduce(r[i] - t);
			r[i     ] = barrett_reduce(r[i] + t);
		}
	}

	for(int step = 32; step >= 8; step >>= 1)
	{
		for(int start = 0; start < NTRUPLUS_N; start += (step << 1))
		{
			zeta = zetas[k++];

			for(int i = start; i < start + step; i++)
			{
				t = fqmul(zeta, r[i + step]);
				
				r[i + step] = (r[i] - t);
				r[i       ] = (r[i] + t);
			}
		}
	}

	for(int start = 0; start < NTRUPLUS_N; start += 8)
	{
		zeta = zetas[k++];

		for(int i = start; i < start + 4; i++)
		{
			t = fqmul(zeta, r[i + 4]);
			
			r[i + 4] = barrett_reduce(r[i] - t);
			r[i    ] = barrett_reduce(r[i] + t);
		}
	}
}
#endif


/*************************************************
* Name:        invntt
*
* Description: inverse number-theoretic transform in Rq and
*              multiplication by Montgomery factor R = 2^16.
*
* Arguments:   - int16_t b[NTRUPLUS_N]: pointer to output vector of elements of Zq
*              - int16_t a[NTRUPLUS_N]: pointer to input vector of elements of Zq
**************************************************/
#ifdef PARAM_1
void invntt(int16_t r[NTRUPLUS_N])
{
	int16_t t;
	int16_t zeta;

	int k = 127;

	for(int step = 4; step <= NTRUPLUS_N/2; step <<= 1)	
	{
		for(int start = 0; start < NTRUPLUS_N; start += (step << 1))
		{
			zeta = zetas[k--];

			for(int i = start; i < start + step; i++)
			{
				t = r[i + step];

				r[i + step] = fqmul(zeta,  t - r[i]);
				r[i       ] = barrett_reduce(r[i] + t);
			}
		}
	}

	for(int i = 0; i < NTRUPLUS_N; i++)
	{
		r[i] = fqmul(512, r[i]);	//256 = R/(2^(10-2)) = 2^16/2^8
	}
}
#endif

#ifdef PARAM_2
void invntt(int16_t r[NTRUPLUS_N])
{
	int16_t t;
	int16_t zeta;

	int k = 255;

	for(int step = 4; step <= NTRUPLUS_N/2; step <<= 1)	
	{
		for(int start = 0; start < NTRUPLUS_N; start += (step << 1))
		{
			zeta = zetas[k--];

			for(int i = start; i < start + step; i++)
			{
				t = r[i + step];

				r[i + step] = fqmul(zeta,  t - r[i]);
				r[i       ] = barrett_reduce(r[i] + t);
			}
		}
	}

	for(int i = 0; i < NTRUPLUS_N; i++)
	{
		r[i] = fqmul(256, r[i]);	//256 = R/(2^(10-2)) = 2^16/2^8
	}
}
#endif
/*************************************************
* Name:        basemul
*
* Description: Multiplication of polynomials in Zq[X]/(X^4-zeta)
*              used for multiplication of elements in Rq in NTT domain
*
* Arguments:   - int16_t c[4]: pointer to the output polynomial
*              - const int16_t a[4]: pointer to the first factor
*              - const int16_t b[4]: pointer to the second factor
*              - int16_t zeta: integer defining the reduction polynomial
**************************************************/
#ifdef PARAM_1
void basemul(int16_t c[4], const int16_t a[4], const int16_t b[4], int16_t zeta)
{
	c[0] = montgomery_reduce(a[1]*b[3]+a[2]*b[2]+a[3]*b[1]);
	c[1] = montgomery_reduce(a[2]*b[3]+a[3]*b[2]);
	c[2] = montgomery_reduce(a[3]*b[3]);
	c[3] = montgomery_reduce(a[0]*b[3]+a[1]*b[2]+a[2]*b[1]+a[3]*b[0]);

	c[0] = montgomery_reduce(a[0]*b[0]+c[0]*zeta);
	c[1] = montgomery_reduce(a[0]*b[1]+a[1]*b[0]+c[1]*zeta);
	c[2] = montgomery_reduce(a[0]*b[2]+a[1]*b[1]+a[2]*b[0]+c[2]*zeta);
}
#endif

#ifdef PARAM_2
void basemul(int16_t c[4], const int16_t a[4], const int16_t b[4], int16_t zeta)
{
	c[0] = montgomery_reduce(a[1]*b[3]+a[2]*b[2]+a[3]*b[1]);
	c[1] = montgomery_reduce(a[2]*b[3]+a[3]*b[2]);
	c[2] = montgomery_reduce(a[3]*b[3]);
	c[3] = montgomery_reduce(a[0]*b[3]+a[1]*b[2]+a[2]*b[1]+a[3]*b[0]);

	c[0] = montgomery_reduce(a[0]*b[0]+c[0]*zeta);
	c[1] = montgomery_reduce(a[0]*b[1]+a[1]*b[0]+c[1]*zeta);
	c[2] = montgomery_reduce(a[0]*b[2]+a[1]*b[1]+a[2]*b[0]+c[2]*zeta);
}
#endif

/*************************************************
* Name:        baseinv
*
* Description: Inversion of polynomial in Zq[X]/(X^4-zeta)
*              used for inversion of element in Rq in NTT domain
*
* Arguments:   - int16_t b[8]: pointer to the output polynomial
*              - const int16_t a[8]: pointer to the input polynomial
*              - int16_t zeta: integer defining the reduction polynomial
**************************************************/
int baseinv(int16_t b[8], const int16_t a[8], int16_t zeta)
{
	int r1, r2;
	int16_t t0, t1, t2;
	int16_t s0, s1, s2;
	int32_t T0, T1, T2, T3;
	int32_t S0, S1, S2, S3;

	t0 = montgomery_reduce(a[2]*a[2] - 2*a[1]*a[3]);
	s0 = montgomery_reduce(a[6]*a[6] - 2*a[5]*a[7]);
	t0 = montgomery_reduce(a[0]*a[0] + t0*zeta);
	s0 = montgomery_reduce(a[4]*a[4] - s0*zeta);
	t1 = montgomery_reduce(a[3]*a[3]);
	s1 = montgomery_reduce(a[7]*a[7]);
	t1 = montgomery_reduce(2*a[0]*a[2] - a[1]*a[1] - t1*zeta);
	s1 = montgomery_reduce(2*a[4]*a[6] - a[5]*a[5] + s1*zeta);

	t2 = montgomery_reduce(t1*t1);
	s2 = montgomery_reduce(s1*s1);
	t2 = montgomery_reduce(t0*t0 - t2*zeta);
	s2 = montgomery_reduce(s0*s0 + s2*zeta);

	t2 = fqinv(t2);
	s2 = fqinv(s2);

	r1 = (uint16_t)t2;
	r2 = (uint16_t)s2;
	r1 = (uint32_t)(-r1) >> 31;
	r2 = (uint32_t)(-r2) >> 31;

	t0 = fqmul(t0,t2);
	s0 = fqmul(s0,s2);
	t1 = fqmul(t1,t2);
	s1 = fqmul(s1,s2);
	t2 = fqmul(t1,zeta);
	s2 = fqmul(s1,zeta);

	T0 = a[0]*t0 - a[2]*t2;
	S0 = a[4]*s0 + a[6]*s2;
	T1 = a[1]*t0 - a[3]*t2;
	S1 = a[5]*s0 + a[7]*s2;
	T2 = a[2]*t0 - a[0]*t1;
	S2 = a[6]*s0 - a[4]*s1;
	T3 = a[3]*t0 - a[1]*t1;
	S3 = a[7]*s0 - a[5]*s1;

	b[0] =  montgomery_reduce(T0);
	b[4] =  montgomery_reduce(S0);
	b[1] = -montgomery_reduce(T1);
	b[5] = -montgomery_reduce(S1);
	b[2] =  montgomery_reduce(T2);
	b[6] =  montgomery_reduce(S2);
	b[3] = -montgomery_reduce(T3);
	b[7] = -montgomery_reduce(S3);

	return (r1 & r2) - 1;
}