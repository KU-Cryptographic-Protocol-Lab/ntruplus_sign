#include "gaussian.h"
#include "fips202.h"
#include "symmetric.h"

#include <stdint.h>
#include "params.h"
#include "poly.h"
#include <stdlib.h>

#ifdef PARAM_1
static int approxBaseGplus(const unsigned char *buf)
{
    int result = 0;

	__int128 r;

	// 85 bits random
	r =  (__int128)buf[0];
	r ^= (__int128)buf[1] <<  8;
	r ^= (__int128)buf[2] << 16;
	r ^= (__int128)buf[3] << 24;
	r ^= (__int128)buf[4] << 32;
	r ^= (__int128)buf[5] << 40;
	r ^= (__int128)buf[6] << 48;
	r ^= (__int128)buf[7] << 56;
	r ^= (__int128)buf[8] << 64;
	r ^= (__int128)buf[9] << 72;
	r ^= (__int128)(buf[10] & 0x1f) << 80;
	//printf(" buf[0] = %d \n ", buf[0]);
	// comparing and counting
	for (int i = 0; i < w1; i++) {
		result += (r > CDT[i]) ? 1 : 0;
	}

	return result;
}
#endif

#ifdef PARAM_2
static int approxBaseGplus(const unsigned char *buf)
{
    int result = 0;

	__int128 r;

	// 86 bits random (all parameters use 86 bits random)
	r = (__int128)buf[0];
	r ^= (__int128)buf[1] << 8;
	r ^= (__int128)buf[2] << 16;
	r ^= (__int128)buf[3] << 24;
	r ^= (__int128)buf[4] << 32;
	r ^= (__int128)buf[5] << 40;
	r ^= (__int128)buf[6] << 48;
	r ^= (__int128)buf[7] << 56;
	r ^= (__int128)buf[8] << 64;
	r ^= (__int128)buf[9] << 72;
	r ^= (__int128)(buf[10] & 0x3f) << 80;
	//printf(" buf[0] = %d \n ", buf[0]);
	// comparing and counting
	for (int i = 0; i < w1; i++) {
		result += (r > CDT[i]) ? 1 : 0;
	}

	return result;
}
#endif

static inline int approxGplus(int16_t *a, const unsigned char *buf)
{
	int y1, y0, x;

	y1 = approxBaseGplus(buf); //10
	y0 = buf[11] & KMASK; //11
	x = y0 * (y0 + (y1 << 1) * K);

	if (SampleBernExpSimple(x, buf+12) == 1) { //19
		*a = y1 * K + y0;
		return 1;
	}
	
	return 0;
	
}

#ifdef PARAM_1
#define POLY_Y_NBLOCKS ((13000 + STREAM256_BLOCKBYTES - 1)/STREAM256_BLOCKBYTES)
void approxG(poly *a, const unsigned char seed[CRHBYTES], uint16_t nonce)
{
	int b, check;
	int pos = 0;
	int ctr = 0;
	int i = 0;
	int off = 0;
    int16_t tmp = 0;
    uint8_t buf[POLY_Y_NBLOCKS*STREAM256_BLOCKBYTES];
    int16_t buflen = POLY_Y_NBLOCKS*STREAM256_BLOCKBYTES;
    stream256_state state;
    
    stream256_init(&state, seed, nonce);
	stream256_squeezeblocks(buf, POLY_Y_NBLOCKS, &state);
 	
  	while(pos + 20 < buflen) {
    	
    	check = approxGplus(&tmp, &buf[pos]);
    	pos += 19;

    	if (pos >= buflen) break;
    	b = buf[pos] & 0x1;

    	while (tmp == 0 && b == 0) {
			pos++;
			if (pos+20 >= buflen) break;
			check = approxGplus(&tmp, &buf[pos]);
			pos += 19;
			if (pos >= buflen) break;
			b = buf[pos] & 0x1;		
		}
    	
    	if (pos >= buflen) break;
    	b = (buf[pos] >> 1) & 0x1;
    	a->coeffs[ctr] = ((b << 1) - 1) * tmp;
    	ctr += check;
    	pos ++;


    	if (ctr == NTRUPLUS_N){
    		break;
    	}

	}
	while(ctr < NTRUPLUS_N) {
		off = buflen % 20;
		for (i = 0; i < off; ++i)
            buf[i] = buf[buflen - off + i];

		stream256_squeezeblocks(buf, 1, &state);
		buflen = STREAM256_BLOCKBYTES + off;

		
		pos = 0;
		while(pos + 20 < buflen) {
			
			check = approxGplus(&tmp, &buf[pos]);
	    	pos += 19;
	    	if (pos >= buflen) break;
	    	b = buf[pos] & 0x1;
	    	
	    	while (tmp == 0 && b == 0) {
	    		pos++;
	    		if (pos+20 >= buflen) break;
				check = approxGplus(&tmp, &buf[pos]);
				pos += 19;
				if (pos >= buflen) break;
				b = buf[pos] & 0x1;
			}
			if (pos >= buflen) break;
    		b = (buf[pos] >> 1) & 0x1;
    		a->coeffs[ctr] = ((b << 1) - 1) * tmp;
    		ctr += check;
    		pos++;

    		if (ctr == NTRUPLUS_N){
    		break;
    	}


		}
		
		if (ctr == NTRUPLUS_N){
    		break;
    	}
	}
	
}
#endif

#ifdef PARAM_2
#define POLY_Y_NBLOCKS2 ((24252 + STREAM256_BLOCKBYTES - 1)/STREAM256_BLOCKBYTES)
void approxG(poly *a, const unsigned char seed[CRHBYTES], uint16_t nonce)
{
	int b, check;
	int pos = 0;
	int ctr = 0;
	int i = 0;
	int off = 0;
    int16_t tmp = 0;
    uint8_t buf[POLY_Y_NBLOCKS2*STREAM256_BLOCKBYTES];
    int16_t buflen = POLY_Y_NBLOCKS2*STREAM256_BLOCKBYTES;
    stream256_state state;
    
    stream256_init(&state, seed, nonce);
	stream256_squeezeblocks(buf, POLY_Y_NBLOCKS2, &state);
 	
  	while(pos + 20 < buflen) {
    	
    	check = approxGplus(&tmp, &buf[pos]);
    	pos += 19;

    	if (pos >= buflen) break;
    	b = buf[pos] & 0x1;

    	while (tmp == 0 && b == 0) {
			pos++;
			if (pos+20 >= buflen) break;
			check = approxGplus(&tmp, &buf[pos]);
			pos += 19;
			if (pos >= buflen) break;
			b = buf[pos] & 0x1;		
		}
    	
    	if (pos >= buflen) break;
    	b = (buf[pos] >> 1) & 0x1;
    	a->coeffs[ctr] = ((b << 1) - 1) * tmp;
    	ctr += check;
    	pos ++;


    	if (ctr == NTRUPLUS_N){
    		break;
    	}

	}
	while(ctr < NTRUPLUS_N) {
		off = buflen % 20;
		for (i = 0; i < off; ++i)
            buf[i] = buf[buflen - off + i];

		stream256_squeezeblocks(buf, 1, &state);
		buflen = STREAM256_BLOCKBYTES + off;

		
		pos = 0;
		while(pos + 20 < buflen) {
			
			check = approxGplus(&tmp, &buf[pos]);
	    	pos += 19;
	    	if (pos >= buflen) break;
	    	b = buf[pos] & 0x1;
	    	
	    	while (tmp == 0 && b == 0) {
	    		pos++;
	    		if (pos+20 >= buflen) break;
				check = approxGplus(&tmp, &buf[pos]);
				pos += 19;
				if (pos >= buflen) break;
				b = buf[pos] & 0x1;
			}
			if (pos >= buflen) break;
    		b = (buf[pos] >> 1) & 0x1;
    		a->coeffs[ctr] = ((b << 1) - 1) * tmp;
    		ctr += check;
    		pos++;

    		if (ctr == NTRUPLUS_N){
    		break;
    	}


		}
		
		if (ctr == NTRUPLUS_N){
    		break;
    	}
	}
	
}
#endif


void sample_y(poly *a1, poly *a2, const unsigned char seed[CRHBYTES], uint16_t nonce)
{
	
	approxG(a1, seed, nonce++);
	approxG(a2, seed, nonce++);
		
}


