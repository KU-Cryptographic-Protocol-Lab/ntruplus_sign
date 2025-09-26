#include "sampler.h"
#include "symmetric.h"
#include <stdint.h>

void hash_H(uint8_t *c_bin, const uint8_t highbits[POLY_HIGHBITS_PACKEDBYTES], const uint8_t mu[SEEDBYTES])
{ // Hash-based function H to generate c'
  uint8_t t[POLY_HIGHBITS_PACKEDBYTES + SEEDBYTES];
  

  for (int i=0; i<POLY_HIGHBITS_PACKEDBYTES; i++) { 
    t[i] = highbits[i];      
  }  

  memcpy(&t[POLY_HIGHBITS_PACKEDBYTES], mu, SEEDBYTES);
  shake256(c_bin, CRYPTO_C_BYTES, t, POLY_HIGHBITS_PACKEDBYTES+SEEDBYTES);
}


void encode_c(poly *c, const uint8_t *c_bin)
{ // Encoding of c' by mapping the output of the hash function H to an N-element vector with entries {-1,0,1}  
    unsigned int i, b, pos = 0;
    uint8_t buf[XOF256_BLOCKBYTES];
    xof256_state state;

#ifndef PARAM_2
	int check;
#endif
  
  // Use the hash value as key to generate some randomness
  	shake256_init(&state);
    shake256_absorb_once(&state, c_bin, CRYPTO_RANDOMBYTES);
	shake256_squeezeblocks(buf,1,&state);

  // Use rejection sampling to determine positions to be set in the new vector
  for (i = 0; i < NTRUPLUS_N; ++i)
        c->coeffs[i] = 0;

    for (i = NTRUPLUS_N - NTRUPLUS_TAU; i < NTRUPLUS_N; ++i) {
        do {
            if (pos >= XOF256_BLOCKBYTES) {
				shake256_squeezeblocks(buf,1,&state);
                pos = 0;
            }

            b = (buf[pos] << 8)  | (buf[pos+1]); //need to check
            pos += 2;

#ifdef PARAM_1
		check = (b >= NTRUPLUS_N*85);
		if(check) continue;
#endif
#ifdef PARAM_3
		check = (b >= NTRUPLUS_N*50);
		if(check) continue;
#endif
		b = b % NTRUPLUS_N;

        } while (b > i);

        c->coeffs[i] = c->coeffs[b];
        c->coeffs[b] = 1;
    }
}


