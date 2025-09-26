#ifndef PARAMS_H
#define PARAMS_H

#define PARAM_1
#define NTT_s

#define SEEDBYTES 32
#define CRHBYTES 64
#define CRYPTO_C_BYTES 32
#define CRYPTO_RANDOMBYTES 32

#ifdef PARAM_1
#define CRYPTO_ALGNAME "NTRU+SIGN512"
#define NTRUPLUS_N 512
#define NTRUPLUS_Q 3329
#define NTRUPLUS_logQ 12
#define NTRUPLUS_T 1
#define NTRUPLUS_TAU 20
#define NTRUPLUS_GAMMA 37.77 //40.50(0.9)
#define NTRUPLUS_d 7
#define B_sc 169 //181
#define B2bound 16000000
#define Binftybound 766

#define qhat 1665
#define Rqhat -522
#define Rmodq -1044
#define R2modq -1976
#define ntt_level 64

#define BASE_ENC_HB_Z1 61
#define BASE_ENC_H 61
#define ENC_bound 269
#define CRYPTO_BYTES 751

#define NTRUPLUS_M (NTRUPLUS_N / NTRUPLUS_TAU)
#define NTRUPLUS_R (NTRUPLUS_N % NTRUPLUS_TAU)

#endif

#ifdef PARAM_2
#define CRYPTO_ALGNAME "NTRU+SIGN1024"
#define NTRUPLUS_N 1024
#define NTRUPLUS_Q 7681
#define NTRUPLUS_logQ 13
#define NTRUPLUS_T 1
#define NTRUPLUS_TAU 36
#define NTRUPLUS_GAMMA 56.71 // 60.15
#define NTRUPLUS_d 8
#define B_sc 341 //360
#define B2bound 100000000
#define Binftybound 1790

#define qhat 3841
#define Rqhat 2044
#define Rmodq -3593
#define R2modq -2112
#define ntt_level 128

#define BASE_ENC_HB_Z1 150
#define BASE_ENC_H 150
#define ENC_bound 493
#define CRYPTO_BYTES 1551 // 32 + 8*1024/8 + 493 + 2

#define NTRUPLUS_M 28 //(NTRUPLUS_N / NTRUPLUS_TAU)
#define NTRUPLUS_R 16 //(NTRUPLUS_N % NTRUPLUS_TAU)

#endif


#define NTRUPLUS_P ((NTRUPLUS_Q - (NTRUPLUS_T)) / ALPHA_HINT) 

#define NTRUPLUS_GAMMA2N (NTRUPLUS_GAMMA * NTRUPLUS_GAMMA * NTRUPLUS_N -1)
#define halfN NTRUPLUS_N/2
#define halfp NTRUPLUS_P/2
#define ALPHA_HINT (1 << NTRUPLUS_d)
#define HALF_ALPHA_HINT (ALPHA_HINT >> 1) // ALPHA / 2
#define B_scsquare (B_sc*B_sc)

#define POLYC_PACKEDBYTES (NTRUPLUS_N >> 3)       // 1bit * N / 8bits
#define POLYQ_PACKEDBYTES ((NTRUPLUS_logQ * NTRUPLUS_N) >> 3) 
#define POLY_HIGHBITS_PACKEDBYTES (NTRUPLUS_N * 3 / 4)
#define CRYPTO_PUBLICKEYBYTES (POLYQ_PACKEDBYTES)                                      // seed + b
#define CRYPTO_SECRETKEYBYTES (3 * POLYQ_PACKEDBYTES + SEEDBYTES)  // pk + s + K

#endif