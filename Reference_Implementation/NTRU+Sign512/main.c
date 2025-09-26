#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "poly.h"
#include "sign.h"
#include "packing.h"
#include "sampler.h"
#include "gaussian.h"
#include "ntt.h"
#include "symmetric.h"
#include "randombytes.h"
#include "decompose.h"
#include "rejection.h"

#define TEST_LOOP 1000

uint64_t t[TEST_LOOP];
extern unsigned long long rejwctr;
extern unsigned long long rejyzctr;
extern unsigned long long ctr_keygen;
extern unsigned long long ctr_sign;
extern int64_t maxg;
int total_count[21] = {0};

static int64_t cpucycles(void)
{
  unsigned int hi, lo;

    __asm__ __volatile__ ("rdtsc\n\t" : "=a" (lo), "=d"(hi));

    return ((int64_t)lo) | (((int64_t)hi) << 32);
}

static int test_Correctness()
{
	double rejctr=.0;

	uint8_t pk[CRYPTO_PUBLICKEYBYTES] = {0};
	uint8_t sk[CRYPTO_SECRETKEYBYTES] = {0};

	unsigned long long siglen = 0;
	uint8_t sig[CRYPTO_BYTES] = {0};
	uint8_t msg[SEEDBYTES] = {0};
	randombytes(msg, SEEDBYTES);

	int result = 0;
	unsigned long long kcycles, ecycles, dcycles;
	unsigned long long cycles1, cycles2, cycles3, cycles4, cycles5, cycles6;

	kcycles=0;
	ecycles=0;
	dcycles=0;
	for(int i = 0; i < TEST_LOOP; i++)
	{
		cycles1 = cpucycles();
		crypto_sign_keypair(pk, sk);
		cycles2 = cpucycles();
		kcycles += cycles2-cycles1;

		cycles3 = cpucycles();
		crypto_sign(sig, &siglen, msg, SEEDBYTES, sk);
		cycles4 = cpucycles();
		ecycles += cycles4-cycles3;
        
		rejctr+=ctr_sign;

		cycles5 = cpucycles();
		result += crypto_sign_verify(msg, SEEDBYTES, sig, siglen, pk);
		cycles6 = cpucycles();
		dcycles += cycles6-cycles5;

	}  

	printf("  KeyGen runs in ................. %8lld cycles", kcycles/TEST_LOOP);
	printf("\n"); 

	printf("  Sign runs in ................. %8lld cycles", ecycles/TEST_LOOP);
	printf("\n"); 

	printf("  Verify runs in ................. %8lld cycles", dcycles/TEST_LOOP);
	printf("\n"); 

	printf("Acceptance rate of Signing: %.2f\n",(double)TEST_LOOP/rejctr);
	printf("\n");   

	if(result == 0)    {
		printf("Verification Success!\n");
	}

	else{
		printf("Verification Fail\n");
	}

	return 0;
}


static int test_speed(void)
{    
    unsigned long long i;

    uint8_t msg[SEEDBYTES] = {0};
    randombytes(msg, SEEDBYTES);
    uint8_t seedbuf[CRHBYTES] = {0};

    unsigned long long kcycles;
    unsigned long long cycles1, cycles2;
    
  
    int nonce = 0;
    poly y1, y2;
    poly cs1, cs2, z1, z2;
    int b1, b2, b3;
    int64_t res;
    int64_t val;

    /*printf("Gaussian_rand\n");
    kcycles=0;
    cycles1 = cpucycles();
    for(i = 0; i < TEST_LOOP; ++i)    
      poly_sample_y2(&y1, &y2);
    cycles2 = cpucycles();
    kcycles = cycles2-cycles1;
    printf("  gaussian_rand runs in ................. %8lld cycles", kcycles/TEST_LOOP);
    printf("\n"); 
    
    printf("Gaussian_shake\n");
    kcycles=0;
    cycles1 = cpucycles();
    for(i = 0; i < TEST_LOOP; ++i){    
      approxG_test(&y3, seedbuf, nonce++);
    }
    cycles2 = cpucycles();
    kcycles = cycles2-cycles1;
    printf("  gaussian_shake runs in ................. %8lld cycles", kcycles/TEST_LOOP);
    printf("\n"); */
 
    


    printf("sample_y_Gaussian_buf\n");
    kcycles=0;
    cycles1 = cpucycles();
    for(i = 0; i < TEST_LOOP; ++i){    
      sample_y(&y1, &y2, seedbuf, nonce++);
    }
    cycles2 = cpucycles();
    kcycles = cycles2-cycles1;
    printf("  sample_y_Gaussian_buf runs in ................. %8lld cycles", kcycles/TEST_LOOP);
    printf("\n"); 

    poly_cbd1(&cs1, seedbuf, 1);
    poly_cbd1(&cs2, seedbuf, 2);
    
    
    poly_add(&z1, &y1, &cs1);
    poly_add(&z2, &y2, &cs2);
    
    printf("poly_sqnorm2\n");
    kcycles=0;
    cycles1 = cpucycles();
    for(i = 0; i < TEST_LOOP; ++i){    
	    val = B_scsquare - poly_sqnorm2(&cs1, &cs2);
    }

    cycles2 = cpucycles();
    kcycles = cycles2-cycles1;
    printf("  poly_sqnorm2 runs in ................. %8lld cycles", kcycles/TEST_LOOP);
    printf("\n"); 


    printf("SampleBernExp\n");
    kcycles=0;
    cycles1 = cpucycles();
    for(i = 0; i < TEST_LOOP; ++i){    
	    b1 = SampleBernExp(val, seedbuf);
    }
    cycles2 = cpucycles();
    kcycles = cycles2-cycles1;
    printf("  SampleBernExp runs in ................. %8lld cycles", kcycles/TEST_LOOP);
    printf("\n"); 


    printf("innerproduct\n");
    kcycles=0;
    cycles1 = cpucycles();
    for(i = 0; i < TEST_LOOP; ++i){    
	    res = innerproduct(&z1, &cs1) + innerproduct(&z2, &cs2);
    }
    cycles2 = cpucycles();
    kcycles = cycles2-cycles1;
    printf("  innerproduct runs in ................. %8lld cycles", kcycles/TEST_LOOP);
    printf("\n"); 


    printf("SampleBernCosh\n");
    kcycles=0;
    cycles1 = cpucycles();
    for(i = 0; i < TEST_LOOP; ++i){    
	    b2 = SampleBernCosh(res, seedbuf);
    }
    cycles2 = cpucycles();
    kcycles = cycles2-cycles1;
    printf("  SampleBernCosh runs in ................. %8lld cycles", kcycles/TEST_LOOP);
    printf("\n"); 

    printf("rejection sampling\n");
    kcycles=0;
    cycles1 = cpucycles();
    for(i = 0; i < TEST_LOOP; ++i){    
	    val = B_scsquare - poly_sqnorm2(&cs1, &cs2);
	    //printf("B_scsquare = %d, poly_sqnorm2 = %ld \n", B_scsquare, poly_sqnorm2(&cs1, &cs2));
	    //printf("sqnorm2 = %d \n", poly_sqnorm2(&cs1, &cs2));
	    b1 = SampleBernExp(val, seedbuf);
	    //printf("b1 = %d \n", b1);
	    res = innerproduct(&z1, &cs1) + innerproduct(&z2, &cs2);
	    res <<= 1;
	    b2 = SampleBernCosh(res, seedbuf);
	    //printf("b2 = %d \n", b2);
	    b3 = b1*b2;
    }
    cycles2 = cpucycles();
    kcycles = cycles2-cycles1;
    printf("  rejection sampling runs in ................. %8lld cycles", kcycles/TEST_LOOP);
    printf("\n"); 

    
    //compute ||Sc||^2 for rejection sampling
    
    /*
    printf("poly_baseinv \n");
    kcycles=0;
    cycles1 = cpucycles();
    for(i = 0; i < TEST_LOOP; ++i) 
      r |= poly_baseinv(&y1, &y1); 
    cycles2 = cpucycles();
    kcycles = cycles2-cycles1;
    printf("  poly_baseinv runs in ................. %8lld cycles", kcycles/TEST_LOOP);
    printf("\n");*/

  return 0;
}


int main(void)
{
   printf("=================== PARAMETERS ===================\n");
   printf("ALGORITHM_NAME  : %s\n", CRYPTO_ALGNAME);
   printf("PUBLICKEYBYTES  : %d\n", CRYPTO_PUBLICKEYBYTES);
   printf("SECRETKEYBYTES  : %d\n", CRYPTO_SECRETKEYBYTES);
   printf("SIGNATUREBYTES  : %d\n", CRYPTO_BYTES);
   printf("\n");

    test_Correctness();
    //test_speed();
    return 0;
}
