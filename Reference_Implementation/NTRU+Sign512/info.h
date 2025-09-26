#ifndef INFO_H
#define INFO_H

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <math.h>
#include "params.h"

int64_t __rdtsc(void);



#ifdef PARAM_1
// discrete Gaussian
#define K       64
#define KMASK   0x3f
#define w1		15
static const __int128 CDT[w1 + 1] = {
	((__int128)0x181d03 << 64) ^ (__int128)0x2ab96eaa07be703b,
	((__int128)0x2c78d7 << 64) ^ (__int128)0x41a6d3e37c638887,
	((__int128)0x38b983 << 64) ^ (__int128)0x1e6c291aec404f67,
	((__int128)0x3dfb2a << 64) ^ (__int128)0xe519ac899c537c1d,
	((__int128)0x3f96ae << 64) ^ (__int128)0x7ae99431c810af82,
	((__int128)0x3ff063 << 64) ^ (__int128)0x350b00544e285ff7,
	((__int128)0x3ffe53 << 64) ^ (__int128)0xb14ef40243bddad2,
	((__int128)0x3fffde << 64) ^ (__int128)0xf61cc08e1bc70f26,
	((__int128)0x3ffffe << 64) ^ (__int128)0x2bbf7238a2f11dad,
	((__int128)0x3fffff << 64) ^ (__int128)0xed6e8af0bbb2e9d0,
	((__int128)0x3fffff << 64) ^ (__int128)0xff792dc8822b344a,
	((__int128)0x3fffff << 64) ^ (__int128)0xfffd448548aa0aed,
	((__int128)0x3fffff << 64) ^ (__int128)0xfffff5e0b46461ef,
	((__int128)0x3fffff << 64) ^ (__int128)0xffffffe53d01a7a1,
	((__int128)0x3fffff << 64) ^ (__int128)0xffffffffcdc85638,
	((__int128)0x3fffff << 64) ^ (__int128)0xffffffffffffffff
};
#define DegG 11
static const uint64_t Coef_G[DegG + 1] = {
	0x001c6fa92d13b80d,
	0x0015d439399cba51,
	0x000ae8cb6ed537d5,
	0x001257783c8d5633,
	0x000d901db61f08df,
	0x001187e3481e5965,
	0x0009b5d470da1fd1,
	0x00023db3eb673c53,
	0x001a7b1b052c3f31,
	0x00075573fdb8e9a5,
	0x0015aa2fe645e647,
	0x0400000000000000 // 2^{58}
};
static const int32_t expDiff_G[DegG] = { 19, 19, 17, 18, 17, 18, 19, 13, 18, 14, 9 };

// rejection sampling
#define DegR 12
static const uint64_t Coef_R[DegR + 1] = {
	0x03c0d5a3b08576e7,
	0x02cf6d287cbb320b,
	0x02fa96e4a6d72f49,
	0x016193581fbaafdb,
	0x024bbaa976f2f685,
	0x03641f58f980f92f,
	0x00461fab123ba2ed,
	0x004daea50fb3554f,
	0x023db3ebcc4039ab,
	0x000d3d8d82a519bb,
	0x003aab9fedc863b7,
	0x02b545fcc8bce3cf,
	0x03ffffffffffffff
};
static const int32_t expDiff_R[DegR] = { 37, 36, 37, 35, 35, 39, 35, 32, 40, 32, 30, 32};

#define vartheta2 58
#define BitGamma 7 // gamma = ceiling[beta2 / (sigma * alpha * ln2)] = 81 < 2^{7}
#define ThetaChat 18
#define ValChat 0x00010618a5a7 // chat = ValChat / 2^{ThetaChat}
#define ThetaEps 52
#define EpsVal 0xe5e9 // epsilon = EpsVal / 2^{ThetaEps}
#endif

#ifdef PARAM_2
// discrete Gaussian
#define K       128
#define KMASK   0x7f
#define w1		14
// 아래 CDT랑 Gaussian용 exp 마무리
static const __int128 CDT[w1 + 1] = {
	((__int128)0x341184 << 64) ^ (__int128)0xc0987fc4fdd127b2,
	((__int128)0x5e7e8c << 64) ^ (__int128)0x0d63c69fe0922494,
	((__int128)0x7571fc << 64) ^ (__int128)0xca9fb44368810d3d,
	((__int128)0x7db02f << 64) ^ (__int128)0x10ef7e9c42c08807,
	((__int128)0x7fa75b << 64) ^ (__int128)0xa4cebd2e6c7ad473,
	((__int128)0x7ff703 << 64) ^ (__int128)0xf7426e922b7effbb,
	((__int128)0x7fff63 << 64) ^ (__int128)0x4a83a3b748e406b6,
	((__int128)0x7ffff8 << 64) ^ (__int128)0xda7d82666f46329d,
	((__int128)0x7fffff << 64) ^ (__int128)0xc84fa7300274a719,
	((__int128)0x7fffff << 64) ^ (__int128)0xfedee03412eb11b3,
	((__int128)0x7fffff << 64) ^ (__int128)0xfffc19181c04826d,
	((__int128)0x7fffff << 64) ^ (__int128)0xfffff70930a4b3e7,
	((__int128)0x7fffff << 64) ^ (__int128)0xfffffff250597c67,
	((__int128)0x7fffff << 64) ^ (__int128)0xfffffffff2276863,
	((__int128)0x7fffff << 64) ^ (__int128)0xffffffffffffffff
};
#define DegG 12
static const uint64_t Coef_G[DegG + 1] = {
	0x0007e963cab3fcf9,
	0x15e02799d5b93bff,
	0x0009efe611cba877,
	0x001eaf3a9f882bd7,
	0x1518b28dd2b1e5ff,
	0x000ce0ac5252f065,
	0x000dc1401d9a0827,
	0x00064be6aca20e61,
	0x000266e4374e20e9,
	0x001774ccac3c6a5f,
	0x0015798ee23089a1,
	0x001a36e2eb1c4327,
	0x2000000000000000 // 2^{61}
};
static const int32_t expDiff_G[DegG] = { 11, 29, 18, 12, 28, 19, 20, 20, 15, 18, 17, 8};

// rejection sampling
#define DegR 13
static const uint64_t Coef_R[DegR + 1] = {
	0x0951ffe75ebb4053,
	0x063f020eba76060d,
	0x02fc1acf540be4c1,
	0x0a10a8f4c4428809,
	0x01eba0f42873041d,
	0x0a8ca6b483445151,
	0x06705920b001672d,
	0x0370501b43a632f1,
	0x0c97cd5c08a72679,
	0x099b90dd482a5b21,
	0x0bba66561e9bfdd3,
	0x055e63b88c230e6b,
	0x068db8bac710cb29,
	0x0fffffffffffffff // 2^{61}-1
};
static const int32_t expDiff_R[DegR] = { 42, 42, 39, 43, 38, 41, 41, 38, 40, 39, 40, 38, 36, };

#define vartheta2 60
#define BitGamma 7 // gamma = ceiling[beta2 / (sigma * alpha * ln2)] = 123 < 2^{7}
#define ThetaChat 21
#define ValChat 0x001b1378c841 // chat = ValChat / 2^{ThetaChat}
#define ThetaEps 53
#define EpsVal 0x72eb // epsilon = EpsVal / 2^{ThetaEps}
#endif


#endif