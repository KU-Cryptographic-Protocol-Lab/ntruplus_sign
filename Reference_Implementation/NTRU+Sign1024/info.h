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
	((__int128)0x0c0e81 << 64) ^ (__int128)0x955cb75503df3821,
	((__int128)0x163c6b << 64) ^ (__int128)0xa0d369f1be31c447,
	((__int128)0x1c5cc1 << 64) ^ (__int128)0x8f36148d762027b7,
	((__int128)0x1efd95 << 64) ^ (__int128)0x728cd644ce29be12,
	((__int128)0x1fcb57 << 64) ^ (__int128)0x3d74ca18e40857c4,
	((__int128)0x1ff831 << 64) ^ (__int128)0x9a85802a27142ffe,
	((__int128)0x1fff29 << 64) ^ (__int128)0xd8a77a0121deed6b,
	((__int128)0x1fffef << 64) ^ (__int128)0x7b0e60470de38795,
	((__int128)0x1fffff << 64) ^ (__int128)0x15dfb91c51788ed8,
	((__int128)0x1fffff << 64) ^ (__int128)0xf6b745785dd974e9,
	((__int128)0x1fffff << 64) ^ (__int128)0xffbc96e441159a26,
	((__int128)0x1fffff << 64) ^ (__int128)0xfffea242a4550577,
	((__int128)0x1fffff << 64) ^ (__int128)0xfffffaf05a3230f8,
	((__int128)0x1fffff << 64) ^ (__int128)0xfffffff29e80d3d1,
	((__int128)0x1fffff << 64) ^ (__int128)0xffffffffe6e42b1c,
	((__int128)0x1fffff << 64) ^ (__int128)0xffffffffffffffff
};
#define DegG 11
static const uint64_t Coef_G[DegG + 1] = {
	0x071bea4b44ee0340,
	0x05750e4e672e9440,
	0x057465b76a9bea80,
	0x0495de0f23558cc0,
	0x06c80edb0f846f80,
	0x0461f8d207965940,
	0x04daea386d0fe880,
	0x047b67d6ce78a600,
	0x069ec6c14b0fcc40,
	0x075573fdb8e9a500,
	0x056a8bf9917991c0,
	0x0800000000000000 // 2^{59}
};

static const int32_t expDiff_G[DegG] = { 19, 18, 18, 17, 18, 17, 17, 16, 16, 16, 14};

// rejection sampling
#define DegR 12
static const uint64_t Coef_R[DegR + 1] = {
	0x0ea261465b064dbf,
	0x0b2ca5ef52ad1a7f,
	0x0be7cecaae82fdc0,
	0x0b0c650d441ac2ff,
	0x092ee808bf89d600,
	0x0d907d3e46d35900,
	0x08c3f561cc379140,
	0x09b5d4a1f80ed740,
	0x08f6cfaf31372400,
	0x0d3d8d82a51ca600,
	0x0eaae7fb7218fb80,
	0x0ad517f322f38f40,
	0x1000000000000000 // 2^{60}
};
static const int32_t expDiff_R[DegR] = { 37, 36, 36, 36, 35, 36, 35, 35, 34, 34, 34, 32 };

#define vartheta2 60
#define BitGamma 7 // gamma = ceiling[beta2 / (sigma * alpha * ln2)] = 78 < 2^{7}
#define ThetaChat 18
#define ValChat 0x000000010618a5a7 // chat = ValChat / 2^{ThetaChat} (Note that ValChat < 2^{33})
#define ThetaEps 52
#define EpsVal 58857 // epsilon = EpsVal / 2^{ThetaEps} (Note that EpsVal < 2^{16})
#endif

#ifdef PARAM_2
// discrete Gaussian
#define K       128
#define KMASK   0x7f
#define w1		14
// 아래 CDT랑 Gaussian용 exp 마무리
static const __int128 CDT[w1 + 1] = {
	((__int128)0x1a08c2 << 64) ^ (__int128)0x604c3fe27ee893da,
	((__int128)0x2f3f46 << 64) ^ (__int128)0x06b1e34ff049124b,
	((__int128)0x3ab8fe << 64) ^ (__int128)0x654fda21b440869f,
	((__int128)0x3ed817 << 64) ^ (__int128)0x8877bf4e21604404,
	((__int128)0x3fd3ad << 64) ^ (__int128)0xd2675e97363d6a3a,
	((__int128)0x3ffb81 << 64) ^ (__int128)0xfba1374915bf7fde,
	((__int128)0x3fffb1 << 64) ^ (__int128)0xa541d1dba472035b,
	((__int128)0x3ffffc << 64) ^ (__int128)0x6d3ec13337a3194e,
	((__int128)0x3fffff << 64) ^ (__int128)0xe427d398013a538c,
	((__int128)0x3fffff << 64) ^ (__int128)0xff6f701a097588d9,
	((__int128)0x3fffff << 64) ^ (__int128)0xfffe0c8c0e024136,
	((__int128)0x3fffff << 64) ^ (__int128)0xfffffb84985259f3,
	((__int128)0x3fffff << 64) ^ (__int128)0xfffffff9282cbe33,
	((__int128)0x3fffff << 64) ^ (__int128)0xfffffffff913b431,
	((__int128)0x3fffff << 64) ^ (__int128)0xffffffffffffffff
};
#define DegG 11
static const uint64_t Coef_G[DegG + 1] = {
	0x00039de6de607651,
	0x00024b6a774a383b,
	0x0001e53d8b8ed795,
	0x001511e36d2e994f,
	0x0019c0ad6a565127,
	0x0000dc13a86bfe81,
	0x000325f3470bd8ff,
	0x000266e436e5a273,
	0x0002ee99958455a1,
	0x00015798ee230265,
	0x000346dc5d638845,
	0x1000000000000000 // 2^{60}
};
static const int32_t expDiff_G[DegG] = { 21, 20, 16, 19, 24, 17, 19, 18, 19, 16, 6 };

// rejection sampling
#define DegR 12
static const uint64_t Coef_R[DegR + 1] = {
	0x024f7f697f05484f,
	0x02dc0a849923ed87,
	0x005029480401f763,
	0x07add9d45d019a15,
	0x05464cbe608c776b,
	0x067058c6b8791005,
	0x01b8280cd17c59f1,
	0x00325f356f827305,
	0x0266e43751e3ac37,
	0x05dd332b0f4ca0a1,
	0x02af31dc461185b3,
	0x01a36e2eb1c432c9,
	0x07ffffffffffffff
};
static const int32_t expDiff_R[DegR] = { 40, 43, 35, 40, 39, 41, 42, 35, 37, 39, 38, 34 };

#define vartheta2 59
#define BitGamma 7 // gamma = ceiling[beta2 / (sigma * alpha * ln2)] = 101 < 2^{7}
#define ThetaChat 20
#define ValChat 0x0000000d89bc6421 // chat = ValChat / 2^{ThetaChat} (Note that ValChat < 2^{36})
#define ThetaEps 52
#define EpsVal 41553 // epsilon = EpsVal / 2^{ThetaEps} (Note that EpsVal < 2^{28})
#endif


#endif