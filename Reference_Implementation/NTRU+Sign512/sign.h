#ifndef SIGN_H
#define SIGN_H

#include "params.h"
extern int ct;
extern int count[21];
extern int count2[21];
int crypto_sign_keypair(unsigned char *pk, unsigned char *sk);

int crypto_sign(unsigned char *sm, unsigned long long *smlen,
            const unsigned char *m, unsigned long long mlen,
            const unsigned char *sk);
int crypto_sign_verify(unsigned char *m, unsigned long long mlen,
                 const unsigned char *sm, unsigned long long smlen,
                 const unsigned char *pk);

#endif