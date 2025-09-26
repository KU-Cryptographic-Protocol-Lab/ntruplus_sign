#include "rans_byte.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define SCALE_BITS 10
#define SCALE (1u << SCALE_BITS)

#define PARAM_1

/*#ifdef PARAM_1 //2d
#define M_H 7 // Alphabet size for HighBits(z1)
static uint32_t f_h[M_H]         = {1, 26, 239, 492, 239, 26, 1};
#define M_HB_Z1 7 // Alphabet size for h
static uint32_t f_hb_z1[M_HB_Z1] = {1, 19, 230, 524, 230, 19, 1};
#endif*/
#ifdef PARAM_1 //2d
#define M_H 13 // Alphabet size for HighBits(z1)
static uint32_t f_h[M_H]         = {1, 1, 7, 37, 112, 218, 272, 218, 112, 37, 7, 1, 1};
#define M_HB_Z1 13 // Alphabet size for h
static uint32_t f_hb_z1[M_HB_Z1] = {1, 1, 6, 35, 110, 220, 278, 220, 110, 35, 6, 1, 1};
#endif

#ifdef PARAM_2 //2d
#define M_HB_Z1 9 // Alphabet size for h
static uint32_t f_hb_z1[M_HB_Z1]  = {1, 1, 27, 239, 488, 239, 27, 1, 1};
#define M_H 9
static uint32_t f_h[M_H] = {1, 2, 35, 244, 460, 244, 35, 2, 1};
#endif

#ifdef PARAM_3 //d=9
#define M_H 7 // Alphabet size for h
static uint32_t f_h[M_H]         = {1, 6, 206, 598, 206, 6, 1};
#define M_HB_Z1 7
static uint32_t f_hb_z1[M_HB_Z1] = {1, 2, 177, 664, 177, 2, 1};
#endif

static RansEncSymbol esyms_h[M_H];
static uint16_t symbol_h[SCALE] = {0};
static RansDecSymbol dsyms_h[M_H];

static RansEncSymbol esyms_hb_z1[M_HB_Z1];
static uint16_t symbol_hb_z1[SCALE] = {0};
static RansDecSymbol dsyms_hb_z1[M_HB_Z1];

void print_array(uint32_t* array, size_t length) {
    for(size_t i=0; i< length; i++) {
        printf("%u, ",array[i]);
    }
}

void compute_freqs(uint32_t *freqs, const uint8_t *samples, int sample_size, int m) {
    double tmp;
    int i_max = 0;
    uint32_t freq_max = 0;
    int sum = 0;
    int diff;
    
    // init freqs with 0:
    for(int i=0; i < m; i++)
        freqs[i] = 0;
    // count symbols:
    for(int i=0; i < sample_size; i++) {
        if(samples[i] >= m)
            printf("ERROR: symbol out of range!\n");
        freqs[samples[i]]++;
    }
    // normalize counts:
    for(int i=0; i < m; i++) {
        tmp = freqs[i] / (double) sample_size;
        tmp = tmp * (double) SCALE;
        freqs[i] = (uint32_t) tmp;
        // freq must be at least one:
        if(freqs[i] == 0)
            freqs[i] = 1;
        // track most frequent symbol:
        if(freqs[i] > freq_max) {
            i_max = i;
            freq_max = freqs[i];
        }        
        // compute sum:
        sum += freqs[i];
    }
    // ensure sum == SCALE
    diff = sum - SCALE;
    freqs[i_max] -= diff;

    printf("Freq table with %u symbols, scale %u, diff was %i: {", m, SCALE, diff);
    print_array(freqs, m);
    printf("};\n");
}

void symbol_table(uint16_t *symbol, const uint32_t *freq, size_t alphabet_size) 
{
    int pos = 0;
    for (size_t sym=0; sym < alphabet_size; sym++) {
        for (uint32_t i=0; i < freq[sym]; i++)
            symbol[pos++] = sym;
    }
}

void cum_freq_table(uint32_t *cum_freq, const uint32_t *freq, size_t alphabet_size) {
    cum_freq[0] = 0;
    for (size_t i = 1; i < alphabet_size; i++) {
      cum_freq[i] = cum_freq[i-1] + freq[i-1];
    }
}

void encode_symbols(RansEncSymbol *esyms, const uint32_t *freq, size_t alphabet_size) {
    uint32_t cum_freq[alphabet_size];
    cum_freq_table(cum_freq, freq, alphabet_size);

    for (size_t i=0; i < alphabet_size; i++) 
    {
        RansEncSymbolInit(&esyms[i], cum_freq[i], freq[i], SCALE_BITS);
    }
}

void decode_symbols(RansDecSymbol *dsyms, uint16_t *symbol, const uint32_t *freq, size_t alphabet_size) {
    uint32_t cum_freq[alphabet_size];
    cum_freq_table(cum_freq, freq, alphabet_size);

    symbol_table(symbol, freq, alphabet_size);
    
    for (size_t i=0; i < alphabet_size; i++) 
    {
        RansDecSymbolInit(&dsyms[i], cum_freq[i], freq[i]);
    }
}

void precomputations_rans() {
    //compute_freqs(f_h, samples_h, SAMPLE_SIZE_H, M_H);
    encode_symbols(esyms_h, f_h, M_H);
    decode_symbols(dsyms_h, symbol_h, f_h, M_H);

    //compute_freqs(f_hb_z1, samples_z, SAMPLE_SIZE_Z, M_HB_Z1);
    encode_symbols(esyms_hb_z1, f_hb_z1, M_HB_Z1);
    decode_symbols(dsyms_hb_z1, symbol_hb_z1, f_hb_z1, M_HB_Z1);
}

void print_esym(RansEncSymbol x){
    printf("{%i, %i, %i, %i, %i}", x.x_max, x.rcp_freq, x.bias, x.cmpl_freq, x.rcp_shift);
}

void print_dsym(RansDecSymbol x){
    printf("{%i, %i}", x.start, x.freq);
}

int main() {
    precomputations_rans();

    printf("\n");

    // encoding h
    printf("static RansEncSymbol esyms_h[M_H] = {");
    for(int i=0; i<M_H; i++) {
        print_esym(esyms_h[i]);
        printf(", ");
    }
    printf("};\n");

    printf("static RansDecSymbol dsyms_h[M_H] = {");
    for(int i=0; i<M_H; i++) {
        print_dsym(dsyms_h[i]);
        printf(", ");
    }
    printf("};\n");

    printf("static uint16_t symbol_h[SCALE] = {");
    for(int i=0; i< (int)SCALE; i++) {
        printf("%i, ", symbol_h[i]);
    }
    printf("};\n");

    // encoding HB(z1)
    printf("static RansEncSymbol esyms_hb_z1[M_HB_Z1] = {");
    for(int i=0; i<M_HB_Z1; i++) {
        print_esym(esyms_hb_z1[i]);
        printf(", ");
    }
    printf("};\n");

    printf("static RansDecSymbol dsyms_hb_z1[M_HB_Z1] = {");
    for(int i=0; i<M_HB_Z1; i++) {
        print_dsym(dsyms_hb_z1[i]);
        printf(", ");
    }
    printf("};\n");

    printf("static uint16_t symbol_hb_z1[SCALE] = {");
    for(int i=0; i< (int)SCALE; i++) {
        printf("%i, ", symbol_hb_z1[i]);
    }
    printf("};\n");

    printf("\n");
}
