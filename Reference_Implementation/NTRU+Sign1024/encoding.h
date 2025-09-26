#ifndef HAETAE_ENCODING_H
#define HAETAE_ENCODING_H

#include "params.h"
#include <stdint.h>
#include <stddef.h>
#include "poly.h"

uint16_t encode_hb_z1(uint8_t *buf, const poly *a);
uint16_t decode_hb_z1(poly *hb_z1, const uint8_t *buf, uint16_t size_in);
uint16_t encode_h(uint8_t *buf, const poly *h);
uint16_t decode_h(poly *h, const uint8_t *buf, uint16_t size_in);
#endif
