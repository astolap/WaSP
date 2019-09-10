/* ycbcr.hh */
/* Author: Pekka Astola */
/* <pekka.astola@tuni.fi>*/

#ifndef YCBCR_HH
#define YCBCR_HH

#include <cstdint>

using std::int32_t;
using std::uint32_t;

using std::int16_t;
using std::uint16_t;

using std::int8_t;
using std::uint8_t;

#define YUV_422 false

void RGB2YUV422(
    uint16_t *rgb, 
    uint16_t **yy, 
    uint16_t **cbb,
    uint16_t **crr, 
    const int32_t NR,
    const int32_t NC,
    const int32_t NCOMP, 
    const int32_t N);

void RGB2YCbCr(
    const uint16_t *rgb, 
    uint16_t *ycbcr, 
    const int32_t nr,
    const int32_t nc, 
    const int32_t N);

void YCbCr2RGB(
    const uint16_t *ycbcr,
    uint16_t *rgb, 
    const int32_t nr,
    const int32_t nc, 
    const int32_t N);

#endif
