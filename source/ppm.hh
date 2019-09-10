/* ppm.hh */
/* Author: Pekka Astola */
/* <pekka.astola@tuni.fi>*/

#ifndef PPM_HH
#define PPM_HH

#include <cstdint>

using std::int32_t;
using std::uint32_t;

using std::int16_t;
using std::uint16_t;

using std::int8_t;
using std::uint8_t;

bool aux_read16PGMPPM(
    const char* filename, 
    int32_t &width, 
    int32_t &height, 
    int32_t &ncomp,
    uint16_t *&img);

bool aux_write16PGMPPM(
    const char* filename, 
    const int32_t width, 
    const int32_t height,
    const int32_t ncomp, 
    const uint16_t *img);

#endif
