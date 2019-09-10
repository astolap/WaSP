/* fastols.hh */
/* Author: Pekka Astola */
/* <pekka.astola@tuni.fi>*/

#ifndef FASTOLS_HH
#define FASTOLS_HH

#include <cstdint>

using std::int32_t;
using std::uint32_t;

using std::int16_t;
using std::uint16_t;

using std::int8_t;
using std::uint8_t;

int32_t FastOLS_new(
    double **AAA, 
    double **Ydd, 
    int32_t *PredRegr0, 
    double *PredTheta0,
    const int32_t Ms, 
    const int32_t MT, 
    const int32_t MPHI, 
    const int32_t N);

#endif
