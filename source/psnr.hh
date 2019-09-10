/* psnr.hh */
/* Author: Pekka Astola */
/* <pekka.astola@tuni.fi>*/

#ifndef PSNR_HH
#define PSNR_HH

#include <cstdio>


#include <cstdint>

using std::int32_t;
using std::uint32_t;

using std::int16_t;
using std::uint16_t;

using std::int8_t;
using std::uint8_t;

double PSNR(
    uint16_t *im0, 
    uint16_t* im1, 
    const int32_t NR,
    const int32_t NC, 
    const int32_t NCOMP, 
    double maxval);

double PSNR(
    uint16_t *im0, 
    uint16_t* im1,
    const int32_t NR,
    const int32_t NC, 
    const int32_t NCOMP);

double getYCbCr_444_PSNR(
    uint16_t *im0,
    uint16_t* im1,
    const int32_t NR,
    const int32_t NC,
    const int32_t NCOMP,
    const int32_t N);

double getYCbCr_422_PSNR(
    uint16_t *im0, 
    uint16_t* im1, 
    const int32_t NR,
    const int32_t NC, 
    const int32_t NCOMP, 
    const int32_t N);

#endif
