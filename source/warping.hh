/* warping.hh */
/* Author: Pekka Astola */
/* <pekka.astola@tuni.fi>*/

#ifndef WARPING_HH
#define WARPING_HH

#include <cstdint>

using std::int32_t;
using std::uint32_t;

using std::int16_t;
using std::uint16_t;

using std::int8_t;
using std::uint8_t;

#include "view.hh"

#define INIT_DISPARITY_VALUE -10000.0 /* Bug fix from VM1.0. Introduced because -1.0 is no longer a good initial value since we can have negative disparity as well.*/

void warpView0_to_View1(
    view *view0, 
    view *view1, 
    uint16_t *warpedColor,
    uint16_t *warpedDepth, 
    float *DispTarg);

#endif
