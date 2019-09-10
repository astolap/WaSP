/* warping.cpp */
/* Author: Pekka Astola */
/* <pekka.astola@tuni.fi>*/

#include "warping.hh"
#include "bitdepth.hh"

#include <cstdlib>
#include <cmath>
#include <cstring>

void warpView0_to_View1(
    view *view0, 
    view *view1, 
    uint16_t *warpedColor,
    uint16_t *warpedDepth, 
    float *DispTarg) {

  /*this function forward warps from view0 to view1 for both color and depth*/

  float ddy = view0->y - view1->y;
  float ddx = view0->x - view1->x;

  uint16_t *AA1 = view0->color;
  uint16_t *DD1 = view0->depth;

  memset(warpedColor, 0, sizeof(uint16_t)*view0->nr*view0->nc * 3);
  memset(warpedDepth, 0, sizeof(uint16_t)*view0->nr*view0->nc);
  //memset(DispTarg, 0, sizeof(float)*view0->nr*view0->nc);

  for (int32_t ij = 0; ij < view0->nr * view0->nc; ij++) {
    DispTarg[ij] = INIT_DISPARITY_VALUE;
  }

  for (int32_t ij = 0; ij < view0->nr * view0->nc; ij++) {

    float disp = 
        (static_cast<float>(DD1[ij]) - static_cast<float>(view0->min_inv_d))
        / static_cast<float>(1 << D_DEPTH);

    float DM_COL = disp * ddx;
    float DM_ROW = disp * ddy;

    int32_t iy = ij % view0->nr;  //row
    int32_t ix = (ij - iy) / view0->nr;  //col

    int32_t ixnew = ix + static_cast<int32_t>( floor(DM_COL + 0.5f) );
    int32_t iynew = iy + static_cast<int32_t>( floor(DM_ROW + 0.5f) );

    if (
        iynew >= 0 &&  
        ixnew >= 0 && 
        ixnew < view0->nc && 
        iynew < view0->nr) {

      int32_t indnew = iynew + ixnew * view0->nr;

      if (DispTarg[indnew] < disp) {

        DispTarg[indnew] = disp;
        warpedDepth[indnew] = DD1[ij];

        for (int32_t icomp = 0; icomp < view0->ncomp; icomp++) {

            warpedColor[indnew + view0->nr * view0->nc*icomp] = 
                AA1[ij+ view0->nr * view0->nc*icomp];

        }
      }
    }
  }
}
