/*BSD 2-Clause License
* Copyright(c) 2019, Pekka Astola
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met :
*
* 1. Redistributions of source code must retain the above copyright notice, this
* list of conditions and the following disclaimer.
*
* 2. Redistributions in binary form must reproduce the above copyright notice,
* this list of conditions and the following disclaimer in the documentation
* and/or other materials provided with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED.IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
*     SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*     CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
*     OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
*     OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

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
