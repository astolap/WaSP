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

#ifndef MERGING_HH
#define MERGING_HH

#include "view.hh"

void init_warping_arrays(
    const int32_t N,
    uint16_t **&warped_texture_views,
    uint16_t **&warped_depth_views,
    float **&DispTargs,
    const int32_t nr,
    const int32_t nc,
    const int32_t ncomp);

void clean_warping_arrays(
    const int32_t N,
    uint16_t **warped_texture_views,
    uint16_t **warped_depth_views,
    float **DispTargs);

void setBMask(view *view0);

void initSegVp(
    view *view0, 
    float **DispTargs);

void initViewW(
    view *view0, 
    float **DispTargs);

void mergeMedian_N(
    uint16_t **warpedColorViews, 
    float **DispTargs,  
    view *view0, 
    const int32_t ncomponents);

void mergeWarped_N_icomp(
    uint16_t **warpedColorViews,
    view *view0,
    const int32_t icomp);

void getViewMergingLSWeights_icomp(
    view *view0,
    uint16_t **warpedColorViews,
    const uint16_t *original_color_view,
    const int32_t icomp);

void getGeomWeight_icomp(
    view *view0,
    view *LF,
    const int32_t icomp);


#endif
