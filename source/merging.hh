/* merging.hh */
/* Author: Pekka Astola */
/* <pekka.astola@tuni.fi>*/

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
