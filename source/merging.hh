#ifndef MERGING_HH
#define MERGING_HH

#include "view.hh"

void setBMask(view *view0);

void getViewMergingLSWeights_N(view *view0, unsigned short **warpedColorViews, float **DispTargs, const unsigned short *original_color_view);

void initSegVp(view *view0, float **DispTargs);
void initViewW(view *view0, float **DispTargs);

void getGeomWeight(view *view0, view *LF);
void mergeMedian_N(unsigned short **warpedColorViews, float **DispTargs, view *view0, const int ncomponents);
void mergeWarped_N(unsigned short **warpedColorViews, float **DispTargs, view *view0, const int ncomponents);

#endif