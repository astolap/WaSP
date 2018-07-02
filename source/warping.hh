#ifndef WARPING_HH
#define WARPING_HH

#include "view.hh"

void warpView0_to_View1(view *view0, view *view1, unsigned short *&warpedColor, unsigned short *&warpedDepth, float *&DispTarg);

#endif