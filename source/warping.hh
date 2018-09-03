#ifndef WARPING_HH
#define WARPING_HH

#include "view.hh"

#define INIT_DISPARITY_VALUE -10000.0 /* Bug fix from VM1.0. Introduced because -1.0 is no longer a good initial value since we can have negative disparity as well.*/

void warpView0_to_View1(view *view0, view *view1, unsigned short *&warpedColor, unsigned short *&warpedDepth, float *&DispTarg);

#endif