#ifndef SPARSEFILTER_HH
#define SPARSEFILTER_HH

#include "view.hh"

void applyGlobalSparseFilter(view *view0);
void getGlobalSparseFilter(view *view0, unsigned short *original_color_view);

#endif