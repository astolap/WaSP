#ifndef MINCONF_HH
#define MINCONF_HH

#include "view.hh"

struct minimal_config { /* this goes to bitstream */

	unsigned short r , c ; // SAI subscript
	float y, x; // camera displacement

	unsigned short min_inv_d; // needed only if inverse depth has negative values, [0,max]-mind = [-mind,max-mind]

	unsigned char n_references, n_depth_references;

	unsigned char use_std;

	unsigned char NNt, Ms; //for global sparse, NNt defines the neighborhood size [ -NNt:NNt,-NNt:NNt ], Ms is the filter order

	unsigned char use_median; //use median merging or not

	unsigned char yuv_transform;

};

minimal_config makeMinimalConfig(view *view0);
void setup_form_minimal_config(minimal_config *mconf, view *view0);


#endif