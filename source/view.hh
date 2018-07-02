#ifndef VIEW_HH
#define VIEW_HH

#include <stdint.h>

#define MEDFILT_DEPTH false

struct view{

	unsigned short *color;
	unsigned short *depth;

	unsigned short *segmentation;

	int r, c; // SAI subscript

	int nr, nc; // image height, width

	float y, x; // camera displacement

	int min_inv_d; // needed only if inverse depth has negative values, [0,max]-mind = [-mind,max-mind]

	int n_references, n_depth_references;

	int *references, *depth_references; /* depth references not necessarily the same as color references
													  we can have, for example, depth warping only from the externally obtained depth but we still
													  warp color from neighbors that don't have depth provided. We don't want to propagate depth errors from
													  badly warped depth views, thus we restrict depth warping to some high quality subset (usually meaning the 
													  externally obtained high quality depth maps)*/

	signed short *merge_weights;
	int32_t *sparse_weights;

	unsigned char *sparse_mask;

	float *merge_weights_float;

	int *number_of_pixels_per_region;

	bool *bmask; /* view mask for merging weights */
	unsigned short *seg_vp; /* class segmentation, used for view merging weights */
	int NB;

	float residual_rate_color;
	float residual_rate_depth;

	float stdd;

	int NNt, Ms; //for global sparse, NNt defines the neighborhood size [ -NNt:NNt,-NNt:NNt ], Ms is the filter order

	int has_segmentation;
	int maxL; // number of regions in segmentation

	int ****region_displacements; /* region displacement vectors [iar][iac][iR][xy], e.g., [13][13][25][2], for 13x13 angular views with 25 regions for segmentation */

	char path_input_pgm[1024], path_input_ppm[1024], path_input_seg[1024];
	char path_out_pgm[1024], path_out_ppm[1024];

	//char path_input_Y_pgm[1024], path_out_Y_pgm[1024];
	//char path_input_Cb_pgm[1024], path_out_Cb_pgm[1024];
	//char path_input_Cr_pgm[1024], path_out_Cr_pgm[1024];

	float *DM_ROW, *DM_COL; /* for lenslet with region displacement vectors */

	int i_order; /* view position in encoding configuration */

	bool use_median; //use median merging or not

	bool yuv_transform;

};

void initView(view* view);


#endif