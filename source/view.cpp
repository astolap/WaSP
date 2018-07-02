#include "view.hh"

void initView(view* view)
{

	view->color = NULL;
	view->depth = NULL;
	view->segmentation = NULL;

	view->r = 0;
	view->c = 0;

	view->min_inv_d = 0;
	view->n_references = 0;

	view->references = NULL;
	view->depth_references = NULL;

	view->merge_weights = NULL;
	view->sparse_weights = NULL;
	view->sparse_mask = NULL;
	view->merge_weights_float = NULL;
	view->number_of_pixels_per_region = NULL;

	view->bmask = NULL;
	view->seg_vp = NULL;

	view->NNt = 0;
	view->Ms = 0;

	view->has_segmentation = 0;
	view->maxL = 0;

	view->DM_ROW = NULL;
	view->DM_COL = NULL;

	view->i_order = 0;

	view->stdd = 0.0;

	view->use_median = false;

	view->yuv_transform = true;

}