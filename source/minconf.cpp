#include "minconf.hh"

minimal_config makeMinimalConfig(view *view0) 
{

	minimal_config min_conf;

	min_conf.r = (unsigned short)view0->r;
	min_conf.c = (unsigned short)view0->c;

	min_conf.x = view0->x;
	min_conf.y = view0->y;

	min_conf.min_inv_d = (unsigned short)view0->min_inv_d;

	min_conf.n_references = (unsigned char)view0->n_references;
	min_conf.n_depth_references = (unsigned char)view0->n_depth_references;

	min_conf.NNt = (unsigned char)view0->NNt;
	min_conf.Ms = (unsigned char)view0->Ms;

	//min_conf.use_median = (unsigned char)( view0->use_median ? 1 : 0 );
	//min_conf.use_std = view0->stdd > 0 ? 1 : 0;
	//min_conf.yuv_transform = (unsigned char)(view0->yuv_transform ? 1 : 0);

	min_conf.encoding_flags = 0;

	min_conf.encoding_flags = view0->use_median ? min_conf.encoding_flags | (1 << 0) : min_conf.encoding_flags;
	min_conf.encoding_flags = view0->stdd>0 ? min_conf.encoding_flags | (1 << 1) : min_conf.encoding_flags;
	min_conf.encoding_flags = view0->yuv_transform ? min_conf.encoding_flags | (1 << 2) : min_conf.encoding_flags;

	min_conf.encoding_flags = view0->has_color_residual ? min_conf.encoding_flags | (1 << 3) : min_conf.encoding_flags;
	min_conf.encoding_flags = view0->has_depth_residual ? min_conf.encoding_flags | (1 << 4) : min_conf.encoding_flags;

	return min_conf;

}

void setup_form_minimal_config(minimal_config *mconf, view *view0) {

	view0->r = (int)mconf->r;
	view0->c = (int)mconf->c;

	view0->x = mconf->x;
	view0->y = mconf->y;

	view0->min_inv_d = (int)mconf->min_inv_d;

	view0->n_references = (int)mconf->n_references;
	view0->n_depth_references = (int)mconf->n_depth_references;

	view0->NNt = (int)mconf->NNt;
	view0->Ms = (int)mconf->Ms;

	view0->stdd = mconf->encoding_flags & (1 << 0) ? (float)1.0 : (float)0.0;
	view0->use_median = mconf->encoding_flags & (1 << 1) ? 1 : 0;
	view0->yuv_transform = mconf->encoding_flags & (1 << 2) ? 1 : 0;

	view0->has_color_residual = mconf->encoding_flags & (1 << 3) ? 1 : 0;
	view0->has_depth_residual = mconf->encoding_flags & (1 << 4) ? 1 : 0;

	//view0->stdd = (float)mconf->use_std;
	//view0->use_median = mconf->use_median > 0 ? true : false;
	//view0->yuv_transform = mconf->yuv_transform > 0 ? true : false;

}
