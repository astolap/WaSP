#include "minconf.hh"

minimal_config makeMinimalConfig(view *view0) 
{

	minimal_config min_conf;

	min_conf.r = (unsigned char)view0->r;
	min_conf.c = (unsigned char)view0->c;

	min_conf.encoding_flags = 0;

	min_conf.encoding_flags = view0->use_median ? min_conf.encoding_flags | (1 << 0) : min_conf.encoding_flags;
	min_conf.encoding_flags = view0->stdd>0.001 ? min_conf.encoding_flags | (1 << 1) : min_conf.encoding_flags;
	min_conf.encoding_flags = view0->yuv_transform ? min_conf.encoding_flags | (1 << 2) : min_conf.encoding_flags;

	min_conf.encoding_flags = view0->has_color_residual ? min_conf.encoding_flags | (1 << 3) : min_conf.encoding_flags;
	min_conf.encoding_flags = view0->has_depth_residual ? min_conf.encoding_flags | (1 << 4) : min_conf.encoding_flags;

	min_conf.encoding_flags = view0->use_global_sparse ? min_conf.encoding_flags | (1 << 5) : min_conf.encoding_flags;

	min_conf.encoding_flags = view0->has_color_references ? min_conf.encoding_flags | (1 << 6) : min_conf.encoding_flags;
	min_conf.encoding_flags = view0->has_depth_references ? min_conf.encoding_flags | (1 << 7) : min_conf.encoding_flags;

	min_conf.encoding_flags = view0->has_x_displacement ? min_conf.encoding_flags | (1 << 8) : min_conf.encoding_flags;
	min_conf.encoding_flags = view0->has_y_displacement ? min_conf.encoding_flags | (1 << 9) : min_conf.encoding_flags;

	//min_conf.encoding_flags = view0->has_min_inv_depth ? min_conf.encoding_flags | (1 << 8) : min_conf.encoding_flags;

	return min_conf;

}

void setup_form_minimal_config(minimal_config *mconf, view *view0) {

	view0->r = (int)mconf->r;
	view0->c = (int)mconf->c;

	view0->stdd = ( mconf->encoding_flags & (1 << 0) )>0 ? (float)1.0 : (float)0.0;
	view0->use_median = (mconf->encoding_flags & (1 << 1))>0 ? 1 : 0;
	view0->yuv_transform = (mconf->encoding_flags & (1 << 2))>0 ? 1 : 0;

	view0->has_color_residual = (mconf->encoding_flags & (1 << 3))>0 ? 1 : 0;
	view0->has_depth_residual = (mconf->encoding_flags & (1 << 4))>0 ? 1 : 0;

	view0->use_global_sparse = (mconf->encoding_flags & (1 << 5))>0 ? 1 : 0;

	view0->has_color_references = (mconf->encoding_flags & (1 << 6))>0 ? 1 : 0;
	view0->has_depth_references = (mconf->encoding_flags & (1 << 7))>0 ? 1 : 0;

	view0->has_x_displacement = (mconf->encoding_flags & (1 << 8))>0 ? 1 : 0;
	view0->has_y_displacement = (mconf->encoding_flags & (1 << 9))>0 ? 1 : 0;

	//view0->has_min_inv_depth = (mconf->encoding_flags & (1 << 8))>0 ? 1 : 0;
}
