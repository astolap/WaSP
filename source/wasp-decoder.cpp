#include <iostream>
#include <stdio.h>
#include <vector>
#include <cstring>
#include <cstdint>
#include <ctime>
#include <algorithm>

#include "bitdepth.hh"
#include "merging.hh"
#include "minconf.hh"
#include "sparsefilter.hh"
#include "ycbcr.hh"
#include "view.hh"
#include "warping.hh"
#include "residualjp2.hh"
#include "ppm.hh"
#include "fileaux.hh"
#include "medianfilter.hh"
#include "psnr.hh"
#include "inpainting.hh"
#include "predictdepth.hh"
#include "codestream.hh"

#define SAVE_PARTIAL_WARPED_VIEWS false

int main(int argc, char** argv) {

	const char *input_file = argv[1];
	const char *output_dir = argv[2];
	const char *kakadu_dir = argv[3];

	char kdu_expand_path[1024];

	sprintf(kdu_expand_path, "%s%s", kakadu_dir, "/kdu_expand");

	FILE *input_LF;
	input_LF = fopen(input_file, "rb");

	int n_views_total;

	int n_bytes_prediction = 0, n_bytes_residual = 0;

	n_bytes_prediction += (int)fread(&n_views_total, sizeof(int), 1,input_LF)* sizeof(int);

	int _NR, _NC;

	n_bytes_prediction += (int)fread(&_NR, sizeof(int), 1, input_LF)* sizeof(int);
	n_bytes_prediction += (int)fread(&_NC, sizeof(int), 1, input_LF)* sizeof(int);

	bool YUV_TRANSFORM = false;

	int yuv_transform_s;
	n_bytes_prediction += (int)fread(&yuv_transform_s, sizeof(int), 1, input_LF)* sizeof(int);

	YUV_TRANSFORM = yuv_transform_s > 0 ? true : false;

	unsigned short MINIMUM_DEPTH = 0;
	n_bytes_prediction += (int)fread(&MINIMUM_DEPTH, sizeof(unsigned short), 1, input_LF) * sizeof(unsigned short);


	const bool RESIDUAL_16BIT_bool = RESIDUAL_16BIT ? 1 : 0;

	std::vector<std::vector<unsigned char>> JP2_dict;

	view *LF = new view[n_views_total]();

	int ii = 0; /*view index*/

	while ( ii < n_views_total ) {

		view *SAI = LF+ii; 
		ii++;

		initView(SAI);

		SAI->nr = _NR;
		SAI->nc = _NC;

		if (MINIMUM_DEPTH > 0) {
			SAI->min_inv_d = (int)MINIMUM_DEPTH;
		}

		minimal_config mconf;

		codestreamToViewHeader(n_bytes_prediction, SAI, input_LF, mconf);

		if ( feof( input_LF ) ) {
			printf("File reading error. Terminating\t...\n");
			exit(0);
		}

		SAI->color = new unsigned short[SAI->nr*SAI->nc * 3]();
		SAI->depth = new unsigned short[SAI->nr*SAI->nc]();

		predictDepth(SAI, LF);

		/* forward warp color */
		if (SAI->has_color_references) {

			/* holds partial warped views for ii */
			unsigned short **warped_color_views = new unsigned short*[SAI->n_references]();
			unsigned short **warped_depth_views = new unsigned short*[SAI->n_references]();
			float **DispTargs = new float*[SAI->n_references]();

			for (int ij = 0; ij < SAI->n_references; ij++)
			{

				view *ref_view = LF + SAI->references[ij];

				int tmp_w, tmp_r, tmp_ncomp;

				aux_read16PGMPPM(ref_view->path_out_pgm, tmp_w, tmp_r, tmp_ncomp, ref_view->depth);
				aux_read16PGMPPM(ref_view->path_out_ppm, tmp_w, tmp_r, tmp_ncomp, ref_view->color);

				/* FORWARD warp color AND depth */
				warpView0_to_View1(ref_view, SAI, warped_color_views[ij], warped_depth_views[ij], DispTargs[ij]);

				delete[](ref_view->depth);
				delete[](ref_view->color);

				ref_view->depth = NULL;
				ref_view->color = NULL;

				char tmp_str[1024];

				if (SAVE_PARTIAL_WARPED_VIEWS) {

					sprintf(tmp_str, "%s/%03d_%03d%s%03d_%03d%s", output_dir, (ref_view)->c, (ref_view)->r, "_warped_to_", SAI->c, SAI->r, ".ppm");
					aux_write16PGMPPM(tmp_str, SAI->nc, SAI->nr, 3, warped_color_views[ij]);

					sprintf(tmp_str, "%s/%03d_%03d%s%03d_%03d%s", output_dir, (ref_view)->c, (ref_view)->r, "_warped_to_", SAI->c, SAI->r, ".pgm");
					aux_write16PGMPPM(tmp_str, SAI->nc, SAI->nr, 1, warped_depth_views[ij]);

				}

			}

			initViewW(SAI, DispTargs);

			/* Bug fix from VM1.0. The logic for choosing median merging over fixed weight merging was faulty. */
			if (!SAI->use_median) {
				if (SAI->stdd < 0.001) {
					/* merge color with prediction */
					mergeWarped_N(warped_color_views, DispTargs, SAI, 3);
					/* hole filling for color*/
					holefilling(SAI->color, 3, SAI->nr, SAI->nc, 0);
				}
				else {
					/* we don't use LS weights but something derived on geometric distance in view array*/
					getGeomWeight(SAI, LF);
					/* merge color with prediction */
					mergeWarped_N(warped_color_views, DispTargs, SAI, 3);
					/* hole filling for color*/
					holefilling(SAI->color, 3, SAI->nr, SAI->nc, 0);
				}
			}
			else {
				int startt = clock();
				mergeMedian_N(warped_color_views, DispTargs, SAI, 3);
				std::cout << "time elapsed in color median merging\t" << (int)clock() - startt << "\n";
				holefilling(SAI->color, 3, SAI->nr, SAI->nc, 0);
			}
			
			/* clean */
			for (int ij = 0; ij < SAI->n_references; ij++)
			{
				delete[](warped_color_views[ij]);
				delete[](warped_depth_views[ij]);
				delete[](DispTargs[ij]);
			}

			delete[](warped_color_views);
			delete[](warped_depth_views);
			delete[](DispTargs);
		}

		if ( SAI->use_global_sparse )
		{
			applyGlobalSparseFilter(SAI);
		}

		/* prediction part OVER, move on to residual */

		int n_bytes_color_residual = 0, n_bytes_depth_residual = 0;

		/* get residual */
		if ( SAI->has_color_residual )
		{
			//n_bytes_residual += (int)fread(&n_bytes_color_residual, sizeof(int), 1, input_LF)* sizeof(int);
			if (SAI->yuv_transform && YUV_TRANSFORM) {

				char pgm_residual_Y_path[1024];
				char jp2_residual_Y_path_jp2[1024];
				char pgm_residual_Cb_path[1024];
				char jp2_residual_Cb_path_jp2[1024];
				char pgm_residual_Cr_path[1024];
				char jp2_residual_Cr_path_jp2[1024];

				char *ycbcr_pgm_names[3];
				char *ycbcr_jp2_names[3];

				sprintf(pgm_residual_Y_path, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_Y_residual.pgm");
				sprintf(jp2_residual_Y_path_jp2, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_Y_residual.jp2");

				sprintf(pgm_residual_Cb_path, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_Cb_residual.pgm");
				sprintf(jp2_residual_Cb_path_jp2, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_Cb_residual.jp2");

				sprintf(pgm_residual_Cr_path, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_Cr_residual.pgm");
				sprintf(jp2_residual_Cr_path_jp2, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_Cr_residual.jp2");

				ycbcr_pgm_names[0] = pgm_residual_Y_path;
				ycbcr_pgm_names[1] = pgm_residual_Cb_path;
				ycbcr_pgm_names[2] = pgm_residual_Cr_path;

				ycbcr_jp2_names[0] = jp2_residual_Y_path_jp2;
				ycbcr_jp2_names[1] = jp2_residual_Cb_path_jp2;
				ycbcr_jp2_names[2] = jp2_residual_Cr_path_jp2;

				for (int icomp = 0; icomp < 3; icomp++) {

					readResidualFromDisk(ycbcr_jp2_names[icomp], n_bytes_residual, input_LF, JP2_dict);

				}

				int offset_v = 0;
				if (RESIDUAL_16BIT) {
					offset_v = (1 << 15) - 1;// pow(2, 15) - 1;
				}
				else {
					offset_v = (1 << 10) - 1;// pow(2, 10) - 1;
				}

				decodeResidualJP2_YUV(SAI->color, kdu_expand_path, ycbcr_jp2_names, ycbcr_pgm_names, 3, offset_v, (1 << 10) - 1, RESIDUAL_16BIT_bool);

			}
			else {

				int ncomp1 = 0; /* temporary to hold the number of components */

				char ppm_residual_path[1024];

				char jp2_residual_path_jp2[1024];

				sprintf(ppm_residual_path, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_residual.ppm");

				sprintf(jp2_residual_path_jp2, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_residual.jp2");

				readResidualFromDisk(jp2_residual_path_jp2, n_bytes_residual, input_LF, JP2_dict);

				decodeResidualJP2(SAI->color, kdu_expand_path, jp2_residual_path_jp2, ppm_residual_path, ncomp1, (1 << BIT_DEPTH) - 1, (1 << BIT_DEPTH) - 1, RESIDUAL_16BIT_bool);
			}
			
		}


		if (SAI->has_depth_residual) { /* residual depth if needed */

			//n_bytes_residual =+ (int)fread(&n_bytes_depth_residual, sizeof(int), 1, input_LF)* sizeof(int);

			int ncomp1 = 0; /* temporary to hold the number of components */

			char pgm_residual_depth_path[1024];

			char jp2_residual_depth_path_jp2[1024];

			sprintf(pgm_residual_depth_path, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_depth_residual.pgm");

			sprintf(jp2_residual_depth_path_jp2, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_depth_residual.jp2");

			readResidualFromDisk(jp2_residual_depth_path_jp2, n_bytes_residual, input_LF, JP2_dict);

			decodeResidualJP2(SAI->depth, kdu_expand_path, jp2_residual_depth_path_jp2, pgm_residual_depth_path, ncomp1, 0, (1 << 16) - 1,1);

		}

		/* median filter depth */
		if (MEDFILT_DEPTH) {
			unsigned short *tmp_depth = new unsigned short[SAI->nr*SAI->nc]();
			int startt = clock();
			medfilt2D(SAI->depth, tmp_depth, 3, SAI->nr, SAI->nc);
			std::cout << "time elapsed in depth median filtering\t" << (int)clock() - startt << "\n";
			memcpy(SAI->depth, tmp_depth, sizeof(unsigned short)*SAI->nr*SAI->nc);
			delete[](tmp_depth);
		}

		sprintf(SAI->path_out_ppm, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, ".ppm");
		sprintf(SAI->path_out_pgm, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, ".pgm");

		aux_write16PGMPPM(SAI->path_out_ppm, SAI->nc, SAI->nr, 3, SAI->color);
		aux_write16PGMPPM(SAI->path_out_pgm, SAI->nc, SAI->nr, 1, SAI->depth);

		if (SAI->color != NULL) {
			delete[](SAI->color);
			SAI->color = NULL;
		}

		if (SAI->depth != NULL) {
			delete[](SAI->depth);
			SAI->depth = NULL;
		}
		
		if (SAI->seg_vp != NULL) {
			delete[](SAI->seg_vp);
			SAI->seg_vp = NULL;
		}

	}

	fclose(input_LF);

	for (int ii = 0; ii < n_views_total; ii++)
	{

		//printf("ii=%d\n", ii);

		view *SAI = LF + ii;

		if (SAI->color != NULL)
			delete[](SAI->color);
		if (SAI->depth != NULL)
			delete[](SAI->depth);
		if (SAI->references != NULL)
			delete[](SAI->references);
		if (SAI->depth_references != NULL)
			delete[](SAI->depth_references);
		if (SAI->merge_weights != NULL)
			delete[](SAI->merge_weights);
		if (SAI->sparse_weights != NULL)
			delete[](SAI->sparse_weights);
		if (SAI->bmask != NULL)
			delete[](SAI->bmask);
		if (SAI->seg_vp != NULL)
			delete[](SAI->seg_vp);
		if (SAI->sparse_mask != NULL)
			delete[](SAI->sparse_mask);

	}

	delete[](LF);

	
	exit(0);
}