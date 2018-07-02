#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <algorithm>

#include <cstdint>

#include <ctime>

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

	size_t f_status;

	f_status = fread(&n_views_total, sizeof(int), 1,input_LF);

	int _NR, _NC;

	f_status = fread(&_NR, sizeof(int), 1, input_LF);
	f_status = fread(&_NC, sizeof(int), 1, input_LF);

	bool YUV_TRANSFORM = false;

	int yuv_transform_s;
	f_status = fread(&yuv_transform_s, sizeof(int), 1, input_LF);

	YUV_TRANSFORM = yuv_transform_s > 0 ? true : false;

	const bool RESIDUAL_16BIT_bool = RESIDUAL_16BIT ? 1 : 0;

	view *LF = new view[n_views_total]();

	int ii = 0; /*view index*/

	while ( ii < n_views_total ) {

		view *SAI = LF+ii; 
		ii++;

		initView(SAI);

		SAI->nr = _NR;
		SAI->nc = _NC;

		minimal_config mconf;
		f_status = fread(&mconf, sizeof(minimal_config), 1, input_LF);

		printf("size of minimal_config %i bytes\n", (int)sizeof(minimal_config));

		setup_form_minimal_config(&mconf, SAI);

		if (SAI->n_references > 0) {
			SAI->references = new int[SAI->n_references]();
			for (int ij = 0; ij < SAI->n_references; ij++) {
				unsigned short nid;
				f_status = fread(&nid, sizeof(unsigned short), 1, input_LF);
				*(SAI->references + ij) = (int)nid;
			}
		}

		if (SAI->n_depth_references > 0) {
			SAI->depth_references = new int[SAI->n_depth_references]();
			for (int ij = 0; ij < SAI->n_depth_references; ij++) {
				unsigned short nid;
				f_status = fread(&nid, sizeof(unsigned short), 1, input_LF);
				*(SAI->depth_references + ij) = (int)nid;
			}
		}

		if (SAI->Ms > 0 && SAI->NNt > 0) {
			SAI->sparse_mask = new unsigned char[SAI->Ms]();
			f_status = fread(SAI->sparse_mask, sizeof(unsigned char), SAI->Ms, input_LF);
			SAI->sparse_weights = new int32_t[SAI->Ms]();
			f_status = fread(SAI->sparse_weights, sizeof(int32_t), SAI->Ms, input_LF);
		}

		SAI->NB = (1 << SAI->n_references)*SAI->n_references;// (pow(2, SAI->n_references)*SAI->n_references);

		if ( !SAI->use_median ) {
			if (SAI->stdd < 0.001) {
				if (SAI->n_references > 0) {
					SAI->merge_weights = new signed short[SAI->NB / 2]();
					f_status = fread(SAI->merge_weights, sizeof(signed short), SAI->NB / 2, input_LF);
				}
			}
			else {
				f_status = fread(&SAI->stdd, sizeof(float), 1, input_LF);
			}
		}

		if (!(f_status > 0)) {
			printf("File reading error. Terminating\t...\n");
			exit(0);
		}

		SAI->color = new unsigned short[SAI->nr*SAI->nc * 3]();
		SAI->depth = new unsigned short[SAI->nr*SAI->nc]();

		predictDepth(SAI, LF);

		/* forward warp color */
		if (SAI->n_references > 0) {

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

				//int uu = SAI->references[ij];

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

			if (SAI->stdd < 0.01) {
				/* merge color with prediction */
				mergeWarped_N(warped_color_views, DispTargs, SAI, 3);
				/* hole filling for color*/
				holefilling(SAI->color, 3, SAI->nr, SAI->nc, 0);
			}
			else {
				if (SAI->use_median) {
					int startt = clock();
					mergeMedian_N(warped_color_views, DispTargs, SAI, 3);
					std::cout << "time elapsed in color median merging\t" << (int)clock() - startt << "\n";
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

		if (SAI->NNt > 0 && SAI->Ms > 0)
		{
			applyGlobalSparseFilter(SAI);
		}

		/* prediction part OVER, move on to residual */

		int n_bytes_color_residual = 0, n_bytes_depth_residual = 0;

		f_status = fread(&n_bytes_color_residual, sizeof(int), 1,  input_LF);

		/* get residual */
		if (n_bytes_color_residual > 0)
		{

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

					if (icomp > 0) {
						f_status = fread(&n_bytes_color_residual, sizeof(int), 1, input_LF);
					}

					unsigned char *jp2_residual = new unsigned char[n_bytes_color_residual]();
					f_status = fread(jp2_residual, sizeof(unsigned char), n_bytes_color_residual, input_LF);

					FILE *jp2_res_file;
					jp2_res_file = fopen(ycbcr_jp2_names[icomp], "wb");
					fwrite(jp2_residual, sizeof(unsigned char), n_bytes_color_residual, jp2_res_file);
					fclose(jp2_res_file);

					delete[](jp2_residual);

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

				unsigned char *jp2_residual = new unsigned char[n_bytes_color_residual]();
				f_status = fread(jp2_residual, sizeof(unsigned char), n_bytes_color_residual, input_LF);

				FILE *jp2_res_file;
				jp2_res_file = fopen(jp2_residual_path_jp2, "wb");
				fwrite(jp2_residual, sizeof(unsigned char), n_bytes_color_residual, jp2_res_file);
				fclose(jp2_res_file);

				delete[](jp2_residual);

				decodeResidualJP2(SAI->color, kdu_expand_path, jp2_residual_path_jp2, ppm_residual_path, ncomp1, (1 << BIT_DEPTH) - 1, (1 << BIT_DEPTH) - 1, RESIDUAL_16BIT_bool);
			}
			
		}

		f_status = fread(&n_bytes_depth_residual, sizeof(int), 1, input_LF);

		if (n_bytes_depth_residual>0) { /* residual depth if needed */

			int ncomp1 = 0; /* temporary to hold the number of components */

			char pgm_residual_depth_path[1024];

			char jp2_residual_depth_path_jp2[1024];

			sprintf(pgm_residual_depth_path, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_depth_residual.pgm");

			sprintf(jp2_residual_depth_path_jp2, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_depth_residual.jp2");

			unsigned char *jp2_depth_residual = new unsigned char[n_bytes_depth_residual]();
			f_status = fread(jp2_depth_residual, sizeof(unsigned char), n_bytes_depth_residual, input_LF);

			FILE *jp2_depth_res_file;
			jp2_depth_res_file = fopen(jp2_residual_depth_path_jp2, "wb");
			fwrite(jp2_depth_residual, sizeof(unsigned char), n_bytes_depth_residual, jp2_depth_res_file);
			fclose(jp2_depth_res_file);

			delete[](jp2_depth_residual);

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