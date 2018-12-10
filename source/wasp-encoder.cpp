#include <cstdio>
#include <cstdint>
#include <ctime>
#include <cstring>
#include <iostream>
#include <vector>
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


#define USE_difftest_ng false

#define YUV_SEARCH_LOW 6.80f
#define YUV_SEARCH_HIGH 8.0f
#define YUV_SEARCH_STEP 0.10f
#define YUV_RATIO_DEFAULT 7.20f
#define MAX_Y_RATIO 0.99

#define STD_SEARCH_LOW 10
#define STD_SEARCH_HIGH 250
#define STD_SEARCH_STEP 10

#define SAVE_PARTIAL_WARPED_VIEWS false

#ifndef FLT_MAX
#define FLT_MAX          3.402823466e+38F        // max value
#endif

int main(int argc, char** argv) {

	const char *input_dir = argv[1];
	const char *output_dir = argv[2];
	const char *kakadu_dir = argv[3];
	const char *config_file = argv[4];

	char kdu_compress_path[1024];
	char kdu_expand_path[1024];

	sprintf(kdu_compress_path, "%s%s", kakadu_dir, "/kdu_compress");
	sprintf(kdu_expand_path, "%s%s", kakadu_dir, "/kdu_expand");

	//const char *difftest_call = "C:/Local/astolap/Data/JPEG_PLENO/RIO_INPUT/ScriptJan2018/ScriptSolution/difftest_ng.exe --toycbcr --psnr ";
	//const char *difftest_call_pgm = "C:/Local/astolap/Data/JPEG_PLENO/RIO_INPUT/ScriptJan2018/ScriptSolution/difftest_ng.exe --psnr ";

	//const char *difftest_call = "D:/JPEG_VM_01/DevelopmentC/Pekka/difftest_ng.exe --toycbcr --psnr ";
	//const char *difftest_call_pgm = "D:/JPEG_VM_01/DevelopmentC/Pekka/difftest_ng.exe --psnr ";

	FILE *filept;

	filept = fopen(config_file, "rb");

	int n_views_total;
	fread(&n_views_total, sizeof(int), 1, filept); /*reading*/

	view *LF = new view[n_views_total](); /* one dimensional view vector */

	int maxR = 0, maxC = 0;

	bool YUV_TRANSFORM = false;
	bool YUV_RATIO_SEARCH = false;
	bool STD_SEARCH = false;

	const bool RESIDUAL_16BIT_bool = RESIDUAL_16BIT ? 1 : 0;

	int yuv_transform_s,yuv_ratio_search_s,std_search_s;
	fread(&yuv_transform_s, sizeof(int), 1, filept); /*reading*/
	fread(&yuv_ratio_search_s, sizeof(int), 1, filept); /*reading*/
	fread(&std_search_s, sizeof(int), 1, filept); /*reading*/

	YUV_TRANSFORM = yuv_transform_s > 0 ? true : false;
	YUV_RATIO_SEARCH = yuv_ratio_search_s > 0 ? true : false;
	STD_SEARCH = std_search_s > 0 ? true : false;

	unsigned short MINIMUM_DEPTH = 0;

	for (int ii = 0; ii < n_views_total; ii++) {

		view *SAI = LF + ii;

		initView(SAI);

		SAI->i_order = ii;

		fread(&(SAI->r), sizeof(int), 1, filept); /*reading*/
		fread(&(SAI->c), sizeof(int), 1, filept); /*reading*/

		/* find number of rows and cols */
		maxR = (SAI->r+1 > maxR) ? SAI->r + 1 : maxR;
		maxC = (SAI->c+1 > maxC) ? SAI->c + 1 : maxC;

		int xx = 0, yy = 0;

		fread(&xx, sizeof(int), 1, filept); /*reading*/
		fread(&yy, sizeof(int), 1, filept); /*reading*/

		if ( abs(xx) > 0) {
			SAI->has_x_displacement = true;
			SAI->x = float(xx) / 100000;
		}

		if ( abs(yy) > 0) {
			SAI->has_y_displacement = true;
			SAI->y = float(yy) / 100000;
		}

		int rate_color, rate_depth;

		fread(&rate_color, sizeof(int), 1, filept); /*reading*/
		fread(&rate_depth, sizeof(int), 1, filept); /*reading*/

		SAI->residual_rate_color = ((float)rate_color) / 100000;
		SAI->residual_rate_depth = ((float)rate_depth) / 100000;

		fread(&SAI->Ms, sizeof(int), 1, filept); /*reading*/
		fread(&SAI->NNt, sizeof(int), 1, filept); /*reading*/

		int stdd = 0;

		fread(&stdd, sizeof(int), 1, filept); /*reading*/
		SAI->stdd = ((float)stdd) / 100000;

		unsigned short tmpminv;

		fread(&tmpminv, sizeof(int), 1, filept); /*reading, if we have negative inverse depth,
												   for example in lenslet, we need to subtract min_inv_d
												   from the inverse depth maps*/
		if (ii == 0) {
			if (tmpminv > 0) {
				MINIMUM_DEPTH = tmpminv;
			}
		}

		if (MINIMUM_DEPTH > 0) {
			//SAI->has_min_inv_depth = true;
			SAI->min_inv_d = (int)MINIMUM_DEPTH;
		}

		fread(&(SAI->n_references), sizeof(int), 1, filept); /*reading*/

		if (SAI->n_references > 0) {

			SAI->has_color_references = true;

			SAI->references = new int[SAI->n_references]();

			fread(SAI->references, sizeof(int), SAI->n_references, filept); /*reading*/

		}

		fread(&(SAI->n_depth_references), sizeof(int), 1, filept); /*reading*/

		if (SAI->n_depth_references > 0) {

			SAI->has_depth_references = true;

			SAI->depth_references = new int[SAI->n_depth_references]();

			fread(SAI->depth_references, sizeof(int), SAI->n_depth_references, filept); /*reading*/

		}

		fread(&SAI->has_segmentation, sizeof(int), 1, filept); // new,13.06.18 /*reading*/

		sprintf(SAI->path_input_ppm, "%s%c%03d_%03d%s", input_dir, '/', SAI->c, SAI->r, ".ppm");
		sprintf(SAI->path_input_pgm, "%s%c%03d_%03d%s", input_dir, '/', SAI->c, SAI->r, ".pgm");

		sprintf(SAI->path_input_seg, "%s%c%03d_%03d%s", input_dir, '/', SAI->c, SAI->r, "_segmentation.pgm");

		sprintf(SAI->path_out_ppm, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, ".ppm");
		sprintf(SAI->path_out_pgm, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, ".pgm");

	}
	fclose(filept);

	char path_out_LF_data[1024];
	sprintf(path_out_LF_data, "%s%c%s", output_dir, '/', "output.LF");

	/* debugging codestream overhead reduction */
	char path_codestream[1024];
	sprintf(path_codestream, "%s%c%s", output_dir, '/', "LF.codestream");
	FILE *tmp_codestream;
	tmp_codestream = fopen(path_codestream, "wb");
	fclose(tmp_codestream);

	/* our bitstream starts here */
	FILE *output_LF_file;
	output_LF_file = fopen(path_out_LF_data, "wb");
	fwrite(&n_views_total, sizeof(int), 1, output_LF_file);
	fclose(output_LF_file);

	bool global_header_written = false;

	FILE *output_results_file;
	char output_results_filename[1024];
	sprintf(output_results_filename, "%s/%s", output_dir, "results.txt");

	char output_results[2048];
	int output_buffer_length = 0;

	//output_buffer_length += sprintf(output_results + output_buffer_length, "%s",
	//	"ROW\tCOL\tPSNR1\tPSNR2\tPSNR3\tPSNR4\tSTD\tYUVRATIO\tPSNR5\t\tbytes_prediction\t\tbytes_residual");
	output_results_file = fopen(output_results_filename, "w");
	//fprintf(output_results_file, "%s\n", output_results);
	fclose(output_results_file);


	/* to get effiency from multiple JP2 files, we remove parts of the files
	which are repetative over all files. For this we have a minimalistic 
	dictionary method. */

	std::vector<std::vector<unsigned char>> JP2_dict;

	double psnr_yuv_mean = 0;

	for (int ii = 0; ii < n_views_total; ii++) {

		view *SAI = LF + ii;

		printf("Encoding view %03d_%03d\t", SAI->c, SAI->r);

		//char output_results[1024];
		memset(output_results, 0x00, sizeof(char) * sizeof(output_results)/sizeof(char));
		output_buffer_length = 0;
		output_buffer_length += sprintf(output_results + output_buffer_length, "%03d\t%03d", SAI->r, SAI->c);

		unsigned short *original_color_view = NULL;
		unsigned short *original_depth_view = NULL;

		int nc1, nr1, ncomp1;
		aux_read16PGMPPM(SAI->path_input_ppm, SAI->nc, SAI->nr, ncomp1, original_color_view);

		///* debug yuv */
		//unsigned short *ycbcr = new unsigned short[SAI->nr*SAI->nc * 3]();
		//RGB2YCbCr(original_color_view, ycbcr, SAI->nr, SAI->nc, 10);
		//aux_write16PGMPPM("C:/Local/astolap/Data/JPEG_PLENO/TUT-HDCA_tmp_output_lenslet/YCbCr.ppm",
		//	SAI->nc, SAI->nr, 3, ycbcr);
		//unsigned short *rgb = new unsigned short[SAI->nr*SAI->nc * 3]();
		//YCbCr2RGB(ycbcr, rgb, SAI->nr, SAI->nc, 10);
		//aux_write16PGMPPM("C:/Local/astolap/Data/JPEG_PLENO/TUT-HDCA_tmp_output_lenslet/RGB.ppm",
		//	SAI->nc, SAI->nr, 3, rgb);
		////aux_write16PGMPPM("C:/Local/astolap/Data/JPEG_PLENO/TUT-HDCA_tmp_output_lenslet/original.ppm",
		////	SAI->nc, SAI->nr, 3, original_color_view);

		//exit(0);

		bool depth_file_exist = false;
		
		if (SAI->residual_rate_depth > 0) {
			depth_file_exist = aux_read16PGMPPM(SAI->path_input_pgm, nc1, nr1, ncomp1, original_depth_view);
		}

		SAI->color = new unsigned short[SAI->nr*SAI->nc * 3]();
		SAI->depth = new unsigned short[SAI->nr*SAI->nc]();

		predictDepth(SAI,LF);

		/* color prediction */
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

				//FILE *tmpf;
				//sprintf(tmp_str, "%s%03d_%03d%s%03d_%03d%s", output_dir, (ref_view)->c, (ref_view)->r, "_warped_to_", SAI->c, SAI->r, "_DispTarg.float");
				//tmpf = fopen(tmp_str, "wb");
				//fwrite(DispTargs[ij], sizeof(float), SAI->nr * SAI->nc, tmpf);
				//fclose(tmpf);

			}

			initViewW(SAI, DispTargs);

			/* get LS weights */
			if (SAI->stdd < 0.0001) {
				getViewMergingLSWeights_N(SAI, warped_color_views, DispTargs, original_color_view);
				/* merge color with prediction */
				mergeWarped_N(warped_color_views, DispTargs, SAI, 3);
				/* hole filling for color*/
				holefilling(SAI->color, 3, SAI->nr, SAI->nc, 0);
			}
			else {
				/* get baseline with median, we then study whether weighting improves */
				/* merge color with median */
				int startt = clock();
				mergeMedian_N(warped_color_views, DispTargs, SAI, 3);
				std::cout << "time elapsed in color median merging\t" << (float)( (int)clock() - startt ) / CLOCKS_PER_SEC << "\n";
				holefilling(SAI->color, 3, SAI->nr, SAI->nc, 0);

				//double psnr_med = getYCbCr_422_PSNR(SAI->color, original_color_view, SAI->nr, SAI->nc, 3, 10);
				double psnr_med = getYCbCr_422_PSNR(SAI->color, original_color_view, SAI->nr, SAI->nc, 3, 10);

				unsigned short *tmp_m = new unsigned short[SAI->nr*SAI->nc*3]();
				memcpy(tmp_m, SAI->color, sizeof(unsigned short)*SAI->nr*SAI->nc * 3);

				double psnr_w = 0;

				if (STD_SEARCH) {

					float stdi = 0;

					for (float stds = STD_SEARCH_LOW; stds < STD_SEARCH_HIGH; stds += STD_SEARCH_STEP) {
						SAI->stdd = stds;
						/* we don't use LS weights but something derived on geometric distance in view array*/
						getGeomWeight(SAI, LF);
						/* merge color with prediction */
						mergeWarped_N(warped_color_views, DispTargs, SAI, 3);
						/* hole filling for color*/
						holefilling(SAI->color, 3, SAI->nr, SAI->nc, 0);
						double tpsnr = getYCbCr_422_PSNR(SAI->color, original_color_view, SAI->nr, SAI->nc, 3, BIT_DEPTH);
						//double tpsnr = getYCbCr_422_PSNR(SAI->color, original_color_view, SAI->nr, SAI->nc, 3, 10);

						if (tpsnr > psnr_w) {
							psnr_w = tpsnr;
							stdi = SAI->stdd;
						}

					}

					SAI->stdd = stdi;

					printf("PSNR RGB median:%f\tPSNR RGB weights:%f\t SAI->std: %f\n", psnr_med, psnr_w, SAI->stdd);

				}

				/* we don't use LS weights but something derived on geometric distance in view array*/
				getGeomWeight(SAI, LF);
				/* merge color with prediction */
				mergeWarped_N(warped_color_views, DispTargs, SAI, 3);
				/* hole filling for color*/
				holefilling(SAI->color, 3, SAI->nr, SAI->nc, 0);
				psnr_w = getYCbCr_422_PSNR(SAI->color, original_color_view, SAI->nr, SAI->nc, 3, BIT_DEPTH);

				if (psnr_w < psnr_med) {
					delete[](SAI->color);
					SAI->color = tmp_m;
					SAI->use_median = true;
					SAI->stdd = 0.0;
				}
				else {
					delete[](tmp_m);
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


		double psnr_without_sparse = 0;

		if (SAI->n_references > 0) {

			//aux_write16PGMPPM(SAI->path_out_ppm, SAI->nc, SAI->nr, 3, SAI->color);

			//psnr_result = USE_difftest_ng ? getPSNR(NULL, SAI->path_out_ppm, SAI->path_input_ppm, difftest_call) : 0;
		
			psnr_without_sparse = getYCbCr_422_PSNR(SAI->color, original_color_view, SAI->nr, SAI->nc, 3, BIT_DEPTH);

			output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", psnr_without_sparse);

		}
		else {
			output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", 0.0);
		}

		unsigned short *colorview_temp = new unsigned short[SAI->nr*SAI->nc * 3]();
		memcpy(colorview_temp, SAI->color, sizeof(unsigned short)*SAI->nr*SAI->nc * 3);

		double psnr_with_sparse = 0;

		if (SAI->NNt > 0 && SAI->Ms > 0)
		{

			SAI->use_global_sparse = true;

			int startt = clock();

			getGlobalSparseFilter(SAI, original_color_view);

			//std::cout << "time elapsed in getGlobalSparseFilter()\t" << (float)( (int)clock() - startt ) / CLOCKS_PER_SEC << "\n";

			applyGlobalSparseFilter(SAI);

			psnr_with_sparse = getYCbCr_422_PSNR(SAI->color, original_color_view, SAI->nr, SAI->nc, 3, BIT_DEPTH);

			//aux_write16PGMPPM(SAI->path_out_ppm, SAI->nc, SAI->nr, 3, SAI->color);
			//psnr_result = USE_difftest_ng ? getPSNR(NULL, SAI->path_out_ppm, SAI->path_input_ppm, difftest_call) : 0;
			
		}

		if (SAI->use_global_sparse) { /* check validity of sparse filter */
			if ( psnr_with_sparse<psnr_without_sparse ) //<0.1
			{
				SAI->use_global_sparse = false;
				memcpy(SAI->color, colorview_temp, sizeof(unsigned short)*SAI->nr*SAI->nc * 3);
			}
		}

		if (SAI->use_global_sparse) {
			output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", psnr_with_sparse);
		}
		else {
			output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", 0.0);
		}

		delete[](colorview_temp);

		char ppm_residual_path[1024];

		char jp2_residual_path_jp2[1024];

		char pgm_residual_Y_path[1024];
		char jp2_residual_Y_path_jp2[1024];
		char pgm_residual_Cb_path[1024];
		char jp2_residual_Cb_path_jp2[1024];
		char pgm_residual_Cr_path[1024];
		char jp2_residual_Cr_path_jp2[1024];

		char pgm_residual_depth_path[1024];

		char jp2_residual_depth_path_jp2[1024];

		char *ycbcr_pgm_names[3];
		char *ycbcr_jp2_names[3];

		float rate_a1 = (float)YUV_RATIO_DEFAULT;

		/* get residual */
		if (SAI->residual_rate_color > 0)
		{

			/* COLOR residual here, lets try YUV */

			sprintf(pgm_residual_Y_path, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_Y_residual.pgm");
			sprintf(jp2_residual_Y_path_jp2, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_Y_residual.jp2");

			sprintf(pgm_residual_Cb_path, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_Cb_residual.pgm");
			sprintf(jp2_residual_Cb_path_jp2, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_Cb_residual.jp2");

			sprintf(pgm_residual_Cr_path, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_Cr_residual.pgm");
			sprintf(jp2_residual_Cr_path_jp2, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_Cr_residual.jp2");

			sprintf(ppm_residual_path, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_residual.ppm");
			sprintf(jp2_residual_path_jp2, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_residual.jp2");

			ycbcr_pgm_names[0] = pgm_residual_Y_path;
			ycbcr_pgm_names[1] = pgm_residual_Cb_path;
			ycbcr_pgm_names[2] = pgm_residual_Cr_path;

			ycbcr_jp2_names[0] = jp2_residual_Y_path_jp2;
			ycbcr_jp2_names[1] = jp2_residual_Cb_path_jp2;
			ycbcr_jp2_names[2] = jp2_residual_Cr_path_jp2;
			

			if (YUV_TRANSFORM) {

				unsigned short *tmpim = new unsigned short[SAI->nr*SAI->nc * 3]();
				memcpy(tmpim, SAI->color, sizeof(unsigned short)*SAI->nr*SAI->nc * 3);

				int offset_v = 0;

				if (RESIDUAL_16BIT) {
					offset_v = (1<<15) - 1;
				}
				else {
					offset_v = (1<<BIT_DEPTH) - 1;
				}

				if (YUV_RATIO_SEARCH) {

					float highest_psnr = 0;

					unsigned short *tmp_im = new unsigned short[SAI->nr*SAI->nc * 3]();

					for (float rate_a = YUV_SEARCH_LOW; rate_a <= YUV_SEARCH_HIGH; rate_a += YUV_SEARCH_STEP) {

						memcpy(tmp_im, SAI->color, sizeof(unsigned short)*SAI->nr*SAI->nc * 3);

						int ncomp_r = 3;

						if (rate_a / 8.0 < MAX_Y_RATIO) {
							SAI->has_chrominance = true;
						}
						else {
							SAI->has_chrominance = false;
							ncomp_r = 1;
						}

						encodeResidualJP2_YUV(SAI->nr, SAI->nc, original_color_view, tmp_im, ycbcr_pgm_names,
							kdu_compress_path, ycbcr_jp2_names, SAI->residual_rate_color, ncomp_r, offset_v, rate_a / 8.0f, RESIDUAL_16BIT_bool);

						decodeResidualJP2_YUV(tmp_im, kdu_expand_path, ycbcr_jp2_names, ycbcr_pgm_names, ncomp_r, offset_v, (1<<BIT_DEPTH) - 1, RESIDUAL_16BIT_bool);

						double psnr_result_yuv = getYCbCr_422_PSNR(tmp_im, original_color_view, SAI->nr, SAI->nc, 3, 10);

						if (psnr_result_yuv > highest_psnr) {
							highest_psnr = (float)psnr_result_yuv;
							rate_a1 = rate_a;
						}

						printf("PSNR YUV:\t%f\t%f\tratio\t%1.2f/%1.2f\n", highest_psnr, psnr_result_yuv, rate_a, 8.0);

					}

					delete[](tmp_im);
				}

				int ncomp_r = 3;

				if (rate_a1 / 8.0 < MAX_Y_RATIO) {
					SAI->has_chrominance = true;
				}
				else {
					SAI->has_chrominance = false;
					ncomp_r = 1;
				}

				encodeResidualJP2_YUV(SAI->nr, SAI->nc, original_color_view, SAI->color, ycbcr_pgm_names,
					kdu_compress_path, ycbcr_jp2_names, SAI->residual_rate_color, ncomp_r, offset_v, 
					rate_a1 / (float)8.0, RESIDUAL_16BIT_bool);

				decodeResidualJP2_YUV(SAI->color, kdu_expand_path, ycbcr_jp2_names, 
					ycbcr_pgm_names, ncomp_r, offset_v, (1<<BIT_DEPTH) - 1, RESIDUAL_16BIT_bool);

				/* also compete against no yuv transformation */

				double psnr_result_yuv_w_trans = getYCbCr_422_PSNR(SAI->color, original_color_view, SAI->nr, SAI->nc, 3, 10);

				offset_v = (1<<BIT_DEPTH) - 1;

				encodeResidualJP2(SAI->nr, SAI->nc, original_color_view, tmpim, ppm_residual_path,
					kdu_compress_path, jp2_residual_path_jp2, SAI->residual_rate_color, 3, offset_v, 
					RESIDUAL_16BIT_bool);

				decodeResidualJP2(tmpim, kdu_expand_path, jp2_residual_path_jp2, 
					ppm_residual_path, ncomp1, offset_v, offset_v, RESIDUAL_16BIT_bool);

				double psnr_result_yuv_wo_trans = getYCbCr_422_PSNR(tmpim, original_color_view, SAI->nr, SAI->nc, 3, 10);

				if (psnr_result_yuv_wo_trans > psnr_result_yuv_w_trans) {
					memcpy(SAI->color, tmpim, sizeof(unsigned short)*SAI->nr*SAI->nc * 3);
					SAI->yuv_transform = false;
					rate_a1 = 0.0;
				}

				delete[](tmpim);

			}
			else {

				int offset_v = (1<<BIT_DEPTH) - 1;

				encodeResidualJP2(SAI->nr, SAI->nc, original_color_view, SAI->color, ppm_residual_path,
					kdu_compress_path, jp2_residual_path_jp2, SAI->residual_rate_color, 3, offset_v, RESIDUAL_16BIT_bool);

				decodeResidualJP2(SAI->color, kdu_expand_path, jp2_residual_path_jp2, ppm_residual_path, ncomp1, offset_v, offset_v, RESIDUAL_16BIT_bool);

			}

			SAI->has_color_residual = true;
		}

		if (SAI->residual_rate_depth > 0 && depth_file_exist) { /* residual depth if needed */

			sprintf(pgm_residual_depth_path, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_depth_residual.pgm");

			sprintf(jp2_residual_depth_path_jp2, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_depth_residual.jp2");

			encodeResidualJP2(SAI->nr, SAI->nc, original_depth_view, SAI->depth, pgm_residual_depth_path,
				kdu_compress_path, jp2_residual_depth_path_jp2, SAI->residual_rate_depth, 1, 0, 1);

			decodeResidualJP2(SAI->depth, kdu_expand_path, jp2_residual_depth_path_jp2, pgm_residual_depth_path, ncomp1, 0, (1<<16) - 1, 1);

			SAI->has_depth_residual = true;
		}

		/* median filter depth */
		if (MEDFILT_DEPTH) {
			unsigned short *tmp_depth = new unsigned short[SAI->nr*SAI->nc]();
			int startt = clock();
			medfilt2D(SAI->depth, tmp_depth, 3, SAI->nr, SAI->nc);
			std::cout << "time elapsed in depth median filtering\t" << (float)( (int)clock() - startt )/CLOCKS_PER_SEC << "\n";
			memcpy(SAI->depth, tmp_depth, sizeof(unsigned short)*SAI->nr*SAI->nc);
			delete[](tmp_depth);
		}

		aux_write16PGMPPM(SAI->path_out_ppm, SAI->nc, SAI->nr, 3, SAI->color);
		aux_write16PGMPPM(SAI->path_out_pgm, SAI->nc, SAI->nr, 1, SAI->depth);

		if (SAI->residual_rate_depth > 0 && depth_file_exist) 
		{
			//psnr_result = USE_difftest_ng ? getPSNR(NULL, SAI->path_out_pgm, SAI->path_input_pgm, difftest_call_pgm) : 0;
			output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", PSNR(SAI->depth, original_depth_view, SAI->nr, SAI->nc, 1) );
		}
		else {
			output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", 0.0);
		}


		//psnr_result = USE_difftest_ng ? getPSNR(NULL, SAI->path_out_ppm, SAI->path_input_ppm, difftest_call) : 0;

		double final_psnr = getYCbCr_422_PSNR(SAI->color, original_color_view, SAI->nr, SAI->nc, 3, BIT_DEPTH);
		psnr_yuv_mean += final_psnr;

		output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", final_psnr );

		output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", rate_a1); 
		output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", SAI->stdd);
		

		delete[](original_color_view);
		delete[](original_depth_view);

		/* write view configuration data to bitstream */

		int n_bytes_prediction = 0, tmp_pred_bytes = 0;
		int n_bytes_residual = 0;

		output_LF_file = fopen(path_out_LF_data, "ab");

		if (!global_header_written) { //these are global
			n_bytes_prediction += (int)fwrite(&SAI->nr, sizeof(int), 1, output_LF_file) * sizeof(int); // needed only once per LF
			n_bytes_prediction += (int)fwrite(&SAI->nc, sizeof(int), 1, output_LF_file) * sizeof(int); // 
			n_bytes_prediction += (int)fwrite(&yuv_transform_s, sizeof(int), 1, output_LF_file) * sizeof(int); 
			n_bytes_prediction += (int)fwrite(&MINIMUM_DEPTH, sizeof(unsigned short), 1, output_LF_file) * sizeof(unsigned short);
			global_header_written = true;
		}

		viewHeaderToCodestream(n_bytes_prediction, SAI, output_LF_file, yuv_transform_s);

		/* debugging */
		tmp_codestream = fopen(path_codestream, "ab");
		viewHeaderToCodestream(tmp_pred_bytes, SAI, tmp_codestream, yuv_transform_s);
		fclose(tmp_codestream);

		if (SAI->residual_rate_color > 0) {

			if (SAI->yuv_transform && YUV_TRANSFORM) {

				int ncomp_r = SAI->has_chrominance ? 3 : 1;

				for (int icomp = 0; icomp < ncomp_r; icomp++) {

					writeResidualToDisk(ycbcr_jp2_names[icomp], output_LF_file, n_bytes_residual, JP2_dict);

				}
			}
			else {

				writeResidualToDisk(jp2_residual_path_jp2, output_LF_file, n_bytes_residual, JP2_dict);
			}
		}

		if (SAI->residual_rate_depth > 0 && depth_file_exist) {

			writeResidualToDisk(jp2_residual_depth_path_jp2, output_LF_file, n_bytes_residual, JP2_dict);

		}

		fclose(output_LF_file);

		printf("encoded: %i kilobytes\t\tPSNR YUV (%03d_%03d): %2.3f\tPSNR YUV mean: %2.3f\n", aux_GetFileSize(path_out_LF_data) / 1000, SAI->r,SAI->c, final_psnr, psnr_yuv_mean/(ii+1));

		output_buffer_length += sprintf(output_results + output_buffer_length, "\t%i", n_bytes_prediction);
		output_buffer_length += sprintf(output_results + output_buffer_length, "\t%i", n_bytes_residual);

		output_results_file = fopen(output_results_filename, "a");
		fprintf(output_results_file, "%s\n", output_results);
		fclose(output_results_file);

		/* to reduce memory usage */
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

		if (SAI->segmentation != NULL) {
			delete[](SAI->segmentation);
			SAI->segmentation = NULL;
		}

	}

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
		if (SAI->segmentation != NULL)
			delete[](SAI->segmentation);

	}

	delete[](LF);

	exit(0);
}