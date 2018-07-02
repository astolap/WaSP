#include <cstdio>
#include <cstdint>
#include <ctime>
#include <iostream>
#include <vector>
#include <string>
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


#define USE_difftest_ng false

#define YUV_SEARCH_LOW 6.80f
#define YUV_SEARCH_HIGH 7.80f
#define YUV_SEARCH_STEP 0.20f
#define YUV_RATIO_DEFAULT 7.20f

#define STD_SEARCH_LOW 10
#define STD_SEARCH_HIGH 250
#define STD_SEARCH_STEP 10

#define SAVE_PARTIAL_WARPED_VIEWS false

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

		SAI->x = float(xx) / 100000;
		SAI->y = float(yy) / 100000;

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

		fread(&SAI->min_inv_d, sizeof(int), 1, filept); /*reading, if we have negative inverse depth,
												   for example in lenslet, we need to subtract min_inv_d
												   from the inverse depth maps*/

		fread(&(SAI->n_references), sizeof(int), 1, filept); /*reading*/

		if (SAI->n_references > 0) {

			SAI->references = new int[SAI->n_references]();

			fread(SAI->references, sizeof(int), SAI->n_references, filept); /*reading*/

		}

		fread(&(SAI->n_depth_references), sizeof(int), 1, filept); /*reading*/
		if (SAI->n_depth_references > 0) {

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


	/* 2D array format for views, useful in some cases */
	view ***LF_mat = new view**[maxR]();
	for (int ii = 0; ii < maxR; ii++) {
		LF_mat[ii] = new view*[maxC]();
	}

	for (int r = 0; r < maxR; r++) {
		for (int c = 0; c < maxC; c++) {
			LF_mat[r][c] = NULL;
		}
	}

	for (int ii = 0; ii < n_views_total; ii++) {
		view *SAI = LF + ii;
		LF_mat[SAI->r][SAI->c] = SAI;
	}

	/* get motion vectors/displacements/warping tables for segmentations, NOT USED IN THIS VERSION */
	for (int ii = 0; ii < n_views_total; ii++)
	{
		view *SAI = LF + ii;

		if (SAI->has_segmentation > 0) { /* displacement/disparity search is only applied to views which have a segmentation */

			unsigned short *temp_segmentation;
			int ncomp1; //nc_temp and nr_temp should equal SAI->nr,SAI->nc

			if (!aux_read16PGMPPM(SAI->path_input_seg, SAI->nc, SAI->nr, ncomp1, temp_segmentation)) {
				printf("Error reading segmentation %s\nExiting...\n", SAI->path_input_seg);
				exit(0);
			}

			SAI->segmentation = new unsigned short[SAI->nr*SAI->nc]();

			int maxL = 0;

			for (int jj = 0; jj < SAI->nr*SAI->nc; jj++) {
				*(SAI->segmentation + jj) = (unsigned short)*(temp_segmentation + jj);

				maxL = *(SAI->segmentation + jj) > maxL ? *(SAI->segmentation + jj) : maxL;

			}

			SAI->maxL = maxL;

			/* initialize region displacement table */
			SAI->region_displacements = new int***[maxR]();
			for (int iar = 0; iar < maxR; iar++) {
				SAI->region_displacements[iar] = new int**[maxC]();
				for (int iac = 0; iac < maxC; iac++) {
					SAI->region_displacements[iar][iac] = new int*[maxL + 1]();
					for (int iR = 0; iR <= maxL; iR++) {
						SAI->region_displacements[iar][iac][iR] = new int[2]();
					}
				}
			}

			delete[](temp_segmentation);

			/* motion vectors for all regions against all views */

			unsigned short *im0;

			aux_read16PGMPPM(SAI->path_input_ppm, SAI->nc, SAI->nr, ncomp1, im0);

			for (int jj = 0; jj < n_views_total; jj++) {

				view *SAI1 = LF + jj;

				if (jj != ii) {

					unsigned short *im1;
					aux_read16PGMPPM(SAI1->path_input_ppm, SAI1->nc, SAI1->nr, ncomp1, im1);

					int search_radius = 10; // search -search_radius:search_radius in both x,y

					float ***match_score = new float**[maxL + 1]();
					for (int ik = 0; ik <= maxL; ik++) {
						match_score[ik] = new float*[search_radius * 2 + 1]();
						for (int isr = 0; isr < search_radius * 2 + 1; isr++) {
							match_score[ik][isr] = new float[search_radius * 2 + 1]();
						}
					}

					int ***counts = new int**[maxL + 1]();
					for (int ik = 0; ik <= maxL; ik++) {
						counts[ik] = new int*[search_radius * 2 + 1]();
						for (int isr = 0; isr < search_radius * 2 + 1; isr++) {
							counts[ik][isr] = new int[search_radius * 2 + 1]();
						}
					}

					unsigned short *segp = SAI->segmentation;

					for (int ijk = 0; ijk < SAI->nr*SAI->nc; ijk++) {
						for (int ik = 0; ik <= maxL; ik++) {
							if (*(segp + ijk) == ik) {
#pragma omp parallel for
								for (int isr = -search_radius; isr <= search_radius; isr++) {
									for (int isc = -search_radius; isc <= search_radius; isc++) {

										int iy = ijk % SAI->nr; //row
										int ix = (ijk - iy) / SAI->nr; //col

										int iy1 = iy + isr;
										int ix1 = ix + isc;

										if (iy1 >= 0 && iy1 < SAI->nr && ix1 >= 0 && ix1 < SAI->nc) {

											int ijk1 = ix1*SAI->nr + iy1;

											for (int ic = 0; ic < 3; ic++) {
												int offc = ic*SAI->nr*SAI->nc;
												match_score[ik][isr + search_radius][isc + search_radius] += abs(((float)*(im0 + ijk + offc) - (float)*(im1 + ijk1 + offc)));
											}

											//printf("%f\n", match_score[ik][isr + search_radius][isc + search_radius]);

											counts[ik][isr + search_radius][isc + search_radius]++;

											//printf("%d\n", counts[ik][isr + search_radius][isc + search_radius]);

										}
									}
								}
							}
						}
					}



					delete[](im1);

					/* find best matching displacement */
					for (int ik = 0; ik <= maxL; ik++) {
						float lowest_score = FLT_MAX;
						for (int isr = 0; isr < search_radius * 2 + 1; isr++) {
							for (int isc = 0; isc < search_radius * 2 + 1; isc++) {

								int nk = counts[ik][isr][isc];

								if (nk > 0) {

									float clow = match_score[ik][isr][isc] / (float)nk;

									//printf("%f\n", clow);

									if (clow < lowest_score) {
										SAI->region_displacements[SAI1->r][SAI1->c][ik][0] = isr - search_radius;
										SAI->region_displacements[SAI1->r][SAI1->c][ik][1] = isc - search_radius;
										lowest_score = clow;
									}

								}
							}
						}
					}

					/* Fill in missing depth values. This step needs to be done sequentially since we propagate the missing values from side-view to side-view. */
					/*% find as reference the already solved view closest to the center that
					% is a neighbor of the current view*/

					int closest_jj = 0;
					float smallest_distance = FLT_MAX;

					for (int dy = -1; dy <= 1; dy++) {
						for (int dx = -1; dx <= 1; dx++) {

							view *SAI_f = LF_mat[SAI->r + dy][SAI->c + dx];

							if (SAI_f->i_order < jj) {

								float dist = sqrt((float)((SAI_f->r - SAI->r)*(SAI_f->r - SAI->r) + 
									(SAI_f->c- SAI->c)*(SAI_f->c - SAI->c)));

								if (dist < smallest_distance) {
									smallest_distance = dist;
									closest_jj = SAI_f->i_order;
								}

							}

						}
					}


					for (int iR = 0; iR <= maxL; iR++) {
						printf("SAI->region_displacements[%d][%d][%d][0] = %d\t", SAI1->r, SAI1->c, iR, SAI->region_displacements[SAI1->r][SAI1->c][iR][0]);
						printf("SAI->region_displacements[%d][%d][%d][1] = %d\n", SAI1->r, SAI1->c, iR, SAI->region_displacements[SAI1->r][SAI1->c][iR][1]);
					}


					unsigned short *seg_warped;
					unsigned short *color_seg_warped;

					seg_warped = new unsigned short[SAI->nr*SAI->nc]();
					color_seg_warped = new unsigned short[SAI->nr*SAI->nc * 3]();

					for (int ijk = 0; ijk < SAI->nr*SAI->nc; ijk++) {
						for (int ik = 0; ik <= maxL; ik++) {
							if (*(segp + ijk) == ik) {
								int iy = ijk % SAI->nr; //row
								int ix = (ijk - iy) / SAI->nr; //col

								int iy1 = iy + SAI->region_displacements[SAI1->r][SAI1->c][ik][0];
								int ix1 = ix + SAI->region_displacements[SAI1->r][SAI1->c][ik][1];

								if (iy1 >= 0 && iy1 < SAI->nr && ix1 >= 0 && ix1 < SAI->nc) {

									int ijk1 = ix1*SAI->nr + iy1;

									if (*(seg_warped + ijk1) == 0) {

										for (int ic = 0; ic < 3; ic++) {
											int offc = ic*SAI->nr*SAI->nc;
											if (ic < 1) {
												*(seg_warped + ijk1 + offc) = ik;
											}
											*(color_seg_warped + ijk1 + offc) = *(im0 + ijk + offc);
										}

									}
								}
							}
						}
					}

					/* now use the chosen view closest_jj and use 3x3 pixel neighborhood to fill in missing depth with the maximum */
					if (closest_jj > 0) {
						view *SAI_f = LF + closest_jj;

						unsigned short *seg_f = SAI_f->segmentation;

						for (int ijk = 0; ijk < SAI->nr*SAI->nc; ijk++) {

							if (*(seg_warped + ijk) == 0) {

								int iy = ijk % SAI->nr; //row
								int ix = (ijk - iy) / SAI->nr; //col

								unsigned short largest_depth = 0;

								for (int dy = -1; dy <= 1; dy++) {
									for (int dx = -1; dx <= 1; dx++) {

										int iy1 = iy + dy;
										int ix1 = ix + dx;

										if (iy1 >= 0 && iy1 < SAI->nr && ix1 >= 0 && ix1 < SAI->nc) {

											unsigned short dval = *(seg_f + ijk + dx*SAI->nr + dy);

											if (dval > largest_depth) {
												largest_depth = dval;
											}

										}
									}
								}

								*(seg_warped + ijk) = largest_depth;

							}
						}
					}

					unsigned short *DM_ROW_tmp = new unsigned short[SAI->nr*SAI->nc]();
					unsigned short *DM_COL_tmp = new unsigned short[SAI->nr*SAI->nc]();

					/* make row and column disparity maps, save to disk */
					for (int ijk = 0; ijk < SAI->nr*SAI->nc; ijk++) {
						int ik = *(seg_warped + ijk);
						DM_ROW_tmp[ijk] = (unsigned short)SAI->region_displacements[SAI1->r][SAI1->c][ik][0] + 20;
						DM_COL_tmp[ijk] = (unsigned short)SAI->region_displacements[SAI1->r][SAI1->c][ik][1] + 20;
					}

					char tempchar[1024];
					sprintf(tempchar, "%s%03d_%03d%s%03d_%03d%s", output_dir, (SAI)->c, (SAI)->r, "_segmentation_warped_to_", SAI1->c, SAI1->r, ".pgm");
					aux_write16PGMPPM(tempchar, SAI->nc, SAI->nr, 1, seg_warped);

					memset(tempchar, 0x00, 1024 * sizeof(char));
					sprintf(tempchar, "%s%03d_%03d%s%03d_%03d%s", output_dir, (SAI)->c, (SAI)->r, "_color_warped_to_", SAI1->c, SAI1->r, ".ppm");
					aux_write16PGMPPM(tempchar, SAI->nc, SAI->nr, 3, color_seg_warped);

					memset(tempchar, 0x00, 1024 * sizeof(char));
					sprintf(tempchar, "%s%03d_%03d%s%03d_%03d%s", output_dir, (SAI)->c, (SAI)->r, "_warped_to_", SAI1->c, SAI1->r, "_DM_ROW.pgm");
					aux_write16PGMPPM(tempchar, SAI->nc, SAI->nr, 1, DM_ROW_tmp);

					memset(tempchar, 0x00, 1024 * sizeof(char));
					sprintf(tempchar, "%s%03d_%03d%s%03d_%03d%s", output_dir, (SAI)->c, (SAI)->r, "_warped_to_", SAI1->c, SAI1->r, "_DM_COL.pgm");
					aux_write16PGMPPM(tempchar, SAI->nc, SAI->nr, 1, DM_COL_tmp);

					SAI1->segmentation = seg_warped;

					delete[](DM_ROW_tmp);
					delete[](DM_COL_tmp);

					//delete[](seg_warped); //!!!!!!!!!!!!!!!!!!!!!!!!!!!!! we need to clean this still
					delete[](color_seg_warped);


					//char dummy;
					//std::cin >> dummy;

					for (int ik = 0; ik <= maxL; ik++) {
						for (int isr = 0; isr < search_radius * 2 + 1; isr++) {
							delete[](match_score[ik][isr]);
							delete[](counts[ik][isr]);
						}
						delete[](match_score[ik]);
						delete[](counts[ik]);
					}
					delete[](match_score);
					delete[](counts);

				}
				else {
					/* nothing to do, for each region displacements are zero*/
				}

			}

			delete[](im0);

		}
	}


	char path_out_LF_data[1024];
	sprintf(path_out_LF_data, "%s%c%s", output_dir, '/', "output.LF");

	/* our bitstream starts here */
	FILE *output_LF_file;
	output_LF_file = fopen(path_out_LF_data, "wb");
	fwrite(&n_views_total, sizeof(int), 1, output_LF_file);
	fclose(output_LF_file);

	bool size_written = false;

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

	for (int ii = 0; ii < n_views_total; ii++) {

		view *SAI = LF + ii;

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
			if (SAI->stdd == 0) {
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
				std::cout << "time elapsed in color median merging\t" << (int)clock() - startt << "\n";
				holefilling(SAI->color, 3, SAI->nr, SAI->nc, 0);

				//double psnr_med = getYCbCr_422_PSNR(SAI->color, original_color_view, SAI->nr, SAI->nc, 3, 10);
				double psnr_med = getYCbCr_422_PSNR(SAI->color, original_color_view, SAI->nr, SAI->nc, 3, 10);

				unsigned short *tmp_m = new unsigned short[SAI->nr*SAI->nc*3]();
				memcpy(tmp_m, SAI->color, sizeof(unsigned short)*SAI->nr*SAI->nc * 3);

				double psnr_w = 0;

				if (STD_SEARCH) {

					float stdi = 0;

					//double stds[] = { 5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85 };
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


		//float psnr_result;

		if (SAI->n_references > 0) {

			//aux_write16PGMPPM(SAI->path_out_ppm, SAI->nc, SAI->nr, 3, SAI->color);

			//psnr_result = USE_difftest_ng ? getPSNR(NULL, SAI->path_out_ppm, SAI->path_input_ppm, difftest_call) : 0;
		
			output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", getYCbCr_422_PSNR(SAI->color, original_color_view, SAI->nr, SAI->nc, 3, BIT_DEPTH) );

		}
		else {
			output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", 0.0);
		}

		if (SAI->NNt > 0 && SAI->Ms > 0)
		{

			int startt = clock();

			getGlobalSparseFilter(SAI, original_color_view);

			std::cout << "time elapsed in getGlobalSparseFilter()\t" << (int)clock() - startt << "\n";

			applyGlobalSparseFilter(SAI);

			//aux_write16PGMPPM(SAI->path_out_ppm, SAI->nc, SAI->nr, 3, SAI->color);
			//psnr_result = USE_difftest_ng ? getPSNR(NULL, SAI->path_out_ppm, SAI->path_input_ppm, difftest_call) : 0;
			output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", getYCbCr_422_PSNR(SAI->color, original_color_view, SAI->nr, SAI->nc, 3, BIT_DEPTH) );

		}
		else
		{
			output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", 0.0);
		}

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

		float rate_a1 = 0;

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

				int offset_v = 0;

				if (RESIDUAL_16BIT) {
					offset_v = (1<<15) - 1;
				}
				else {
					offset_v = (1<<BIT_DEPTH) - 1;
				}

				//float rate_a = 6.5 / 8.0;// 7.2 / 8.0;

				if (YUV_RATIO_SEARCH) {

					float highest_psnr = 0;

					unsigned short *tmp_im = new unsigned short[SAI->nr*SAI->nc * 3]();

					for (float rate_a = YUV_SEARCH_LOW; rate_a <= YUV_SEARCH_HIGH; rate_a += YUV_SEARCH_STEP) {

						memcpy(tmp_im, SAI->color, sizeof(unsigned short)*SAI->nr*SAI->nc * 3);

						encodeResidualJP2_YUV(SAI->nr, SAI->nc, original_color_view, tmp_im, ycbcr_pgm_names,
							kdu_compress_path, ycbcr_jp2_names, SAI->residual_rate_color, 3, offset_v, rate_a / 8.0f, RESIDUAL_16BIT_bool);

						decodeResidualJP2_YUV(tmp_im, kdu_expand_path, ycbcr_jp2_names, ycbcr_pgm_names, 3, offset_v, (1<<BIT_DEPTH) - 1, RESIDUAL_16BIT_bool);

						double psnr_result_yuv = getYCbCr_422_PSNR(tmp_im, original_color_view, SAI->nr, SAI->nc, 3, 10);

						if (psnr_result_yuv > highest_psnr) {
							highest_psnr = (float)psnr_result_yuv;
							rate_a1 = rate_a;
							printf("PSNR YUV:\t%f\tratio\t%f/%f\n", highest_psnr, rate_a1, 8.0);
						}

					}

					delete[](tmp_im);
				}
				else {
					rate_a1 = (float)YUV_RATIO_DEFAULT;
				}

				unsigned short *tmpim = new unsigned short[SAI->nr*SAI->nc * 3]();
				memcpy(tmpim, SAI->color, sizeof(unsigned short)*SAI->nr*SAI->nc * 3);

				encodeResidualJP2_YUV(SAI->nr, SAI->nc, original_color_view, SAI->color, ycbcr_pgm_names,
					kdu_compress_path, ycbcr_jp2_names, SAI->residual_rate_color, 3, offset_v, rate_a1 / (float)8.0, RESIDUAL_16BIT_bool);

				decodeResidualJP2_YUV(SAI->color, kdu_expand_path, ycbcr_jp2_names, ycbcr_pgm_names, 3, offset_v, (1<<BIT_DEPTH) - 1, RESIDUAL_16BIT_bool);

				/* also compete against no yuv transformation */

				double psnr_result_yuv_w_trans = getYCbCr_422_PSNR(SAI->color, original_color_view, SAI->nr, SAI->nc, 3, 10);

				offset_v = (1<<BIT_DEPTH) - 1;

				encodeResidualJP2(SAI->nr, SAI->nc, original_color_view, tmpim, ppm_residual_path,
					kdu_compress_path, jp2_residual_path_jp2, SAI->residual_rate_color, 3, offset_v, RESIDUAL_16BIT_bool);

				decodeResidualJP2(tmpim, kdu_expand_path, jp2_residual_path_jp2, ppm_residual_path, ncomp1, offset_v, offset_v, RESIDUAL_16BIT_bool);

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

		}

		if (SAI->residual_rate_depth > 0 && depth_file_exist) { /* residual depth if needed */

			sprintf(pgm_residual_depth_path, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_depth_residual.pgm");

			sprintf(jp2_residual_depth_path_jp2, "%s%c%03d_%03d%s", output_dir, '/', SAI->c, SAI->r, "_depth_residual.jp2");

			encodeResidualJP2(SAI->nr, SAI->nc, original_depth_view, SAI->depth, pgm_residual_depth_path,
				kdu_compress_path, jp2_residual_depth_path_jp2, SAI->residual_rate_depth, 1, 0, 1);

			decodeResidualJP2(SAI->depth, kdu_expand_path, jp2_residual_depth_path_jp2, pgm_residual_depth_path, ncomp1, 0, (1<<16) - 1, 1);

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
		output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", getYCbCr_422_PSNR(SAI->color, original_color_view, SAI->nr, SAI->nc, 3, BIT_DEPTH) );

		output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", rate_a1); 
		output_buffer_length += sprintf(output_results + output_buffer_length, "\t%f", SAI->stdd);
		

		delete[](original_color_view);
		delete[](original_depth_view);

		/* write view configuration data to bitstream */

		int n_bytes_prediction = 0;
		int n_bytes_residual = 0;

		output_LF_file = fopen(path_out_LF_data, "ab");

		if ( !size_written ) {
			n_bytes_prediction += (int)fwrite(&SAI->nr, sizeof(int), 1, output_LF_file)* sizeof(int); // needed only once per LF
			n_bytes_prediction += (int)fwrite(&SAI->nc, sizeof(int), 1, output_LF_file)* sizeof(int); // 
			n_bytes_prediction += (int)fwrite(&yuv_transform_s, sizeof(int), 1, output_LF_file)* sizeof(int);
			size_written = true;
		}

		minimal_config mconf = makeMinimalConfig(SAI);

		printf("size of minimal_config %i bytes\n", (int)sizeof(minimal_config));

		n_bytes_prediction += (int)fwrite(&mconf, sizeof(minimal_config), 1, output_LF_file)* sizeof(minimal_config);

		/* lets see what else needs to be written to bitstream */

		if (mconf.n_references > 0) {
			for (int ij = 0; ij < mconf.n_references; ij++) {
				unsigned short nid = (unsigned short) *(SAI->references + ij);
				n_bytes_prediction += (int)fwrite(&nid, sizeof(unsigned short), 1, output_LF_file)* sizeof(unsigned short);
			}
		}

		if (mconf.n_depth_references > 0) {
			for (int ij = 0; ij < mconf.n_depth_references; ij++) {
				unsigned short nid = (unsigned short) *(SAI->depth_references + ij);
				n_bytes_prediction += (int)fwrite(&nid, sizeof(unsigned short), 1, output_LF_file) * sizeof(unsigned short);
			}
		}

		if (mconf.Ms > 0 && mconf.NNt > 0) {
			n_bytes_prediction += (int)fwrite(SAI->sparse_mask, sizeof(unsigned char), SAI->Ms, output_LF_file)* sizeof(unsigned char);
			n_bytes_prediction += (int)fwrite(SAI->sparse_weights, sizeof(int32_t), SAI->Ms, output_LF_file)* sizeof(int32_t);
		}

		if (mconf.use_median < 1) {
			if (mconf.use_std < 1) {
				if (mconf.n_references>0) {
					/* use LS merging weights */
					n_bytes_prediction += (int)fwrite(SAI->merge_weights, sizeof(signed short), SAI->NB / 2, output_LF_file)* sizeof(signed short);
				}
			}
			else {
				/* use standard deviation */
				n_bytes_prediction += (int)fwrite(&SAI->stdd, sizeof(float), 1, output_LF_file) * sizeof(signed short);
			}
		}

		if (SAI->residual_rate_color > 0) {

			if (SAI->yuv_transform && YUV_TRANSFORM) {
				for (int icomp = 0; icomp < 3; icomp++) {
					int n_bytes_color_residual = aux_GetFileSize(ycbcr_jp2_names[icomp]);

					unsigned char *jp2_residual = new unsigned char[n_bytes_color_residual]();
					FILE *jp2_color_residual_file = fopen(ycbcr_jp2_names[icomp], "rb");
					fread(jp2_residual, sizeof(unsigned char), n_bytes_color_residual, jp2_color_residual_file);
					fclose(jp2_color_residual_file);

					n_bytes_residual += (int)fwrite(&n_bytes_color_residual, sizeof(int), 1, output_LF_file)* sizeof(int);
					n_bytes_residual += (int)fwrite(jp2_residual, sizeof(unsigned char), n_bytes_color_residual, output_LF_file)* sizeof(unsigned char);

					delete[](jp2_residual);
				}
			}
			else {
				int n_bytes_color_residual = aux_GetFileSize(jp2_residual_path_jp2);

				unsigned char *jp2_residual = new unsigned char[n_bytes_color_residual]();
				FILE *jp2_color_residual_file = fopen(jp2_residual_path_jp2, "rb");
				fread(jp2_residual, sizeof(unsigned char), n_bytes_color_residual, jp2_color_residual_file);
				fclose(jp2_color_residual_file);

				n_bytes_residual += (int)fwrite(&n_bytes_color_residual, sizeof(int), 1, output_LF_file)* sizeof(int);
				n_bytes_residual += (int)fwrite(jp2_residual, sizeof(unsigned char), n_bytes_color_residual, output_LF_file)* sizeof(unsigned char);

				delete[](jp2_residual);
			}
		}
		else {
			int n_bytes_color_residual = 0;
			n_bytes_residual += (int)fwrite(&n_bytes_color_residual, sizeof(int), 1, output_LF_file)* sizeof(int);
		}

		if (SAI->residual_rate_depth > 0 && depth_file_exist) {
			int n_bytes_depth_residual = aux_GetFileSize(jp2_residual_depth_path_jp2);

			unsigned char *jp2_depth_residual = new unsigned char[n_bytes_depth_residual]();
			FILE *jp2_depth_residual_file = fopen(jp2_residual_depth_path_jp2, "rb");
			fread(jp2_depth_residual, sizeof(unsigned char), n_bytes_depth_residual, jp2_depth_residual_file);
			fclose(jp2_depth_residual_file);

			n_bytes_residual += (int)fwrite(&n_bytes_depth_residual, sizeof(int), 1, output_LF_file)* sizeof(int);
			n_bytes_residual += (int)fwrite(jp2_depth_residual, sizeof(unsigned char), n_bytes_depth_residual, output_LF_file)* sizeof(unsigned char);

			delete[](jp2_depth_residual);
		}
		else {
			int n_bytes_depth_residual = 0;
			n_bytes_residual += (int)fwrite(&n_bytes_depth_residual, sizeof(int), 1, output_LF_file)* sizeof(int);
		}

		fclose(output_LF_file);

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

		printf("ii=%d\n", ii);

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

	
	for (int ii = 0; ii < maxR; ii++) {
		delete[](LF_mat[ii]);
	}
	delete[](LF_mat);

	delete[](LF);

	exit(0);
}