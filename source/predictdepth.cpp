#include "predictdepth.hh"
#include "warping.hh"
#include "medianfilter.hh"
#include "ppm.hh"
#include "inpainting.hh"

#include <ctime>
#include <vector>

void predictDepth(view* SAI, view *LF)
{
	/* forward warp depth */
	if (SAI->n_depth_references > 0) {
		/* currently we forward warp the depth from the N (for HDCA N = 5, lenslet maybe 1?) references */
		unsigned short **warped_color_views_0_N = new unsigned short*[SAI->n_depth_references]();
		unsigned short **warped_depth_views_0_N = new unsigned short*[SAI->n_depth_references]();
		float **DispTargs_0_N = new float*[SAI->n_depth_references]();

		for (int ij = 0; ij < SAI->n_depth_references; ij++)
		{
			view *ref_view = LF + SAI->depth_references[ij];

			int tmp_w, tmp_r, tmp_ncomp;

			aux_read16PGMPPM(ref_view->path_out_pgm, tmp_w, tmp_r, tmp_ncomp, ref_view->depth);
			aux_read16PGMPPM(ref_view->path_out_ppm, tmp_w, tmp_r, tmp_ncomp, ref_view->color);

			warpView0_to_View1(ref_view, SAI, warped_color_views_0_N[ij], warped_depth_views_0_N[ij], DispTargs_0_N[ij]);

			delete[](ref_view->depth);
			delete[](ref_view->color);

			ref_view->depth = NULL;
			ref_view->color = NULL;
		}

		/* merge depth using median*/

		int startt = clock();

#pragma omp parallel for
		for (int ij = 0; ij < SAI->nr*SAI->nc; ij++) {
			std::vector<unsigned short> depth_values;
			for (int uu = 0; uu < SAI->n_depth_references; uu++) {
				//for (int uu = 0; uu < SAI->n_references; uu++) {
				unsigned short *pp = warped_depth_views_0_N[uu];
				float *pf = DispTargs_0_N[uu];
				if (*(pf + ij) > -1) {
					depth_values.push_back(*(pp + ij));
				}
			}
			if (depth_values.size() > 0)
				SAI->depth[ij] = getMedian(depth_values);
		}

		//std::cout << "time elapsed in depth merging\t" << (int)clock() - startt << "\n";

		/* hole filling for depth */
		holefilling(SAI->depth, 1, SAI->nr, SAI->nc, 0);

		for (int ij = 0; ij < SAI->n_depth_references; ij++)
		{
			delete[](warped_color_views_0_N[ij]);
			delete[](warped_depth_views_0_N[ij]);
			delete[](DispTargs_0_N[ij]);
		}

		delete[](warped_color_views_0_N);
		delete[](warped_depth_views_0_N);
		delete[](DispTargs_0_N);
	}

}
