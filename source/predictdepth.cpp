/*BSD 2-Clause License
* Copyright(c) 2019, Pekka Astola
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met :
*
* 1. Redistributions of source code must retain the above copyright notice, this
* list of conditions and the following disclaimer.
*
* 2. Redistributions in binary form must reproduce the above copyright notice,
* this list of conditions and the following disclaimer in the documentation
* and/or other materials provided with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED.IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
*     SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*     CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
*     OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
*     OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "predictdepth.hh"
#include "warping.hh"
#include "medianfilter.hh"
#include "ppm.hh"
#include "inpainting.hh"
#include "merging.hh"

#include <ctime>
#include <vector>
#include <cstdint>

using std::int32_t;
using std::uint32_t;

using std::int16_t;
using std::uint16_t;

using std::int8_t;
using std::uint8_t;

void WaSP_predict_depth(view* SAI, view *LF) {
    /* forward warp depth */
    if (SAI->n_depth_references > 0) {

        printf("Predicting normalized disparity for view %03d_%03d\n", SAI->c, SAI->r);

        uint16_t **warped_texture_views_0_N =
            new uint16_t*[SAI->n_depth_references]();
        uint16_t **warped_depth_views_0_N =
            new uint16_t*[SAI->n_depth_references]();
        float **DispTargs_0_N =
            new float*[SAI->n_depth_references]();

        init_warping_arrays(
            SAI->n_depth_references,
            warped_texture_views_0_N,
            warped_depth_views_0_N,
            DispTargs_0_N,
            SAI->nr,
            SAI->nc,
            SAI->ncomp);

        for (int32_t ij = 0; ij < SAI->n_depth_references; ij++) {
            view *ref_view = LF + SAI->depth_references[ij];

            int32_t tmp_w, tmp_r, tmp_ncomp;

            aux_read16PGMPPM(
                ref_view->path_out_pgm,
                tmp_w,
                tmp_r,
                tmp_ncomp,
                ref_view->depth);

            ref_view->color = new uint16_t[ref_view->nr*ref_view->nc * 3]();

            //aux_read16PGMPPM(ref_view->path_out_ppm, tmp_w, tmp_r, tmp_ncomp,
            //                 ref_view->color);

            warpView0_to_View1(
                ref_view,
                SAI,
                warped_texture_views_0_N[ij],
                warped_depth_views_0_N[ij],
                DispTargs_0_N[ij]);

            delete[](ref_view->depth);
            delete[](ref_view->color);

            ref_view->depth = nullptr;
            ref_view->color = nullptr;
        }

        /* merge depth using median*/

        //int32_t startt = clock();

        double *hole_mask = new double[SAI->nr*SAI->nc]();

        for (int32_t ij = 0; ij < SAI->nr * SAI->nc; ij++) {

            hole_mask[ij] = INIT_DISPARITY_VALUE;

            std::vector<uint16_t> depth_values;
            for (int32_t uu = 0; uu < SAI->n_depth_references; uu++) {
                uint16_t *pp = warped_depth_views_0_N[uu];
                float *pf = DispTargs_0_N[uu];
                if (*(pf + ij) > INIT_DISPARITY_VALUE) {
                    depth_values.push_back(*(pp + ij));
                }
            }
            if (depth_values.size() > 0) {
                SAI->depth[ij] = getMedian(depth_values);
                hole_mask[ij] = 1.0f;
            }
        }

        uint32_t nholes = holefilling(
            SAI->depth,
            SAI->nr,
            SAI->nc,
            INIT_DISPARITY_VALUE,
            hole_mask);

        delete[](hole_mask);

        clean_warping_arrays(
            SAI->n_depth_references,
            warped_texture_views_0_N,
            warped_depth_views_0_N,
            DispTargs_0_N);

    }

}
