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


#include <cstdio>
#include <cstring>
#include <cmath>

#include "decoder.hh"
#include "view.hh"
#include "minconf.hh"
#include "codestream.hh"
#include "predictdepth.hh"
#include "residualjp2.hh"
#include "medianfilter.hh"
#include "ppm.hh"
#include "warping.hh"
#include "merging.hh"
#include "inpainting.hh"
#include "sparsefilter.hh"
#include "WaSPConf.hh"

#define SAVE_PARTIAL_WARPED_VIEWS false


decoder::decoder(const WaSPsetup decoder_setup)
{
    setup = decoder_setup;

    input_LF = fopen(setup.input_directory.c_str(), "rb");

}

decoder::~decoder() {
    if (LF != nullptr) {
        dealloc();
    }
    fclose(input_LF);
}

void decoder::decode() {
    decode_header();
    decode_views();
}

void decoder::decode_header() {

    n_bytes_prediction += (int32_t)fread(
        &number_of_views,
        sizeof(int32_t),
        1,
        input_LF)
        * sizeof(int32_t);

    n_bytes_prediction += (int32_t)fread(
        &number_of_rows,
        sizeof(int32_t),
        1,
        input_LF)
        * sizeof(int32_t);

    n_bytes_prediction += (int32_t)fread(
        &number_of_columns,
        sizeof(int32_t),
        1,
        input_LF) * sizeof(int32_t);

    n_bytes_prediction += (int32_t)fread(
        &minimum_depth,
        sizeof(uint16_t),
        1,
        input_LF) * sizeof(uint16_t);

    uint8_t colorspace_enumerator;

    n_bytes_prediction += (int32_t)fread(
        &colorspace_enumerator,
        sizeof(uint8_t),
        1,
        input_LF) * sizeof(uint8_t);

    if (colorspace_enumerator == 0) {
        colorspace_LF = "RGB";
    }
    if (colorspace_enumerator == 1) {
        colorspace_LF = "YCbCr";
    }

    n_bytes_prediction += (int32_t)fread(
        &maxh,
        sizeof(int32_t),
        1,
        input_LF) * sizeof(uint8_t);

}

void decoder::forward_warp_texture_references(
    view *LF,
    view *SAI,
    uint16_t **warped_texture_views,
    uint16_t **warped_depth_views,
    float **DispTargs) {

    for (int32_t ij = 0; ij < SAI->n_references; ij++) {

        view *ref_view = LF + SAI->references[ij];

        int32_t tmp_w, tmp_r, tmp_ncomp;

        aux_read16PGMPPM(
            ref_view->path_out_pgm,
            tmp_w,
            tmp_r,
            tmp_ncomp,
            ref_view->depth);

        aux_read16PGMPPM(
            ref_view->path_internal_colorspace_out_ppm,
            tmp_w,
            tmp_r,
            tmp_ncomp,
            ref_view->color);

        /* FORWARD warp color */
        warpView0_to_View1(
            ref_view,
            SAI,
            warped_texture_views[ij],
            warped_depth_views[ij],
            DispTargs[ij]);

        delete[](ref_view->depth);
        delete[](ref_view->color);

        ref_view->depth = nullptr;
        ref_view->color = nullptr;

        if (SAVE_PARTIAL_WARPED_VIEWS) {

            char tmp_str[1024];

            sprintf(
                tmp_str,
                "%s/%03d_%03d_warped_to_%03d_%03d.ppm",
                setup.output_directory.c_str(),
                (ref_view)->c,
                (ref_view)->r,
                SAI->c,
                SAI->r);

            aux_write16PGMPPM(
                tmp_str, 
                SAI->nc, 
                SAI->nr, 3, 
                warped_texture_views[ij]);

            sprintf(
                tmp_str, 
                "%s/%03d_%03d_warped_to_%03d_%03d.pgm",
                setup.output_directory.c_str(), 
                (ref_view)->c, 
                (ref_view)->r,
                SAI->c, 
                SAI->r);

            aux_write16PGMPPM(
                tmp_str, 
                SAI->nc, 
                SAI->nr, 
                1, 
                warped_depth_views[ij]);

        }

    }


}


void decoder::merge_texture_views(
    view *SAI,
    view *LF,
    uint16_t **warped_texture_views,
    float **DispTargs) {

    initViewW(SAI, DispTargs);

    if (SAI->mmode == 0) {

        for (int32_t icomp = 0; icomp < SAI->ncomp; icomp++) {
            mergeWarped_N_icomp(
                warped_texture_views,
                SAI,
                icomp);
        }

        /* hole filling for texture*/
        for (int32_t icomp = 0; icomp < SAI->ncomp; icomp++) {
            uint32_t nholes = holefilling(
                SAI->color + icomp*SAI->nr*SAI->nc,
                SAI->nr,
                SAI->nc,
                (uint16_t)0,
                SAI->seg_vp);
        }
    }

    if (SAI->mmode == 1) {
        /* we don't use LS weights but something derived on geometric distance in view array*/
        for (int32_t icomp = 0; icomp < SAI->ncomp; icomp++) {
            getGeomWeight_icomp(
                SAI,
                LF,
                icomp);
        }

        /* merge color with prediction */
        for (int32_t icomp = 0; icomp < SAI->ncomp; icomp++) {
            mergeWarped_N_icomp(
                warped_texture_views,
                SAI,
                icomp);
        }

        /* hole filling for texture*/
        for (int32_t icomp = 0; icomp < 3; icomp++) {
            uint32_t nholes = holefilling(
                SAI->color + icomp*SAI->nr*SAI->nc,
                SAI->nr,
                SAI->nc,
                (uint16_t)0,
                SAI->seg_vp);
        }
    }

    if (SAI->mmode == 2) {

        mergeMedian_N(warped_texture_views, DispTargs, SAI, 3);

        /* hole filling for texture*/
        for (int32_t icomp = 0; icomp < 3; icomp++) {
            uint32_t nholes = holefilling(
                SAI->color + icomp*SAI->nr*SAI->nc,
                SAI->nr,
                SAI->nc,
                (uint16_t)0,
                SAI->seg_vp);
        }

    }

    delete[](SAI->seg_vp);
    SAI->seg_vp = nullptr;
    delete[](SAI->bmask);
    SAI->bmask = nullptr;
}

void decoder::predict_texture_view(view* SAI) {

    if (SAI->n_references > 0) {

        printf("Predicting texture for view %03d_%03d\n", SAI->c, SAI->r);

        uint16_t **warped_texture_views = nullptr;
        uint16_t **warped_depth_views = nullptr;
        float **DispTargs = nullptr;

        init_warping_arrays(
            SAI->n_references,
            warped_texture_views,
            warped_depth_views,
            DispTargs,
            SAI->nr,
            SAI->nc,
            SAI->ncomp);

        forward_warp_texture_references(
            LF,
            SAI,
            warped_texture_views,
            warped_depth_views,
            DispTargs);

        merge_texture_views(
            SAI,
            LF,
            warped_texture_views,
            DispTargs);

        clean_warping_arrays(
            SAI->n_references,
            warped_texture_views,
            warped_depth_views,
            DispTargs);

        if (SAI->use_global_sparse) {

            uint16_t *sp_filtered_image_padded =
                new uint16_t[(SAI->nr + 2 * SAI->NNt)*(SAI->nc + 2 * SAI->NNt)*SAI->ncomp]();

            for (int32_t icomp = 0; icomp < SAI->ncomp; icomp++) {

                dequantize_and_reorder_spfilter(
                    SAI->sparse_filters.at(icomp));

                uint16_t *padded_icomp_sai =
                    padArrayUint16_t(SAI->color + SAI->nr*SAI->nc*icomp,
                        SAI->nr,
                        SAI->nc,
                        SAI->NNt);

                std::vector<double> filtered_icomp = applyGlobalSparseFilter(
                    padded_icomp_sai,
                    SAI->nr + 2 * SAI->NNt,
                    SAI->nc + 2 * SAI->NNt,
                    SAI->Ms,
                    SAI->NNt,
                    SPARSE_BIAS_TERM,
                    SAI->sparse_filters.at(icomp).filter_coefficients);

                delete[](padded_icomp_sai);

                for (int32_t iij = 0; iij < (SAI->nr + 2 * SAI->NNt)*(SAI->nc + 2 * SAI->NNt); iij++) {

                    double mmax = static_cast<double>((1 << BIT_DEPTH) - 1);

                    filtered_icomp[iij] =
                        clip(filtered_icomp[iij], 0.0, mmax);

                    sp_filtered_image_padded[iij + (SAI->nr + 2 * SAI->NNt)*(SAI->nc + 2 * SAI->NNt)*icomp] =
                        static_cast<uint16_t>(floor(filtered_icomp[iij] + 0.5));

                }

                uint16_t *cropped_icomp =
                    cropImage(sp_filtered_image_padded + (SAI->nr + 2 * SAI->NNt)*(SAI->nc + 2 * SAI->NNt)*icomp,
                    (SAI->nr + 2 * SAI->NNt),
                        (SAI->nc + 2 * SAI->NNt),
                        SAI->NNt);

                memcpy(
                    SAI->color + SAI->nr*SAI->nc*icomp,
                    cropped_icomp,
                    sizeof(uint16_t)*SAI->nr*SAI->nc);

                delete[](cropped_icomp);

            }

            delete[](sp_filtered_image_padded);
        }
    }
}

void decoder::decode_views() {

    LF = new view[number_of_views]();

    int32_t ii = 0; /*view index*/

    while (ii < number_of_views) {

        view *SAI = LF + ii;
        ii++;

        initView(SAI);

        SAI->nr = number_of_rows;
        SAI->nc = number_of_columns;

        SAI->colorspace = colorspace_LF;

        if (minimum_depth > 0) {
            SAI->min_inv_d = static_cast<int32_t>(minimum_depth);
        }

        minimal_config mconf;

        codestreamToViewHeader(
            n_bytes_prediction,
            SAI,
            input_LF,
            mconf);

        setPaths(
            SAI,
            "",
            setup.output_directory.c_str());

        if (feof(input_LF)) {
            printf("File reading error. Terminating\t...\n");
            exit(0);
        }

        printf("Decoding view %03d_%03d\n", SAI->c, SAI->r);

        SAI->color = new uint16_t[SAI->nr * SAI->nc * 3]();
        SAI->depth = new uint16_t[SAI->nr * SAI->nc]();

        /*main texture prediction here*/
        predict_texture_view(SAI);

        /*extract residuals from codestream*/
        if (SAI->has_color_residual) {

            printf("Decoding texture residual for view %03d_%03d\n", SAI->c, SAI->r);

            readResidualFromDisk(
                SAI->jp2_residual_path_jp2,
                n_bytes_residual,
                input_LF,
                JP2_dict);
        }

        if (SAI->has_depth_residual) {

            printf("Decoding normalized disparity residual for view %03d_%03d\n", SAI->c, SAI->r);

            readResidualFromDisk(
                SAI->jp2_residual_depth_path_jp2,
                n_bytes_residual,
                input_LF,
                JP2_dict);
        }

        /* apply texture residual */
        if (SAI->has_color_residual) {

            int32_t Q = 1;
            int32_t offset = 0;

            if (SAI->level > 1) {
                Q = 2;
                const int32_t bpc = 10;
                offset = (1 << bpc) - 1; /* 10bit images currently */
            }

            uint16_t *decoded_residual_image = decode_residual_JP2(
                SAI->path_raw_texture_residual_at_decoder_ppm,
                (setup.wasp_kakadu_directory + "/kdu_expand").c_str(),
                SAI->jp2_residual_path_jp2);

            aux_write16PGMPPM(
                SAI->path_raw_texture_residual_at_decoder_ppm,
                SAI->nc,
                SAI->nr,
                3,
                decoded_residual_image);

            double *residual = dequantize_residual(
                decoded_residual_image,
                SAI->nr,
                SAI->nc,
                SAI->ncomp,
                10,
                Q,
                offset);

            uint16_t *corrected = apply_residual(
                SAI->color,
                residual,
                SAI->nr,
                SAI->nc,
                SAI->ncomp,
                10);

            /* update SAI->color to contain
            corrected (i.e., prediction + residual) version*/
            memcpy(
                SAI->color,
                corrected,
                sizeof(uint16_t)*SAI->nr*SAI->nc*SAI->ncomp);

            delete[](residual);
            delete[](corrected);
            delete[](decoded_residual_image);

        }

        if (SAI->has_depth_residual) {

            delete[](SAI->depth);

            /* has JP2 encoded depth */

            decodeKakadu(
                SAI->path_out_pgm,
                (setup.wasp_kakadu_directory + "/kdu_expand").c_str(),
                SAI->jp2_residual_depth_path_jp2);

            int32_t nr1, nc1, ncomp1;

            aux_read16PGMPPM(
                SAI->path_out_pgm,
                nc1,
                nr1,
                ncomp1,
                SAI->depth);

        }
        else {

            /*inverse depth prediction*/

            if (SAI->level < maxh) {
                WaSP_predict_depth(SAI, LF);
            }

        }

        if (MEDFILT_DEPTH) {

            uint16_t *filtered_depth = medfilt2D(
                SAI->depth,
                3,
                SAI->nr,
                SAI->nc);

            memcpy(
                SAI->depth,
                filtered_depth,
                sizeof(uint16_t) * SAI->nr * SAI->nc);

            delete[](filtered_depth);

        }

        /*internal colorspace version*/

        aux_write16PGMPPM(
            SAI->path_internal_colorspace_out_ppm,
            SAI->nc,
            SAI->nr,
            SAI->ncomp,
            SAI->color);

        /*colorspace transformation back to input*/

        write_output_ppm(
            SAI->color,
            SAI->path_out_ppm,
            SAI->nr,
            SAI->nc,
            SAI->ncomp,
            10,
            SAI->colorspace);

        /*write inverse depth .pgm*/
        if (SAI->level < maxh) {

            aux_write16PGMPPM(
                SAI->path_out_pgm,
                SAI->nc,
                SAI->nr,
                1,
                SAI->depth);
        }

        if (SAI->color != nullptr) {
            delete[](SAI->color);
            SAI->color = nullptr;
        }

        if (SAI->depth != nullptr) {
            delete[](SAI->depth);
            SAI->depth = nullptr;
        }

        if (SAI->seg_vp != nullptr) {
            delete[](SAI->seg_vp);
            SAI->seg_vp = nullptr;
        }

    }
}

void decoder::dealloc() {

    for (int32_t ii = 0; ii < number_of_views; ii++) {
        view* SAI = LF + ii;
        if (SAI->color != nullptr)
            delete[](SAI->color);

        if (SAI->depth != nullptr)
            delete[](SAI->depth);

        if (SAI->references != nullptr)
            delete[](SAI->references);

        if (SAI->depth_references != nullptr)
            delete[](SAI->depth_references);

        if (SAI->merge_weights != nullptr)
            delete[](SAI->merge_weights);

        if (SAI->sparse_weights != nullptr)
            delete[](SAI->sparse_weights);

        if (SAI->bmask != nullptr)
            delete[](SAI->bmask);

        if (SAI->seg_vp != nullptr)
            delete[](SAI->seg_vp);

        if (SAI->sparse_mask != nullptr)
            delete[](SAI->sparse_mask);
    }
    delete[](LF);
}
