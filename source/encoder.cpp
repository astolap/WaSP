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

#include <fstream>
#include <iomanip>

#include "encoder.hh"
#include "ppm.hh"
#include "fileaux.hh"
#include "codestream.hh"
#include "residualjp2.hh"
#include "clip.hh"
#include "json.hh"
#include "view.hh"
#include "merging.hh"
#include "medianfilter.hh"
#include "inpainting.hh"
#include "predictdepth.hh"
#include "ycbcr.hh"
#include "psnr.hh"
#include "warping.hh"
#include "bitdepth.hh"

encoder::encoder(const WaSPsetup encoder_setup)
{

    setup = encoder_setup;

    load_config_json(setup.config_file);

}

encoder::~encoder() {
}

void encoder::encode() {

    generate_normalized_disparity();
    generate_texture();
    write_bitstream();

}

void encoder::write_config(string config_json_file_out) {

    nlohmann::json conf_out;

    conf_out["type"] = "WaSP";

    conf_out["number_of_views_to_be_encoded"] = n_views_total;
    conf_out["enable_fixed_weight_parameter_search"] = std_search_s;

    conf_out["image_height"] = nr;
    conf_out["image_width"] = nc;

    conf_out["colorspace"] = colorspace_LF;

    vector<nlohmann::json::object_t> views;

    for (int32_t ii = 0; ii < n_views_total; ii++) {

        view *SAI = LF + ii;

        nlohmann::json view_configuration;

        view_configuration["id"] = SAI->i_order;

        view_configuration["level"] = SAI->level;

        view_configuration["row_index"] = SAI->r;
        view_configuration["column_index"] = SAI->c;

        view_configuration["cweight_search"] = SAI->cweight_search;

        view_configuration["rate_texture"] = SAI->residual_rate_color,
            view_configuration["rate_inverse_depth"] = SAI->residual_rate_depth;

        view_configuration["view_merging_mode"] = SAI->mmode;

        view_configuration["horizontal_camera_center_position"] = SAI->x;
        view_configuration["vertical_camera_center_position"] = SAI->y;

        view_configuration["sparse_filter_order"] = SAI->Ms;
        view_configuration["sparse_filter_neighborhood_size"] = SAI->NNt;

        view_configuration["fixed_merging_weight_parameter"] = SAI->stdd;

        view_configuration["minimum_inverse_depth"] = SAI->min_inv_d;

        view_configuration["number_of_texture_references"] = SAI->n_references;
        view_configuration["number_of_inverse_depth_references"] =
            SAI->n_depth_references;

        std::vector<int32_t> sai_texture_references;
        for (int32_t ij = 0; ij < SAI->n_references; ij++) {
            sai_texture_references.push_back(SAI->references[ij]);
        }
        view_configuration["texture_reference_indices"] =
            sai_texture_references;

        std::vector<int32_t> sai_inverse_depth_references;
        for (int32_t ij = 0; ij < SAI->n_depth_references; ij++) {
            sai_inverse_depth_references.push_back(SAI->depth_references[ij]);
        }
        view_configuration["inverse_depth_reference_indices"] =
            sai_inverse_depth_references;


        views.push_back(view_configuration);

    }

    conf_out["views"] = views;

    std::ofstream file(config_json_file_out);
    file << std::setw(2) << conf_out << std::endl;
}

void encoder::load_config_json(string config_json_file) {

    ifstream ifs(config_json_file);
    nlohmann::json conf = nlohmann::json::parse(ifs);

    n_views_total = conf["number_of_views_to_be_encoded"].get<int32_t>();
    std_search_s = conf["enable_fixed_weight_parameter_search"].get<int32_t>();

    STD_SEARCH = std_search_s > 0 ? true : false;

    nr = conf["image_height"].get<int32_t>();
    nc = conf["image_width"].get<int32_t>();

    vector<nlohmann::json::object_t> conf_views =
        conf["views"].get<vector<nlohmann::json::object_t>>();

    colorspace_LF = conf["colorspace"].get<std::string>();

    LF = new view[n_views_total]();
    for (nlohmann::json::object_t view_configuration : conf_views) {

        int32_t ii = view_configuration["id"].get<int32_t>();
        view *SAI = LF + ii;

        SAI->nr = nr;
        SAI->nc = nc;

        SAI->colorspace = colorspace_LF;

        initView(SAI);
        SAI->i_order = ii;

        SAI->ncomp = 3;

        SAI->r = view_configuration["row_index"].get<int32_t>();
        SAI->c = view_configuration["column_index"].get<int32_t>();

        SAI->cweight_search =
            view_configuration["cweight_search"].get<int32_t>() > 0 ? true : false;

        SAI->residual_rate_color =
            view_configuration["rate_texture"].get<float>();
        SAI->residual_rate_depth =
            view_configuration["rate_inverse_depth"].get<float>();

        SAI->mmode = view_configuration["view_merging_mode"].get<unsigned char>();

        SAI->x = view_configuration["horizontal_camera_center_position"].get<float>();
        SAI->y = view_configuration["vertical_camera_center_position"].get<float>();

        SAI->Ms = view_configuration["sparse_filter_order"].get<int32_t>();

        SAI->NNt = view_configuration["sparse_filter_neighborhood_size"].get<int32_t>();

        SAI->stdd = view_configuration["fixed_merging_weight_parameter"].get<float>();

        SAI->min_inv_d = view_configuration["minimum_inverse_depth"].get<int32_t>();
        SAI->n_references =
            view_configuration["number_of_texture_references"].get<int32_t>();
        SAI->n_depth_references =
            view_configuration["number_of_inverse_depth_references"].get<int32_t>();

        SAI->level = view_configuration["level"].get<int32_t>();

        if (abs(SAI->x) > 0.0001) {
            SAI->has_x_displacement = true;
        }

        if (abs(SAI->y) > 0.0001) {
            SAI->has_y_displacement = true;
        }

        if (SAI->n_references > 0) {

            SAI->has_color_references = true;

            SAI->references = new int32_t[SAI->n_references]();

            vector<int32_t> texture_references =
                view_configuration["texture_reference_indices"].get<vector<int32_t>>();

            memcpy(
                SAI->references,
                &texture_references[0],
                sizeof(int32_t)*SAI->n_references);

        }

        if (SAI->n_depth_references > 0) {

            SAI->has_depth_references = true;

            SAI->depth_references = new int32_t[SAI->n_depth_references]();

            vector<int32_t> depth_references =
                view_configuration["inverse_depth_reference_indices"].get<vector<int32_t>>();

            memcpy(
                SAI->depth_references,
                &depth_references[0],
                sizeof(int32_t)*SAI->n_depth_references);

        }

        setPaths(
            SAI,
            setup.input_directory.c_str(),
            setup.output_directory.c_str());

    }
}

void encoder::forward_warp_texture_references(
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

            sprintf(tmp_str, "%s/%03d_%03d_warped_to_%03d_%03d.ppm",
                setup.output_directory.c_str(), (ref_view)->c, (ref_view)->r,
                SAI->c, SAI->r);

            aux_write16PGMPPM(tmp_str, SAI->nc, SAI->nr, 3, warped_texture_views[ij]);

            sprintf(tmp_str, "%s/%03d_%03d_warped_to_%03d_%03d.pgm",
                setup.output_directory.c_str(), (ref_view)->c, (ref_view)->r,
                SAI->c, SAI->r);

            aux_write16PGMPPM(tmp_str, SAI->nc, SAI->nr, 1, warped_depth_views[ij]);

        }

    }


}

void encoder::merge_texture_views(
    view *SAI,
    view *LF,
    uint16_t **warped_texture_views,
    float **DispTargs) {

    initViewW(SAI, DispTargs);

    if (SAI->mmode == 0) {

        int32_t nr1, nc1, ncomp1;

        uint16_t *original_color_view = read_input_ppm(
            SAI->path_input_ppm,
            nr1,
            nc1,
            ncomp1,
            10,
            SAI->colorspace);

        for (int32_t icomp = 0; icomp < SAI->ncomp; icomp++) {
            getViewMergingLSWeights_icomp(
                SAI,
                warped_texture_views,
                original_color_view,
                icomp);
        }

        delete[](original_color_view);

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

        /*merge with median operator*/
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

void encoder::generate_normalized_disparity() {

    /*the term inverse depth is used interchangeably with normalized disparity */

    int32_t maxh = get_highest_level(LF, n_views_total);

    /*we dont have to do the last level*/
    for (int32_t hlevel = 1; hlevel < maxh; hlevel++) {

        std::vector< int32_t > view_indices;

        for (int32_t ii = 0; ii < n_views_total; ii++) {
            if ((LF + ii)->level == hlevel) {
                view_indices.push_back(ii);
            }
        }

        /*ascending order of view index at a particular level*/
        sort(view_indices.begin(), view_indices.end());

        for (int32_t ii = 0; ii < view_indices.size(); ii++) {

            view *SAI = LF + view_indices.at(ii);

            SAI->depth = new uint16_t[SAI->nr * SAI->nc]();

            if (hlevel == 1) { /*intra coding of inverse depth*/

                bool depth_file_exist = false;

                int32_t nc1, nr1, ncomp1;

                delete[](SAI->depth);
                SAI->depth = nullptr;

                SAI->depth_file_exist = aux_read16PGMPPM(
                    SAI->path_input_pgm,
                    nc1,
                    nr1,
                    ncomp1,
                    SAI->depth);

                if (SAI->depth_file_exist && SAI->residual_rate_depth > 0) {

                    /* ------------------------------
                    INVERSE DEPTH ENCODING STARTS
                    -------------------------------*/

                    printf("Encoding normalized disparity for view %03d_%03d\n", SAI->c, SAI->r);

                    /*write inverse depth to .pgm (path SAI->path_out_pgm) */
                    aux_write16PGMPPM(
                        SAI->path_out_pgm,
                        SAI->nc,
                        SAI->nr,
                        1,
                        SAI->depth);

                    delete[](SAI->depth);
                    SAI->depth = nullptr;

                    char *oparams = kakadu_oparams(
                        SAI->residual_rate_depth,
                        "YCbCr"); /*cycc shouldn't matter for single-channel*/

                    char *encoding_parameters = new char[65535]();
                    sprintf(
                        encoding_parameters,
                        "%s",
                        oparams);

                    encodeKakadu(
                        SAI->path_out_pgm,
                        (setup.wasp_kakadu_directory + "/kdu_compress").c_str(),
                        SAI->jp2_residual_depth_path_jp2,
                        encoding_parameters,
                        SAI->residual_rate_depth);

                    delete[](encoding_parameters);
                    /* ------------------------------
                    INVERSE DEPTH ENCODING ENDS
                    ------------------------------*/

                    /* ------------------------------
                    INVERSE DEPTH DECODING STARTS
                    ------------------------------*/

                    printf("Decoding normalized disparity for view %03d_%03d\n", SAI->c, SAI->r);

                    decodeKakadu(
                        SAI->path_out_pgm,
                        (setup.wasp_kakadu_directory + "/kdu_expand").c_str(),
                        SAI->jp2_residual_depth_path_jp2);

                    /*------------------------------
                    INVERSE DEPTH DECODING ENDS
                    ------------------------------*/

                    SAI->has_depth_residual = true;

                    aux_read16PGMPPM(
                        SAI->path_out_pgm,
                        nc1,
                        nr1,
                        ncomp1,
                        SAI->depth);

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

                    //aux_write16PGMPPM(
                    //    SAI->path_out_pgm,
                    //    nc1,
                    //    nr1,
                    //    ncomp1,
                    //    SAI->depth);

                }


            }
            else { /*prediction only*/

                printf("Predicting normalized disparity for view %03d_%03d\n", SAI->c, SAI->r);

                /*color array needed in warping function ...*/
                SAI->color = new uint16_t[SAI->nr * SAI->nc * 3]();
                /*SAI->depth = new uint16_t[SAI->nr * SAI->nc]();*/

                WaSP_predict_depth(SAI, LF);

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

                delete[](SAI->color);
                SAI->color = nullptr;

            }


            aux_write16PGMPPM(
                SAI->path_out_pgm,
                SAI->nc,
                SAI->nr,
                1,
                SAI->depth);

            delete[](SAI->depth);
            SAI->depth = nullptr;

        }
    }
}

void encoder::generate_texture() {

    //FILE *tmp;
    //tmp = fopen("C:/Temp/coeffs.data", "wb");
    //int32_t qvalc = (1<<BIT_DEPTH_SPARSE);
    //fwrite(&qvalc, sizeof(int32_t), 1, tmp);
    //fclose(tmp);

    maxh = get_highest_level(LF, n_views_total);

    for (int32_t hlevel = 1; hlevel <= maxh; hlevel++) {

        printf("\n\tProcessing of hierarchical level: %d\n\n", hlevel);

        const int32_t bpc = 10;

        int32_t Q = 1;
        int32_t offset = 0;

        if (hlevel > 1) {
            Q = 2;

            offset = (1 << bpc) - 1; /* 10bit images currently */
        }

        std::vector< int32_t > view_indices;

        for (int32_t ii = 0; ii < n_views_total; ii++) {
            if ((LF + ii)->level == hlevel) {
                view_indices.push_back(ii);
            }
        }

        /*ascending order of view index at level=hlevel*/
        sort(view_indices.begin(), view_indices.end());

        /* predict (i.e., warp and merge) all views at level=hlevel*/
        for (int32_t ii = 0; ii < view_indices.size(); ii++) {

            view *SAI = LF + view_indices.at(ii);

            printf("Encoding view %03d_%03d\n", SAI->c, SAI->r);

            SAI->color = new uint16_t[SAI->nr * SAI->nc * 3]();

            if (SAI->n_references > 0) {

                printf("View prediction for view %03d_%03d\n", SAI->c, SAI->r);

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

                if (SAI->Ms > 0 && SAI->NNt > 0) {

                    printf("Sparse filtering stage for view %03d_%03d\n", SAI->c, SAI->r);

                    uint16_t *original_color_view = read_input_ppm(
                        SAI->path_input_ppm,
                        SAI->nr,
                        SAI->nc,
                        SAI->ncomp,
                        bpc,
                        SAI->colorspace);

                    SAI->sparse_filters.clear();

                    for (int32_t icomp = 0; icomp < SAI->ncomp; icomp++) {

                        uint16_t *padded_icomp_sai =
                            padArrayUint16_t(SAI->color + SAI->nr*SAI->nc*icomp,
                                SAI->nr,
                                SAI->nc,
                                SAI->NNt);

                        uint16_t *padded_icomp_orig =
                            padArrayUint16_t(original_color_view + SAI->nr*SAI->nc*icomp,
                                SAI->nr,
                                SAI->nc,
                                SAI->NNt);

                        //aux_write16PGMPPM(
                        //    "C:/Temp/padded.pgm",
                        //    SAI->nc + 2 * SAI->NNt,
                        //    SAI->nr + 2 * SAI->NNt,
                        //    1,
                        //    padded_icomp_sai);

                        SAI->sparse_filters.push_back(getGlobalSparseFilter(
                            padded_icomp_orig,
                            padded_icomp_sai,
                            SAI->nr + 2 * SAI->NNt,
                            SAI->nc + 2 * SAI->NNt,
                            SAI->NNt,
                            SAI->Ms,
                            SPARSE_BIAS_TERM,
                            setup.sparse_subsampling));

                        //SAI->sparse_filters.push_back(getSP_FILTER_EIGEN(
                        //    padded_icomp_orig,
                        //    padded_icomp_sai,
                        //    SAI->nr + 2 * SAI->NNt,
                        //    SAI->nc + 2 * SAI->NNt,
                        //    SAI->NNt,
                        //    SAI->Ms,
                        //    SPARSE_BIAS_TERM,
                        //    setup.sparse_subsampling));

                        //FILE *tmp;
                        //tmp = fopen("C:/Temp/coeffs.data", "ab");

                        //fwrite(&SAI->sparse_filters.at(icomp).filter_coefficients[icomp],
                        //    sizeof(double),
                        //    SAI->sparse_filters.at(icomp).filter_coefficients.size(),
                        //    tmp);

                        //fclose(tmp);

                      /*  aux_write16PGMPPM("C:/Temp/padded_sai.pgm", SAI->nc + 2 * SAI->NNt, SAI->nr + 2 * SAI->NNt, 1, padded_icomp_sai);
                        aux_write16PGMPPM("C:/Temp/padded_orig.pgm", SAI->nc + 2 * SAI->NNt, SAI->nr + 2 * SAI->NNt, 1, padded_icomp_orig);*/

                        delete[](padded_icomp_sai);
                        delete[](padded_icomp_orig);

                        /* exit(0);*/
                    }

                    uint16_t *sp_filtered_image_padded =
                        new uint16_t[(SAI->nr + 2 * SAI->NNt)*(SAI->nc + 2 * SAI->NNt)*SAI->ncomp]();

                    uint16_t *sp_filtered_image =
                        new uint16_t[SAI->nr*SAI->nc*SAI->ncomp]();

                    for (int32_t icomp = 0; icomp < SAI->ncomp; icomp++) {

                        quantize_and_reorder_spfilter(
                            SAI->sparse_filters.at(icomp));

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
                            sp_filtered_image + SAI->nr*SAI->nc*icomp,
                            cropped_icomp,
                            sizeof(uint16_t)*SAI->nr*SAI->nc);

                        delete[](cropped_icomp);

                    }

                    delete[](sp_filtered_image_padded);

                    double psnr_without_sparse = PSNR(
                        original_color_view,
                        SAI->color,
                        SAI->nr,
                        SAI->nc,
                        SAI->ncomp,
                        (1 << bpc) - 1);

                    double psnr_with_sparse = PSNR(
                        original_color_view,
                        sp_filtered_image,
                        SAI->nr,
                        SAI->nc,
                        SAI->ncomp,
                        (1 << bpc) - 1);

                    if(1){//psnr_with_sparse > psnr_without_sparse) {

                        memcpy(
                            SAI->color,
                            sp_filtered_image,
                            sizeof(uint16_t)*SAI->nr*SAI->nc*SAI->ncomp);

                        SAI->use_global_sparse = true;

                    }
                    else {
                        SAI->use_global_sparse = false; /*should never happen. implies something wrong with filter design or filter coefficient quantization*/
                    }

                    delete[](sp_filtered_image);
                    delete[](original_color_view);

                }

            }

            /* write raw prediction to .ppm */

            aux_write16PGMPPM(
                SAI->path_raw_prediction_at_encoder_ppm,
                SAI->nc,
                SAI->nr,
                3,
                SAI->color);

            delete[](SAI->color);
            SAI->color = nullptr;

        }

        /* get residue for all views at level=hlevel,
        AFTER THIS LOOP,
        YOU CAN FIND ALL RESIDUAL IMAGES IN directories "outputdir/residual/RAW/<level>"
        */
        for (int32_t ii = 0; ii < view_indices.size(); ii++) {

            view *SAI = LF + view_indices.at(ii);

            if (SAI->residual_rate_color > 0) {

                printf("Obtaining texture residual for view %03d_%03d\n", SAI->c, SAI->r);

                aux_read16PGMPPM(
                    SAI->path_raw_prediction_at_encoder_ppm,
                    SAI->nc,
                    SAI->nr,
                    SAI->ncomp,
                    SAI->color);

                uint16_t *original_color_view = read_input_ppm(
                    SAI->path_input_ppm,
                    SAI->nr,
                    SAI->nc,
                    SAI->ncomp,
                    bpc,
                    SAI->colorspace);

                double *residual_image_double = get_residual(
                    original_color_view,
                    SAI->color,
                    SAI->nr,
                    SAI->nc,
                    3);

                delete[](original_color_view);

                SAI->residual_image = quantize_residual(
                    residual_image_double,
                    SAI->nr,
                    SAI->nc,
                    SAI->ncomp,
                    bpc,
                    Q,
                    offset);

                delete[](residual_image_double);

                /* write raw quantized residual to .ppm */

                aux_write16PGMPPM(
                    SAI->path_raw_texture_residual_at_encoder_ppm,
                    SAI->nc,
                    SAI->nr,
                    SAI->ncomp,
                    SAI->residual_image);

                delete[](SAI->color);
                SAI->color = nullptr;
                delete[](SAI->residual_image);
                SAI->residual_image = nullptr;

            }

        }

        /* encode residual images using JP2 for all views at level=hlevel
        here we can substitute JP2 with MuLE etc*/
        for (int32_t ii = 0; ii < view_indices.size(); ii++) {

            view *SAI = LF + view_indices.at(ii);

            if (SAI->residual_rate_color > 0) {

                /* ------------------------------
                TEXTURE RESIDUAL ENCODING STARTS
                ------------------------------*/

                printf("Encoding texture residual for view %03d_%03d\n", SAI->c, SAI->r);

                std::vector< std::pair<double, double> > psnr_cweight;

                uint16_t *original_color_view = read_input_ppm(
                    SAI->path_input_ppm,
                    SAI->nr,
                    SAI->nc,
                    SAI->ncomp,
                    bpc,
                    SAI->colorspace);

                aux_read16PGMPPM(
                    SAI->path_raw_prediction_at_encoder_ppm,
                    SAI->nc,
                    SAI->nr,
                    SAI->ncomp,
                    SAI->color);

                aux_read16PGMPPM(
                    SAI->path_raw_texture_residual_at_encoder_ppm,
                    SAI->nc,
                    SAI->nr,
                    SAI->ncomp,
                    SAI->residual_image);

                /* look for best allocation of rate component wise,
                optimizes for the weighted PSNR */
                if (SAI->cweight_search) {
                    for (double cweight = 1; cweight < 20; cweight += 2) {

                        std::vector<double> cweights = { cweight,1.0,1.0 };

                        char *cparams = kakadu_cparams(
                            &cweights[0],
                            SAI->ncomp);

                        char *oparams = kakadu_oparams(
                            SAI->residual_rate_color,
                            colorspace_LF);

                        char *encoding_parameters = new char[65535]();
                        sprintf(
                            encoding_parameters,
                            "%s %s",
                            cparams,
                            oparams);

                        encode_residual_JP2(
                            SAI->path_raw_texture_residual_at_encoder_ppm,
                            (setup.wasp_kakadu_directory + "/kdu_compress").c_str(),
                            SAI->jp2_residual_path_jp2,
                            encoding_parameters,
                            SAI->residual_rate_color);

                        uint16_t *decoded_residual_image = decode_residual_JP2(
                            SAI->path_raw_texture_residual_at_decoder_ppm,
                            (setup.wasp_kakadu_directory + "/kdu_expand").c_str(),
                            SAI->jp2_residual_path_jp2);

                        double *residual = dequantize_residual(
                            decoded_residual_image,
                            SAI->nr,
                            SAI->nc,
                            SAI->ncomp,
                            bpc,
                            Q,
                            offset);

                        uint16_t *corrected = apply_residual(
                            SAI->color,
                            residual,
                            SAI->nr,
                            SAI->nc,
                            SAI->ncomp,
                            bpc);

                        double psnr_c = getYCbCr_444_PSNR(
                            corrected,
                            original_color_view,
                            SAI->nr,
                            SAI->nc,
                            SAI->ncomp,
                            bpc);

                        delete[](residual);
                        delete[](corrected);
                        delete[](decoded_residual_image);

                        delete[](encoding_parameters);
                        delete[](oparams);
                        delete[](cparams);

                        psnr_cweight.push_back(
                            std::pair<double, double>(psnr_c, cweight));

                    }
                }
                else {
                    psnr_cweight.push_back(
                        std::pair<double, double>(0.0, 1.0));
                }

                delete[](original_color_view);

                sort(psnr_cweight.begin(), psnr_cweight.end());

                double best_cweight =
                    psnr_cweight.at(psnr_cweight.size() - 1).second;

                /* final encoding with best cweight */
                std::vector<double> cweights = { best_cweight,1.0,1.0 };

                char *cparams = kakadu_cparams(&cweights[0], 3);
                char *oparams = kakadu_oparams(
                    SAI->residual_rate_color,
                    colorspace_LF);

                char *encoding_parameters = new char[65535]();
                sprintf(
                    encoding_parameters,
                    "%s%s",
                    cparams,
                    oparams);

                encode_residual_JP2(
                    SAI->path_raw_texture_residual_at_encoder_ppm,
                    (setup.wasp_kakadu_directory + "/kdu_compress").c_str(),
                    SAI->jp2_residual_path_jp2,
                    encoding_parameters,
                    SAI->residual_rate_color);

                delete[](encoding_parameters);
                delete[](oparams);
                delete[](cparams);

                SAI->has_color_residual = true;

                delete[](SAI->color);
                SAI->color = nullptr;
                delete[](SAI->residual_image);
                SAI->residual_image = nullptr;

                /* ------------------------------
                TEXTURE RESIDUAL ENCODING ENDS
                ------------------------------*/

            }

        }

        /* now decode residual images for all views at level=hlevel,
        and write the result to disk*/
        for (int32_t ii = 0; ii < view_indices.size(); ii++) {

            view *SAI = LF + view_indices.at(ii);

            aux_read16PGMPPM(
                SAI->path_raw_prediction_at_encoder_ppm,
                SAI->nc,
                SAI->nr,
                SAI->ncomp,
                SAI->color);

            if (SAI->has_color_residual) {

                /* ------------------------------
                TEXTURE RESIDUAL DECODING STARTS
                ------------------------------*/

                printf("Decoding texture residual for view %03d_%03d\n", SAI->c, SAI->r);

                uint16_t *decoded_residual_image = decode_residual_JP2(
                    SAI->path_raw_texture_residual_at_decoder_ppm,
                    (setup.wasp_kakadu_directory + "/kdu_expand").c_str(),
                    SAI->jp2_residual_path_jp2);

                /* ------------------------------
                TEXTURE RESIDUAL DECODING ENDS
                ------------------------------*/

                double *residual = dequantize_residual(
                    decoded_residual_image,
                    SAI->nr,
                    SAI->nc,
                    SAI->ncomp,
                    bpc,
                    Q,
                    offset);

                uint16_t *corrected = apply_residual(
                    SAI->color,
                    residual,
                    SAI->nr,
                    SAI->nc,
                    SAI->ncomp,
                    bpc);

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

            /*WRITE result to disk*/

            /*writing .ppm in internal colorspace*/
            aux_write16PGMPPM(
                SAI->path_internal_colorspace_out_ppm,
                SAI->nc,
                SAI->nr,
                SAI->ncomp,
                SAI->color);

            /*writing .ppm in output colorspace*/
            write_output_ppm(
                SAI->color,
                SAI->path_out_ppm,
                SAI->nr,
                SAI->nc,
                SAI->ncomp,
                bpc,
                SAI->colorspace);

            delete[](SAI->color);
            SAI->color = nullptr;

        }

    }

}

void encoder::write_bitstream() {

    printf("Writing header information to codestream\n");

    uint8_t colorspace_enumerator;

    if (colorspace_LF.compare("RGB") == 0) {
        colorspace_enumerator = 0;
    }
    if (colorspace_LF.compare("YCbCr") == 0) {
        colorspace_enumerator = 1;
    }

    char path_out_LF_data[1024];
    sprintf(
        path_out_LF_data,
        "%s/%s",
        setup.output_directory.c_str(),
        "output.LF");

    FILE* output_LF_file = fopen(path_out_LF_data, "wb");

    int32_t n_bytes_prediction = 0;
    int32_t n_bytes_residual = 0;

    n_bytes_prediction += (int32_t)fwrite(
        &n_views_total,
        sizeof(int32_t),
        1,
        output_LF_file);
    n_bytes_prediction += (int32_t)fwrite(
        &LF->nr,
        sizeof(int32_t),
        1,
        output_LF_file)
        * sizeof(int32_t);  // needed only once per LF
    n_bytes_prediction += (int32_t)fwrite(
        &LF->nc, sizeof(int32_t),
        1,
        output_LF_file)
        * sizeof(int32_t);  //
    n_bytes_prediction += (int32_t)fwrite(
        &LF->min_inv_d,
        sizeof(uint16_t),
        1,
        output_LF_file) * sizeof(uint16_t);
    n_bytes_prediction += (int32_t)fwrite(
        &colorspace_enumerator,
        sizeof(uint8_t),
        1,
        output_LF_file) * sizeof(uint8_t);
    n_bytes_prediction += (int32_t)fwrite(
        &maxh,
        sizeof(int32_t),
        1,
        output_LF_file) * sizeof(int32_t);

    for (int32_t ii = 0; ii < n_views_total; ii++) {

        view *SAI = LF + ii;

        printf("Writing codestream for view %03d_%03d\n", SAI->c, SAI->r);

        viewHeaderToCodestream(
            n_bytes_prediction,
            SAI,
            output_LF_file);

        if (SAI->has_color_residual) {
            writeResidualToDisk(
                SAI->jp2_residual_path_jp2,
                output_LF_file,
                n_bytes_residual,
                JP2_dict);
        }

        if (SAI->has_depth_residual) {
            writeResidualToDisk(
                SAI->jp2_residual_depth_path_jp2,
                output_LF_file,
                n_bytes_residual,
                JP2_dict);
        }

    }

    fclose(output_LF_file);
}
