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

#include <string.h>

#include "view.hh"
#include "ppm.hh"
#include "ycbcr.hh"
#include "sparsefilter.hh"

void initView(view* view) {

  view->color = nullptr;
  view->depth = nullptr;
  view->segmentation = nullptr;

  view->r = 0;
  view->c = 0;

  view->ncomp = 3;

  view->mmode = 0;

  view->min_inv_d = 0;
  view->n_references = 0;

  view->references = nullptr;
  view->depth_references = nullptr;

  view->merge_weights = nullptr;
  view->merge_weights_double = nullptr;
  view->sparse_weights = nullptr;
  view->sparse_mask = nullptr;
  view->merge_weights_float = nullptr;
  view->number_of_pixels_per_region = nullptr;

  view->bmask = nullptr;
  view->seg_vp = nullptr;

  view->NNt = 0;
  view->Ms = 0;

  view->has_segmentation = 0;
  view->maxL = 0;

  view->i_order = 0;

  view->stdd = 0.0;

  view->yuv_transform = false;

  view->has_color_residual = false;
  view->has_depth_residual = false;
  view->use_global_sparse = false;

  view->has_color_references = false;
  view->has_depth_references = false;

  view->has_x_displacement = false;
  view->has_y_displacement = false;

  view->level = 0;

  view->depth_file_exist = false;

  view->residual_image = nullptr;

  view->cweight_search = false;

}

int32_t get_highest_level(
    view *LF, 
    const int32_t n_views_total) {

    int32_t maxh = 0;

    for (int32_t ii = 0; ii < n_views_total; ii++) {
        maxh = (LF + ii)->level > maxh ? (LF + ii)->level : maxh;
    }

    return maxh;
}

int32_t* determine_levels(
    view *LF, 
    const int32_t n_views_total) {
    /* Here we obtain the hierchical level of the view.
    We obtain the level by inspecting which views are used as references for a given view.
    Since only lower levels can be used as references, we can infer the level of the current view.*/

    int32_t maxR = 0;
    int32_t maxC = 0;

    for (int32_t ii = 0; ii < n_views_total; ii++) {
        view *SAI = (LF + ii);
        /* find number of rows and cols */
        maxR = (SAI->r + 1 > maxR) ? SAI->r + 1 : maxR;
        maxC = (SAI->c + 1 > maxC) ? SAI->c + 1 : maxC;
    }

    int32_t *hierarchy_matrix = new int32_t[maxR * maxC]();
    for (int32_t ii = 0; ii < n_views_total; ii++) {

        view *SAI = (LF + ii);

        /* here we collect the maximum level appearing in the set of references of SAI*/
        int32_t max_level = 0;

        for (int32_t jj = 0; jj < SAI->n_references; jj++) {
            view *ref_SAI = (LF + SAI->references[jj]);
            int32_t hval = hierarchy_matrix[ref_SAI->r + maxR * ref_SAI->c];
            max_level = hval > max_level ? hval : max_level;
        }

        hierarchy_matrix[SAI->r + maxR * SAI->c] = max_level + 1;
        SAI->level = hierarchy_matrix[SAI->r + maxR * SAI->c];
    }

    return hierarchy_matrix;

    /* print hierarchy matrix for debug */
    //  printf("\n");
    //  for (int32_t rr = 0; rr < maxR; rr++) {
    //    for (int32_t cc = 0; cc < maxC; cc++) {
    //      printf("%i\t", hierarchy_matrix[rr + maxR*cc]);
    //    }
    //    printf("\n");
    //  }
    //  for (int32_t ii = 0; ii < n_views_total; ii++) {
    //    printf("%i\n", (LF+ii)->level);
    //  }
    //  exit(0);
}

void setPaths(
    view *SAI,
    const char *input_dir,
    const char *output_dir){

    /*INPUT*/
    sprintf(
        SAI->path_input_ppm,
        "%s/%03d_%03d.ppm",
        input_dir,
        SAI->c,
        SAI->r);
    /*INPUT*/
    sprintf(
        SAI->path_input_pgm,
        "%s/%03d_%03d.pgm",
        input_dir,
        SAI->c,
        SAI->r);

    /*FINAL OUTPUT*/
    sprintf(
        SAI->path_out_ppm,
        "%s/PPM/%03d_%03d.ppm",
        output_dir,
        SAI->c,
        SAI->r);

    /*FINAL OUTPUT*/
    sprintf(
        SAI->path_out_pgm,
        "%s/PGM/%03d_%03d.pgm",
        output_dir,
        SAI->c,
        SAI->r);

    /*INTERNAL OUTPUT*/
    sprintf(
        SAI->path_internal_colorspace_out_ppm,
        "%s/internal_colorspace/PPM/%01d/%03d_%03d.ppm",
        output_dir,
        SAI->level,
        SAI->c,
        SAI->r);

    sprintf(
        SAI->jp2_residual_path_jp2,
        "%s/residual/JP2/%01d/%03d_%03d.jp2",
        output_dir,
        SAI->level,
        SAI->c,
        SAI->r);

    sprintf(
        SAI->pgm_residual_depth_path,
        "%s/residual/depth/PGM/%01d/%03d_%03d.pgm",
        output_dir,
        SAI->level,
        SAI->c,
        SAI->r);

    sprintf(
        SAI->jp2_residual_depth_path_jp2,
        "%s/residual/depth/JP2/%01d/%03d_%03d.jp2",
        output_dir,
        SAI->level,
        SAI->c,
        SAI->r);

    sprintf(
        SAI->path_raw_texture_residual_at_encoder_ppm,
        "%s/residual/RAW/%01d/%03d_%03d.ppm",
        output_dir,
        SAI->level,
        SAI->c,
        SAI->r);

    sprintf(
        SAI->path_raw_prediction_at_encoder_ppm,
        "%s/prediction/RAW/%01d/%03d_%03d.ppm",
        output_dir,
        SAI->level,
        SAI->c,
        SAI->r);

    //sprintf(
    //    SAI->path_raw_texture_residual_at_decoder_ppm,
    //    "%s/residual/JP2_decoded/%01d/%03d_%03d.ppm",
    //    output_dir,
    //    SAI->level,
    //    SAI->c,
    //    SAI->r);

    /*to save some disk space*/
    sprintf(
        SAI->path_raw_texture_residual_at_decoder_ppm,
        "%s",
        SAI->path_internal_colorspace_out_ppm);

    //sprintf(
    //    SAI->inverse_depth_raw_pgm,
    //    "%s/residual/depth/RAW/%01d/%03d_%03d.pgm",
    //    output_dir,
    //    SAI->level,
    //    SAI->c,
    //    SAI->r);
}

uint16_t *read_input_ppm(
    const char *input_ppm_path,
    int32_t &nr,
    int32_t &nc,
    int32_t &ncomp,
    const int32_t bpc,
    const std::string colorspace) {

    uint16_t *texture_in_input_colorspace = nullptr;

    aux_read16PGMPPM(
        input_ppm_path,
        nc,
        nr,
        ncomp,
        texture_in_input_colorspace);

    /* color transformation here. we now apply YCbCr by default,
    input must be rgb. in the future we can define arbitrary input-output
    color transformations based on the user requirements. */

    uint16_t *texture_in_encoder_colorspace =
        new uint16_t[nr*nc*ncomp]();

    /* we can encode/decode in any colorspace.
    currently we have two options: RGB and YCbCr*/

    if (colorspace.compare("RGB")==0) {
        memcpy(
            texture_in_encoder_colorspace,
            texture_in_input_colorspace,
            sizeof(uint16_t)*nr*nc*ncomp);
    }

    if (colorspace.compare("YCbCr")==0) {
        RGB2YCbCr(
            texture_in_input_colorspace,
            texture_in_encoder_colorspace,
            nr,
            nc,
            bpc);
    }

    delete[](texture_in_input_colorspace);

    return texture_in_encoder_colorspace;
}

void write_output_ppm(
    const uint16_t *texture_view_in_encoder_colorspace,
    const char *output_ppm_path,
    const int32_t nr,
    const int32_t nc,
    const int32_t ncomp,
    const int32_t bpc,
    const std::string colorspace) {

    /* transform back to original input colorspace */
    uint16_t *texture_view_in_output_colorspace =
        new uint16_t[nr*nc*ncomp]();

    /* we can encode/decode in any colorspace.
    currently we have two options: RGB and YCbCr*/

    if (colorspace.compare("RGB")==0) {
        memcpy(
            texture_view_in_output_colorspace,
            texture_view_in_encoder_colorspace,
            sizeof(uint16_t)*nr*nc*ncomp);
    }

    if (colorspace.compare("YCbCr")==0) {
        YCbCr2RGB(
            texture_view_in_encoder_colorspace,
            texture_view_in_output_colorspace,
            nr,
            nc,
            bpc);
    }

    aux_write16PGMPPM(
        output_ppm_path,
        nc,
        nr,
        ncomp,
        texture_view_in_output_colorspace);

    delete[](texture_view_in_output_colorspace);

}