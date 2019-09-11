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

#ifndef VIEW_HH
#define VIEW_HH

#include "sparsefilter.hh"

#include <cstdint>
#include <string>
#include <vector>

using std::int32_t;
using std::uint32_t;

using std::int16_t;
using std::uint16_t;

using std::int8_t;
using std::uint8_t;

#define MEDFILT_DEPTH false

struct view {

  uint16_t *color;
  uint16_t *depth;
  uint16_t *residual_image;

  uint16_t *segmentation;

  uint8_t mmode; /* view merging mode, 0=LS, 1=geometric weight, 2=median*/

  int32_t r, c;  // SAI subscript

  int32_t nr, nc;  // image height, width

  int32_t ncomp; //number of components

  float y, x;  // camera displacement

  int32_t min_inv_d;  // needed only if inverse depth has negative values, [0,max]-mind = [-mind,max-mind]

  int32_t n_references, n_depth_references;

  int32_t *references, *depth_references; /* depth references not necessarily the same as color references
   we can have, for example, depth warping only from the externally obtained depth but we still
   warp color from neighbors that don't have depth provided. We don't want to propagate depth errors from
   badly warped depth views, thus we restrict depth warping to some high quality subset (usually meaning the
   externally obtained high quality depth maps)*/

  int16_t *merge_weights;
  int32_t *sparse_weights;

  uint8_t *sparse_mask;

  float *merge_weights_float;

  double *merge_weights_double;

  int32_t *number_of_pixels_per_region;

  bool *bmask; /* view mask for merging weights */
  uint16_t *seg_vp; /* class segmentation, used for view merging weights */
  int32_t NB;

  bool cweight_search;

  float residual_rate_color;
  float residual_rate_depth;

  float stdd;

  int32_t NNt, Ms;  //for global sparse, NNt defines the neighborhood size [ -NNt:NNt,-NNt:NNt ], Ms is the filter order

  int32_t has_segmentation;
  int32_t maxL;  // number of regions in segmentation

  int32_t ****region_displacements; /* region displacement vectors [iar][iac][iR][xy], e.g., [13][13][25][2], for 13x13 angular views with 25 regions for segmentation */

  char path_input_pgm[1024], path_input_ppm[1024], path_input_seg[1024];
  char path_out_pgm[1024], path_out_ppm[1024];

  //char path_input_Y_pgm[1024], path_out_Y_pgm[1024];
  //char path_input_Cb_pgm[1024], path_out_Cb_pgm[1024];
  //char path_input_Cr_pgm[1024], path_out_Cr_pgm[1024];

  float *DM_ROW, *DM_COL; /* for lenslet with region displacement vectors */

  int32_t i_order; /* view position in encoding configuration */

  bool use_median;  //use median merging or not

  bool yuv_transform;

  bool has_color_residual, has_depth_residual, use_global_sparse;
  bool has_color_references, has_depth_references;
  //bool has_min_inv_depth;

  bool has_x_displacement, has_y_displacement;

  bool has_chrominance;

  int32_t level;  // hierarchical level

  std::string colorspace;

  float yuv_rates[3];
  float rgb_rate;

  bool depth_file_exist;

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

  char path_raw_texture_residual_at_encoder_ppm[1024];
  char path_raw_prediction_at_encoder_ppm[1024];
  char path_raw_texture_residual_at_decoder_ppm[1024];

  char path_internal_colorspace_out_ppm[1024];

  char inverse_depth_raw_pgm[1024];

  std::vector<spfilter> sparse_filters;

};

void initView(view* view);

uint16_t *read_input_ppm(
    const char *input_ppm_path,
    int32_t &nr,
    int32_t &nc,
    int32_t &ncomp,
    const int32_t bpc,
    const std::string colorspace);

void write_output_ppm(
    const uint16_t *texture_view_in_encoder_colorspace,
    const char *output_ppm_path,
    const int32_t nr,
    const int32_t nc,
    const int32_t ncomp,
    const int32_t bpc,
    const std::string colorspace);

void setPaths(
    view *SAI,
    const char *input_dir,
    const char *output_dir);

int32_t get_highest_level(
    view *LF,
    const int32_t n_views_total);

int32_t* determine_levels(
    view *LF,
    const int32_t n_views_total);

#endif
