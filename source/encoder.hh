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

#include <cstdint>
#include <vector>

#include "WaSPConf.hh"
#include "view.hh"

using namespace std;

using std::int32_t;
using std::uint32_t;

using std::int16_t;
using std::uint16_t;

using std::int8_t;
using std::uint8_t;

#ifndef ENCODER_HH
#define ENCODER_HH

#define USE_difftest_ng false

#define STD_SEARCH_LOW 10
#define STD_SEARCH_HIGH 250
#define STD_SEARCH_STEP 10

#define SAVE_PARTIAL_WARPED_VIEWS false

#ifndef FLT_MAX
#define FLT_MAX 3.402823466e+38F        // max value
#endif

class encoder {
 private:

  int32_t n_views_total = 0;
  int32_t std_search_s;

  view* LF;
  view ***LF_mat;

  //FILE *output_results_file;

  int32_t nr, nc; /*number of rows in SAI, number of columns in SAI, i.e., V and U*/
  int32_t MINIMUM_DEPTH = 0;
  int32_t maxh; /*maximum hierarchy level*/

  bool STD_SEARCH = false;

  std::string colorspace_LF;

  /* to get efficiency from multiple JP2 files, we remove parts of the files
   which are repetitive over all files. For this we have a minimalistic
   dictionary method. */

  vector<vector<uint8_t>> JP2_dict;

  char path_out_LF_data[1024];

  WaSPsetup setup;

 protected:

  void load_config_json(string config_json_file);
  void write_config(string config_json_file_out); /*debug reasons*/

  void write_bitstream();

  void merge_texture_views(
      view *SAI,
      view *LF,
      uint16_t **warped_texture_views,
      float **DispTargs);
  
  void forward_warp_texture_references(
      view *LF,
      view *SAI,
      uint16_t **warped_texture_views,
      uint16_t **warped_depth_views,
      float **DispTargs);

  void generate_normalized_disparity();
  void generate_texture();

 public:
  virtual ~encoder();

  encoder(WaSPsetup encoder_setup);

  void encode();
};

#endif
