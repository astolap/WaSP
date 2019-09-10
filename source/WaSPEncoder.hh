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

#ifndef SOURCE_WASPENCODER_HH_
#define SOURCE_WASPENCODER_HH_

#define USE_difftest_ng false

#define STD_SEARCH_LOW 10
#define STD_SEARCH_HIGH 250
#define STD_SEARCH_STEP 10

#define SAVE_PARTIAL_WARPED_VIEWS false

#ifndef FLT_MAX
#define FLT_MAX 3.402823466e+38F        // max value
#endif

class WaSPEncoder {
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

  void generate_WaSP_bitstream();

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

  void generate_texture_residual_level_wise();

  void generate_inverse_depth_levelwise();

 public:
  virtual ~WaSPEncoder();

  WaSPEncoder(WaSPsetup encoder_setup);

  void encode();
};

#endif /* SOURCE_WASPENCODER_HH_ */
