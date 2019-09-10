/* WaSPDecoder.hh */
/* Author: Pekka Astola */
/* <pekka.astola@tuni.fi>*/


#include <vector>
#include <cstdint>

#include "view.hh"
#include "bitdepth.hh"
#include "WaSPConf.hh"

using std::int32_t;
using std::uint32_t;

using std::int16_t;
using std::uint16_t;

using std::int8_t;
using std::uint8_t;


#ifndef SOURCE_WASPDECODER_HH_
#define SOURCE_WASPDECODER_HH_

class WaSPDecoder {
private:

    int32_t n_bytes_prediction = 0;
    int32_t n_bytes_residual = 0;
    int32_t number_of_views = 0;
    int32_t number_of_rows = 0;
    int32_t number_of_columns = 0;
    int32_t maxh = 0; /*maximum hierarchy level*/

    std::string colorspace_LF;

    uint16_t minimum_depth = 0;

    bool use_color_transform = false;

    WaSPsetup setup;

    view* LF = nullptr;
    FILE* input_LF = nullptr;
    std::vector<std::vector<uint8_t>> JP2_dict;

    void dealloc();

protected:
    void WaSP_decode_header();
    void WaSP_decode_views();

    void predict_texture_view(view* SAI);

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

    template<class T>
    T clip(T value, T min, T max) {
        value = value < min ? min : value;
        value = value > max ? max : value;
        return value;
    }

public:
    WaSPDecoder(const WaSPsetup decoder_setup);
    virtual ~WaSPDecoder();
    void decode();
};

#endif
