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


#ifndef DECODER_HH
#define DECODER_HH

class decoder {
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
    void decode_header();
    void decode_views();

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
    decoder(const WaSPsetup decoder_setup);
    virtual ~decoder();
    void decode();
};

#endif
