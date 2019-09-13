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

#ifndef SPARSEFILTER_HH
#define SPARSEFILTER_HH

#include <vector>
#include <cstdint>

using std::int32_t;
using std::uint32_t;

using std::int16_t;
using std::uint16_t;

using std::int8_t;
using std::uint8_t;

#define SPARSE_BIAS_TERM 0.5

struct spfilter {

    std::vector<double> filter_coefficients;
    std::vector<int32_t> regressor_indexes;

    std::vector<int16_t> quantized_filter_coefficients;

    int32_t Ms;
    int32_t NNt;
    int32_t MT;

    double bias_term_value;

};

uint16_t *padArrayUint16_t(
    uint16_t *input_image,
    const uint32_t nr,
    const uint32_t nc,
    const uint32_t NNt);

std::vector<uint16_t> padArrayUint16_t_vec(
    uint16_t *input_image,
    const uint32_t nr,
    const uint32_t nc,
    const uint32_t NNt);

uint16_t *cropImage(
    uint16_t *input_image,
    const uint32_t nr,
    const uint32_t nc,
    const uint32_t NNt);

void quantize_and_reorder_spfilter(
    spfilter &sparse_filter);

void dequantize_and_reorder_spfilter(
    spfilter &sparse_filter);

spfilter getGlobalSparseFilter_vec(
    const uint16_t *original_image,
    const std::vector<std::vector<uint16_t>> &input_images,
    const int32_t nr,
    const int32_t nc,
    const int32_t NNt,
    const int32_t Ms,
    const double bias_term_value);

spfilter getGlobalSparseFilter(
    const uint16_t *original_image,
    const uint16_t *input_image,
    const int32_t nr,
    const int32_t nc,
    const int32_t NNt,
    const int32_t Ms,
    const double bias_term_value);

std::vector<double> applyGlobalSparseFilter_vec(
    const std::vector<std::vector<uint16_t>> &input_images,
    const int32_t nr,
    const int32_t nc,
    const int32_t Ms,
    const int32_t NNt,
    const double bias_term_value,
    const std::vector<double> filter_coeffs);

std::vector<double> applyGlobalSparseFilter(
    const uint16_t *input_image,
    const int32_t nr,
    const int32_t nc,
    const int32_t Ms,
    const int32_t NNt,
    const double bias_term_value,
    const std::vector<double> filter_coeffs);

spfilter getGlobalSparseFilter_vec_reg(
    const uint16_t *original_image,
    const std::vector<std::vector<uint16_t>> &input_images,
    const std::vector<int32_t> seg,
    const int32_t regi,
    const int32_t nr,
    const int32_t nc,
    const int32_t NNt,
    const int32_t Ms,
    const double bias_term_value);

std::vector<double> applyGlobalSparseFilter_vec_reg(
    const std::vector<std::vector<uint16_t>> &input_images,
    const std::vector<int32_t> seg,
    const int32_t regi,
    const int32_t nr,
    const int32_t nc,
    const int32_t Ms,
    const int32_t NNt,
    const double bias_term_value,
    const std::vector<double> filter_coeffs);

#endif
