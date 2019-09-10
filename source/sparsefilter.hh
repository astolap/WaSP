/* sparsefilter.hh */
/* Author: Pekka Astola */
/* <pekka.astola@tuni.fi>*/

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
    int32_t coeff_bit_precision;

    double bias_term_value;

};

uint16_t *padArrayUint16_t(
    const uint16_t *input_image,
    const uint32_t nr,
    const uint32_t nc,
    const uint32_t NNt);

uint16_t *cropImage(
    const uint16_t *input_image,
    const uint32_t nr,
    const uint32_t nc,
    const uint32_t NNt);

void quantize_and_reorder_spfilter(
    spfilter &sparse_filter);

void dequantize_and_reorder_spfilter(
    spfilter &sparse_filter);

spfilter getGlobalSparseFilter(
    const uint16_t *original_image,
    const uint16_t *input_image,
    const int32_t nr,
    const int32_t nc,
    const int32_t NNt,
    const int32_t Ms,
    const double bias_term_value);

std::vector<double> applyGlobalSparseFilter(
    const uint16_t *input_image,
    const int32_t nr,
    const int32_t nc,
    const int32_t Ms,
    const int32_t NNt,
    const double bias_term_value,
    const std::vector<double> filter_coeffs);

#endif
