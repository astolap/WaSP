/* residualjp2.hh */
/* Author: Pekka Astola */
/* <pekka.astola@tuni.fi>*/

#ifndef RESIDUALJP2_HH
#define RESIDUALJP2_HH

#include <cstdint>

using std::int32_t;
using std::uint32_t;

using std::int16_t;
using std::uint16_t;

using std::int8_t;
using std::uint8_t;

#include "view.hh"

void getJP2Header(
    unsigned char *JP2, 
    unsigned char *&header, 
    int32_t JP2Size,
    int32_t &headerSize);

int32_t getJP2DictionaryIndex(
    unsigned char *JP2header, 
    int32_t headerSize,
    std::vector<std::vector<unsigned char>> JP2_dict);

void readResidualFromDisk(
    const char *jp2_residual_path_jp2,
    int32_t &n_bytes_residual,
    FILE *input_LF,
    std::vector<std::vector<unsigned char>> &JP2_dict);

void updateJP2Dictionary(
    std::vector<std::vector<unsigned char>> &JP2_dict,
    unsigned char *header,
    int32_t headerSize);

void writeResidualToDisk(
    const char *jp2_residual_path_jp2,
    FILE *output_LF_file, 
    int32_t &n_bytes_residual,
    std::vector<std::vector<unsigned char>> &JP2_dict);

char *kakadu_oparams(
    const double rate,
    const std::string colorspace);

char *kakadu_cparams(
    const double *cweights,
    const int32_t ncomp);

int32_t encodeKakadu(
    const char *ppm_pgm_input_path,
    const char *kdu_compress_path,
    const char *jp2_output_path,
    const char *encoding_parameters,
    const double rate);

int32_t decodeKakadu(
    const char *ppm_pgm_output_path,
    const char *kdu_expand_path,
    const char *jp2_input_path);

double* get_residual(
    const uint16_t *original,
    const uint16_t *prediction,
    const int32_t nr,
    const int32_t nc,
    const int32_t ncomp);

uint16_t* quantize_residual(
    const double *residual,
    const int32_t nr,
    const int32_t nc,
    const int32_t ncomp,
    const int32_t bpc,
    const int32_t Q_i,
    const int32_t offset_i);

double* dequantize_residual(
    const uint16_t *qresidual,
    const int32_t nr,
    const int32_t nc,
    const int32_t ncomp,
    const int32_t bpc,
    const int32_t Q_i,
    const int32_t offset_i);

uint16_t *apply_residual(
    const uint16_t *prediction,
    const double *residual,
    const int32_t nr,
    const int32_t nc,
    const int32_t ncomp,
    const int32_t bpc);

uint16_t* decode_residual_JP2(
    const char *ppm_pgm_output_path,
    const char *kdu_expand_path,
    const char *jp2_input_path);

void encode_residual_JP2(
    const char *ppm_pgm_input_path,
    const char *kdu_compress_path,
    const char *jp2_output_path,
    const char *encoding_parameters,
    const double rate);

#endif
