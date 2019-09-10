/* codestream.hh */
/* Author: Pekka Astola */
/* <pekka.astola@tuni.fi>*/

#ifndef CODESTREAM_HH
#define CODESTREAM_HH

#include "view.hh"
#include "minconf.hh"
#include <iostream>

#include <cstdint>

using std::int32_t;
using std::uint32_t;

using std::int16_t;
using std::uint16_t;

using std::int8_t;
using std::uint8_t;

void viewHeaderToCodestream(
    int32_t &n_bytes_prediction, 
    view *SAI,
    FILE *output_LF_file);

void codestreamToViewHeader(
    int32_t &n_bytes_prediction, 
    view *SAI, 
    FILE *input_LF,
    minimal_config &mconf);

#endif
