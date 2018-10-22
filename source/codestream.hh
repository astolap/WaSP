#ifndef CODESTREAM_HH
#define CODESTREAM_HH

#include "view.hh"
#include "minconf.hh"
#include <iostream>

void viewHeaderToCodestream(int &n_bytes_prediction, view *SAI, FILE *output_LF_file, const int yuv_transform_s);
void codestreamToViewHeader(int &n_bytes_prediction, view *SAI, FILE *input_LF, minimal_config &mconf);

#endif