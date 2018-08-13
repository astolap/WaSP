#ifndef CODESTREAM_HH
#define CODESTREAM_HH

#include "view.hh"
#include <iostream>

void viewHeaderToCodestream(bool &size_written, int &n_bytes_prediction, view *SAI, FILE *output_LF_file, const int yuv_transform_s);

#endif