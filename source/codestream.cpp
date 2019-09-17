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

#include "codestream.hh"
#include "view.hh"
#include "minconf.hh"

#include <iostream>
#include <vector>


void viewHeaderToCodestream(
    int32_t &n_bytes_prediction, 
    view *SAI,
    FILE *output_LF_file) {

  minimal_config mconf = makeMinimalConfig(SAI);

  //printf("size of minimal_config %i bytes\n", (int32_t)sizeof(minimal_config));

  n_bytes_prediction += (int32_t) fwrite(
      &mconf, 
      sizeof(minimal_config), 
      1,
      output_LF_file) * sizeof(minimal_config);

  /* lets see what else needs to be written to bitstream */

  if (SAI->has_x_displacement) {
    n_bytes_prediction += (int32_t) fwrite(&SAI->x, sizeof(float), 1,
                                       output_LF_file) * sizeof(float);
  }

  if (SAI->has_y_displacement) {
    n_bytes_prediction += (int32_t) fwrite(&SAI->y, sizeof(float), 1,
                                       output_LF_file) * sizeof(float);
  }

  if (SAI->has_color_references) {
    unsigned char tmpNREF = (unsigned char) SAI->n_references;
    n_bytes_prediction += (int32_t) fwrite(&tmpNREF, sizeof(unsigned char), 1,
                                       output_LF_file) * sizeof(unsigned char);
    for (int32_t ij = 0; ij < SAI->n_references; ij++) {
      uint16_t nid = (uint16_t) *(SAI->references + ij);
      n_bytes_prediction += (int32_t) fwrite(&nid, sizeof(uint16_t), 1,
                                         output_LF_file)
          * sizeof(uint16_t);
    }
  }

  if (SAI->has_depth_references) {
    unsigned char tmpNDREF = (unsigned char) SAI->n_depth_references;
    n_bytes_prediction += (int32_t) fwrite(&tmpNDREF, sizeof(unsigned char), 1,
                                       output_LF_file) * sizeof(unsigned char);
    for (int32_t ij = 0; ij < SAI->n_depth_references; ij++) {
      uint16_t nid = (uint16_t) *(SAI->depth_references + ij);
      n_bytes_prediction += (int32_t) fwrite(&nid, sizeof(uint16_t), 1,
                                         output_LF_file)
          * sizeof(uint16_t);
    }
  }

  int32_t MMM = (1 << SAI->n_references);
  int32_t N_LS = SAI->ncomp*((MMM*SAI->n_references) / 2);

  if (SAI->has_color_references) {
      if (SAI->mmode == 0) {
          /* use LS merging weights */
          n_bytes_prediction += (int32_t)fwrite(
              SAI->merge_weights,
              sizeof(int16_t),
              N_LS,
              output_LF_file)
              * sizeof(int16_t);
      }
      if (SAI->mmode == 1) {
          /* use standard deviation */
          n_bytes_prediction += (int32_t)fwrite(
              &SAI->stdd,
              sizeof(float),
              1,
              output_LF_file) * sizeof(float);
      }
  }

  if (SAI->use_global_sparse) {

      unsigned char tmpNNt = (unsigned char)SAI->NNt;
      unsigned char tmpMs = (unsigned char)SAI->Ms;

      n_bytes_prediction += (int32_t)fwrite(
          &tmpNNt,
          sizeof(unsigned char),
          1,
          output_LF_file) * sizeof(unsigned char);

      n_bytes_prediction += (int32_t)fwrite(
          &tmpMs,
          sizeof(unsigned char),
          1,
          output_LF_file) * sizeof(unsigned char);

      for (int32_t icomp = 0; icomp < SAI->ncomp; icomp++) {
          n_bytes_prediction += (int32_t)fwrite(
              &SAI->sparse_filters.at(icomp).quantized_filter_coefficients[0],
              sizeof(int16_t),
              SAI->Ms,
              output_LF_file)
              * sizeof(int16_t);
      }

      int32_t Nsp = (SAI->NNt * 2 + 1)* (SAI->NNt * 2 + 1) + 1;
      int32_t sp_mask_nbytes = (Nsp % 8) > 0 ? Nsp / 8 + 1 : Nsp / 8;

      uint8_t *sparsemask = new uint8_t[sp_mask_nbytes*SAI->ncomp]();

      for (int32_t icomp = 0; icomp < SAI->ncomp; icomp++) {

          for (int32_t ij = 0; ij < SAI->Ms; ij++) {

              uint32_t regr_indx = 
                  SAI->sparse_filters.at(icomp).regressor_indexes.at(ij);

              uint32_t q = regr_indx / 8;

              uint8_t *sparse_mask_byte = &sparsemask[q+sp_mask_nbytes*icomp];

              *sparse_mask_byte = *sparse_mask_byte
                  | (1 << (regr_indx - q * 8));

          }
      }

      n_bytes_prediction += (int32_t)fwrite(
          sparsemask,
          sizeof(uint8_t),
          sp_mask_nbytes*SAI->ncomp,
          output_LF_file)
          * sizeof(uint8_t);

      delete[](sparsemask);

  }

  return;

}

void codestreamToViewHeader(
    int32_t &n_bytes_prediction, 
    view *SAI, 
    FILE *input_LF,
    minimal_config &mconf) {

  n_bytes_prediction += (int32_t) fread(
      &mconf, 
      sizeof(minimal_config),
      1, 
      input_LF)
      * sizeof(minimal_config);

  setup_form_minimal_config(&mconf, SAI);

  if (SAI->has_x_displacement) {
    n_bytes_prediction += (int32_t) fread(&SAI->x, sizeof(float), 1, input_LF)
        * sizeof(float);
  }

  if (SAI->has_y_displacement) {
    n_bytes_prediction += (int32_t) fread(&SAI->y, sizeof(float), 1, input_LF)
        * sizeof(float);
  }

  if (SAI->has_color_references) {

    unsigned char tmpNREF = 0;

    n_bytes_prediction += (int32_t) fread(&tmpNREF, sizeof(unsigned char), 1,
                                      input_LF) * sizeof(unsigned char);

    SAI->n_references = tmpNREF;

    SAI->references = new int32_t[SAI->n_references]();
    for (int32_t ij = 0; ij < SAI->n_references; ij++) {
      uint16_t nid;
      n_bytes_prediction += (int32_t) fread(&nid, sizeof(uint16_t), 1,
                                        input_LF) * sizeof(uint16_t);
      *(SAI->references + ij) = (int32_t) nid;
    }
  }

  if (SAI->has_depth_references) {

    unsigned char tmpNDREF = 0;

    n_bytes_prediction += (int32_t) fread(&tmpNDREF, sizeof(unsigned char), 1,
                                      input_LF) * sizeof(unsigned char);

    SAI->n_depth_references = tmpNDREF;

    SAI->depth_references = new int32_t[SAI->n_depth_references]();
    for (int32_t ij = 0; ij < SAI->n_depth_references; ij++) {
      uint16_t nid;
      n_bytes_prediction += (int32_t) fread(&nid, sizeof(uint16_t), 1,
                                        input_LF) * sizeof(uint16_t);
      *(SAI->depth_references + ij) = (int32_t) nid;
    }
  }

  SAI->NB = (1 << SAI->n_references) * SAI->n_references;

  int32_t MMM = (1 << SAI->n_references);
  int32_t N_LS = SAI->ncomp*((MMM*SAI->n_references) / 2);

  if (SAI->has_color_references) {
      if (SAI->mmode == 0) {
          SAI->merge_weights = new int16_t[N_LS]();
          /* use LS merging weights */
          n_bytes_prediction += (int32_t)fread(
              SAI->merge_weights,
              sizeof(int16_t),
              N_LS,
              input_LF)
              * sizeof(int16_t);
      }
      if (SAI->mmode == 1) {
          /* use standard deviation */
          n_bytes_prediction += (int32_t)fread(
              &SAI->stdd,
              sizeof(float),
              1,
              input_LF) * sizeof(float);
      }
  }

  if (SAI->use_global_sparse) {

      unsigned char tmpNNt = 0;
      unsigned char tmpMs = 0;

      n_bytes_prediction += (int32_t)fread(
          &tmpNNt,
          sizeof(unsigned char),
          1,
          input_LF) * sizeof(unsigned char);

      n_bytes_prediction += (int32_t)fread(
          &tmpMs,
          sizeof(unsigned char),
          1,
          input_LF) * sizeof(unsigned char);

      SAI->NNt = (int32_t)tmpNNt;
      SAI->Ms = (int32_t)tmpMs;

      SAI->sparse_filters.clear();

      for (int32_t icomp = 0; icomp < SAI->ncomp; icomp++) {

          spfilter tmpsp;
          SAI->sparse_filters.push_back(tmpsp);

          SAI->sparse_filters.at(icomp).Ms = SAI->Ms;
          SAI->sparse_filters.at(icomp).NNt = SAI->NNt;

          SAI->sparse_filters.at(icomp).quantized_filter_coefficients.clear();

          for (int32_t ii = 0; ii < SAI->Ms; ii++) {
              SAI->sparse_filters.at(icomp).quantized_filter_coefficients.push_back(0);
          }

          n_bytes_prediction += (int32_t)fread(
              &SAI->sparse_filters.at(icomp).quantized_filter_coefficients[0],
              sizeof(int16_t),
              SAI->Ms,
              input_LF)
              * sizeof(int16_t);
      }


      int32_t Nsp = (SAI->NNt * 2 + 1)* (SAI->NNt * 2 + 1) + 1;
      int32_t sp_mask_nbytes = (Nsp % 8) > 0 ? Nsp / 8 + 1 : Nsp / 8;

      uint8_t *sparsemask = new uint8_t[sp_mask_nbytes*SAI->ncomp]();

      n_bytes_prediction += (int32_t)fread(
          sparsemask,
          sizeof(uint8_t),
          sp_mask_nbytes*SAI->ncomp,
          input_LF)
          * sizeof(uint8_t);

      for (int32_t icomp = 0; icomp < SAI->ncomp; icomp++) {

          SAI->sparse_filters.at(icomp).regressor_indexes.clear();

          for (int32_t ii = 0; ii < Nsp; ii++) {
              SAI->sparse_filters.at(icomp).regressor_indexes.push_back(0);
          }

          uint32_t ik = 0;

          for (int32_t ij = 0; ij < Nsp; ij++) {

              uint32_t q = ij / 8;

              uint8_t *sparse_mask_byte = &sparsemask[sp_mask_nbytes*icomp + q];

              if (*sparse_mask_byte & (1 << (ij - q * 8))) {
                  SAI->sparse_filters.at(icomp).regressor_indexes.at(ik) = ij;
                  ik++;
              }

          }
      }

      delete[](sparsemask);

  }

  return;

}
