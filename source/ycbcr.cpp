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

#include "ycbcr.hh"
#include "clip.hh"

#include <cstring>
#include <cstdint>

using std::int32_t;
using std::uint32_t;

using std::int16_t;
using std::uint16_t;

using std::int8_t;
using std::uint8_t;

void RGB2YUV422(uint16_t *rgb, uint16_t **yy, uint16_t **cbb,
                uint16_t **crr, const int32_t NR, const int32_t NC,
                const int32_t NCOMP, const int32_t N) {

  uint16_t *ycbcr;

  ycbcr = new uint16_t[NR * NC * NCOMP]();

  RGB2YCbCr(rgb, ycbcr, NR, NC, N);

  uint16_t *y, *cb, *cr;

  *yy = new uint16_t[NR * NC * NCOMP]();
  *cbb = new uint16_t[NR * NC / 2 * NCOMP]();
  *crr = new uint16_t[NR * NC / 2 * NCOMP]();

  y = *yy;
  cb = *cbb;
  cr = *crr;

  memcpy(y, ycbcr, sizeof(uint16_t) * NR * NC);

  for (int32_t cc = 0; cc < NC; cc += 2) {
    memcpy(cb + (cc / 2) * NR, ycbcr + cc * NR + NR * NC,
           sizeof(uint16_t) * NR);
    memcpy(cr + (cc / 2) * NR, ycbcr + cc * NR + NR * NC * 2,
           sizeof(uint16_t) * NR);
  }

  delete[] (ycbcr);
}

void RGB2YCbCr(
    const uint16_t *rgb, 
    uint16_t *ycbcr, 
    const int32_t nr,
    const int32_t nc, 
    const int32_t N) {

  /* N-bit RGB 444 -> YCbCr 444 conversion */

  double M[] = { 0.212600000000000, -0.114572000000000, 0.500000000000000,
      0.715200000000000, -0.385428000000000, -0.454153000000000,
      0.072200000000000, 0.500000000000000, -0.045847000000000, };

  double *rgbD = new double[nr * nc * 3]();
  double *ycbcrD = new double[nr * nc * 3]();

  double nd = (double) (1 << (N - 8));  // pow(2, (double)N - 8);

  double clipval = (double) (1 << N) - 1;  // pow(2, N) - 1;

  for (int32_t ii = 0; ii < nr * nc * 3; ii++) {

    *(rgbD + ii) = (double) *(rgb + ii);
    *(rgbD + ii) = *(rgbD + ii) / clipval;

  }

  for (int32_t ii = 0; ii < nr * nc; ii++) {
    for (int32_t icomp = 0; icomp < 3; icomp++) {

      *(ycbcrD + ii + icomp * nr * nc) = *(rgbD + ii) * M[icomp + 0]
          + *(rgbD + ii + 1 * nr * nc) * M[icomp + 3]
          + *(rgbD + ii + 2 * nr * nc) * M[icomp + 6];

      //printf("rgb\t%f\tycbcr\t%f\n", *(rgbD + ii + icomp*nr*nc), *(ycbcrD + ii + icomp*nr*nc));

      if (icomp < 1) {
        *(ycbcrD + ii + icomp * nr * nc) = (219
            * (*(ycbcrD + ii + icomp * nr * nc)) + 16) * nd;
      } else {
        *(ycbcrD + ii + icomp * nr * nc) = (224
            * (*(ycbcrD + ii + icomp * nr * nc)) + 128) * nd;
      }

      *(ycbcr + ii + icomp * nr * nc) = (uint16_t) *(ycbcrD + ii
          + icomp * nr * nc);
    }

  }

  delete[] (rgbD);
  delete[] (ycbcrD);

}

void YCbCr2RGB(
    const uint16_t *ycbcr, 
    uint16_t *rgb, 
    const int32_t nr,
    const int32_t nc, 
    const int32_t N) {

  double M[] = { 1.000000000000000, 1.000000000000000, 1.000000000000000, 0,
      -0.187330000000000, 1.855630000000000, 1.574800000000000,
      -0.468130000000000, 0 };

  double *rgbD = new double[nr * nc * 3]();
  double *ycbcrD = new double[nr * nc * 3]();

  double nd = (double) (1 << (N - 8));  // pow(2, (double)N - 8);

  uint16_t clipval = (uint16_t) (1 << N) - 1;  // pow(2, N) - 1;

  double sval1 = 16 * nd;
  double sval2 = 219 * nd;
  double sval3 = 128 * nd;
  double sval4 = 224 * nd;

  for (int32_t ii = 0; ii < nr * nc; ii++) {

    for (int32_t icomp = 0; icomp < 3; icomp++) {

      *(ycbcrD + ii + icomp * nr * nc) =
          (double) *(ycbcr + ii + icomp * nr * nc);

      if (icomp < 1) {
        *(ycbcrD + ii + icomp * nr * nc) = clip(
            (*(ycbcrD + ii + icomp * nr * nc) - sval1) / sval2, 0.0, 1.0);
      } else {
        *(ycbcrD + ii + icomp * nr * nc) = clip(
            (*(ycbcrD + ii + icomp * nr * nc) - sval3) / sval4, -0.5, 0.5);
      }

    }

    for (int32_t icomp = 0; icomp < 3; icomp++) {

      *(rgbD + ii + icomp * nr * nc) = *(ycbcrD + ii) * M[icomp + 0]
          + *(ycbcrD + ii + 1 * nr * nc) * M[icomp + 3]
          + *(ycbcrD + ii + 2 * nr * nc) * M[icomp + 6];

      *(rgb + ii + icomp * nr * nc) = (uint16_t) clip(
          (*(rgbD + ii + icomp * nr * nc) * clipval), 0.0, (double) clipval);
    }

  }

  delete[] (rgbD);
  delete[] (ycbcrD);

}
