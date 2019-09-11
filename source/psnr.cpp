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

#include "psnr.hh"
#include "ycbcr.hh"
#include "fileaux.hh"

#include <stdio.h>
#include <cmath>

double getYCbCr_444_PSNR(
    uint16_t *im0,
    uint16_t* im1,
    const int32_t NR,
    const int32_t NC,
    const int32_t NCOMP,
    const int32_t N)
{

    uint16_t *ycbcr0, *ycbcr1;

    ycbcr0 = im0;
    ycbcr1 = im1;

    double nd = (double)(1 << N) - 1;

    double PSNR_Y = PSNR(ycbcr0, ycbcr1, NR, NC, 1, nd);
    double PSNR_Cb = PSNR(ycbcr0 + NR*NC, ycbcr1 + NR*NC, NR, NC, 1, nd);
    double PSNR_Cr = PSNR(ycbcr0 + 2 * NR*NC, ycbcr1 + 2 * NR*NC, NR, NC, 1, nd);

    return (6 * PSNR_Y + PSNR_Cb + PSNR_Cr) / 8;

}


double getYCbCr_422_PSNR(uint16_t *im0, uint16_t* im1, const int32_t NR,
                         const int32_t NC, const int32_t NCOMP, const int32_t N) {

  uint16_t *im0_y, *im1_y, *im0_cb, *im1_cb, *im0_cr, *im1_cr;

  RGB2YUV422(im0, &im0_y, &im0_cb, &im0_cr, NR, NC, NCOMP, N);
  RGB2YUV422(im1, &im1_y, &im1_cb, &im1_cr, NR, NC, NCOMP, N);

  double nd = (double) (1 << N) - 1;  // pow(2, N) - 1;

  double PSNR_Y = PSNR(im0_y, im1_y, NR, NC, 1, nd);
  double PSNR_Cb = PSNR(im0_cb, im1_cb, NR, NC / 2, 1, nd);
  double PSNR_Cr = PSNR(im0_cb, im1_cb, NR, NC / 2, 1, nd);

  delete[] (im0_y);
  delete[] (im1_y);
  delete[] (im0_cb);
  delete[] (im1_cb);
  delete[] (im0_cr);
  delete[] (im1_cr);

  return (6 * PSNR_Y + PSNR_Cb + PSNR_Cr) / 8;

}

double PSNR(
    uint16_t *im0, 
    uint16_t* im1, 
    const int32_t NR,
    const int32_t NC, 
    const int32_t NCOMP, 
    double maxval) {

  double se = 0;

  for (int32_t ii = 0; ii < NR * NC * NCOMP; ii++)

  {
    double dx = (double) (*(im0 + ii)) - (double) (*(im1 + ii));

    se += dx * dx;
  }

  double mse = se / NR / NC / NCOMP;

  return 10 * log10((maxval * maxval) / mse);

}

double PSNR(
    uint16_t *im0, 
    uint16_t* im1, 
    const int32_t NR,
    const int32_t NC, 
    const int32_t NCOMP) {

  double se = 0;

  double maxval = 0;

  for (int32_t ii = 0; ii < NR * NC * NCOMP; ii++)

  {
    double dx = (double) (*(im0 + ii)) - (double) (*(im1 + ii));

    maxval = *(im1 + ii) > maxval ? *(im1 + ii) : maxval;

    se += dx * dx;
  }

  double mse = se / NR / NC / NCOMP;

  return 10 * log10((maxval * maxval) / mse);

}

