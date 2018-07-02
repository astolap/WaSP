#ifndef PSNR_HH
#define PSNR_HH

#include <cstdio>

float getPSNR(FILE *fileout, const char *path_out_ppm, const char *path_input_ppm, const char *difftest_call);
double PSNR(unsigned short *im0, unsigned short* im1, const int NR, const int NC, const int NCOMP, double maxval);
double PSNR(unsigned short *im0, unsigned short* im1, const int NR, const int NC, const int NCOMP);
double getYCbCr_422_PSNR(unsigned short *im0, unsigned short* im1, const int NR, const int NC, const int NCOMP, const int N);

#endif