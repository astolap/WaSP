#ifndef YCBCR_HH
#define YCBCR_HH

#define YUV_422 false

void RGB2YUV422(unsigned short *rgb, unsigned short **yy, unsigned short **cbb, unsigned short **crr,
	const int NR, const int NC, const int NCOMP, const int N);

void RGB2YCbCr(unsigned short *rgb, unsigned short *ycbcr, const int nr, const int nc,const int N);

void YCbCr2RGB(unsigned short *ycbcr, unsigned short *rgb, const int nr, const int nc, const int N);


#endif