#include "psnr.hh"
#include "ycbcr.hh"

#include <stdio.h>
#include <cmath>

float getPSNR(FILE *fileout, const char *path_out_ppm, const char *path_input_ppm, const char *difftest_call)
{

	if (fileout == NULL)
		fileout = stdout;

	/* run psnr here */

	char psnr_call[1024];
	sprintf(psnr_call, "%s%s%s%s", difftest_call, path_out_ppm, " ", path_input_ppm);

	FILE *pfile;
	pfile = _popen(psnr_call, "r");

	char psnr_buffer[1024];
	while (fgets(psnr_buffer, sizeof(psnr_buffer), pfile) != 0) {
		/*...*/
	}
	_pclose(pfile);

	char tmp_char[1024];
	float psnr_value = 0;

	sscanf(psnr_buffer, "%s\t%f", tmp_char, &psnr_value);

	fprintf(fileout, "%s\n", psnr_buffer);

	return psnr_value;
}


double getYCbCr_422_PSNR(unsigned short *im0, unsigned short* im1, const int NR, const int NC, const int NCOMP, const int N)
{

	unsigned short 
		*im0_y, *im1_y, 
		*im0_cb, *im1_cb, 
		*im0_cr, *im1_cr;

	RGB2YUV422(im0, &im0_y, &im0_cb, &im0_cr, NR, NC, NCOMP, N);
	RGB2YUV422(im1, &im1_y, &im1_cb, &im1_cr, NR, NC, NCOMP, N);

	double nd = (double)(1 << N) - 1; // pow(2, N) - 1;

	double PSNR_Y = PSNR(im0_y, im1_y, NR, NC, 1, nd);
	double PSNR_Cb = PSNR(im0_cb, im1_cb, NR, NC/2, 1, nd);
	double PSNR_Cr = PSNR(im0_cb, im1_cb, NR, NC/2, 1, nd);

	delete[](im0_y);
	delete[](im1_y);
	delete[](im0_cb);
	delete[](im1_cb);
	delete[](im0_cr);
	delete[](im1_cr);

	return (6 * PSNR_Y + PSNR_Cb + PSNR_Cr) / 8;

}

double PSNR(unsigned short *im0, unsigned short* im1, const int NR, const int NC, const int NCOMP, double maxval)
{

	double se = 0;

	for (int ii = 0; ii < NR*NC*NCOMP; ii++)

	{
		double dx = (double)(*(im0 + ii)) - (double)(*(im1 + ii));

		se += dx*dx;
	}

	double mse = se / NR / NC / NCOMP;

	return 10 * log10((maxval*maxval) / mse);

}

double PSNR(unsigned short *im0, unsigned short* im1, const int NR, const int NC, const int NCOMP)
{

	double se = 0;

	double maxval = 0;

	for (int ii = 0; ii < NR*NC*NCOMP; ii++)

	{
		double dx = (double)(*(im0 + ii)) - (double)(*(im1 + ii));

		maxval = *(im1 + ii) > maxval ? *(im1 + ii) : maxval;

		se += dx*dx;
	}

	double mse = se / NR / NC / NCOMP;

	return 10 * log10((maxval*maxval) / mse);

}

