#include "ycbcr.hh"
#include "clip.hh"

#include <cstring>

void RGB2YUV422(unsigned short *rgb, unsigned short **yy, unsigned short **cbb, unsigned short **crr,
	const int NR, const int NC, const int NCOMP, const int N) 
{

	unsigned short *ycbcr;

	ycbcr = new unsigned short[NR*NC*NCOMP]();

	RGB2YCbCr(rgb, ycbcr, NR, NC, N);

	unsigned short *y, *cb, *cr;

	*yy = new unsigned short[NR*NC*NCOMP]();
	*cbb = new unsigned short[NR*NC / 2 * NCOMP]();
	*crr = new unsigned short[NR*NC / 2 * NCOMP]();

	y = *yy;
	cb = *cbb;
	cr = *crr;

	memcpy(y, ycbcr, sizeof(unsigned short)*NR*NC);

	for (int cc = 0; cc < NC; cc += 2) {
		memcpy(cb + (cc / 2)*NR, ycbcr + cc*NR + NR*NC, sizeof(unsigned short)*NR);
		memcpy(cr + (cc / 2)*NR, ycbcr + cc*NR + NR*NC * 2, sizeof(unsigned short)*NR);
	}

	delete[](ycbcr);
}


void RGB2YCbCr(unsigned short *rgb, unsigned short *ycbcr, const int nr, const int nc,const int N) 
{

	/* N-bit RGB 444 -> YCbCr 444 conversion */

	double M[] = { 0.212600000000000, -0.114572000000000,   0.500000000000000,
		0.715200000000000, -0.385428000000000, -0.454153000000000,
		0.072200000000000,   0.500000000000000, -0.045847000000000, };

	double *rgbD = new double[nr*nc * 3]();
	double *ycbcrD = new double[nr*nc * 3]();

	double nd = (double)(1 << (N-8));// pow(2, (double)N - 8);

	double clipval = (double)(1 << N) - 1; // pow(2, N) - 1;

	for (int ii = 0; ii < nr*nc*3; ii++) {

		*(rgbD + ii) = (double) *(rgb + ii);
		*(rgbD + ii) = *(rgbD + ii) / clipval;

	}

	for (int ii = 0; ii < nr*nc; ii++) {
		for (int icomp = 0; icomp < 3; icomp++) {

			*(ycbcrD + ii + icomp*nr*nc) = *(rgbD + ii)*M[icomp + 0] + 
				*(rgbD + ii + 1*nr*nc)*M[icomp + 3] + 
				*(rgbD + ii + 2*nr*nc)*M[icomp + 6];

			//printf("rgb\t%f\tycbcr\t%f\n", *(rgbD + ii + icomp*nr*nc), *(ycbcrD + ii + icomp*nr*nc));

			if (icomp < 1) {
				*(ycbcrD + ii + icomp*nr*nc) = (219 * (*(ycbcrD + ii + icomp*nr*nc)) + 16)*nd;
			}
			else {
				*(ycbcrD + ii + icomp*nr*nc) = (224 * (*(ycbcrD + ii + icomp*nr*nc)) + 128)*nd;
			}

			*(ycbcr + ii + icomp*nr*nc) = (unsigned short)*(ycbcrD + ii + icomp*nr*nc);
		}

	}

	delete[](rgbD);
	delete[](ycbcrD);

}

void YCbCr2RGB(unsigned short *ycbcr, unsigned short *rgb, const int nr, const int nc, const int N)
{

	double M[] = {		1.000000000000000,		1.000000000000000,				1.000000000000000,
						0,						-0.187330000000000,				1.855630000000000,
						1.574800000000000,		-0.468130000000000,                   0 };

	double *rgbD = new double[nr*nc * 3]();
	double *ycbcrD = new double[nr*nc * 3]();

	double nd = (double)(1 << (N - 8));// pow(2, (double)N - 8);

	unsigned short clipval = (unsigned short)(1 << N) - 1;// pow(2, N) - 1;

	double sval1 = 16 * nd;
	double sval2 = 219 * nd;
	double sval3 = 128 * nd;
	double sval4 = 224 * nd;

	for (int ii = 0; ii < nr*nc; ii++) {

		for (int icomp = 0; icomp < 3; icomp++) {

			*(ycbcrD + ii + icomp*nr*nc) = (double) *(ycbcr + ii + icomp*nr*nc);
			
			if (icomp < 1) {
				*(ycbcrD + ii + icomp*nr*nc) = clip( (*(ycbcrD + ii + icomp*nr*nc) - sval1) / sval2, 0.0, 1.0);
			}
			else {
				*(ycbcrD + ii + icomp*nr*nc) = clip( (*(ycbcrD + ii + icomp*nr*nc) - sval3) / sval4, -0.5, 0.5);
			}

		}

		for (int icomp = 0; icomp < 3; icomp++) {

			*(rgbD + ii + icomp*nr*nc) = *(ycbcrD + ii)*M[icomp + 0] +
				*(ycbcrD + ii + 1 * nr*nc)*M[icomp + 3] +
				*(ycbcrD + ii + 2 * nr*nc)*M[icomp + 6];

				*(rgb + ii + icomp*nr*nc) = (unsigned short)clip(( *(rgbD + ii + icomp*nr*nc)*clipval ), 0.0, (double)clipval);
		}

	}

	delete[](rgbD);
	delete[](ycbcrD);

}
