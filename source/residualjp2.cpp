#include "residualjp2.hh"
#include "ycbcr.hh"
#include "ppm.hh"
#include "fileaux.hh"
#include "clip.hh"
#include "medianfilter.hh"

#include <cstdlib>
#include <cstdio>
#include <cstring>

#define CLEVELS 6

void decodeResidualJP2(unsigned short *ps, const char *kdu_expand_path, const char *jp2_residual_path_jp2, const char *ppm_residual_path, int ncomp, const int offset, const int maxvali,
	const bool RESIDUAL_16BIT_bool)
{
	/* decode residual with kakadu */
	char kdu_expand_s[1024];
	sprintf(kdu_expand_s, "\"%s\"%s%s%s%s", kdu_expand_path, " -i ", jp2_residual_path_jp2, " -o ", ppm_residual_path);

	//std::cout << kdu_expand_s << "\n";

	int status = system_1(kdu_expand_s);

	signed int dv = RESIDUAL_16BIT_bool ? 1 : 2;
	signed int BP = RESIDUAL_16BIT_bool ? 16 : 10;
	signed int maxval = (1 << BP) - 1;// pow(2, BP) - 1;

	/* apply residual */

	unsigned short* jp2_residual;

	int nc1, nr1;

	if (aux_read16PGMPPM(ppm_residual_path, nc1, nr1, ncomp, jp2_residual))
	{

		for (int iir = 0; iir < nc1*nr1 * ncomp; iir++)
		{
			signed int val = (signed int)*(ps + iir) + (signed int)(jp2_residual[iir] * dv) - offset; // we assume that for 10bit case we have offset as 2^10-1, so go from 2^11 range to 2^10 and lose 1 bit of precision
			val = clip(val, 0, maxvali);
			*(ps + iir) = (unsigned short)(val);
		}

		delete[](jp2_residual);
	}
}

void decodeResidualJP2_YUV(unsigned short *ps, const char *kdu_expand_path, char *ycbcr_jp2_names[], char *ycbcr_pgm_names[], const int ncomp, const int offset, const int maxvali,
	const bool RESIDUAL_16BIT_bool)
{
	/* decode residual with kakadu */
	char kdu_expand_s[1024];

	for (int icomp = 0; icomp < ncomp; icomp++) {
		sprintf(kdu_expand_s, "\"%s\"%s%s%s%s", kdu_expand_path, " -i ", ycbcr_jp2_names[icomp], " -o ", ycbcr_pgm_names[icomp]);
		int status = system_1(kdu_expand_s);
		if (status < 0) {
			printf("KAKADU ERROR\nTERMINATING ... \n");
			exit(0);
		}
	}

	unsigned short *ycbcr = NULL;

	unsigned short *jp2_residual;

	int nc1, nr1,ncomp1;

	for (int icomp = 0; icomp < ncomp; icomp++) {
		if (aux_read16PGMPPM(ycbcr_pgm_names[icomp], nc1, nr1, ncomp1, jp2_residual))
		{
			if (ycbcr == NULL) {
				ycbcr = new unsigned short[nc1*nr1*ncomp]();
			}

			if (YUV_422) {
				if (icomp < 1) {
					memcpy(ycbcr + icomp*nr1*nc1, jp2_residual, sizeof(unsigned short)*nr1*nc1);
				}
				else {
					for (int cc = 0; cc < nc1; cc++) {
						memcpy(ycbcr + cc * 2 * nr1 + icomp*nr1*nc1*2, jp2_residual + cc*nr1, sizeof(unsigned short)*nr1);
						memcpy(ycbcr + (cc*2+1)*nr1 + icomp*nr1*nc1*2, jp2_residual + cc*nr1, sizeof(unsigned short)*nr1);
					}
					nc1 = nc1 * 2; // since we keep using this variable...
				}
			}
			else {
				memcpy(ycbcr + icomp*nr1*nc1, jp2_residual, sizeof(unsigned short)*nr1*nc1);
			}

			

			delete[](jp2_residual);
		}
	}

	signed int dv = RESIDUAL_16BIT_bool ? 1 : 2;
	signed int BP = RESIDUAL_16BIT_bool ? 16 : 10;
	signed int maxval = (1 << BP) - 1;// pow(2, BP) - 1;

	for (int ii = 0; ii < nr1*nc1*ncomp; ii++) {
		*(ycbcr + ii) = clip(*(ycbcr + ii), (unsigned short)0, (unsigned short)maxval);
	}

	unsigned short *rgb = new unsigned short[nr1*nc1*ncomp]();


	if (RESIDUAL_16BIT_bool) {
		YCbCr2RGB(ycbcr, rgb, nr1, nc1, 16);
	}
	else {
		YCbCr2RGB(ycbcr, rgb, nr1, nc1, 10);
	}

	/* apply residual */

	for (int iir = 0; iir < nc1*nr1 * ncomp; iir++)
	{
		signed int val = (signed int)*(ps + iir) + (signed int)(rgb[iir]*dv) - offset; // we assume that for 10bit case we have offset as 2^10-1, so go from 2^11 range to 2^10 and lose 1 bit of precision
		val = clip(val, 0, maxvali);
		*(ps + iir) = (unsigned short)(val);
	}

	delete[](ycbcr);
	delete[](rgb);
}


void encodeResidualJP2_YUV(const int nr, const int nc, unsigned short *original_intermediate_view, unsigned short *ps, char *ycbcr_pgm_names[],
	const char *kdu_compress_path, char *ycbcr_jp2_names[], const float residual_rate, const int ncomp, const int offset, float rate_a,
	const bool RESIDUAL_16BIT_bool)
{
	/*establish residual*/
	unsigned short *residual_image = new unsigned short[nr*nc * ncomp]();

	signed int dv = RESIDUAL_16BIT_bool ? 1 : 2;
	signed int BP = RESIDUAL_16BIT_bool ? 16 : 10;
	signed int maxval = (1 << BP) - 1;// pow(2, BP) - 1;

	for (int iir = 0; iir < nr*nc*ncomp; iir++) {
		signed int res_val = ( (((signed int)*(original_intermediate_view + iir)) - ((signed int)*(ps + iir)) + offset) )/dv;
		res_val = clip(res_val, 0, maxval);
		*(residual_image + iir) = (unsigned short)(res_val);
	}

	unsigned short *ycbcr = new unsigned short[nr*nc*ncomp]();

	if (RESIDUAL_16BIT_bool) {
		RGB2YCbCr(residual_image, ycbcr, nr, nc, 16);
	}
	else {
		RGB2YCbCr(residual_image, ycbcr, nr, nc, 10);
	}

	unsigned short *tmp_im = new unsigned short[nr*nc]();

	

	for (int icomp = 0; icomp < ncomp; icomp++) {

		float rateR = residual_rate;

		if (icomp < 1) {
			rateR = rate_a*rateR;
		}
		else {
			rateR = (1-rate_a)*rateR/2;
		}

		if (YUV_422) {
			if (icomp < 1) {
				memcpy(tmp_im, ycbcr + icomp*nr*nc, sizeof(unsigned short)*nr*nc);
				aux_write16PGMPPM(ycbcr_pgm_names[icomp], nc, nr, 1, tmp_im);
			}
			else {
				for (int cc = 0; cc < nc; cc += 2) {
					memcpy(tmp_im + (cc / 2)*nr, ycbcr + cc*nr + nr*nc*icomp, sizeof(unsigned short)*nr);
				}
				aux_write16PGMPPM(ycbcr_pgm_names[icomp], nc/2, nr, 1, tmp_im);
			}
		}
		else {
			memcpy(tmp_im, ycbcr + icomp*nr*nc, sizeof(unsigned short)*nr*nc);
			aux_write16PGMPPM(ycbcr_pgm_names[icomp], nc, nr, 1, tmp_im);
		}

		char kdu_compress_s[1024];
	
		sprintf(kdu_compress_s, "\"%s\"%s%s%s%s%s%f%s%d", kdu_compress_path, " -i ", ycbcr_pgm_names[icomp], " -o ", ycbcr_jp2_names[icomp], " -no_weights -no_info -precise -rate ", rateR,
			" Clevels=",CLEVELS);

		int status = system_1(kdu_compress_s);

		if (status < 0) {
			printf("KAKADU ERROR\nTERMINATING ... \n");
			exit(0);
		}


	}

	delete[](tmp_im);
	delete[](ycbcr);
	delete[](residual_image);


}

void encodeResidualJP2(const int nr, const int nc, unsigned short *original_intermediate_view, unsigned short *ps, const char *ppm_residual_path,
	const char *kdu_compress_path, const char *jp2_residual_path_jp2, const float residual_rate, const int ncomp, const int offset, const bool RESIDUAL_16BIT_bool)
{
	/*establish residual*/
	unsigned short *residual_image = new unsigned short[nr*nc * ncomp]();

	signed int dv = RESIDUAL_16BIT_bool ? 1 : 2;
	signed int BP = RESIDUAL_16BIT_bool ? 16 : 10;
	signed int maxval = (1 << BP) - 1;// pow(2, BP) - 1;

	for (int iir = 0; iir < nr*nc*ncomp; iir++) {
		signed int res_val = ((((signed int)*(original_intermediate_view + iir)) - ((signed int)*(ps + iir)) + offset)) / dv;
		res_val = clip(res_val, 0, maxval);
		*(residual_image + iir) = (unsigned short)(res_val);
	}

	aux_write16PGMPPM(ppm_residual_path, nc, nr, ncomp, residual_image);

	delete[](residual_image);

	/* here encode residual with kakadu */

	char kdu_compress_s[1024]; // tolerance 0 ?
	sprintf(kdu_compress_s, "\"%s\"%s%s%s%s%s%f%s%d", kdu_compress_path, " -i ", ppm_residual_path, " -o ", jp2_residual_path_jp2, " -no_weights -no_info -precise -rate ", residual_rate,
		" Clevels=", CLEVELS);

	//std::cout << kdu_compress_s << "\n";

	int status = system_1(kdu_compress_s);
}
