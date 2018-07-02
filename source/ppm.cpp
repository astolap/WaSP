#include "ppm.hh"

#include <cstdio>
#include <cstring>

#define IO_V false

bool aux_read16PGMPPM(const char* filename, int &width, int &height, int &ncomp, unsigned short *&img)
{
	if (IO_V)
		printf("Reading %s\n", filename);

	int  max, x, y;
	int red, green, blue;
	char dummy[100];

	unsigned short *Image16bit = NULL;

	FILE *filept;

	filept = fopen(filename, "rb");

	if (filept == NULL) {
		printf("%s does not exist\n", filename);
		return false;
	}


	fscanf(filept, "%s", dummy);
	fscanf(filept, "%d %d\n", &width, &height);
	fscanf(filept, "%d", &max);
	fgetc(filept);

	if (!strncmp(dummy, "P6", 2)) {
		ncomp = 3;
	}
	else if (!strncmp(dummy, "P5", 2)) {
		ncomp = 1;
	}
	else { printf("ERROR NOT PGM OR PPM\n"); return false; }


	//std::cout << width << "\t" << height << "\n";

	img = new unsigned short[width*height * ncomp]();

	Image16bit = new unsigned short[width*height * ncomp]();

	/*--< Read 16bit ppm image from filept >--*/
	int nread = (int)fread(Image16bit, sizeof(unsigned short), width*height * ncomp, filept);

	if (nread != width*height * ncomp)
	{
		fprintf(stderr, "READ ERROR aux_read16ppm() %s\n", filename);
		delete[](img);
		delete[](Image16bit);
		return false;
	}

	fclose(filept);

	int i = 0;

	for (x = 0; x < width; x++) {
		for (y = 0; y < height; y++) {

			red = Image16bit[(x + y*width) * ncomp];
			if (ncomp == 3) {
				green = Image16bit[(x + y*width) * ncomp + 1];
				blue = Image16bit[(x + y*width) * ncomp + 2];
			}

			// Exhange upper 8bit and lower 8bit for Intel x86
			red = ((red & 0x00ff) << 8) | ((red & 0xff00) >> 8);
			if (ncomp == 3) {
				green = ((green & 0x00ff) << 8) | ((green & 0xff00) >> 8);
				blue = ((blue & 0x00ff) << 8) | ((blue & 0xff00) >> 8);
			}

			img[i] = red;
			if (ncomp == 3) {
				img[i + height*width] = green;
				img[i + 2 * height*width] = blue;
			}

			i++;

		}
	}

	delete[](Image16bit);

	return true;

}

bool aux_write16PGMPPM(const char* filename, const int width, const int height, const int ncomp, unsigned short *img)
{

	if (IO_V)
		printf("Writing %s\n", filename);

	unsigned char *p;
	int i, tmp, j;

	unsigned short maxi = 0;

	FILE *filept;

	filept = fopen(filename, "wb");

	if (filept == NULL) {
		printf("Cannot open %s\n", filename);
		return false;
	}

	unsigned short* img16bit = new unsigned short[height*width * ncomp]();

	int lin_ind = 0;
	for (j = 0; j < height; j++) {
		for (i = 0; i < width; i++) {
			img16bit[lin_ind] = img[j + i*height];
			if (ncomp == 3) {
				img16bit[lin_ind + 1] = img[j + i*height + width*height];
				img16bit[lin_ind + 2] = img[j + i*height + 2 * width*height];
				lin_ind = lin_ind + 3;
			}
			else {
				lin_ind = lin_ind + 1;
			}

		}
	}

	for (i = 0; i < width*height * ncomp; i++)
	{
		if (*(img16bit + i) > maxi)
			maxi = *(img16bit + i);
	}

	p = (unsigned char *)img16bit;
	for (i = 0; i < width*height; i++) {
		tmp = *p; *p = *(p + 1); *(p + 1) = tmp; p += 2;
		if (ncomp == 3) {
			tmp = *p; *p = *(p + 1); *(p + 1) = tmp; p += 2;
			tmp = *p; *p = *(p + 1); *(p + 1) = tmp; p += 2;
		}
	}

	if (maxi > 1023)
		maxi = 65535;
	else
		maxi = 1023;

	if (ncomp == 3) {
		fprintf(filept, "P6\n%d %d\n%d\n", width, height, maxi);
	}
	else {
		fprintf(filept, "P5\n%d %d\n%d\n", width, height, maxi);
	}

	fwrite(img16bit, sizeof(unsigned short), width*height * ncomp, filept);

	fclose(filept);

	delete[](img16bit);

	return true;
}
