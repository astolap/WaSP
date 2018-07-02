#ifndef PPM_HH
#define PPM_HH

bool aux_read16PGMPPM(const char* filename, int &width, int &height, int &ncomp, unsigned short *&img);
bool aux_write16PGMPPM(const char* filename, const int width, const int height, const int ncomp, unsigned short *img);

#endif