#ifndef RESIDUALJP2_HH
#define RESIDUALJP2_HH

#include <cstdio>
#include <vector>

void getJP2Header(unsigned char *JP2, unsigned char *&header, int JP2Size, int &headerSize);

int getJP2DictionaryIndex(unsigned char *JP2header, int headerSize,
	std::vector< std::vector<unsigned char>> JP2_dict);

void readResidualFromDisk(const char *jp2_residual_path_jp2, int &n_bytes_residual, FILE *input_LF,
	std::vector< std::vector<unsigned char>> &JP2_dict);

void updateJP2Dictionary(std::vector< std::vector<unsigned char>> &JP2_dict, unsigned char *header, int headerSize);

void writeResidualToDisk(const char *jp2_residual_path_jp2, FILE *output_LF_file, int &n_bytes_residual,
	std::vector< std::vector<unsigned char>> &JP2_dict);

void decodeResidualJP2(unsigned short *ps, const char *kdu_expand_path, const char *jp2_residual_path_jp2, const char *ppm_residual_path, int ncomp, const int offset, const int maxvali,
	const bool RESIDUAL_16BIT_bool);
	
void decodeResidualJP2_YUV(unsigned short *ps, const char *kdu_expand_path, char *ycbcr_jp2_names[], char *ycbcr_pgm_names[], const int ncomp, const int offset, const int maxvali,
	const bool RESIDUAL_16BIT_bool);

void encodeResidualJP2_YUV(const int nr, const int nc, unsigned short *original_intermediate_view, unsigned short *ps, char *ycbcr_pgm_names[],
	const char *kdu_compress_path, char *ycbcr_jp2_names[], const float residual_rate, const int ncomp, const int offset, float rate_a,
	const bool RESIDUAL_16BIT_bool);

void encodeResidualJP2(const int nr, const int nc, unsigned short *original_intermediate_view, unsigned short *ps, const char *ppm_residual_path,
	const char *kdu_compress_path, const char *jp2_residual_path_jp2, const float residual_rate, const int ncomp, const int offset, const bool RESIDUAL_16BIT_bool);


#endif