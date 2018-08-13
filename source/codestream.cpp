#include "codestream.hh"
#include "view.hh"
#include "minconf.hh"

#include <iostream>
#include <algorithm>
#include <vector>

void viewHeaderToCodestream(bool &size_written, int &n_bytes_prediction, view *SAI, FILE *output_LF_file, const int yuv_transform_s) {

	if (!size_written) {
		n_bytes_prediction += (int)fwrite(&SAI->nr, sizeof(int), 1, output_LF_file) * sizeof(int); // needed only once per LF
		n_bytes_prediction += (int)fwrite(&SAI->nc, sizeof(int), 1, output_LF_file) * sizeof(int); // 
		n_bytes_prediction += (int)fwrite(&yuv_transform_s, sizeof(int), 1, output_LF_file) * sizeof(int);
		size_written = true;
	}

	minimal_config mconf = makeMinimalConfig(SAI);

	printf("size of minimal_config %i bytes\n", (int)sizeof(minimal_config));

	n_bytes_prediction += (int)fwrite(&mconf, sizeof(minimal_config), 1, output_LF_file) * sizeof(minimal_config);

	/* lets see what else needs to be written to bitstream */

	if (mconf.n_references > 0) {
		for (int ij = 0; ij < mconf.n_references; ij++) {
			unsigned short nid = (unsigned short) *(SAI->references + ij);
			n_bytes_prediction += (int)fwrite(&nid, sizeof(unsigned short), 1, output_LF_file) * sizeof(unsigned short);
		}
	}

	if (mconf.n_depth_references > 0) {
		for (int ij = 0; ij < mconf.n_depth_references; ij++) {
			unsigned short nid = (unsigned short) *(SAI->depth_references + ij);
			n_bytes_prediction += (int)fwrite(&nid, sizeof(unsigned short), 1, output_LF_file) * sizeof(unsigned short);
		}
	}

	if (mconf.Ms > 0 && mconf.NNt > 0) {


		std::vector<std::pair<int, int32_t>> sparsefilter;

		for (int ij = 0; ij < SAI->Ms; ij++) {
			std::pair<int, int32_t> tmp_sp;
			tmp_sp.first = SAI->sparse_mask[ij];
			tmp_sp.second = SAI->sparse_weights[ij];
			sparsefilter.push_back(tmp_sp);
		}

		sort(sparsefilter.begin(), sparsefilter.end());

		int32_t *tmp_sw = new int32_t[SAI->Ms]();

		int32_t sparse_mask_binary_p1 = 0;
		int32_t sparse_mask_binary_p2 = 0;

		for (int ij = 0; ij < SAI->Ms; ij++) {
			
			tmp_sw[ij] = sparsefilter.at(ij).second;

			if (sparsefilter.at(ij).first > 0) {
				
				if ( sparsefilter.at(ij).first < 32 ) {
					sparse_mask_binary_p1 = sparse_mask_binary_p1 | 1<<(int32_t)sparsefilter.at(ij).first; // note: regressor indexing starts from 1
				}
				else {
					sparse_mask_binary_p2 = sparse_mask_binary_p2 | 1<<( (int32_t)sparsefilter.at(ij).first-32 );
				}

			}

		}

		n_bytes_prediction += (int)fwrite(tmp_sw, sizeof(int32_t), SAI->Ms, output_LF_file) * sizeof(int32_t);
		n_bytes_prediction += (int)fwrite(&sparse_mask_binary_p1, sizeof(int32_t), 1, output_LF_file) * sizeof(int32_t);
		n_bytes_prediction += (int)fwrite(&sparse_mask_binary_p2, sizeof(int32_t), 1, output_LF_file) * sizeof(int32_t);

		delete[](tmp_sw);

		//n_bytes_prediction += (int)fwrite(SAI->sparse_mask, sizeof(unsigned char), SAI->Ms, output_LF_file) * sizeof(unsigned char);
		//n_bytes_prediction += (int)fwrite(SAI->sparse_weights, sizeof(int32_t), SAI->Ms, output_LF_file) * sizeof(int32_t);
	}

	if (!SAI->use_median) {
		if (SAI->stdd < 0.001) {
			if (mconf.n_references > 0) {
				/* use LS merging weights */
				n_bytes_prediction += (int)fwrite(SAI->merge_weights, sizeof(signed short), SAI->NB / 2, output_LF_file) * sizeof(signed short);
			}
		}
		else {
			/* use standard deviation */
			n_bytes_prediction += (int)fwrite(&SAI->stdd, sizeof(float), 1, output_LF_file) * sizeof(float);
		}
	}

	return;

}