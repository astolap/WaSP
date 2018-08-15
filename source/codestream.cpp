#include "codestream.hh"
#include "view.hh"
#include "minconf.hh"

#include <iostream>

void viewHeaderToCodestream(int &n_bytes_prediction, view *SAI, FILE *output_LF_file, const int yuv_transform_s) {

	minimal_config mconf = makeMinimalConfig(SAI);

	printf("size of minimal_config %i bytes\n", (int)sizeof(minimal_config));

	n_bytes_prediction += (int)fwrite(&mconf, sizeof(minimal_config), 1, output_LF_file) * sizeof(minimal_config);

	/* lets see what else needs to be written to bitstream */

	if (SAI->has_color_references) {
		unsigned char tmpNREF = (unsigned char)SAI->n_references;
		n_bytes_prediction += (int)fwrite(&tmpNREF, sizeof(unsigned char), 1, output_LF_file) * sizeof(unsigned char);
		for (int ij = 0; ij < SAI->n_references; ij++) {
			unsigned short nid = (unsigned short) *(SAI->references + ij);
			n_bytes_prediction += (int)fwrite(&nid, sizeof(unsigned short), 1, output_LF_file) * sizeof(unsigned short);
		}
	}

	if (SAI->has_depth_references) {
		unsigned char tmpNDREF = (unsigned char)SAI->n_depth_references;
		n_bytes_prediction += (int)fwrite(&tmpNDREF, sizeof(unsigned char), 1, output_LF_file) * sizeof(unsigned char);
		for (int ij = 0; ij < SAI->n_depth_references; ij++) {
			unsigned short nid = (unsigned short) *(SAI->depth_references + ij);
			n_bytes_prediction += (int)fwrite(&nid, sizeof(unsigned short), 1, output_LF_file) * sizeof(unsigned short);
		}
	}

	if (SAI->has_min_inv_depth) {
		unsigned short tmpminv = (unsigned short)SAI->min_inv_d;
		n_bytes_prediction += (int)fwrite(&tmpminv, sizeof(unsigned short), 1, output_LF_file) * sizeof(unsigned short);
	}

	if (!SAI->use_median) {
		if (SAI->stdd < 0.001) {
			if (SAI->has_color_references) {
				/* use LS merging weights */
				n_bytes_prediction += (int)fwrite(SAI->merge_weights, sizeof(signed short), SAI->NB / 2, output_LF_file) * sizeof(signed short);
			}
		}
		else {
			/* use standard deviation */
			n_bytes_prediction += (int)fwrite(&SAI->stdd, sizeof(float), 1, output_LF_file) * sizeof(float);
		}
	}

	if ( SAI->use_global_sparse ) {

		int32_t sparse_mask_binary_p1 = 0;
		int32_t sparse_mask_binary_p2 = 0;

		for (int ij = 0; ij < SAI->Ms; ij++) {

			if (SAI->sparse_mask[ij] > 0) {

				if (SAI->sparse_mask[ij]  <= 30) {
					sparse_mask_binary_p1 = sparse_mask_binary_p1 | 1 << ((int32_t)SAI->sparse_mask[ij]); // note: regressor indexing starts from 1
				}
				else {
					sparse_mask_binary_p2 = sparse_mask_binary_p2 | 1 << ((int32_t)SAI->sparse_mask[ij] - 31);
				}

			}

		}

		unsigned char tmpNNt = (unsigned char)SAI->NNt;
		unsigned char tmpMs = (unsigned char)SAI->Ms;

		n_bytes_prediction += (int)fwrite(&tmpNNt, sizeof(unsigned char), 1, output_LF_file) * sizeof(unsigned char);
		n_bytes_prediction += (int)fwrite(&tmpMs, sizeof(unsigned char), 1, output_LF_file) * sizeof(unsigned char);

		n_bytes_prediction += (int)fwrite(SAI->sparse_weights, sizeof(int32_t), SAI->Ms, output_LF_file) * sizeof(int32_t);
		n_bytes_prediction += (int)fwrite(&sparse_mask_binary_p1, sizeof(int32_t), 1, output_LF_file) * sizeof(int32_t);
		n_bytes_prediction += (int)fwrite(&sparse_mask_binary_p2, sizeof(int32_t), 1, output_LF_file) * sizeof(int32_t);

	}



	return;

}

void codestreamToViewHeader( int &n_bytes_prediction, view *SAI, FILE *input_LF, minimal_config &mconf ) {

	n_bytes_prediction += (int)fread(&mconf, sizeof(minimal_config), 1, input_LF)* sizeof(minimal_config);

	printf("size of minimal_config %i bytes\n", (int)sizeof(minimal_config));

	setup_form_minimal_config(&mconf, SAI);

	if (SAI->has_color_references) {

		unsigned char tmpNREF = 0;

		n_bytes_prediction += (int)fread(&tmpNREF, sizeof(unsigned char), 1, input_LF) * sizeof(unsigned char);

		SAI->n_references = tmpNREF;

		SAI->references = new int[SAI->n_references]();
		for (int ij = 0; ij < SAI->n_references; ij++) {
			unsigned short nid;
			n_bytes_prediction += (int)fread(&nid, sizeof(unsigned short), 1, input_LF)* sizeof(unsigned short);
			*(SAI->references + ij) = (int)nid;
		}
	}

	if (SAI->has_depth_references) {

		unsigned char tmpNDREF = 0;

		n_bytes_prediction += (int)fread(&tmpNDREF, sizeof(unsigned char), 1, input_LF) * sizeof(unsigned char);

		SAI->n_depth_references = tmpNDREF;


		SAI->depth_references = new int[SAI->n_depth_references]();
		for (int ij = 0; ij < SAI->n_depth_references; ij++) {
			unsigned short nid;
			n_bytes_prediction += (int)fread(&nid, sizeof(unsigned short), 1, input_LF)* sizeof(unsigned short);
			*(SAI->depth_references + ij) = (int)nid;
		}
	}

	if (SAI->has_min_inv_depth) {
		unsigned short tmpminv = (unsigned short)SAI->min_inv_d;
		n_bytes_prediction += (int)fread(&tmpminv, sizeof(unsigned short), 1, input_LF) * sizeof(unsigned short);
		SAI->min_inv_d = tmpminv;
	}

	SAI->NB = (1 << SAI->n_references)*SAI->n_references;

	if (!SAI->use_median) {
		if (SAI->stdd < 0.001) {
			if (SAI->n_references > 0) {
				SAI->merge_weights = new signed short[SAI->NB / 2]();
				n_bytes_prediction += (int)fread(SAI->merge_weights, sizeof(signed short), SAI->NB / 2, input_LF) * sizeof(signed short);
			}
		}
		else {
			n_bytes_prediction += (int)fread(&SAI->stdd, sizeof(float), 1, input_LF) * sizeof(float);
		}
	}

	int32_t sparse_mask_binary_p1 = 0;
	int32_t sparse_mask_binary_p2 = 0;

	if ( SAI->use_global_sparse ) {

		unsigned char tmpNNt = (unsigned char)SAI->NNt;
		unsigned char tmpMs = (unsigned char)SAI->Ms;

		n_bytes_prediction += (int)fread(&tmpNNt, sizeof(unsigned char), 1, input_LF) * sizeof(unsigned char);
		n_bytes_prediction += (int)fread(&tmpMs, sizeof(unsigned char), 1, input_LF) * sizeof(unsigned char);

		SAI->NNt = tmpNNt;
		SAI->Ms = tmpMs;

		SAI->sparse_weights = new int32_t[SAI->Ms]();
		n_bytes_prediction += (int)fread(SAI->sparse_weights, sizeof(int32_t), SAI->Ms, input_LF)* sizeof(int32_t);

		SAI->sparse_mask = new unsigned char[SAI->Ms]();

		n_bytes_prediction += (int)fread(&sparse_mask_binary_p1, sizeof(int32_t), 1, input_LF)* sizeof(int32_t);
		n_bytes_prediction += (int)fread(&sparse_mask_binary_p2, sizeof(int32_t), 1, input_LF)* sizeof(int32_t);
		
		int ik = 0;

		for (int ij = 0; ij < 64; ij++) {
			if (ij <= 30) {
				if ( ( sparse_mask_binary_p1 & (1 << ij) )>0 ) {
					SAI->sparse_mask[ik] = ij; ik++;
				}
			}
			else {
				if ( (sparse_mask_binary_p2 & (1 << (ij-31)))>0 ) {
					SAI->sparse_mask[ik] = ij; ik++;
				}
			}
		}
		
	}

	return;

}