/*BSD 2-Clause License
* Copyright(c) 2019, Pekka Astola
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met :
* 
* 1. Redistributions of source code must retain the above copyright notice, this
* list of conditions and the following disclaimer.
*
* 2. Redistributions in binary form must reproduce the above copyright notice,
* this list of conditions and the following disclaimer in the documentation
* and/or other materials provided with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED.IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
*     SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*     CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
*     OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
*     OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <cstring>

#include "sparsefilter.hh"
#include "fastols.hh"
#include "bitdepth.hh"
#include "Eigen\Dense" /*only needed if you plan to use getSP_FILTER_EIGEN()
                            instead of FastOLS.*/

uint16_t *cropImage(
    const uint16_t *input_image,
    const uint32_t nr,
    const uint32_t nc,
    const uint32_t NNt) {

    int32_t nr_cropped = nr - 2 * NNt;
    int32_t nc_cropped = nc - 2 * NNt;

    uint16_t *input_image_cropped = new uint16_t[nr_cropped*nc_cropped]();

    for (uint32_t ir = NNt; ir < nr - NNt; ir++) {
        for (uint32_t ic = NNt; ic < nc - NNt; ic++) {
            int32_t offset_padded = ir + nr * ic;
            int32_t offset_cropped = (ir - NNt) + nr_cropped * (ic - NNt);

            *(input_image_cropped + offset_cropped) =
                *(input_image + offset_padded);

        }
    }

    return input_image_cropped;

}

uint16_t *padArrayUint16_t(
    const uint16_t *input_image, 
    const uint32_t nr,
    const uint32_t nc,
    const uint32_t NNt) {

    int32_t nr_padded = nr + 2 * NNt;
    int32_t nc_padded = nc + 2 * NNt;

    uint16_t *input_image_padded = new uint16_t[nr_padded*nc_padded]();

    for (uint32_t ir = NNt; ir < nr_padded - NNt; ir++) {
        for (uint32_t ic = NNt; ic < nc_padded - NNt; ic++) {

            int32_t offset_padded = ir + nr_padded * ic;
            int32_t offset = (ir - NNt) + nr * (ic - NNt);

            *(input_image_padded + offset_padded) =
                *(input_image + offset);

        }
    }

    /*top borders*/
    for (uint32_t ir = 0; ir < NNt; ir++) {
        for (uint32_t ic = NNt; ic < nc_padded - NNt; ic++) {

            int32_t offset_from = NNt + nr_padded * ic;
            int32_t offset_to = ir + nr_padded * ic;

            *(input_image_padded + offset_to) =
                *(input_image_padded + offset_from);

        }
    }

    /*bottom borders*/
    for (uint32_t ir = nr_padded - NNt; ir < nr_padded; ir++) {
        for (uint32_t ic = NNt; ic < nc_padded - NNt; ic++) {

            int32_t offset_from = nr_padded - NNt - 1 + nr_padded * ic;
            int32_t offset_to = ir + nr_padded * ic;

            *(input_image_padded + offset_to) =
                *(input_image_padded + offset_from);

        }
    }

    /*left borders*/
    for (uint32_t ir = 0; ir < nr_padded; ir++) {
        for (uint32_t ic = 0; ic < NNt; ic++) {

            int32_t offset_from = ir + nr_padded * NNt;
            int32_t offset_to = ir + nr_padded * ic;

            *(input_image_padded + offset_to) =
                *(input_image_padded + offset_from);

        }
    }

    /*right borders*/
    for (uint32_t ir = 0; ir < nr_padded; ir++) {
        for (uint32_t ic = nc_padded - NNt; ic < nc_padded; ic++) {

            int32_t offset_from = ir + nr_padded * (nc_padded - NNt - 1);
            int32_t offset_to = ir + nr_padded * ic;

            *(input_image_padded + offset_to) =
                *(input_image_padded + offset_from);

        }
    }

    return input_image_padded;

}

spfilter getGlobalSparseFilter(
    const uint16_t *original_image,
    const uint16_t *input_image,
    const int32_t nr,
    const int32_t nc,
    const int32_t NNt,
    const int32_t Ms,
    const double bias_term_value,
    const int32_t sub_sampling_factor) {

    int32_t MT = (NNt * 2 + 1) * (NNt * 2 + 1) + 1; /* number of regressors */

    //int32_t totalP = (nr - NNt * 2) * (nc - NNt * 2);

    //const int32_t skipv = sub_sampling_factor;
    //const int32_t skiph = sub_sampling_factor;

    //uint32_t Npp = ((nr - NNt * 2) / skipv + 1)*((nc - NNt * 2) / skipv + 1);

    int32_t Npp0 = (nr - NNt * 2) * (nc - NNt * 2);
    int32_t Npp = static_cast<int32_t>( ceil( static_cast<float>( Npp0 / sub_sampling_factor) ));

    double *AA0 = new double[Npp0 * MT]();
    double *Yd0 = new double[Npp0]();

    /*init bias column*/
    for (int32_t ii = 0; ii < Npp0; ii++) {
        *(AA0 + ii + (NNt * 2 + 1) * (NNt * 2 + 1) * Npp0) = bias_term_value;
    }

    int32_t iiu = 0;

    double Q = static_cast<double>( (1 << BIT_DEPTH) - 1);

    for (uint32_t ir = NNt; ir < nr - NNt; ir++) {
        for (uint32_t ic = NNt; ic < nc - NNt; ic++) {
            int32_t ai = 0;
            for (int32_t dy = -NNt; dy <= NNt; dy++) {
                for (int32_t dx = -NNt; dx <= NNt; dx++) {

                    int32_t offset = ir + dy + nr * (ic + dx);

                    /* get the desired Yd*/
                    if (dy == 0 && dx == 0) {
                        *(Yd0 + iiu) +=
                            ((double) *(original_image + offset)) / Q;	
                    }

                    /* get the regressors */
                    *(AA0 + iiu + ai * Npp0) += 
                        ((double) *(input_image + offset)) / Q;	

                    ai++;
                }
            }

            iiu++;
        }
    }

    double *AA;
    double *Yd;

    if (sub_sampling_factor > 1) {

        /* subsampled stored here */
        AA = new double[Npp * MT]();
        Yd = new double[Npp]();

        iiu = 0; /*index for subsampled array*/

        for (int32_t ii = 0; ii < Npp0; ii = ii + sub_sampling_factor) {

            *(Yd + iiu) = *(Yd0 + ii);

            for (int32_t ai = 0; ai < MT; ai++) {
                *(AA + iiu + ai * Npp) = *(AA0 + ii + ai * Npp0);
            }

            iiu++;

        }

        delete[](AA0);
        delete[](Yd0);
    }
    else { /*do nothing*/

        AA = AA0;
        Yd = Yd0;

    }

    int32_t *PredRegr0 = new int32_t[MT]();
    double *PredTheta0 = new double[MT]();

    int32_t Mtrue = FastOLS_new(
        &AA,
        &Yd,
        PredRegr0,
        PredTheta0,
        Ms,
        MT,
        MT,
        Npp);

    if (AA != nullptr) {
        delete[](AA);
    }
    if (Yd != nullptr) {
        delete[](Yd);
    }

    spfilter sparse_filter;

    for (int32_t ii = 0; ii < MT; ii++) {
        sparse_filter.regressor_indexes.push_back(PredRegr0[ii]);
        sparse_filter.filter_coefficients.push_back(PredTheta0[ii]);
    }

    for (int iu = 0; iu < MT; iu++) {
        printf("\n\t%d", sparse_filter.regressor_indexes.at(iu));
        printf("\n\t%f", sparse_filter.filter_coefficients.at(iu));
    }

    sparse_filter.Ms = Ms;
    sparse_filter.NNt = NNt;
    sparse_filter.bias_term_value = bias_term_value;

    return sparse_filter;
}

void quantize_and_reorder_spfilter(
    spfilter &sparse_filter) {

    std::vector<std::pair<int32_t, double>> sparsefilter_pair;

    for (int32_t ij = 0; ij < sparse_filter.Ms; ij++) {
        std::pair<int32_t, double> tmp_sp;
        tmp_sp.first = sparse_filter.regressor_indexes.at(ij);
        tmp_sp.second = sparse_filter.filter_coefficients.at(ij);
        sparsefilter_pair.push_back(tmp_sp);
    }

    /*ascending sort based on regressor index 
    (e.g., integer between 1 and 50 since we have 50 regressors */
    sort(sparsefilter_pair.begin(),
        sparsefilter_pair.end());

    sparse_filter.quantized_filter_coefficients.clear();
    sparse_filter.regressor_indexes.clear();

    double Q = static_cast<double>(1 << BIT_DEPTH_SPARSE);

    for (int32_t ii = 0; ii < sparse_filter.Ms; ii++) {

        sparse_filter.regressor_indexes.push_back(sparsefilter_pair.at(ii).first);

        int32_t quantized_coeff = static_cast<int32_t>(floor(sparsefilter_pair.at(ii).second * Q + 0.5));

        quantized_coeff = quantized_coeff > (1 << 15)-1 ? (1 << 15)-1 : quantized_coeff;
        quantized_coeff = quantized_coeff < -(1 << 15) ? -(1 << 15) : quantized_coeff;

        sparse_filter.quantized_filter_coefficients.push_back(static_cast<int16_t>(quantized_coeff));

    }

}

void dequantize_and_reorder_spfilter(
    spfilter &sparse_filter) {

    sparse_filter.filter_coefficients.clear();

    int32_t MT = (sparse_filter.NNt * 2 + 1) * (sparse_filter.NNt * 2 + 1) + 1;

    for (int32_t ii = 0; ii < MT; ii++) {
        sparse_filter.filter_coefficients.push_back(0.0);
    }

    double Q = static_cast<double>(1 << BIT_DEPTH_SPARSE);

    for (int32_t ii = 0; ii < sparse_filter.Ms; ii++) {

        double theta0 =
            static_cast<double>(sparse_filter.quantized_filter_coefficients.at(ii));

        sparse_filter.filter_coefficients.at(sparse_filter.regressor_indexes.at(ii)) = theta0 / Q;

    }
}

std::vector<double> applyGlobalSparseFilter(
    const uint16_t *input_image,
    const int32_t nr,
    const int32_t nc,
    const int32_t Ms,
    const int32_t NNt,
    const double bias_term_value,
    const std::vector<double> filter_coeffs) {

    std::vector<double> final_view;

    double Q = ((double)(1 << BIT_DEPTH) - 1);

    for (int32_t ii = 0; ii < nr * nc; ii++) {
        final_view.push_back(static_cast<double>(input_image[ii]) / Q);
    }

    for (int32_t rr = NNt; rr < nr - NNt; rr++) {
        for (int32_t cc = NNt; cc < nc - NNt; cc++) {

            final_view.at(rr + cc * nr) = 0.0;

            int32_t ee = 0;

            for (int32_t dy = -NNt; dy <= NNt; dy++) {
                for (int32_t dx = -NNt; dx <= NNt; dx++) {

                    double regr_value =
                        static_cast<double>(input_image[rr + dy + (cc + dx) * nr]) / Q;

                    final_view.at(rr + cc * nr) += filter_coeffs.at(ee)*regr_value;

                    ee++;
                }
            }

            /* bias term */
            final_view.at(rr + cc * nr) +=
                bias_term_value*filter_coeffs.at((2 * NNt + 1)* (2 * NNt + 1));

        }
    }

    for (int32_t ii = 0; ii < nr * nc; ii++) {
        final_view.at(ii) = final_view.at(ii) * Q;
    }

    return final_view;

}

spfilter getSP_FILTER_EIGEN(
    const uint16_t *original_image,
    const uint16_t *input_image,
    const int32_t nr,
    const int32_t nc,
    const int32_t NNt,
    const int32_t Ms,
    const double bias_term_value,
    const int32_t sub_sampling_factor) {

    int32_t MT = (NNt * 2 + 1) * (NNt * 2 + 1) + 1; /* number of regressors */                                      

    const int32_t skipv = sub_sampling_factor;
    const int32_t skiph = sub_sampling_factor;

    uint32_t Npp = ((nr - NNt * 2) / skipv + 1)*((nc - NNt * 2) / skipv + 1);

    Eigen::MatrixXf A = Eigen::MatrixXf::Zero(Npp, MT);
    Eigen::VectorXf b = Eigen::VectorXf::Zero(Npp);

    for (int32_t ii = 0; ii < Npp; ii++) {
        A(ii, MT - 1) = static_cast<float>(bias_term_value);
    }

    int32_t iiu = 0;

    float Q = static_cast<float>(1 << BIT_DEPTH) - 1;

    for (int32_t ic = NNt; ic < nc - NNt; ic++) {
        for (int32_t ir = NNt; ir < nr - NNt; ir++) {

            int32_t ai = 0;
            for (int32_t dy = -NNt; dy <= NNt; dy++) {
                for (int32_t dx = -NNt; dx <= NNt; dx++) {

                    int32_t offset = ir + dy + nr * (ic + dx);

                    /* get the desired Yd*/
                    if (dy == 0 && dx == 0) {
                        b(iiu) = static_cast<float>(*(original_image + offset)) / Q;
                    }

                    /* get the regressors */
                    A(iiu, ai) = static_cast<float>(*(input_image + offset)) / Q;

                    ai++;
                }
            }

            iiu++;
        }
    }

    /* full solution */
    Eigen::VectorXf X1 = A.colPivHouseholderQr().solve(b);

    std::vector<int> sparse_subset = get_SP_SUBSET(X1, Ms);

    /* solve LS for (sparse) subset */
    Eigen::MatrixXf A_sp = Eigen::MatrixXf::Zero(Npp, Ms);
    for (int jj = 0; jj<Ms; jj++) {
        for (int ii = 0; ii < Npp; ii++) {
            A_sp(ii, jj) = A(ii, sparse_subset.at(jj));
        }
    }

    A.resize(0, 0); /*free memory*/

    Eigen::VectorXf X1_sp = A_sp.colPivHouseholderQr().solve(b);

    bool failure = false;

    const int SP_THRESH = (1 << (15 - BIT_DEPTH_SPARSE - 1));
    //const int SP_THRESH = 3;

    for (int jj = 0; jj < Ms; jj++) {
        if (abs(X1_sp(jj)) > SP_THRESH) {
            failure = true;
            break;
        }
    }

    spfilter sparse_filter;

    /*pick only one regressor (the largest) */
    if (failure) {
        std::vector<int> sparse_subset_2 = get_SP_SUBSET(X1_sp, 1);

        /* solve LS for (sparse) subset */
        Eigen::MatrixXf A_sp_2 = Eigen::MatrixXf::Zero(Npp, 1);
        for (int ii = 0; ii < Npp; ii++) {
            A_sp_2(ii, 0) = A_sp(ii, sparse_subset_2.at(0));
        }

        Eigen::VectorXf X2_sp = A_sp_2.colPivHouseholderQr().solve(b);

        sparse_filter.regressor_indexes.push_back(sparse_subset.at(sparse_subset_2.at(0)));
        sparse_filter.filter_coefficients.push_back(X2_sp(0));

        int32_t iki = 0;
        while (iki<Ms - 1) {
            /* PROBLEM  with zeros, what if sparse_subset.at(sparse_subset_2.at(0))==0*/
            if (iki == sparse_filter.regressor_indexes.back()) {
                iki++;
            }
            sparse_filter.regressor_indexes.push_back(iki);
            sparse_filter.filter_coefficients.push_back(0);
        }

    }
    else {
        for (int32_t jj = 0; jj < Ms; jj++) {
            sparse_filter.regressor_indexes.push_back(sparse_subset.at(jj));
            sparse_filter.filter_coefficients.push_back(X1_sp(jj));
        }
    }

    for (int iu = 0; iu < MT; iu++) {
        printf("\n\t%d", sparse_filter.regressor_indexes.at(iu));
        printf("\n\t%f", sparse_filter.filter_coefficients.at(iu));
    }

    sparse_filter.Ms = Ms;
    sparse_filter.NNt = NNt;
    sparse_filter.bias_term_value = bias_term_value;

    return sparse_filter;
}

bool sortinrev(
    const std::pair<float, int> &a,
    const std::pair<float, int> &b)
{
    return (a.first > b.first);
}

std::vector<int> get_SP_SUBSET(
    const Eigen::VectorXf X1,
    const int Ms) {

    std::vector<std::pair<float, int>> filter_coeffs_abs;

    for (int ii = 0; ii < X1.size(); ii++) {

        std::pair<float, int> tmp;
        tmp.first = abs(X1(ii));
        tmp.second = ii;

        filter_coeffs_abs.push_back(tmp);

    }

    std::sort(filter_coeffs_abs.begin(), filter_coeffs_abs.end(), sortinrev);

    std::vector<int> sparse_subset;

    /* pick regressors corresponding to Ms largest coeffs */
    for (int ii = 0; ii < Ms; ii++) {
        sparse_subset.push_back(filter_coeffs_abs.at(ii).second);
    }

    return sparse_subset;

}