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
#include "merging.hh"
#include "medianfilter.hh"
#include "fastols.hh"
#include "bitdepth.hh"
#include "warping.hh"

#include <cstring>
#include <cmath>

void init_warping_arrays(
    const int32_t N,
    uint16_t **&warped_texture_views,
    uint16_t **&warped_depth_views,
    float **&DispTargs,
    const int32_t nr,
    const int32_t nc,
    const int32_t ncomp) {

    /* holds partial warped views for ii */
    warped_texture_views = new uint16_t*[N]();
    warped_depth_views = new uint16_t*[N]();
    DispTargs = new float*[N]();

    for (int32_t ij = 0; ij < N; ij++) {
        warped_texture_views[ij] = new uint16_t[nr*nc*ncomp]();
        warped_depth_views[ij] = new uint16_t[nr*nc]();
        DispTargs[ij] = new float[nr*nc]();
    }
}

void clean_warping_arrays(
    const int32_t N,
    uint16_t **warped_texture_views,
    uint16_t **warped_depth_views,
    float **DispTargs) {

    /* clean */
    for (int32_t ij = 0; ij < N; ij++) {
        delete[](warped_texture_views[ij]);
        delete[](warped_depth_views[ij]);
        delete[](DispTargs[ij]);
    }

    delete[](warped_texture_views);
    delete[](warped_depth_views);
    delete[](DispTargs);
}

void setBMask(view *view0) {
  /* sets the binary mask used to derive view availability in each of the MMM classes,
   size of the binary mask is [MMM x n_references] */
  int32_t MMM = 1 << view0->n_references;  // pow(2, view0->n_references);
  bool *bmask = new bool[MMM * view0->n_references]();
  (view0)->bmask = bmask;
  for (int32_t ij = 0; ij < MMM; ij++) {
    int32_t uu = ij;
    for (int32_t ik = view0->n_references - 1; ik >= 0; ik--) {
      //if (floor(uu / pow(2, ik)) > 0)
      if (floor(uu / (1 << ik)) > 0) {
        //uu = uu - pow(2, ik);
        uu = uu - (1 << ik);
        bmask[ij + ik * MMM] = 1;
      }
    }
  }
}

void initSegVp(view *view0, float **DispTargs) {

  int32_t nr = view0->nr;
  int32_t nc = view0->nc;
  int32_t n_references = view0->n_references;

  uint16_t *seg_vp = new uint16_t[nr * nc]();

  (view0)->seg_vp = seg_vp;

  int32_t MMM = 1 << view0->n_references;  // pow(2, (view0)->n_references);

  int32_t *number_of_pixels_per_region = new int32_t[MMM]();

  for (int32_t ii = 0; ii < nr * nc; ii++) {

    uint16_t ci = 0;

    for (int32_t ik = 0; ik < n_references; ik++) {
      float *pf = DispTargs[ik];
      if (*(pf + ii) > INIT_DISPARITY_VALUE)
        ci = ci + (uint16_t) (1 << ik);  // pow(2, ik);
    }

    seg_vp[ii] = ci;

    number_of_pixels_per_region[ci]++;

  }

  view0->number_of_pixels_per_region = number_of_pixels_per_region;
}

void initViewW(view *view0, float **DispTargs) {

    /* sets some of the parameters for a view in the light view structure */

    int32_t MMM = (1 << view0->n_references);

    (view0)->NB = MMM * view0->n_references;

    if ((view0)->merge_weights == nullptr) {
        int16_t *merge_weights = 
            new int16_t[view0->ncomp*(view0->NB / 2)]();
        (view0)->merge_weights = merge_weights;
    }

    double *merge_weights_double = new double[MMM*view0->n_references]();

    (view0)->merge_weights_double = merge_weights_double;

    setBMask(view0);

    initSegVp(view0, DispTargs);

}

void mergeMedian_N(uint16_t **warpedColorViews, float **DispTargs,
                   view *view0, const int32_t ncomponents) {

  uint16_t *AA2 =
      new uint16_t[view0->nr * view0->nc * ncomponents]();

//#pragma omp parallel for
  for (int32_t ii = 0; ii < view0->nr * view0->nc; ii++) {
    for (int32_t icomp = 0; icomp < ncomponents; icomp++) {
      std::vector<uint16_t> vals;
      for (int32_t ik = 0; ik < view0->n_references; ik++) {

        float *pf = DispTargs[ik];
        uint16_t *ps = warpedColorViews[ik];

        if (*(pf + ii) > INIT_DISPARITY_VALUE) {

          vals.push_back(*(ps + ii + icomp * view0->nr * view0->nc));

        }

      }

      if (vals.size() > 0) {
        *(AA2 + ii + icomp * view0->nr * view0->nc) = getMedian(vals);
      }
    }
  }

  memcpy(view0->color, AA2, sizeof(uint16_t) * view0->nr * view0->nc * 3);
  delete[] (AA2);

}

void printMatrix(
    const double *matrix,
    const int32_t nr,
    const int32_t nc) {

    printf("----------\n");
    for (int32_t r = 0; r < nr; r++) {
        for (int32_t c = 0; c < nc; c++) {
            int32_t lin_ind = c*nr + r;
            printf("%2.4f\t", *(matrix + lin_ind));
        }
        printf("\n");
    }
    printf("----------\n");

}

void mergeWarped_N_icomp(
    uint16_t **warpedColorViews,
    view *view0,
    const int32_t icomp) {

    int32_t MMM = 1 << view0->n_references;  // pow(2, (view0)->n_references);

    int32_t N_LS = MMM*view0->n_references;

    bool *bmask = view0->bmask;

    int32_t uu = icomp*((MMM*view0->n_references) / 2);

    for (int32_t ii = 0; ii < N_LS; ii++) {
        if (bmask[ii]) {
            view0->merge_weights_double[ii] =
                ((double)(view0)->merge_weights[uu++])
                / (double)(1 << BIT_DEPTH_MERGE);
        }
        else {
            (view0)->merge_weights_double[ii] = 0.0;
        }
    }

    /*printMatrix(
        view0->merge_weights_double,
        MMM,
        view0->n_references);*/

    int32_t nr = view0->nr;
    int32_t nc = view0->nc;
    int32_t n_views = view0->n_references;

    double *LSw = view0->merge_weights_double;

    uint16_t *seg_vp = view0->seg_vp;

    double *AA1 = new double[nr * nc]();
    uint16_t *AA2 = new uint16_t[nr * nc]();

    for (int32_t ii = 0; ii < nr * nc; ii++) {

        int32_t ci = seg_vp[ii]; /* occlusion class index and row index in thetas */

        for (int32_t ik = 0; ik < n_views; ik++) {
            uint16_t *ps = warpedColorViews[ik];
            AA1[ii] +=
                LSw[ci + ik * MMM] * ((double)(*(ps + ii + icomp * nr * nc)));
        }

        if (AA1[ii] < 0)
            AA1[ii] = 0;
        if (AA1[ii] > (1 << BIT_DEPTH) - 1)
            AA1[ii] = (1 << BIT_DEPTH) - 1;

        AA2[ii] = (uint16_t)(floor(AA1[ii] + 0.5));

    }

    memcpy(
        view0->color + icomp*nr*nc,
        AA2,
        sizeof(uint16_t) * nr * nc);

    delete[](AA1);
    delete[](AA2);

}

void getViewMergingLSWeights_icomp(
    view *view0,
    uint16_t **warpedColorViews,
    const uint16_t *original_color_view,
    const int32_t icomp) {

    /* This function puts the LS view merging weights into LSw */

    int32_t n_references = (view0)->n_references;
    int32_t nr = (view0)->nr;
    int32_t nc = (view0)->nc;

    int32_t MMM = 1 << n_references;  // pow(2, n_references);

    int16_t *LScoeffs = view0->merge_weights;

    bool *bmask = (view0)->bmask;

    uint16_t *seg_vp = (view0)->seg_vp;

    /* go through all regions, collect regressors from references */
    uint16_t **reference_view_pixels_in_classes = new uint16_t*[MMM
        * n_references]();

    /* also collect desired values */
    uint16_t **original_view_in_classes = new uint16_t*[MMM]();

    int32_t *number_of_pixels_per_region = (view0)->number_of_pixels_per_region;

    for (int32_t ij = 1; ij < MMM; ij++) {  // region

        int32_t NN = number_of_pixels_per_region[ij];

        original_view_in_classes[ij] = new uint16_t[NN]();
        uint16_t *p3s = original_view_in_classes[ij];

        int32_t jj = 0;

        for (int32_t ii = 0; ii < nr * nc; ii++) {
            if (seg_vp[ii] == ij) {
                *(p3s + jj) = *(original_color_view + ii + icomp*nr*nc);
                jj++;
            }
        }

        for (int32_t ik = 0; ik < n_references; ik++) {  // reference view

            if (bmask[ij + ik * MMM]) {

                /* allocate pixels for region */
                reference_view_pixels_in_classes[ij + ik * MMM] =
                    new uint16_t[NN]();

                uint16_t *ps = reference_view_pixels_in_classes[ij + ik * MMM];
                uint16_t *pss = warpedColorViews[ik];

                jj = 0;

                for (int32_t ii = 0; ii < nr * nc; ii++) {
                    if (seg_vp[ii] == ij) {
                        *(ps + jj) = *(pss + ii + icomp*nr*nc);
                        jj++;
                    }
                }
            }
        }
    }

    /* run fastOLS on the classes */
    int16_t *thetas = new int16_t[MMM * n_references]();
    for (int32_t ij = 1; ij < MMM; ij++) {
        /* form A for this class, compute A'*A (phi)
        also compute A'*y (psi), where y is the desired data from the original view */

        int32_t M = 0; /* number of active reference views for class ij */

        for (int32_t ik = 0; ik < n_references; ik++) {
            if (bmask[ij + MMM * ik])
                M++;
        }

        int32_t N = number_of_pixels_per_region[ij];  // number of rows in A

        double *AA = new double[N * M](); 
        double *Yd = new double[N]();

        uint16_t *ps;

        int32_t ikk = 0;

        for (int32_t ik = 0; ik < n_references; ik++) {
            if (bmask[ij + ik * MMM]) {
                ps = reference_view_pixels_in_classes[ij + ik * MMM];
                for (int32_t ii = 0; ii < N; ii++) {
                    *(AA + ii + ikk * N) = ((double) *(ps + ii))
                        / (double)((1 << BIT_DEPTH) - 1);
                }
                ikk++;
            }
        }

        ps = original_view_in_classes[ij];

        for (int32_t ii = 0; ii < N; ii++) {
            *(Yd + ii) = ((double) *(ps + ii))
                / (double)((1 << BIT_DEPTH) - 1);
        }

        /* fastols */

        int32_t *PredRegr0 = new int32_t[M]();
        double *PredTheta0 = new double[M]();

        //int32_t Mtrue = FastOLS(ATA, ATYd, YdTYd, PredRegr0, PredTheta0, M, M, M);

        int32_t Mtrue = FastOLS_new(
            &AA,
            &Yd,
            PredRegr0,
            PredTheta0,
            M, M, M, N);

        if (AA != nullptr) {
            delete[](AA);
        }
        if (Yd != nullptr) {
            delete[](Yd);
        }

        /* establish the subset of reference views available for class */
        int32_t *iks = new int32_t[M]();
        int32_t ee = 0;
        for (int32_t ik = 0; ik < n_references; ik++) {
            if (bmask[ij + ik * MMM]) {
                *(iks + ee) = ik;
                ee++;
            }
        }

        for (int32_t ii = 0; ii < M; ii++) {
            double quant_theta = floor(*(PredTheta0 + ii) * static_cast<double>(1 << BIT_DEPTH_MERGE) + 0.5);
            double max_t_val = static_cast<double>((1 << (16 - 1)));
            quant_theta = quant_theta > max_t_val ? max_t_val : quant_theta;
            quant_theta = quant_theta < -max_t_val ? -max_t_val : quant_theta;
            thetas[ij + MMM * iks[PredRegr0[ii]]] = static_cast<int16_t>(quant_theta);  // pow(2, BIT_DEPTH_MERGE) + 0.5);
        }

        delete[](iks);

        delete[](PredRegr0);
        delete[](PredTheta0);

    }

    /* columnwise collecting of thetas */
    int32_t N_LS = (MMM*n_references) / 2;
    int32_t in = N_LS*icomp;
    for (int32_t ik = 0; ik < n_references; ik++) {
        for (int32_t ij = 0; ij < MMM; ij++) {
            if (bmask[ij + ik * MMM]) {
                LScoeffs[in] = thetas[ij + MMM * ik];
                in++;
            }
        }
    }

    delete[](thetas);

    for (int32_t ij = 0; ij < MMM; ij++) {

        delete[](original_view_in_classes[ij]);

        for (int32_t ik = 0; ik < n_references; ik++) {
            delete[](reference_view_pixels_in_classes[ij + MMM * ik]);
        }

    }

    delete[](original_view_in_classes);
    delete[](reference_view_pixels_in_classes);
    //delete[](seg_vp);
    //delete[](number_of_pixels_per_region);

}

void getGeomWeight_icomp(
    view *view0,
    view *LF,
    const int32_t icomp) {

    float stdd = view0->stdd;

    int32_t MMM = 1 << view0->n_references;

    double *thetas = new double[MMM*(view0->n_references)]();

    bool *bmask = (view0)->bmask;

    for (int32_t ii = 0; ii < MMM; ii++) {
        double sumw = 0;

        for (int32_t ij = 0; ij < (view0)->n_references; ij++) {
            view *view1 = LF + (view0)->references[ij];
            double vdistance = (view0->x - view1->x) * (view0->x - view1->x)
                + (view0->y - view1->y) * (view0->y - view1->y);
            if (bmask[ii + ij * MMM]) {
                thetas[ii + ij * MMM] = exp(-(vdistance) / (2 * stdd * stdd));
                sumw = sumw + thetas[ii + ij * MMM];
            }
        }
        for (int32_t ij = 0; ij < (view0)->n_references; ij++)
            thetas[ii + ij * MMM] = thetas[ii + ij * MMM] / sumw;
    }

    /* columnwise collecting of thetas */
    int16_t *LScoeffs = (view0)->merge_weights;
    int32_t in = icomp*((MMM*view0->n_references) / 2);
    for (int32_t ik = 0; ik < (view0)->n_references; ik++) {
        for (int32_t ij = 0; ij < MMM; ij++) {
            if (bmask[ij + ik * MMM]) {
                LScoeffs[in] = (int16_t)floor(
                    thetas[ij + MMM * ik] * (int16_t)(1 << BIT_DEPTH_MERGE)
                    + 0.5);  // pow(2, BIT_DEPTH_MERGE) + 0.5);
                in++;
            }
        }
    }

    delete[](thetas);

}