#include "merging.hh"
#include "medianfilter.hh"
#include "fastols.hh"
#include "bitdepth.hh"

void setBMask(view *view0) 
{

	/* sets the binary mask used to derive view availability in each of the MMM classes,
	size of the binary mask is [MMM x n_references] */

	int MMM = 1 << view0->n_references;// pow(2, view0->n_references);

	bool *bmask = new bool[MMM * view0->n_references]();

	(view0)->bmask = bmask;

	for (int ij = 0; ij < MMM; ij++) {

		int uu = ij;

		for (int ik = view0->n_references - 1; ik >= 0; ik--) {

			//if (floor(uu / pow(2, ik)) > 0)
			if (floor(uu / (1<<ik) ) > 0)
			{
				//uu = uu - pow(2, ik);
				uu = uu - (1<<ik);
				bmask[ij + ik * MMM] = 1;
			}

		}
	}
}

void getViewMergingLSWeights_N(view *view0, unsigned short **warpedColorViews, float **DispTargs, const unsigned short *original_color_view)
{

	/* This function puts the LS view merging weights into LSw */

	int n_references = (view0)->n_references;
	int nr = (view0)->nr;
	int nc = (view0)->nc;

	signed short *LScoeffs = (view0)->merge_weights;

	int MMM = 1 << n_references;// pow(2, n_references);
	
	bool *bmask = (view0)->bmask;

	unsigned short *seg_vp = (view0)->seg_vp;

	/* go through all regions, collect regressors from references */
	unsigned short **reference_view_pixels_in_classes = new unsigned short*[MMM * n_references]();

	/* also collect desired values */
	unsigned short **original_view_in_classes = new unsigned short*[MMM]();

	int *number_of_pixels_per_region = (view0)->number_of_pixels_per_region;

	for (int ij = 1; ij < MMM; ij++){ // region

		int NN = number_of_pixels_per_region[ij];

		original_view_in_classes[ij] = new unsigned short[NN * 3]();
		unsigned short *p3s = original_view_in_classes[ij];

		int jj = 0;

		for (int ii = 0; ii < nr*nc; ii++){
			if (seg_vp[ii] == ij){
				for (int icomp = 0; icomp < 3; icomp++)
					*(p3s + jj + icomp*NN) = *(original_color_view + ii + icomp*nr*nc);
				jj++;
			}
		}

		for (int ik = 0; ik < n_references; ik++){ // reference view

			if (bmask[ij + ik * MMM]){

				/* allocate pixels for region */
				reference_view_pixels_in_classes[ij + ik * MMM] = new unsigned short[NN * 3]();

				unsigned short *ps = reference_view_pixels_in_classes[ij + ik * MMM];
				unsigned short *pss = warpedColorViews[ik];


				jj = 0;

				for (int ii = 0; ii < nr*nc; ii++){
					if (seg_vp[ii] == ij){
						for (int icomp = 0; icomp < 3; icomp++)
							*(ps + jj + icomp*NN) = *(pss + ii + icomp*nr*nc);
						jj++;
					}
				}

			}

		}

	}


	/* run fastOLS on the classes */

	signed short *thetas = new signed short[MMM * n_references]();

	for (int ij = 1; ij < MMM; ij++){

		/* form A for this class, compute A'*A (phi)
		also compute A'*y (psi), where y is the desired data from the original view */

		int M = 0;

		for (int ik = 0; ik<n_references; ik++)
		{
			if (bmask[ij + MMM * ik])
				M++;
		}

		int N = number_of_pixels_per_region[ij] * 3; // number of rows in A

		double *AA = new double[N*M]();
		double *Yd = new double[N]();

		unsigned short *ps;

		int ikk = 0;

		for (int ik = 0; ik < n_references; ik++){
			if (bmask[ij + ik * MMM]){
				ps = reference_view_pixels_in_classes[ij + ik * MMM];
				for (int ii = 0; ii < N; ii++){
					*(AA + ii + ikk*N) = ((double)*(ps + ii)) / (double)((1 << BIT_DEPTH) - 1);// (pow(2, BIT_DEPTH) - 1);
				}
				ikk++;
			}
		}

		ps = original_view_in_classes[ij];

		for (int ii = 0; ii < N; ii++){
			*(Yd + ii) = ((double)*(ps + ii)) / (double)((1 << BIT_DEPTH) - 1);// (pow(2, BIT_DEPTH) - 1);
		}

		/* fastols */

		int *PredRegr0 = new int[M]();
		double *PredTheta0 = new double[M]();

		//int Mtrue = FastOLS(ATA, ATYd, YdTYd, PredRegr0, PredTheta0, M, M, M);

		int Mtrue = FastOLS_new(&AA, &Yd, PredRegr0, PredTheta0, M, M, M, N);

		if (AA != NULL) {
			delete[](AA);
		}
		if (Yd != NULL) {
			delete[](Yd);
		}

		/* establish the subset of reference views available for class */
		int *iks = new int[M]();
		int ee = 0;
		for (int ik = 0; ik < n_references; ik++){
			if (bmask[ij + ik * MMM]){
				*(iks + ee) = ik;
				ee++;
			}
		}

		for (int ii = 0; ii < M; ii++){
			thetas[ij + MMM * iks[PredRegr0[ii]]] = (signed short)floor(*(PredTheta0 + ii)*(signed short)(1 << BIT_DEPTH_MERGE) + 0.5);// pow(2, BIT_DEPTH_MERGE) + 0.5);
		}

		delete[](iks);

		delete[](PredRegr0);
		delete[](PredTheta0);
		


	}

	/* columnwise collecting of thetas */
	int in = 0;
	for (int ik = 0; ik < n_references; ik++){
		for (int ij = 0; ij < MMM; ij++){
			if (bmask[ij + ik * MMM]){
				LScoeffs[in] = thetas[ij + MMM * ik];
				in++;
			}
		}
	}

	delete[](thetas);

	for (int ij = 0; ij < MMM; ij++){

		delete[](original_view_in_classes[ij]);

		for (int ik = 0; ik < n_references; ik++){
			delete[](reference_view_pixels_in_classes[ij + MMM * ik]);
		}

	}

	delete[](original_view_in_classes);
	delete[](reference_view_pixels_in_classes);
	//delete[](seg_vp);
	//delete[](number_of_pixels_per_region);

}

void initSegVp(view *view0, float **DispTargs) {

	int nr = view0->nr;
	int nc = view0->nc;
	int n_references = view0->n_references;

	unsigned short *seg_vp = new unsigned short[nr*nc]();

	(view0)->seg_vp = seg_vp;

	int MMM = 1 << view0->n_references;// pow(2, (view0)->n_references);

	int *number_of_pixels_per_region = new int[MMM]();

	for (int ii = 0; ii < nr*nc; ii++) {

		unsigned short ci = 0;

		for (int ik = 0; ik < n_references; ik++) {
			float *pf = DispTargs[ik];
			if (*(pf + ii) > -1)
				ci = ci + (unsigned short)( 1 << ik );// pow(2, ik);
		}

		seg_vp[ii] = ci;

		number_of_pixels_per_region[ci]++;

	}

	view0->number_of_pixels_per_region = number_of_pixels_per_region;
}

void initViewW(view *view0, float **DispTargs) {

	/* sets some of the parameters for a view in the light view structure */

	(view0)->NB = (1 << view0->n_references)*view0->n_references; // (pow(2, (view0)->n_references)*(view0)->n_references);
	
	if ((view0)->merge_weights == NULL) {
		signed short *merge_weights = new signed short[(view0)->NB / 2]();
		(view0)->merge_weights = merge_weights;
	}

	float *merge_weights_float = new float[(view0)->NB]();

	(view0)->merge_weights_float = merge_weights_float;

	setBMask(view0);

	initSegVp(view0, DispTargs);

}

void getGeomWeight(view *view0, view *LF) {

	float stdd = view0->stdd;

	int MMM = 1 << view0->n_references;// pow(2, (view0)->n_references);

	double *thetas = new double[view0->NB]();

	bool *bmask = (view0)->bmask;

	for (int ii = 0; ii < MMM; ii++) {
		double sumw = 0;

		for (int ij = 0; ij < (view0)->n_references; ij++) {
			view *view1 = LF + (view0)->references[ij];
			double vdistance = (view0->x - view1->x)*(view0->x - view1->x) + (view0->y - view1->y)*(view0->y - view1->y);
			if (bmask[ii + ij * MMM]) {
				thetas[ii + ij * MMM] = exp(-(vdistance) / (2 * stdd*stdd));
				sumw = sumw + thetas[ii + ij * MMM];
			}
		}
		for (int ij = 0; ij < (view0)->n_references; ij++)
			thetas[ii + ij * MMM] = thetas[ii + ij * MMM] / sumw;
	}

	/* columnwise collecting of thetas */
	signed short *LScoeffs = (view0)->merge_weights;
	int in = 0;
	for (int ik = 0; ik < (view0)->n_references; ik++) {
		for (int ij = 0; ij < MMM; ij++) {
			if (bmask[ij + ik * MMM]) {
				LScoeffs[in] = (signed short)floor(thetas[ij + MMM * ik] * (signed short)(1 << BIT_DEPTH_MERGE) + 0.5);// pow(2, BIT_DEPTH_MERGE) + 0.5);
				in++;
			}
		}
	}

	delete[](thetas);

}

void mergeMedian_N(unsigned short **warpedColorViews, float **DispTargs, view *view0, const int ncomponents) {

	unsigned short *AA2 = new unsigned short[view0->nr*view0->nc*ncomponents]();

#pragma omp parallel for
	for (int ii = 0; ii < view0->nr*view0->nc; ii++) {
		for (int icomp = 0; icomp < ncomponents; icomp++) {
			std::vector< unsigned short > vals;
			for (int ik = 0; ik < view0->n_references; ik++) {

				float *pf = DispTargs[ik];
				unsigned short *ps = warpedColorViews[ik];

				if (*(pf + ii) > -1) {

					vals.push_back(*(ps + ii + icomp*view0->nr*view0->nc));

				}

			}

			if (vals.size() > 0) {
				*(AA2 + ii + icomp*view0->nr*view0->nc) = getMedian(vals);
			}
		}
	}

	memcpy(view0->color, AA2, sizeof(unsigned short)*view0->nr*view0->nc * 3);
	delete[](AA2);

}

void mergeWarped_N(unsigned short **warpedColorViews, float **DispTargs, view *view0, const int ncomponents)
{

	int MMM = 1 << view0->n_references;// pow(2, (view0)->n_references);

	bool *bmask = (view0)->bmask;
	
	int uu = 0;

	for (int ii = 0; ii < MMM * (view0)->n_references; ii++){
		if (bmask[ii]){
			(view0)->merge_weights_float[ii] = ((float)(view0)->merge_weights[uu++]) / (float)(1 << BIT_DEPTH_MERGE);// pow(2, BIT_DEPTH_MERGE);
			//std::cout << (view0)->merge_weights_float[ii] << "\t";
			//printf("%f\t", (view0)->merge_weights_float[ii]);
		}
		else{
			(view0)->merge_weights_float[ii] = 0.0;
		}
	}

	int nr = (view0)->nr;
	int nc = (view0)->nc;
	int n_views = (view0)->n_references;
	float *LSw = (view0)->merge_weights_float;

	unsigned short *seg_vp = (view0)->seg_vp;

	float *AA1 = new float[nr*nc*ncomponents]();
	unsigned short *AA2 = new unsigned short[nr*nc*ncomponents]();

	for (int ii = 0; ii < nr*nc; ii++){

		int ci = seg_vp[ii];

		for (int ik = 0; ik < n_views; ik++){
			unsigned short *ps = warpedColorViews[ik];
			for (int icomp = 0; icomp < 3; icomp++){
				AA1[ii + icomp*nr*nc] = AA1[ii + icomp*nr*nc] + LSw[ci + ik * MMM] * ((float)(*(ps + ii + icomp*nr*nc)));
			}
		}

		for (int icomp = 0; icomp < 3; icomp++){
			if (AA1[ii + icomp*nr*nc] < 0)
				AA1[ii + icomp*nr*nc] = 0;
			if (AA1[ii + icomp*nr*nc] > (1 << BIT_DEPTH) - 1)//(pow(2, BIT_DEPTH) - 1))
				AA1[ii + icomp*nr*nc] = (1 << BIT_DEPTH) - 1;// (pow(2, BIT_DEPTH) - 1);

			AA2[ii + icomp*nr*nc] = (unsigned short)(floor(AA1[ii + icomp*nr*nc]));

		}
	}

	memcpy((view0)->color, AA2, sizeof(unsigned short)*nr*nc*ncomponents);

	//delete[](seg_vp);
	//delete[](bmask);
	delete[](AA1);
	delete[](AA2);

}
