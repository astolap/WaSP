#include "warping.hh"
#include "bitdepth.hh"
#include <cmath>

void warpView0_to_View1(view *view0, view *view1, unsigned short *&warpedColor, unsigned short *&warpedDepth, float *&DispTarg)
{

	/*this function forward warps from view0 to view1 for both color and depth*/

	float ddy = view0->y - view1->y;
	float ddx = view0->x - view1->x;

	unsigned short *AA1 = view0->color;
	unsigned short *DD1 = view0->depth;

	warpedColor = new unsigned short[view0->nr*view0->nc * 3]();
	warpedDepth = new unsigned short[view0->nr*view0->nc]();
	DispTarg = new float[view0->nr*view0->nc]();

	for (int ij = 0; ij < view0->nr*view0->nc; ij++){
		DispTarg[ij] = -1.0;
	}

	for (int ij = 0; ij < view0->nr*view0->nc; ij++)
	{


		float disp = ((float)DD1[ij] - (float)view0->min_inv_d) / (float)(1 << D_DEPTH);// pow(2, D_DEPTH);
		float DM_COL = disp*ddx;
		float DM_ROW = -disp*ddy;
		float disp0 = abs(DM_COL) + abs(DM_ROW);

		//if (view0->DM_ROW != NULL && view0->DM_COL != NULL)
		//{
		//	DM_COL = view0->DM_COL[ij];
		//	DM_ROW = view0->DM_ROW[ij];
		//	disp0 = 1 / (abs(DM_COL) + abs(DM_ROW)); /*for lenslet large disparity means further away ... */
		//}

		int ix = ij % view0->nr; //row
		int iy = (ij - ix) / view0->nr; //col

		int iynew = iy + (int)floor(DM_COL + 0.5);
		int ixnew = ix + (int)floor(DM_ROW + 0.5);

		if (iynew >= 0 && iynew < view0->nc && ixnew >= 0 && ixnew < view0->nr){
			int indnew = ixnew + iynew*view0->nr;
			if (DispTarg[indnew] < disp0){
				DispTarg[indnew] = disp0;
				warpedColor[indnew] = AA1[ij];
				warpedColor[indnew + view0->nr*view0->nc] = AA1[ij + view0->nr*view0->nc];
				warpedColor[indnew + 2 * view0->nr*view0->nc] = AA1[ij + 2 * view0->nr*view0->nc];
				warpedDepth[indnew] = DD1[ij];
			}
		}

	}

}