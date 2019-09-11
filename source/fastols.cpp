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

#include "fastols.hh"

int32_t FastOLS_new(
    double **AAA, 
    double **Ydd, 
    int32_t *PredRegr0,
    double *PredTheta0,
    const int32_t Ms, 
    const int32_t MT, 
    const int32_t MPHI, 
    const int32_t N) {

  int32_t M, iM, iM1;
  double *B, *C, sigerr, *Ag, *g;
  double valm1, temp, crit, sabsval;
  int32_t p, j_p, i, j, k, itemp;

  double *AA = *AAA;
  double *Yd = *Ydd;

  double *PHI = new double[MT * MT]();
  double *PSI = new double[MT]();

  //int32_t startt = clock();

  /* make ATA. this is slow. */
  for (int32_t i1 = 0; i1 < MT; i1++) {
#pragma omp parallel for shared(i1)
    for (int32_t j1 = 0; j1 < MT; j1++) {
      for (int32_t ii = 0; ii < N; ii++) {
        *(PHI + i1 + j1 * MT) += (*(AA + ii + i1 * N)) * (*(AA + ii + j1 * N));
      }
    }
  }

  ///* try ATA with simple SSE optimization, down to around 9sec for 1080p*/
  //double* result = (double*)_aligned_malloc(2 * sizeof(double), 16);

  //__m128d x, y, z;

  //for (int32_t i1 = 0; i1 < MT; i1++) {
  //	for (int32_t j1 = 0; j1 < MT; j1++) {
  //		int32_t ii = 0;
  //		while (ii + 2 < N)
  //		{

  //			x = _mm_set_pd(*(AA + ii + i1*N), *(AA + ii + i1*N + 1));
  //			y = _mm_set_pd(*(AA + ii + j1*N), *(AA + ii + j1*N + 1));
  //			z = _mm_mul_pd(x, y);

  //			_mm_store_pd(result, z);

  //			*(PHI + i1 + j1*MT) += result[0] + result[1];

  //			ii += 2;
  //		}
  //		for (int32_t ee = ii; ee < N; ee++) {
  //			*(PHI + i1 + j1*MT) += (*(AA + ee + i1*N))*(*(AA + ee + j1*N));
  //		}
  //	}
  //}
  //_aligned_free(result);

  //std::cout << "time elapsed in ATA\t" << (int32_t)clock() - startt << "\n";

  //for (int32_t i1 = 0; i1 < MT; i1++){
  //	for (int32_t j1 = 0; j1 < MT; j1++){
  //		std::cout << *(ATA + i1 + j1*MT) << "\n";
  //	}
  //}

  //std::cout << "------------------------------------------------\n";

  /* make ATYd */
//#pragma omp parallel for
  for (int32_t i1 = 0; i1 < MT; i1++) {
    for (int32_t ii = 0; ii < N; ii++) {
      *(PSI + i1) += (*(AA + ii + i1 * N)) * (*(Yd + ii));
    }
  }

  delete[] (*AAA);
  *AAA = nullptr;

  /* YdTYd */
  double yd2 = 0;
//#pragma omp parallel for
  for (int32_t ii = 0; ii < N; ii++) {
      yd2 += (*(Yd + ii)) * (*(Yd + ii));
  }

  delete[] (*Ydd);
  *Ydd = nullptr;

  // Usage example: Ms= 3 says the sparsity (length of final predictor) and MT =42 tells how many regressors are available
  // Finally, MPHI = 63 tells the dimensions of the matrices, for getting linear indices in PHI
  M = MT + 1;
  //B = alocadoubleVector(M*M);// ((M+2)*(M+2));
  //C = alocadoubleVector(M*M);// ((M+2)*(M+2));
  //Ag = alocadoubleVector(M*M);// ((M+2)*(M+2));
  //g = alocadoubleVector(M);

  B = new double[M * M]();
  C = new double[M * M]();
  Ag = new double[M * M]();
  g = new double[M]();

  // Inputs: PHI is MTxMT, PSI is MTx1;
  // Outputs: PredRegr and PredTheta are also MTx1 but only first Ms entries are needed
  // Internal variables: B and C are (MT+1)x(MT+1) i.e. MxM

  B[MT + MT * M] = yd2;  //B[MT,MT] = yd2; // we start from B[0,0]
  for (iM = 0; iM < MT; iM++) {
    PredRegr0[iM] = iM;
    B[iM + MT * M] = PSI[iM];  // B[iM,MT] = PSI[iM];
    B[MT + iM * M] = PSI[iM];  //B[MT,iM] = PSI[iM];
    for (iM1 = 0; iM1 < MT; iM1++) {
      B[iM + iM1 * M] = PHI[iM + iM1 * MPHI];  //B[iM,iM1]=PHI[iM,iM1];
    }
  }
  for (iM = 0; iM < M; iM++)
    for (iM1 = 0; iM1 < M; iM1++)
      C[iM + iM1 * M] = 0;  //C[iM,iM1] = 0
  for (iM = 0; iM < MT; iM++)
    C[iM + iM * M] = 1;  // C[iM,iM] = 1;
  crit = B[MT + MT * M];  //crit = B[MT,MT];
  if (crit < 0.0000001) {
    //printf("crir, yd2 [%f] [%f] ", crit, yd2);
    i = 0;
    return i;
  }

  for (p = 0; p < Ms; p++) {
    valm1 = 0;
    j_p = 0;  // pick the max value in next loop
    for (j = p; j < MT; j++) {
      //if(B[j+j*M] > 0.00000000000000001)
      sigerr = B[j + MT * M] * B[j + MT * M] / B[j + j * M];  //sigerr  = B[j,MT]*B[j,MT]/B[j,j];
      //else
      //	sigerr = 0;
      if (sigerr > valm1) {
        valm1 = sigerr;
        j_p = j;
      }
    }  // j_p is the index of maximum
    crit = crit - valm1;
    itemp = PredRegr0[j_p];
    PredRegr0[j_p] = PredRegr0[p];
    PredRegr0[p] = itemp;
    for (j = p; j < M; j++) {
      //% interchange B(p:end,j_p) with B(p:end,p)
      temp = B[j + j_p * M];  //temp = B[j,j_p];
      B[j + j_p * M] = B[j + p * M];  //B[j,j_p] = B[j,p];
      B[j + p * M] = temp;  // B[j,p] = temp;
    }
    for (j = p; j < M; j++) {
      //% interchange B(j_p,p:end) with B(p,p:end)
      temp = B[j_p + j * M];  // temp = B[j_p,j];
      B[j_p + j * M] = B[p + j * M];  //B[j_p,j] = B[p,j];
      B[p + j * M] = temp;  //B[p,j] = temp;
    }

    //% Fast
    for (j = 0; j <= p - 1; j++) {
      // % interchange C(1:p-1,j_p) with C(1:p-1,p)
      temp = C[j + j_p * M];  // temp = C[j,j_p];
      C[j + j_p * M] = C[j + p * M];  //C[j,j_p] = C[j,p];
      C[j + p * M] = temp;  //C[j,p] = temp;
    }
    for (j = (p + 1); j < M; j++) {
      //if(B[p+p*M] > 0.00000000000000000000001)
      C[p + j * M] = B[p + j * M] / B[p + p * M];  //C[p,j] = B[p,j]/B[p,p];
      //else
      //C[p+j*M] = 0;
    }

    for (j = (p + 1); j < MT; j++)
      for (k = j; k <= MT; k++) {
        B[j + k * M] = B[j + k * M]
            - C[p + j * M] * C[p + k * M] * B[p + p * M];	//B[j,k] = B[j,k]-C[p,j]*C[p,k]*B[p,p];
      }
    for (j = (p + 1); j < MT; j++)
      for (k = j; k <= MT; k++) {
        B[k + j * M] = B[j + k * M];												 //B[k,j] = B[j,k];
      }
    //for j = (p+1):M
    //    for k = j:(M+1)
    //        B(j,k) = B(j,k)-C(p,j)*C(p,k)*B(p,p);
    //        %B(j,k) = B(j,k)-C(p,j)*B(p,k);
    //        B(k,j) = B(j,k);
    //    end
    //end
    //        for( iM=0; iM<Ms; iM++)
    //        {
    //        for( iM1=0; iM1<Ms; iM1++)
    //        	printf("C[%f] ",C[iM+iM1*M]);
    //        printf(" \n" );
    //        }
    // scanf("%d",&i);
  }  //% for( p=0; p<M; p++ )

  // final triangular backsolving
  for (i = 0; i < Ms; i++) {
    g[i] = C[i + MT * M];  //g[i] = C[i,MT];
    for (j = 0; j < Ms; j++)
      Ag[i + j * M] = C[i + j * M];  //Ag[i,j] = C[i,j];
  }
  PredTheta0[Ms - 1] = g[Ms - 1];
  for (i = Ms - 2; i >= 0; i--) {
    PredTheta0[i] = g[i];
    for (j = i + 1; j < Ms; j++)
      PredTheta0[i] = PredTheta0[i] - Ag[i + j * M] * PredTheta0[j];  //PredTheta[i] = PredTheta[i]- Ag[i,j]*PredTheta[j];
  }
  //printf("pred FASTOLS [%f][%f][%f][%f][%f][%f]\n",PredTheta0[0], PredTheta0[1], PredTheta0[2], PredTheta0[3], PredTheta0[4], PredTheta0[5]);
  // printf("pred [%d][%d][%d][%d][%d]\n",PredRegr0[0], PredRegr0[1], PredRegr0[2], PredRegr0[3], PredRegr0[4] );

  if (PredTheta0[0] != PredTheta0[0]) {	// if is nan
                                        //printf("PredTheta0[0]  is NaN\n");
    PredTheta0[0] = 1.0;
    for (i = 1; i < Ms; i++) {
      PredTheta0[i] = 0.0;
    }
  }

  sabsval = 0;
  for (i = 0; i < Ms; i++) {

    if (PredTheta0[i] != PredTheta0[i]) {		// if is nan
      PredTheta0[i] = 0.0;
      //printf("PredTheta0[%i]  is NaN\n", i);
    }

    if (PredTheta0[i] > 0)
      sabsval = sabsval + PredTheta0[i];
    else
      sabsval = sabsval - PredTheta0[i];
  }
  //printf("%f\n", sabsval);
  if (sabsval > 2 * Ms)  // if average coefficients are too high forget about intrpolation
      {
    PredTheta0[0] = C[0 + MT * M];  //g[0] = C[0,MT];

    //PredTheta0[0] = 1;

    if (PredTheta0[0] != PredTheta0[0]) {  // if is nan
                                           //printf("C[0+MT*M]  is NaN\n", i);
      PredTheta0[0] = 1.0;
    }

    // fix ???
    //if (abs(PredTheta0[0]) > 100)
    //PredTheta0[0] = PredTheta0[0] / abs(PredTheta0[0]);

    for (i = 1; i < Ms; i++) {
      PredTheta0[i] = 0.0;
    }

    i = 1;
  } else {
    i = Ms;
  }

  delete[] (B);
  delete[] (C);
  delete[] (Ag);
  delete[] (g);

  delete[] (PSI);
  delete[] (PHI);

  //printf("pred Ms [%d]",Ms);
  return i;

}