#include "mex.h"
#include <cmath>
//#include <cassert>
// Mex function to update beta - this was the slowest of the update functions


void mexFunction(
                 int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{

  
  // ovals - struct, the structure containing the ovals
  if(nlhs != 0 || nrhs != 0)
    mexErrMsgTxt("Error this function takes and returns no arguments");

  const mxArray* pmxBeta = mexGetArrayPtr("BETA", "global");
  double*  beta = mxGetPr(pmxBeta);
  const mxArray* pmxNdata = mexGetArrayPtr("NDATA", "global");
  int nData = (int)*mxGetPr(pmxNdata);
  const mxArray* pmxDataDim = mexGetArrayPtr("DATADIM", "global");
  int dataDim = (int)*mxGetPr(pmxDataDim);
  const mxArray* pmxLatentDim = mexGetArrayPtr("LATENTDIM", "global");
  int latentDim = (int)*mxGetPr(pmxLatentDim);
  const mxArray* pmxX = mexGetArrayPtr("X", "global");
  double* X = mxGetPr(pmxX);
  const mxArray* pmxA = mexGetArrayPtr("A", "global");
  double* A = mxGetPr(pmxA);
  const mxArray* pmxSBar = mexGetArrayPtr("SBAR", "global");
  double* sBar = mxGetPr(pmxSBar);
  const mxArray* pmxSigma_s = mexGetArrayPtr("SIGMA_S", "global");
  double* Sigma_s = mxGetPr(pmxSigma_s);
  const mxArray* pmxFANoise = mexGetArrayPtr("FANOISE", "global");
  int FANoise = (int)*mxGetPr(pmxFANoise);


  if(FANoise !=0) {//  Factor analysis type noise model
    for(int i = 0; i < dataDim; i++) {
      beta[i] = 0;
      for(int n = 0; n < nData; n++) {
        double sbarTA = 0;
        for(int j = 0; j < latentDim; j++) {
          sbarTA += sBar[n+ j*nData]*A[i + j*dataDim];
          double temp = 0;
          for(int j2 = 0; j2 < latentDim; j2++) // this adds ATSigma_sA
            temp += A[i + j2*dataDim]
              *Sigma_s[j2 + j*latentDim + n*latentDim*latentDim];
          beta[i] += temp*A[i + j*dataDim];
              
        }
        double expectedOutputDiff = sbarTA - X[n + i*nData]; 
        beta[i] = beta[i] +  expectedOutputDiff*expectedOutputDiff;
      }
      beta[i] = (double)nData/beta[i];
    }
  }
  else { // PPCA type noise model
    beta[0] = 0;
    for(int i = 0; i < dataDim; i++) {
      for(int n = 0; n < nData; n++) {
        double sbarTA = 0;
        for(int j = 0; j < latentDim; j++) {
          sbarTA += sBar[n+ j*nData]*A[i + j*dataDim];
          for(int j2 = 0; j2 < latentDim; j2++) // this adds ATSigma_sA
            beta[0] += A[i + j2*dataDim]
              *Sigma_s[j2 + j*latentDim + n*latentDim*latentDim]
              *A[i + j*dataDim];
        }
        double expectedOutputDiff = sbarTA - X[n + i*nData]; 
        beta[0] = beta[0] +  expectedOutputDiff*expectedOutputDiff;
      }
    }
    beta[0] = (double)nData*(double)dataDim/beta[0];
  }
}



