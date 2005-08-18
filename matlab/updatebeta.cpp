#include "updatebeta.h"
void mexFunction(
                 int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{

  
  // ovals - struct, the structure containing the ovals
  if(nrhs != 2)
    mexErrMsgTxt("Error this function takes two arguments");

  if(nlhs != 1)
    mexErrMsgTxt("Error this function returns one argument");

  if(mxGetClassID(prhs[0]) != mxSTRUCT_CLASS)
    mexErrMsgTxt("Error model should be a structure");  

  double* A = mxGetPr(mxGetField(prhs[0], 0, "A"));
  double* sBar = mxGetPr(mxGetField(prhs[0], 0, "sBar"));
  double* Sigma_s = mxGetPr(mxGetField(prhs[0], 0, "Sigma_s"));
  int nData = (int)*mxGetPr(mxGetField(prhs[0], 0, "numData"));
  int dataDim = (int)*mxGetPr(mxGetField(prhs[0], 0, "dataDim"));
  int latentDim = (int)*mxGetPr(mxGetField(prhs[0], 0, "latentDim"));
  int FANoise = (int)*mxGetPr(mxGetField(prhs[0], 0, "FANoise"));

  if(mxGetClassID(prhs[1]) != mxDOUBLE_CLASS)
    mexErrMsgTxt("Error X should be DOUBLE");  

  double* X = mxGetPr(prhs[1]);

  int dims[2];
  dims[0] = 1;
  if(FANoise !=0) {
    dims[1] = dataDim;
  }
  else {
    dims[1] = 1;
  }

  plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);

  double* beta = mxGetPr(plhs[0]);

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



