#include "mex.h"
#include <cmath>
#include "lapack.h"
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
  const mxArray* pmxTau = mexGetArrayPtr("TAU", "global");
  double* tau = mxGetPr(pmxTau);


//for n = 1:NDATA
//end
  int lda = latentDim;
  int length = latentDim;
  int info = 0;
  int* ipiv = (int*)mxMalloc(length*sizeof(int));
  int order = latentDim;
  int lwork = order*16;
  double* work = (double*)mxMalloc(lwork*sizeof(double));

  for(int n = 0; n < nData; n++) {
    //   invSigma_s = diag(TAU(n, :)) + ATBA;
    for(int j = 0; j < latentDim; j++) {
      for(int j2 = 0; j2 < latentDim; j2++) {
        if(j2==j) {
          // Add the diagonal term
          Sigma_s[j + j2*latentDim + n*latentDim*latentDim] = tau[n + j*nData];
        }
        else {
          Sigma_s[j + j2*latentDim + n*latentDim*latentDim] = 0;
        }
        double temp = 0;
        if (FANoise != 0){
          for(int i = 0; i < dataDim; i++) {
            temp += A[i + j2*dataDim]*A[i + j*dataDim]*beta[i];
          }
        }
        else {
          for(int i = 0; i < dataDim; i++) {
            temp += A[i + j2*dataDim]*A[i + j*dataDim];
          }
          temp *= beta[0];		
        }
        // This is really inv(Sigma_s) but it is stored here for convenience
		Sigma_s[j + j2*latentDim + n*latentDim*latentDim] += temp;
      }
    }
    
    // It is not being done by cholesky decomposition in the c++ code
    // but it should be
    // C = chol(invSigma_s);
    // Cinv = eye(LATENTDIM)/C;
    // SIGMA_S(:, :, n) = Cinv*Cinv'; 
    
    // create inverse first by lu decomposition of input    
    
    // call lapack
    dgetrf(latentDim, latentDim, 
            Sigma_s+n*latentDim*latentDim, lda, ipiv, info);
    if(info != 0)
      mexErrMsgTxt("Problems in lu factorisation of matrix");
    
    info = 0;
    // peform the matrix inversion.
    dgetri(order, Sigma_s+n*latentDim*latentDim, lda, 
            ipiv, work, lwork, info);
    // check for successfull inverse
    if(info > 0)
      mexErrMsgTxt("Matrix is singular");
    else if(info < 0)
      mexErrMsgTxt("Problem in matrix inverse");
    
    
    // SBAR(n, :) = (X(n, :).*BETA)*A*SIGMA_S(:, :, n);
    
    if(FANoise != 0) {
      for(int j = 0; j < latentDim; j++) {
        sBar[n + j*nData] = 0;
        for(int j2 = 0; j2 < latentDim; j2++) {
          for(int i = 0; i < dataDim; i++) {
            sBar[n + j*nData] += X[n + i*nData]*beta[i]*A[i + j2*dataDim]*Sigma_s[j + j2*latentDim + n*latentDim*latentDim];  
          }
        }
      }
    }
    else {
      for(int j = 0; j < latentDim; j++) {
        sBar[n + j*nData] = 0;
        for(int j2 = 0; j2 < latentDim; j2++) {
          for(int i = 0; i < dataDim; i++) {
            sBar[n + j*nData] += X[n + i*nData]*beta[0]*A[i + j2*dataDim]*Sigma_s[j + j2*latentDim + n*latentDim*latentDim];  
          }
        }
      }
    }
  }
  mxFree(work);
  mxFree(ipiv);

}
    
 
