#include "updatelatent.h"

void mexFunction(
                 int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{

  
  // ovals - struct, the structure containing the ovals
  if(nrhs != 2)
    mexErrMsgTxt("Error this function takes two arguments");

  if(nlhs != 2)
    mexErrMsgTxt("Error this function returns two arguments");

  if(mxGetClassID(prhs[0]) != mxSTRUCT_CLASS)
    mexErrMsgTxt("Error model should be a structure");  

  double* A = mxGetPr(mxGetField(prhs[0], 0, "A"));
  double* beta = mxGetPr(mxGetField(prhs[0], 0, "beta"));
  double* tau = mxGetPr(mxGetField(prhs[0], 0, "tau"));
  int nData = (int)*mxGetPr(mxGetField(prhs[0], 0, "numData"));
  int dataDim = (int)*mxGetPr(mxGetField(prhs[0], 0, "dataDim"));
  int latentDim = (int)*mxGetPr(mxGetField(prhs[0], 0, "latentDim"));
  int FANoise = (int)*mxGetPr(mxGetField(prhs[0], 0, "FANoise"));

  if(mxGetClassID(prhs[1]) != mxDOUBLE_CLASS)
    mexErrMsgTxt("Error X should be DOUBLE");  

  double* X = mxGetPr(prhs[1]);

  int dims[3];
  dims[0] = nData;
  dims[1] = latentDim;
  plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
  double* sBar = mxGetPr(plhs[0]);

  dims[0] = latentDim;
  dims[1] = latentDim;
  dims[2] = nData;
  plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
  double* Sigma_s = mxGetPr(plhs[1]);


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
    
 
