#include "mex.h"
#include <cmath>
#include <cassert>
#include <cstdlib>


//using namespace std;



void randperm(int number, int order[])
{
  for(int i = 0; i<number; i++)
    order[i] = i;
  for(int j = 1; j<number; j++){
    int location = rand()%j;
    int tempint = 0;
    tempint = order[j];
    order[j] = order[location];
    order[location] = tempint;
  }
}


void mfupdatelatent(double XULinv[],
                    double V[],
                    double tau[],
                    double beta,
                    double lambda[],
                    double sbar[],
                    double Sigma_s[],
                    const int ndata,
                    const int latentDim,
                    const int dataDim,
                    const int maxCount,
                    const double tolerance)
                  
{

  double *Bl2 = (double *)mxCalloc(dataDim, sizeof(double));
  double *V2 = (double *)mxCalloc(dataDim*latentDim, sizeof(double));
  double *vTs = (double *)mxCalloc(dataDim, sizeof(double));
  double *Bl2V2 = 
    (double *)mxCalloc(dataDim*latentDim, sizeof(double));
  int *order = (int *)mxCalloc(latentDim, sizeof(int));
  
  for(int i = 0; i < dataDim; i++){
    Bl2[i] = beta*lambda[i]*lambda[i];
    for(int j = 0; j < latentDim; j++){
      V2[i + dataDim*j] = V[j + latentDim*i]*V[j + latentDim*i];
      Bl2V2[i + dataDim*j] = V2[i + dataDim*j]*Bl2[i];
    }
  }
  for(int n = 0; n < ndata; n++){
    float maxsbardiff = 0.0;
    float sbardiff = 0.0;
    int counter = 0;
    do{
      counter ++;
      maxsbardiff = 0.0;
      randperm(latentDim, order);

      int j = 0;
      for(int index = 0; index < latentDim; index++){        
        j = order[index];
        //for(int i2 = 0; i2 < latentDim; i2++)
        //  mexPrintf("Order %d: %d\n", i2, order[i2]);
        if(j>=latentDim)
          mexErrMsgTxt("Error j is greater than latent dimension");
        for(int i = 0; i < dataDim; i++){
          vTs[i] = 0.0;
          
          for(int notj = 0; notj < latentDim; notj++){
            if(notj != j){
              vTs[i] +=  V[notj + latentDim*i]*sbar[n + notj*ndata];
            }
          }
        }
        
        double oldsbar = sbar[n + ndata*j];
        double sigma2 = tau[n + ndata*j];
        double tempDouble = 0.0;
        
        for(i = 0; i < dataDim; i++){
          sigma2 += Bl2V2[i + dataDim*j];
          tempDouble += Bl2[i]*(XULinv[n + ndata*i] 
                                - vTs[i])*V[j + latentDim*i];
        }
        sigma2 = 1/sigma2;
        sbar[n + ndata*j] = sigma2*tempDouble;
        sbardiff = fabs(float(sbar[n + ndata*j] - oldsbar));
        if (sbardiff > maxsbardiff)
          maxsbardiff = sbardiff;
        Sigma_s[j + latentDim*j + latentDim*latentDim*n] = sigma2;
      }
    }
    while (maxsbardiff>tolerance && counter <= maxCount);
    if(counter > maxCount)
      mexPrintf("Max Count exceeded in mfupdatelatent, difference %f\n", (float)maxsbardiff);
  
  }      
  
  mxFree(Bl2);
  mxFree(V2);
  mxFree(vTs);
  mxFree(Bl2V2);
  mxFree(order);
}

void mexFunction(int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  double *XULinv = (double *)mxGetData(prhs[0]);
  const int *dataSetDims;
  dataSetDims = mxGetDimensions(prhs[0]);
  int ndata = dataSetDims[0];
  int dataDim = dataSetDims[1];
  
  double *V = (double *)mxGetData(prhs[1]);
  const int *VDims;
  VDims = mxGetDimensions(prhs[1]);
  int latentDim = VDims[0];
  if(VDims[1] != dataDim)
    mexErrMsgTxt("V matrix data dimensions do not match XULinv");

  double *tau = (double *)mxGetData(prhs[2]);
  VDims = mxGetDimensions(prhs[2]);
  if(VDims[0] != ndata)
    mexErrMsgTxt("Tau matrix does not match data size");
  if(VDims[1] != latentDim)
    mexErrMsgTxt("Tau matrix does not match latent dimension size");
    

  double beta = mxGetScalar(prhs[3]);
  VDims = mxGetDimensions(prhs[3]);
  if(VDims[0] != 1 || VDims[1] != 1)
    mexErrMsgTxt("Beta should be a scalar");
   
  double *lambda = (double *)mxGetData(prhs[4]);
  VDims = mxGetDimensions(prhs[4]);
  if(!(VDims[0] == 1 && VDims[1] == dataDim) && !(VDims[1] == 1 && VDims[0] == dataDim))
    mexErrMsgTxt("Lambda should be a vector of dataDim length");

  
  int sbardims[2];
  sbardims[0] = ndata;
  sbardims[1] = dataDim;
  plhs[0] = mxDuplicateArray(prhs[5]);
  double *sbar = (double *)mxGetData(plhs[0]);
  int maxcount = 100;
  if (nrhs > 6)
    maxcount = (int)mxGetScalar(prhs[6]);
  double tolerance = 1e-5;
  if (nrhs > 7)
    tolerance = mxGetScalar(prhs[7]);

  int dims[3];
  dims[0] = latentDim;
  dims[1] = latentDim;
  dims[2] = ndata;
  plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
  double *Sigma_s = (double *)mxGetData(plhs[1]);

  plhs[2] = mxCreateScalarDouble(0.0);
  double *sbardiff = (double *)mxGetData(plhs[2]);
  mfupdatelatent(XULinv, V, tau,
                 beta, lambda, 
                 sbar, Sigma_s, ndata, latentDim, dataDim, maxcount, tolerance);
}





