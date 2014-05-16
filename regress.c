#include <gsl/gsl_multifit.h>

void regress(size_t nInd, size_t nCov, double *x, double *y, double *rOut){
  int i, j;
  gsl_matrix *xMat, *cov;
  gsl_vector *yVec, *c, *r;
  double chisq;

  xMat = gsl_matrix_alloc(nInd, nCov);
  yVec = gsl_vector_alloc(nInd);
  r = gsl_vector_alloc(nInd);
  c = gsl_vector_alloc(nCov);
  cov = gsl_matrix_alloc(nCov, nCov);
  for(i = 0; i < nInd; i++)
    {
      gsl_vector_set(yVec, i, *(y+i));
      for(j = 0; j < nCov; j++)
	gsl_matrix_set(xMat, i, j, *(x + i * nCov + j));
    }
 
  gsl_multifit_linear_workspace *work =  gsl_multifit_linear_alloc(nInd, nCov);
  gsl_multifit_linear(xMat, yVec, c, cov,
		      &chisq, work);
  gsl_multifit_linear_residuals(xMat, yVec, c, r);
  gsl_multifit_linear_free(work);
 

  for(i = 0; i < nInd; i++)
    rOut[i] = gsl_vector_get(r, i); 
}
