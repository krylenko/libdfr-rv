/********************* random variable generator *************************/
/*                          daniel ford                                  */
/*                            may 2011                                   */
/*************************************************************************/

#ifndef LIBDFR_RV_H
#define LIBDFR_RV_H

#include <iostream>
#include <vector>
#include <time.h>

//#include "libdfr-matrix.h"

#include "C:/Users/ghost/Desktop/robotics/code/dfrobot/matrix/libdfr-matrix.h"

// note: could make this into a static class rather than just a group of methods

using namespace std;

// default seed for uniform random generators
/*
struct timeval tv;
gettimeofday(&tv,NULL);
rvSeed(tv.tv_usec);
*/
void rvSeed(int seedIn=time(NULL));

// quick & dirty random number generator from Numerical Recipes
double rvDirtyUniform(const double& begin=0., const double& end=1.);

// return one sample from a uniform distribution on [0,1)
double rvStdUniform(const double& begin=0., const double& end=1.);

// return one sample from a Gaussian distribution with given parameters
double rvGaussian(const double& mean=0., const double& variance=1.);

// return $SAMPLES samples from a given distribution using inverse CDF sampling
// input pdf must be column vector
Matrix invCDFsample(const int& samples, const double& step,
					const Matrix& pdf);

// sample from desired distribution p(x) using rejection sampling
// M is the height of the enveloping distribution, PDFscale its range
double rvRejectionSample(const double& M=1.0, const double& PDFscale=4.);

// generate N samples from desired distribution p(x) using SIR sampling
Matrix rvSIRBatch(const int& N, const double& PDFscale=4.);

// definition of p(x) for rejection sampling
double p_x(const double& x);

// sort input vectors into given number of bins
// output matrix of bin centers and counts
Matrix rvBin(	const Matrix& dist, const double& begin,
				const double& end, const int& binNum);

// vector input version of histogram function
Matrix rvBin(	const vector<double>& dist, const double& begin,
				const double& end, const int& binNum);

// multivariate Gaussian class        
class Gaussian{
  
  public:

    // std constructor
    Gaussian(const Matrix& mean, const Matrix& covar);

    // public so we can inspect/change them w/o defining new methods
    Matrix mu;     
    Matrix sigma;    
    
    // calculate the probability of a point
    double pX(const Matrix& X);
    //double pX(const vector<double>& X);
    
    // dim getter
    int dim() const {return dim_;}
    
  private:

    int dim_;
    
};
				
#endif