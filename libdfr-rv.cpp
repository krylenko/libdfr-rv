/********************* random variable generator *************************/
/*                          daniel ford                                  */
/*                            may 2011                                   */
/*************************************************************************/

/*
G++/MINGW
g++ rvTest.cpp libdfr-rv.cpp ../matrix/libdfr-matrix.cpp -ffast-math -I../matrix -o rvTest.exe

GNUPLOT
load "histplot.gpi"

ACTION
-M-H algorithm
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <time.h>

#include "libdfr-rv.h"
#include "../libdfr-matrix/libdfr-matrix.h"

#define FASTPOW 4294967296.		// 2^32 for fast division
#define PI      3.14151926535897932

// global variables to track seed status
int rvSeedVal=0;
int seedUpdate=0;
int seedFlag=0;

// ctor for Gaussian class
Gaussian::Gaussian(const Matrix& mean, const Matrix& covar)
: mu(mean), sigma(covar), dim_(mu.rows())
{
  
  if(sigma.rows() != dim_){
    cout << "dimension mismatch: Gaussian::Gaussian" << endl;
    exit(-1);
  }

}
  
// return one sample from the distribution (matrix version)
// expects a column vector X
double Gaussian::pX(const Matrix& X){

  double indim = X.rows();

  if( indim != dim_ ){
    cout << "dimension mismatch: Gaussian::pX" << endl;
    exit(-1);
  }

  double p = 0;
  double norm = 0;
  norm = pow((2*PI),(-dim_/2.))*pow(sigma.det(),(-1./2));
  p = exp((-1./2)*( (X-mu).T()*sigma.inv()*(X-mu) ).sum());
  return norm*p;

}

// seed random number generators
void rvSeed(int seedIn)
{
	srand(seedIn);
	rvSeedVal = seedIn;
	seedUpdate = 0;
	seedFlag = 0;
}

// desired distribution for rejection and SIR sampling, change per application
// currently a mixture of Gaussians
double p_x(const double& x)
{
	double p;
	//p = 0.3*exp(-pow((x-0.3),2)) + 0.7* exp(-pow((x-2.),2)/0.3);	// test mixture from book
    //p = 0.75*exp(-pow((x-2),2)/.22) + 1.5* exp(-pow((x-.1),2)/9);	// second test mixture
    p = (1/sqrt(2*3.14159))*exp(-pow(x,2.)/2);											// simple 0,1 Gaussian
	//p = 1./15.;													// uniform for 15 bins
	return p;
}

// generate one sample from quick & dirty PRNG from Numerical Recipes
double rvDirtyUniform(const double& begin, const double& end)
{
	
	double length, variate, result;

	length = abs(end-begin);
	
	if(seedFlag == 0)
	{
		seedUpdate = 1664525L*rvSeedVal + 1013904223L;
		seedFlag++;
	}
	else
		seedUpdate = 1664525L*seedUpdate + 1013904223L;

	variate = seedUpdate/FASTPOW+0.5;	// convert int to floating-point and scale to [0,1)

	result = (variate*length)+begin;
	
	return result;
}

// generate one sample from a uniform distribution
double rvStdUniform(const double& begin, const double& end)
{
	double length, variate, result;
	
	length = abs(end-begin);
	variate = (double)rand()/RAND_MAX;
	result = (variate*length)+begin;
	
	return result;
}

// generate one sample from a Gaussian with given parameters
double rvGaussian(const double& mean, const double& variance)
{
	double x1,x2,y1,y2,w;
	
	// grab random uniform numbers 
	x1 = rvStdUniform(-1,1);
	x2 = rvStdUniform(-1,1);
	w = x1*x1+x2*x2;
	
	// if they're outside the unit circle, throw them out
	while(w >= 1)
	{
		x1 = rvStdUniform(-1,1);
		x2 = rvStdUniform(-1,1);
		w = x1*x1+x2*x2;
	}
	
	// if they're OK, generate a Gaussian of zero mean, unit variance
	y1 = x1*sqrt(-2.0*log(w)/w);
	
	// scale Gaussian to desired shape
	y1 = mean + sqrt(variance)*y1;
	
	return y1;
}

// generate samples from a distribution using inverse CDF sampling
// pdf input is expected to be an Nx2 matrix with col 1 = values and col 2 = probabilities of those values
// step size scales pdf to equal 1, set equal 
Matrix invCDFsample(const int& samples, const double& step,
					const Matrix& pdf)
{

	Matrix samp(samples,1);
	double tempSamp = 0.;
	int j = 0;

	// integrate to generate CDF
	Matrix CDF(pdf.rows(),1);
	CDF[0][0] = pdf[0][1]*step;
	for(int i=0;i<CDF.rows()-1;i++)
	{
    //cout << "CDF[i] = " << CDF[i][0] << " pdf[i] = " << pdf[i][1] << endl;
    CDF[i+1][0] = CDF[i][0] + pdf[i+1][1];	
    //cout << "CDF[i+1] = " << CDF[i+1][0] << endl;
  }

	// multiply by step size to scale integral correctly
	CDF = CDF*step;
	//CDF.print();
  
	// choose sample from uniform generator according to CDF 
	for(int i=0;i<samples;i++)
	{
		j = 0;
		tempSamp = rvStdUniform();
		
		while( CDF[j][0] < tempSamp )
		{
			j++;
		}
		
		if( j == 0 )
			samp[i][0] = pdf[j][0];
		else
			samp[i][0] = (pdf[j-1][0]+pdf[j][0])/2;	
	}
	
	return samp;
	
}

// generate one sample from a target distribution p(x) using rejection sampling
double rvRejectionSample(const double& M, const double& PDFscale)
{
	double x, u, target, rejectRatio;
	
	// note: PDFscale scales the related density g(x) to ensure
	// it's larger than the distribution we want to sample;
	// otherwise we'll just end up sampling g(x)
	
	// note: machine learning book says to use u < p(x)/Mg(x)
	// as the acceptance criteria, but actually uses u < p(x)
	// also, book says to sample from g(x) Gaussian (for Gaussian desired density)
	// but actually uses uniform density
	
	// the following code is "correct" but doesn't seem
	// to match desired density as well as book's method
	
	// to use book's method, comment out all *^ lines
	// and comment in all ** lines
	
	// generate initial sample
	x = rvGaussian(0,1)*PDFscale; 		// related density, here Gaussian
										// can also just choose a constant value
										// larger than the max of the target distribution
										// *^
	
	//x = rvStdUniform(0,1)*PDFscale;	// **
								
	u = rvStdUniform(0,1)*M;			// enveloping uniform distribution
	//rejectRatio = p_x(x)/u;				// threshold to accept/reject samples *^
	
	// generate sample x until it falls within the target distribution
	//while( u >= rejectRatio)			// *^
	while( u >= p_x(x) )				// **
	{
		x = rvGaussian(0,1)*PDFscale;		// *^
		//x = rvStdUniform(0,1)*PDFscale;	// **
		u = rvStdUniform(0,1)*M;
		//rejectRatio = p_x(x)/u;				// *^	
	}
		
	return x;
}

// generate batch of N samples from a target distribution p(x)
// using sampling-importance-rejection sampling w/ one resampling stage
Matrix rvSIRBatch(const int& N, const double& PDFscale)
{

	int j=0, k=0, m=0;
	
	Matrix Qsamples(1,N);
	Matrix w(1,N);
	Matrix wCDF(1,N);
	Matrix rand(1,N);
	Matrix x(1,N);
	
	// sample from q(x) and calculate weights (and create random array for later)
	for(int i=0;i<N;i++)
	{
		Qsamples[0][i] = rvStdUniform(-1,1)*PDFscale;		// [-1,1] to handle all PDFs
		
		//w[0][i] = p_x(Qsamples[0][i])/Qsamples[0][i];
		// line below works, line above doesn't;
		w[0][i] = p_x(Qsamples[0][i])/PDFscale;			
		
		rand[0][i] = rvStdUniform(0,1);
	}
		
	// normalize weights
	w = w/w.sum();
	
	// create CDF of weights
	wCDF[0][0] = w[0][0];
	for(int j=1;j<N;j++)
		wCDF[0][j] = wCDF[0][j-1] + w[0][j];
	
	//make copies as many times as numbers occur in CDF
	for(k = 0; k < N; ++k)
	{
		m = 0;
		while(wCDF[0][m] < rand[0][k])
		{
			++m;
		}
		x[0][k] = Qsamples[0][m];
	}
	
	// reset weights to uniform values
	double wUniform = 1./N;
	w.fill(wUniform);
	
	return x;

}

// sort input matrix into bins for histogram plotting
Matrix rvBin(	const Matrix& dist, const double& begin,
				const double& end, const int& binNum)
{

	double resolution = abs(end-begin)/binNum;
	double width = resolution/2;
	double offset = width+begin;
	Matrix hist(binNum,2);			// first column is bin centers, second is counts 

	// initialize with equally spaced bin centers
	for(int i=0;i<binNum;i++)
	{
		if(i==0)
			hist[i][0] = offset;
		else 	
		{
			hist[i][0] = hist[i-1][0] + resolution;
		}
	}
	// sort values into bins
	for(int j=0;j<dist.cols();j++)
	{
		for(int i=0;i<binNum;i++)
		{
			if( 	( dist[0][j] > (hist[i][0]-width) ) &&
					( dist[0][j] <= (hist[i][0]+width) )  	)
			{
				hist[i][1] += 1;
			}
		}
	}

	
	// normalize counts to probabilities based on the number of samples
	// actually present in the specified x range
	double sampSum = 0;
	for(int p=0;p<hist.rows();p++)
		sampSum += hist[p][1];
	for(int m=0;m<hist.rows();m++)
	{
		hist[m][1] /= sampSum;
		hist[m][1] /= resolution;		// scale probabilities
										// to reflect the fact
										// that the histogram
										// is an integral
										// (i.e., multiply by dx)
	}
	
	return hist;
}

// sort input matrix into bins for histogram plotting (vector input)
Matrix rvBin(	const vector<double> dist, const double& begin,
				const double& end, const int& binNum)
{

	double resolution = abs(end-begin)/binNum;
	double width = resolution/2;
	double offset = width+begin;
	Matrix hist(binNum,2);			// first column is bin centers, second is counts 

	// initialize with equally spaced bin centers
	for(int i=0;i<binNum;i++)
	{
		if(i==0)
			hist[i][0] = offset;
		else 	
		{
			hist[i][0] = hist[i-1][0] + resolution;
		}
	}
	// sort values into bins
	for(int j=0;j<dist.size();j++)
	{
		for(int i=0;i<binNum;i++)
		{
			if( 	( dist[j] > (hist[i][0]-width) ) &&
					( dist[j] <= (hist[i][0]+width) )  	)
			{
				hist[i][1] += 1;
			}
		}
	}
	
	// normalize counts to probabilities based on the number of samples
	// actually present in the specified x range
	double sampSum = 0;
	for(int p=0;p<hist.rows();p++)
		sampSum += hist[p][1];
	for(int m=0;m<hist.rows();m++)
		hist[m][1] /= sampSum;
		
	return hist;
}
