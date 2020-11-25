#include "libdfr-rv.h"
#include "../libdfr-matrix/libdfr-matrix.h"

int main()
{

	int samples = 10000;
	int bins = 50;
	
	// seed all random generators
	rvSeed();
	
	// generate samples
	/*
	Matrix values(1,samples);
	for(int i=0;i<samples;i++)
		//values[0][i] = rvStdUniform(-1,1);
		//values[0][i] = rvDirtyUniform(-1,1);
		//values[0][i] = rvGaussian();
		values[0][i] = rvRejectionSample(4, 16);
	*/
    Matrix values = rvSIRBatch(samples, 3);

	// generate Pr histogram and print to file
	Matrix hist = rvBin(values, -4, 4, bins);

	// plot desired density for comparison
	Matrix targets(hist);
	
	// calc desired function
	double tSum=0;
	for(int j=0;j<targets.rows();j++)
	{
		targets[j][1] = p_x(targets[j][0]);
		tSum += targets[j][1];
	}
	// normalize to make probabilities
	for(int j=0;j<targets.rows();j++)
		targets[j][1] /= tSum;
	
	// print sampled and desired densities to file	
    hist.printFile("dist.dat");
    targets.printFile("real.dat");
	//hist.print();

}
