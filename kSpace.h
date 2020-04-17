#include "fftw3.h"

void kSpace_line_recon_alg(double* rf_data,
	double* source_locations,
	double c, double fsamp,
	int Nsrc, int Nt, int N,
	double* xaxis, double* zaxis,
	double* image);

void compute_factor(fftw_complex constant, double dist, fftw_complex output);