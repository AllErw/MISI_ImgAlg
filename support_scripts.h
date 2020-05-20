#include "fftw3.h"

double squa(double base);
void Compute_envelope_FFTW(double *image,int Nx,int Ny,int Nz);
void complex_mult(fftw_complex in1, fftw_complex in2, fftw_complex output);
void complex_add(fftw_complex in1, fftw_complex in2, fftw_complex output);
void PCChirp(double* in, double dur, double fmin, double fmax, double fsamp,
			 double Tukpar, int Nt, int Nscan, double* out);