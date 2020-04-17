#include "stdafx.h"
#include "fftw3.h"
#include "math.h"

double squa(double base) { return base*base; }

void Compute_envelope_FFTW(double *image,int Nx,int Ny,int Nz)
{
	fftw_complex	*in, *out, *hilb;
	fftw_plan		p, pinv;
	int cnt,zcnt;
	int NxNy = Nx*Ny;
	double Nzinv = 1.0/Nz;

	in			= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nz);
	out			= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nz);
	hilb		= (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nz);
	for (zcnt=0; zcnt<Nz; zcnt++){ in[zcnt][0] = 0.0;	in[zcnt][1] = 0.0; }
	for (zcnt=0; zcnt<Nz; zcnt++){ hilb[zcnt][0] = 0.0;	hilb[zcnt][1] = 0.0; }

	p		= fftw_plan_dft_1d(Nz, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	pinv	= fftw_plan_dft_1d(Nz, hilb, out, FFTW_BACKWARD, FFTW_ESTIMATE);

	for (cnt=0; cnt<NxNy; cnt++){
		if (Nz>1){
			// Load each z-trace in image:
			for (zcnt=0; zcnt<Nz; zcnt++){ in[zcnt][0] = image[cnt + zcnt*NxNy]; }

			// Envelope computed as in http://cdn.intechopen.com/pdfs-wm/36434.pdf, p.295:
			// compute FFT(image selection):
			fftw_execute(p);

			for (zcnt=1; zcnt<=Nz/2; zcnt++){
				hilb[zcnt][0] = out[zcnt][0];
				hilb[zcnt][1] = out[zcnt][1];
			}
			fftw_execute(pinv);
			for (zcnt=0; zcnt<Nz; zcnt++){
				image[cnt + zcnt*NxNy] = 2.0*Nzinv * sqrt((squa(out[zcnt][0]) + squa(out[zcnt][1])));
			}
		}
		else {
			image[cnt] = fabs(image[cnt]);
		}
	}
	
	fftw_destroy_plan(p);	fftw_destroy_plan(pinv);
	fftw_free(in); fftw_free(out); fftw_free(hilb);
	return;
}

void complex_mult(fftw_complex in1, fftw_complex in2, fftw_complex output) {
	// (a+b*i) * (c+d*i) = ac-bd+(ad+bc)*i
	// NOTE: THIS FUNCTION DOES NOT WORK IN-PLACE!!!!!
	output[0] = in1[0] * in2[0] - in1[1] * in2[1];
	output[1] = in1[0] * in2[1] + in1[1] * in2[0];
	return;
}

void complex_add(fftw_complex in1, fftw_complex in2, fftw_complex output) {
	// (a+b*i) + (c+d*i) = a+c + (b+d)*i
	// NOTE: THIS FUNCTION DOES NOT WORK IN-PLACE!!!!!
	output[0] = in1[0] + in2[0];
	output[1] = in1[1] + in2[1];
	return;
}