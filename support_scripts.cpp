#include "stdafx.h"
#include "fftw3.h"
#include <stdlib.h>             // for malloc and calloc
#include "math.h"
using namespace std;

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

void PCChirp(double* in, double dur, double fmin, double fmax, double fsamp, 
					double Tukpar, int Nt, int Nscan, double* out) {
	int Nchirp = int(ceil(dur * fsamp));
	int tCcnt,scnt,tcnt;
	double tchirp,Tukwin;
	double* chirp;

	// Define chirp:
	chirp = (double*)calloc(Nchirp, sizeof(double));
	for (tCcnt = 0; tCcnt < Nchirp; tCcnt++) {
		tchirp = tCcnt / fsamp;
		chirp[tCcnt] = sin(2.0 * M_PI * ((fmax - fmin) / 2.0 / dur * tchirp + fmin) * tchirp);
	}

	// Apply Tukey window for mismatched filtering:
	if (Tukpar > 0.0 && Tukpar <= 1.0) {
		for (tCcnt = 0; tCcnt < ceil(Nchirp * Tukpar / 2.0); tCcnt++) {
			Tukwin = 0.5 * (1.0 + cos(2.0 * M_PI * (tCcnt / double(Nchirp - 1) / Tukpar - 0.5)));
			chirp[tCcnt] *= Tukwin;
			chirp[Nchirp - 1 - tCcnt] *= Tukwin;
		}
	}

	// Cross-correlation for pulse compression:
	/*for (scnt = 0; scnt < Nscan; scnt++) {
		for (tcnt = 0; tcnt < Nt; tcnt++) {
			for (tCcnt = 0; tCcnt < min(Nchirp, Nt - tcnt); tCcnt++) {
				out[tcnt + scnt*Nt] += in[tcnt+tCcnt + scnt*Nt] * chirp[tCcnt];
			}
		}
	}*/

	// Cross-correlation for pulse compression, frequency domain:
	// 1. Create variables, allocate memory, create FFTW plans:
	fftw_complex *input, *INPUT, *OUTPUT, *CHIRP;
	fftw_plan	  p_in, p_chirp, p_out;
	CHIRP = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (Nt+Nchirp-1));
	input = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (Nt + Nchirp - 1));
	INPUT = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (Nt + Nchirp - 1));
	OUTPUT = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * (Nt + Nchirp - 1));
	for (tcnt = 0; tcnt<Nt+Nchirp-1; tcnt++){
		CHIRP[tcnt][0] = 0.0; CHIRP[tcnt][1] = 0.0;
		input[tcnt][0] = 0.0; input[tcnt][1] = 0.0;
		INPUT[tcnt][0] = 0.0; INPUT[tcnt][1] = 0.0;
		OUTPUT[tcnt][0] = 0.0; OUTPUT[tcnt][1] = 0.0;
	}
	p_chirp = fftw_plan_dft_1d(Nt + Nchirp - 1, CHIRP, CHIRP, FFTW_FORWARD, FFTW_ESTIMATE);
	p_in	= fftw_plan_dft_1d(Nt + Nchirp - 1, input, INPUT, FFTW_FORWARD, FFTW_ESTIMATE);
	p_out   = fftw_plan_dft_1d(Nt + Nchirp - 1, OUTPUT, OUTPUT, FFTW_FORWARD, FFTW_ESTIMATE);

	// 2. compute fft(chirp(end:-1:1)):
	for (tCcnt = 0; tCcnt < Nchirp; tCcnt++) {
		CHIRP[tCcnt][0] = chirp[Nchirp - tCcnt - 1];
	}
	fftw_execute(p_chirp);

	// 3. Perform cross-correlation in freq dom:
	for (scnt = 0; scnt < Nscan; scnt++) {	// loop over all A-scans
		for (tcnt = 0; tcnt < Nt; tcnt++) {	// extract A-scan from "in"
			input[tcnt][0] = in[tcnt + scnt * Nt];
		}
		fftw_execute(p_in);					// compute fft(in)
		for (tcnt = 0; tcnt < Nt + Nchirp - 1; tcnt++) { // compute IN*CHIRP
			complex_mult(INPUT[tcnt], CHIRP[tcnt], OUTPUT[tcnt]);
		}
		fftw_execute(p_out);				// compute ifft(IN*CHIRP)
		for (tcnt = 0; tcnt < Nt; tcnt++) { // extract correct part of cross-correlation
			//out[tcnt + scnt * Nt] = OUTPUT[Nchirp-1+tcnt][0] / (1.0*(Nt+Nchirp-1));
			out[tcnt + scnt * Nt] = OUTPUT[Nt - tcnt][0] / (1.0 * (Nt + Nchirp - 1));
		}
	}

	return;
}