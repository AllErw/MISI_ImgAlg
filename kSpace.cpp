#include "stdafx.h"
#include "fftw3.h"
#include "math.h"
#include <stdlib.h>             // for malloc and calloc
#include <cmath>				// for abs in double precision
#include "support_scripts.h"
#include "kSpace.h"
#include <complex>
using namespace std;

//void kSpace_line_recon_alg(double* rf_data,
//	double* source_locations,
//	double c, double fsamp,
//	int Nsrc, int Nt, int N,
//	double* xaxis, double* zaxis,
//	double* image)
//	/*	k-space reconstruction of initial pressure distribution, computed using
//		the method of [Jaeger, 2007], "Fourier reconstruction in optoacoustic
//		imaging using truncated regularized inverse k - space interpolation". This
//		method is a generalisation of the method implemented in k-Wave's
//		kspacelinerecon function, with the dstinct difference that the interpolation
//		scheme resembles "sinc" instead of linear / cubic / etc.
//		Note: setting N = 0 results in "nearest neighbour" interpolation.	*/
//{
//	double dx, dz, constant, kt, mindiff, diff, dist, J;
//	int kxcnt, kzcnt, wcnt, Ncnt, nearest;
//	double* kxvect, * kzvect, * wvect;
//	fftw_complex* RF, * P0, RFelt, factor, factRFelt, complconst;
//	fftw_plan plan_fwd, plan_inv;
//
//	// Compute x- and z-axes:
//	dx = source_locations[1] - source_locations[0];
//	dz = c / fsamp;
//	for (kxcnt = 0; kxcnt < Nsrc; kxcnt++) { xaxis[kxcnt] = kxcnt * dx; }
//	for (kzcnt = 0; kzcnt < Nt; kzcnt++) { zaxis[kzcnt] = kzcnt * dz; }
//
//
//	// Compute wave vector components:
//	kxvect = (double*)malloc(Nsrc * sizeof(double));
//	constant = 2.0 * M_PI / dx / (double)Nsrc;
//	for (kxcnt = 0; kxcnt < ceil(Nsrc / 2.0); kxcnt++) { kxvect[kxcnt] = (double)kxcnt * constant; }
//	for (kxcnt = (int)ceil(Nsrc / 2.0); kxcnt < Nsrc; kxcnt++) { kxvect[kxcnt] = (double)(kxcnt - Nsrc) * constant; }
//
//	kzvect = (double*)malloc(Nt * sizeof(double));
//	wvect = (double*)malloc(Nt * sizeof(double));
//	constant = 2.0 * M_PI * fsamp / c / (double)Nt;
//	for (kzcnt = 0; kzcnt < ceil(Nt / 2.0); kzcnt++) { kzvect[kzcnt] = (double)kzcnt * constant; }
//	for (kzcnt = (int)ceil(Nt / 2.0); kzcnt < Nt; kzcnt++) { kzvect[kzcnt] = (double)(kzcnt - Nt) * constant; }
//	for (kzcnt = 0; kzcnt < Nt; kzcnt++) { wvect[kzcnt] = c * kzvect[kzcnt]; }
//
//
//	// Preallocate RF and P0 and fill by copying rf_data:
//	RF = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nt * Nsrc);
//	P0 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nt * Nsrc);
//	for (kxcnt = 0; kxcnt < Nt * Nsrc; kxcnt++) {
//		RF[kxcnt][0] = rf_data[kxcnt];
//		RF[kxcnt][1] = 0.0;
//		P0[kxcnt][0] = 0.0;
//		P0[kxcnt][1] = 0.0;
//	}
//
//	// Create plans for in-place forward and inverse DFTs:
//	plan_fwd = fftw_plan_dft_2d(Nsrc, Nt, RF, RF, FFTW_FORWARD, FFTW_ESTIMATE);
//	plan_inv = fftw_plan_dft_2d(Nt, Nsrc, P0, P0, FFTW_BACKWARD, FFTW_ESTIMATE);
//	fftw_execute(plan_fwd);
//
//	complconst[0] = 0.0; complconst[1] = (1.0 * Nt) / fsamp;
//	for (kxcnt = 0; kxcnt < Nsrc; kxcnt++) {
//		for (kzcnt = 0; kzcnt < Nt; kzcnt++) {
//
//			// Compute time-offset kz value using dispersion relation (hence "c * "):
//			if (kzcnt != 0) { kt = copysign(c * sqrt(squa(kxvect[kxcnt]) + squa(kzvect[kzcnt])), kzvect[kzcnt]); }
//			else { kt = 0.0; }
//
//			// find closest match between kt and kzvect:
//			nearest = 0;
//			mindiff = 2.0 * M_PI * fsamp;
//			for (wcnt = 0; wcnt < Nt; wcnt++) {
//				diff = fabs(kt - wvect[wcnt]);
//				if ((diff) < (mindiff)) {
//					mindiff = diff;
//					nearest = wcnt;
//				}
//			}
//
//			// Perform interpolation using -N : N elements centered around nearest:
//			for (Ncnt = -N; Ncnt <= N; Ncnt++) {
//				if (nearest + Ncnt >= 0 && nearest + Ncnt < Nt) {
//					dist = kt - wvect[nearest + Ncnt];
//					if (dist != 0.0) {
//						compute_factor(complconst, dist, factor);
//						RFelt[0] = RF[nearest + Ncnt + kxcnt * Nt][0]; RFelt[1] = RF[nearest + Ncnt + kxcnt * Nt][1];
//						complex_mult(factor, RFelt, factRFelt);
//						P0[kxcnt + kzcnt * Nsrc][0] += factRFelt[0];
//						P0[kxcnt + kzcnt * Nsrc][1] += factRFelt[1];
//					}
//				}
//			}
//
//			// Compute and apply Jacobian:
//			if (kt == 0.0) { J = 1.0; }
//			else { J = kzvect[kzcnt] / kt; }
//			P0[kxcnt + kzcnt * Nsrc][0] *= J;	P0[kxcnt + kzcnt * Nsrc][1] *= J;
//
//		}
//	}
//	//P0[10000][0] += 0.04; // adds a sinusoidal oscillation along lateral distance
//	//P0[5][0] += 0.02; // adds a sinusoidal oscillation along time
//	//P0[10005][0] += 0.02; // adds a sinusoidal oscillation along time and space
//
//	fftw_execute(plan_inv);
//	for (kxcnt = 0; kxcnt < Nt * Nsrc; kxcnt++) { image[kxcnt] = P0[kxcnt][0]; }
//
//	fftw_destroy_plan(plan_fwd); fftw_destroy_plan(plan_inv);
//	fftw_free(RF); fftw_free(P0);
//}

void kSpace_line_recon_alg(double* rf_data,
	double* source_locations,
	double c, double fsamp,
	int Nsrc, int Nt, int N,
	double* xaxis, double* zaxis,
	double* image)
	/*	k-space reconstruction of initial pressure distribution, computed using
		the method of [Jaeger, 2007], "Fourier reconstruction in optoacoustic
		imaging using truncated regularized inverse k - space interpolation". This
		method is a generalisation of the method implemented in k-Wave's
		kspacelinerecon function, with the dstinct difference that the interpolation
		scheme resembles "sinc" instead of linear / cubic / etc.
		Note: setting N = 0 results in "nearest neighbour" interpolation.	*/
{
	double dx, dz, constant, kt, mindiff, diff, dist, J;
	int kxcnt, kzcnt, wcnt, Ncnt, nearest;
	double* kxvect, * kzvect, * wvect;
	fftw_complex* RF, * P0, RFelt, factor, factRFelt, complconst;
	fftw_plan plan_fwd, plan_inv;

	// Compute x- and z-axes:
	dx = source_locations[1] - source_locations[0];
	dz = c / fsamp;
	for (kxcnt = 0; kxcnt < Nsrc; kxcnt++) { xaxis[kxcnt] = kxcnt * dx; }
	for (kzcnt = 0; kzcnt < Nt; kzcnt++) { zaxis[kzcnt] = kzcnt * dz; }


	// Compute wave vector components:
	kxvect = (double*)malloc(Nsrc * sizeof(double));
	constant = 2.0 * M_PI / dx / (double)Nsrc;
	for (kxcnt = 0; kxcnt < ceil(Nsrc / 2.0); kxcnt++) { kxvect[kxcnt] = (double)kxcnt * constant; }
	for (kxcnt = (int)ceil(Nsrc / 2.0); kxcnt < Nsrc; kxcnt++) { kxvect[kxcnt] = (double)(kxcnt - Nsrc) * constant; }

	kzvect = (double*)malloc(Nt * sizeof(double));
	wvect = (double*)malloc(Nt * sizeof(double));
	constant = 2.0 * M_PI * fsamp / c / (double)Nt;
	for (kzcnt = 0; kzcnt < ceil(Nt / 2.0); kzcnt++) { kzvect[kzcnt] = (double)kzcnt * constant; }
	for (kzcnt = (int)ceil(Nt / 2.0); kzcnt < Nt; kzcnt++) { kzvect[kzcnt] = (double)(kzcnt - Nt) * constant; }
	for (kzcnt = 0; kzcnt < Nt; kzcnt++) { wvect[kzcnt] = c * kzvect[kzcnt]; }


	// Preallocate RF and P0 and fill by copying rf_data:
	RF = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nt * Nsrc);
	P0 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nt * Nsrc);
	for (kxcnt = 0; kxcnt < Nt * Nsrc; kxcnt++) {
		RF[kxcnt][0] = rf_data[kxcnt];
		RF[kxcnt][1] = 0.0;
		P0[kxcnt][0] = 0.0;
		P0[kxcnt][1] = 0.0;
	}

	// Create plans for in-place forward and inverse DFTs:
	plan_fwd = fftw_plan_dft_2d(Nsrc, Nt, RF, RF, FFTW_FORWARD, FFTW_ESTIMATE);
	plan_inv = fftw_plan_dft_2d(Nt, Nsrc, P0, P0, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(plan_fwd);

	complconst[0] = 0.0; complconst[1] = (1.0 * Nt) / fsamp;
	for (kxcnt = 0; kxcnt < Nsrc; kxcnt++) {
		for (kzcnt = 0; kzcnt < ceil(Nt / 2.0)+1; kzcnt++) {

			// Compute time-offset kz value using dispersion relation (hence "c * "):
			if (kzcnt != 0) { kt = copysign(c * sqrt(squa(kxvect[kxcnt]) + squa(kzvect[kzcnt])), kzvect[kzcnt]); }
			else { kt = 0.0; }

			// find closest match between kt and kzvect:
			nearest = 0;
			mindiff = 2.0 * M_PI * fsamp;
			for (wcnt = 0; wcnt < Nt; wcnt++) {
				diff = fabs(kt - wvect[wcnt]);
				if ((diff) < (mindiff)) {
					mindiff = diff;
					nearest = wcnt;
				}
			}

			// Perform interpolation using -N : N elements centered around nearest:
			for (Ncnt = -N; Ncnt <= N; Ncnt++) {
				if (nearest + Ncnt >= 0 && nearest + Ncnt < Nt) {
					dist = kt - wvect[nearest + Ncnt];
					if (dist != 0.0) {
						compute_factor(complconst, dist, factor);
						RFelt[0] = RF[nearest + Ncnt + kxcnt * Nt][0]; RFelt[1] = RF[nearest + Ncnt + kxcnt * Nt][1];
						complex_mult(factor, RFelt, factRFelt);
						P0[kxcnt + kzcnt * Nsrc][0] += factRFelt[0];
						P0[kxcnt + kzcnt * Nsrc][1] += factRFelt[1];
					}
				}
			}

			// Compute and apply Jacobian:
			if (kt == 0.0) { J = 1.0; }
			else { J = kzvect[kzcnt] / kt; }
			P0[kxcnt + kzcnt * Nsrc][0] *= J;	P0[kxcnt + kzcnt * Nsrc][1] *= J;

		}
	}
	//P0[10000][0] += 0.04; // adds a sinusoidal oscillation along lateral distance
	//P0[5][0] += 0.02; // adds a sinusoidal oscillation along time
	//P0[10005][0] += 0.02; // adds a sinusoidal oscillation along time and space

	/* Exploit symmetry of Fourier transforms of real-valued signals 
	   to fill second half of P0: */
	for (kzcnt = ceil(Nt / 2.0) + 1; kzcnt < Nt; kzcnt++) {
		for (kxcnt = 0; kxcnt < Nsrc; kxcnt++) {
			//P0[kxcnt + kzcnt * Nsrc][0] = 
			//P0[kxcnt + kzcnt * Nsrc][0] = -1.0 * 
		}
	}

	fftw_execute(plan_inv);
	for (kxcnt = 0; kxcnt < Nt * Nsrc; kxcnt++) { image[kxcnt] = P0[kxcnt][0]; }

	fftw_destroy_plan(plan_fwd); fftw_destroy_plan(plan_inv);
	fftw_free(RF); fftw_free(P0);
}

void compute_factor(fftw_complex constant, double dist, fftw_complex output) {
	std::complex<double> temporary;
	const std::complex<double> i(0, 1.0);

	/* *** The updated line below now results in mostly correct values of "factor", but suffers
	   from a few minus-sign issues!! *** 
	   This is most likely an indexing issue in the search for "nearest", but unconfirmed.
	 */
	temporary = (constant[0], constant[1]*i);
	temporary *= dist;
	temporary = (1.0 - exp(-temporary)) / temporary;

	output[0] = real(temporary);
	output[1] = imag(temporary);
}