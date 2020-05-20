// MISI_ImgAlg.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include <stdio.h> 
#include <stdlib.h>

#include "MISI_ImgAlg.h"
#include "support_scripts.h"
#include "delay_and_sum.h"
#include "DMAS.h"
#include "SLSC.h"
#include "kSpace.h"






//==========================================================================================================================
//
// DELAY AND SUM:
//
//==========================================================================================================================

__declspec(dllexport) void DnS_1rec_fixed_pos(	double *rf_data, double *source_locations,
												double *receiver_location, double *image_coordinates,
												double c, double fsamp,
												int Nsrc, int Nt, int Nimg,
												double *image)
{
    /* call the computational routine */
	DnS_1rec_fixed_pos_alg(rf_data, source_locations, receiver_location, image_coordinates, c, fsamp, Nsrc, Nt, Nimg, image);
	return;
}

__declspec(dllexport) void DnS_1rec_fixed_pos_precomp(	double *source_locations,
														double *receiver_location, double *image_coordinates,
														double c, double fsamp,
														int Nsrc, int Nimg,
														short unsigned int *delays)
{
    /* call the computational routine */
	DnS_1rec_fixed_pos_precomp_delays_alg(source_locations,receiver_location,image_coordinates,c,fsamp,Nsrc,Nimg,delays);
	return;
}

__declspec(dllexport) void DnS_1rec_fixed_pos_from_precomp(	double *rf_data,
															short unsigned int *delays,
															int Nsrc, int Nt, int Nimg,
															double *image)
{
    /* call the computational routine */
	DnS_1rec_fixed_pos_from_precomp_alg(rf_data,delays,Nsrc,Nt,Nimg,image);
	return;
}




__declspec(dllexport) void DnS_1rec_as_src(	double *rf_data, double *source_locations,
											double *image_coordinates,
											double c, double fsamp,
											int Nsrc, int Nt, int Nimg,
											double *image)
{
    /* call the computational routine */
	DnS_1rec_as_src_alg(rf_data, source_locations, image_coordinates, c, fsamp, Nsrc, Nt, Nimg, image);
	return;
}




__declspec(dllexport) void DnS_1rec_free_pos(	double *rf_data, double *source_locations,
												double *receiver_location, double *image_coordinates,
												double c, double fsamp,
												int Nsrc, int Nt, int Nimg,
												double *image)
{
    /* call the computational routine */
	DnS_1rec_free_pos_alg(rf_data, source_locations, receiver_location, image_coordinates, c, fsamp, Nsrc, Nt, Nimg, image);
	return;
}






//==========================================================================================================================
//
// DELAY, MULTIPLY AND SUM:
//
//==========================================================================================================================

__declspec(dllexport) void DMnS_1rec_fixed_pos(	double* rf_data, double* source_locations,
												double* receiver_location, double* image_coordinates,
												double c, double fsamp,
												int Nsrc, int Nt, int Nimg,
												double* image)
{
	/* call the computational routine */
	DMnS_1rec_fixed_pos_alg(rf_data, source_locations, receiver_location, image_coordinates, c, fsamp, Nsrc, Nt, Nimg, image);
	return;
}

__declspec(dllexport) void DMnS_1rec_as_src(double* rf_data, double* source_locations,
											double* image_coordinates,
											double c, double fsamp,
											int Nsrc, int Nt, int Nimg,
											double* image)
{
	/* call the computational routine */
	DMnS_1rec_as_src_alg(rf_data, source_locations, image_coordinates, c, fsamp, Nsrc, Nt, Nimg, image);
	return;
}

__declspec(dllexport) void DMnS_1rec_free_pos(	double* rf_data, double* source_locations,
												double* receiver_locations, double* image_coordinates,
												double c, double fsamp,
												int Nsrc, int Nt, int Nimg,
												double* image)
{
	/* call the computational routine */
	DMnS_1rec_free_pos_alg(rf_data, source_locations, receiver_locations, image_coordinates, c, fsamp, Nsrc, Nt, Nimg, image);
	return;
}




//==========================================================================================================================
//
// Supporting scripts:
//
//==========================================================================================================================

__declspec(dllexport) void Detect_envelope( double *image, int Nx, int Ny, int Nz)
{
	Compute_envelope_FFTW(image,Nx,Ny,Nz);
	return;
}

__declspec(dllexport) void PulseCompChirp(double* in, double dur, double fmin, double fmax, double fsamp,
										  double Tukpar, int Nt, int Nscan, double* out)
{
	PCChirp(in, dur, fmin, fmax, fsamp, Tukpar, Nt, Nscan, out);
	return;
}




//==========================================================================================================================
//
// SLSC:
//
//==========================================================================================================================

__declspec(dllexport) void SLSC_1rec_fixed_pos(double* rf_data, double* source_locations,
	double* receiver_location,
	double* image_coordinates,
	double c, double fsamp,
	int Nsrc, int Nt, int Nimg, int m, int w,
	double* R)
{
	/* call the computational routine */
	SLSC_1rec_fixed_pos_alg(rf_data, source_locations, receiver_location, image_coordinates, c, fsamp, Nsrc, Nt, Nimg, m, w, R);
	return;
}

__declspec(dllexport) void SLSC_1rec_as_src(double* rf_data, double* source_locations,
	double* image_coordinates,
	double c, double fsamp,
	int Nsrc, int Nt, int Nimg, int m, int w,
	double* R)
{
	/* call the computational routine */
	SLSC_1rec_as_src_alg(rf_data, source_locations, image_coordinates, c, fsamp, Nsrc, Nt, Nimg, m, w, R);
	return;
}

__declspec(dllexport) void SLSC_1rec_free_pos(double* rf_data, double* source_locations,
	double* receiver_locations,
	double* image_coordinates,
	double c, double fsamp,
	int Nsrc, int Nt, int Nimg, int m, int w,
	double* R)
{
	/* call the computational routine */
	SLSC_1rec_free_pos_alg(rf_data, source_locations, receiver_locations, image_coordinates, c, fsamp, Nsrc, Nt, Nimg, m, w, R);
	return;
}





//==========================================================================================================================
//
// k space:
//
//==========================================================================================================================

__declspec(dllexport) void kSpace_line_recon(double* rf_data,
	double* source_locations,
	double c, double fsamp,
	int Nsrc, int Nt, int N,
	double* xaxis, double* zaxis,
	double* image)
{
	/* call the computational routine */
	kSpace_line_recon_alg(rf_data, source_locations, c, fsamp, Nsrc, Nt, N, xaxis, zaxis, image);
	return;
}