//==========================================================================================================================
//
// DELAY AND SUM:
//
//==========================================================================================================================

#ifdef __cplusplus
extern "C" {
#endif
__declspec(dllexport) void DnS_1rec_fixed_pos(	double* rf_data, double* source_locations,
												double* receiver_location, double* image_coordinates,
												double c, double fsamp,
												int Nsrc, int Nt, int Nimg,
												double* image);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
__declspec(dllexport) void DnS_1rec_fixed_pos_precomp(	double* source_locations,
														double* receiver_location, double* image_coordinates,
														double c, double fsamp,
														int Nsrc, int Nimg,
														short unsigned int* delays);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
__declspec(dllexport) void DnS_1rec_fixed_pos_from_precomp(	double* rf_data,
															short unsigned int* delays,
															int Nsrc, int Nt, int Nimg,
															double *image);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
__declspec(dllexport) void DnS_1rec_as_src(	double* rf_data, double* source_locations,
											double* image_coordinates,
											double c, double fsamp,
											int Nsrc, int Nt, int Nimg,
											double* image);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
__declspec(dllexport) void DnS_1rec_free_pos(	double* rf_data, double* source_locations,
												double* receiver_location, double* image_coordinates,
												double c, double fsamp,
												int Nsrc, int Nt, int Nimg,
												double* image);
#ifdef __cplusplus
}
#endif





//==========================================================================================================================
//
// DELAY, MULTIPLY AND SUM:
//
//==========================================================================================================================

#ifdef __cplusplus
extern "C" {
#endif
	__declspec(dllexport) void DMnS_1rec_fixed_pos(	double* rf_data, double* source_locations,
													double* receiver_location, double* image_coordinates,
													double c, double fsamp,
													int Nsrc, int Nt, int Nimg,
													double* image);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
	__declspec(dllexport) void DMnS_1rec_as_src(double* rf_data, double* source_locations,
												double* image_coordinates,
												double c, double fsamp,
												int Nsrc, int Nt, int Nimg,
												double* image);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
	__declspec(dllexport) void DMnS_1rec_free_pos(	double* rf_data, double* source_locations,
													double* receiver_locations, double* image_coordinates,
													double c, double fsamp,
													int Nsrc, int Nt, int Nimg,
													double* image);
#ifdef __cplusplus
}
#endif




//==========================================================================================================================
//
// Supporting functions:
//
//==========================================================================================================================

#ifdef __cplusplus
extern "C" {
#endif
__declspec(dllexport) void Detect_envelope( double *image, int Nx, int Ny, int Nz);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
	__declspec(dllexport) void PulseCompChirp(double* in, double dur, double fmin, double fmax, double fsamp,
		double Tukpar, int Nt, int Nscan, double* out);
#ifdef __cplusplus
}
#endif




//==========================================================================================================================
//
// SLSC:
//
//==========================================================================================================================

#ifdef __cplusplus
extern "C" {
#endif
	__declspec(dllexport) void SLSC_1rec_fixed_pos(double* rf_data, double* source_locations,
		double* receiver_location,
		double* image_coordinates,
		double c, double fsamp,
		int Nsrc, int Nt, int Nimg, int m, int w,
		double* R);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
	__declspec(dllexport) void SLSC_1rec_as_src(double* rf_data, double* source_locations,
		double* image_coordinates,
		double c, double fsamp,
		int Nsrc, int Nt, int Nimg, int m, int w,
		double* R);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
	__declspec(dllexport) void SLSC_1rec_free_pos(double* rf_data, double* source_locations,
		double* receiver_locations,
		double* image_coordinates,
		double c, double fsamp,
		int Nsrc, int Nt, int Nimg, int m, int w,
		double* R);
#ifdef __cplusplus
}
#endif





//==========================================================================================================================
//
// k space:
//
//==========================================================================================================================

#ifdef __cplusplus
extern "C" {
#endif
	__declspec(dllexport) void kSpace_line_recon(double* rf_data,
		double* source_locations,
		double c, double fsamp,
		int Nsrc, int Nt, int N,
		double* xaxis, double* zaxis,
		double* image);
#ifdef __cplusplus
}
#endif
