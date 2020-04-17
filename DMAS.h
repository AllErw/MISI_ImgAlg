void DMnS_1rec_fixed_pos_alg(	double* rf_data,
								double* source_locations,
								double* receiver_location,
								double* image_coordinates,
								double c, double fsamp,
								int Nsrc, int Nt, int Nimg,
								double* image);

void DMnS_1rec_as_src_alg(	double* rf_data,
							double* source_locations,
							double* image_coordinates,
							double c, double fsamp,
							int Nsrc, int Nt, int Nimg,
							double* image);

void DMnS_1rec_free_pos_alg(double* rf_data,
							double* source_locations,
							double* receiver_locations,
							double* image_coordinates,
							double c, double fsamp,
							int Nsrc, int Nt, int Nimg,
							double* image);