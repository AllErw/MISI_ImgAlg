void DnS_1rec_fixed_pos_alg(double *rf_data,
							double *source_locations,
							double *receiver_location,
			                double *image_coordinates,
					        double c, double fsamp,
							int Nsrc, int Nt, int Nimg,
							double *image);
void DnS_1rec_fixed_pos_precomp_delays_alg(	double *source_locations,
						                    double *receiver_location,
								            double *image_coordinates,
										    double c, double fsamp, 
									        int Nsrc, int Nimg,
											short unsigned int *delays);
void DnS_1rec_fixed_pos_from_precomp_alg(	double *rf_data,
											short unsigned int *delays,
										    int Nsrc, int Nt, int Nimg,
											double *image);




void DnS_1rec_as_src_alg(	double *rf_data,
								double *source_locations,
								double *image_coordinates,
								double c, double fsamp,
								int Nsrc, int Nt, int Nimg,
								double *image);




void DnS_1rec_free_pos_alg(	double *rf_data,
							double *source_locations,
							double *receiver_location,
			                double *image_coordinates,
					        double c, double fsamp,
							int Nsrc, int Nt, int Nimg,
							double *image);