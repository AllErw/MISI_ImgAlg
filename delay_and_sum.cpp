#include <stdafx.h>
#include <math.h>
#include "support_scripts.h"



//==========================================================================================================================
//
// a single, fixed receiver position interrogating pulse-echo signals originating from multiple source locations:
//
//==========================================================================================================================

/* The computational routine - no TGC */
void DnS_1rec_fixed_pos_alg(double *rf_data,
		                    double *source_locations,
				            double *receiver_location,
						    double *image_coordinates,
							double c, double fsamp,
		                    int Nsrc, int Nt, int Nimg,
				            double *image)
// Compute D&S image in the case of a single receiver at a fixed location using TGC_power = 1
{
    int srccnt,imgcnt,tcnt;
    int Nsrc2 = 2*Nsrc,     Nimg2 = 2*Nimg;
    double feff = fsamp/c;
    double ximg,yimg,zimg;
	double dist_src_refl,dist_refl_hydr;
    
    for (imgcnt=0; imgcnt<Nimg; imgcnt++){
        ximg = image_coordinates[imgcnt];
        yimg = image_coordinates[imgcnt+Nimg];
        zimg = image_coordinates[imgcnt+Nimg2];
        dist_refl_hydr = sqrt( squa(ximg - receiver_location[0]) + 
                               squa(yimg - receiver_location[1]) + 
                               squa(zimg - receiver_location[2]) );
        for (srccnt=0; srccnt<Nsrc; srccnt++){
            dist_src_refl = sqrt( squa(ximg - source_locations[srccnt]) + 
                                  squa(yimg - source_locations[srccnt+Nsrc]) + 
                                  squa(zimg - source_locations[srccnt+Nsrc2]) );
            tcnt = int(round( (dist_refl_hydr + dist_src_refl)*feff ));
            if (tcnt < Nt){
                image[imgcnt] += rf_data[tcnt + srccnt*Nt];
            }
        }
    }
    return;
}

/* Pre-compute time delays for all source locations in all image locations */
void DnS_1rec_fixed_pos_precomp_delays_alg(	double *source_locations,
					                        double *receiver_location,
							                double *image_coordinates,
    									    double c, double fsamp, 
	    							        int Nsrc, int Nimg,
		    								short unsigned int *delays)
// Precompute delays for D&S imaging - single receiver at a fixed location
{
    int srccnt,imgcnt,tcnt,imgcntNsrc;
    int Nsrc2 = 2*Nsrc,     Nimg2 = 2*Nimg;
    double feff = fsamp/c;
    double ximg,yimg,zimg;
	double dist_src_refl,dist_refl_hydr;
    
    for (imgcnt=0; imgcnt<Nimg; imgcnt++){
        ximg = image_coordinates[imgcnt];
        yimg = image_coordinates[imgcnt+Nimg];
        zimg = image_coordinates[imgcnt+Nimg2];
        dist_refl_hydr = sqrt( squa(ximg - receiver_location[0]) + 
                               squa(yimg - receiver_location[1]) + 
                               squa(zimg - receiver_location[2]) );
		imgcntNsrc = imgcnt * Nsrc;
        for (srccnt=0; srccnt<Nsrc; srccnt++){
            dist_src_refl = sqrt( squa(ximg - source_locations[srccnt]) + 
                                  squa(yimg - source_locations[srccnt+Nsrc]) + 
                                  squa(zimg - source_locations[srccnt+Nsrc2]) );
            tcnt = int(round( (dist_refl_hydr + dist_src_refl)*feff ));
			delays[srccnt + imgcntNsrc] = tcnt;
        }
    }
    return;
}

/* The computational routine - no TGC */
void DnS_1rec_fixed_pos_from_precomp_alg(   double *rf_data,
										    short unsigned int *delays,
    									    int Nsrc, int Nt, int Nimg,
	    									double *image)
// Compute D&S image in the case of a single receiver at a fixed location from precomputed delays:
{
    int srccnt,imgcnt,tcnt,imgcntNsrc;

	for (imgcnt=0; imgcnt<Nimg; imgcnt++){
		imgcntNsrc = imgcnt * Nsrc;
        for (srccnt=0; srccnt<Nsrc; srccnt++){
			tcnt = delays[srccnt + imgcntNsrc];
            if (tcnt < Nt){
                image[imgcnt] += rf_data[tcnt + srccnt*Nt];
            }
        }
    }
    return;
}




//==========================================================================================================================
//
// The receiver coincides with the source, which is only approximately true in fibre assemblies:
//
//==========================================================================================================================

/* The computational routine - no TGC */
void DnS_1rec_as_src_alg(   double *rf_data,
	                        double *source_locations,
		                    double *image_coordinates,
    			            double c, double fsamp,
	    			        int Nsrc, int Nt, int Nimg,
		    			    double *image)
{
    int srccnt,imgcnt,tcnt;
    int Nsrc2 = 2*Nsrc,     Nimg2 = 2*Nimg;
    double feff = fsamp/c;
    double ximg,yimg,zimg;
	double dist_src_refl;
    
    for (imgcnt=0; imgcnt<Nimg; imgcnt++){
        ximg = image_coordinates[imgcnt];
        yimg = image_coordinates[imgcnt+Nimg];
        zimg = image_coordinates[imgcnt+Nimg2];
        for (srccnt=0; srccnt<Nsrc; srccnt++){
            dist_src_refl = sqrt( squa(ximg - source_locations[srccnt]) + 
                                  squa(yimg - source_locations[srccnt+Nsrc]) + 
                                  squa(zimg - source_locations[srccnt+Nsrc2]) );
            tcnt = int(round( 2.0*dist_src_refl*feff ));
            if (tcnt < Nt){
                image[imgcnt] += rf_data[tcnt + srccnt*Nt];
            }
        }
    }
    return;
}




//==========================================================================================================================
//
// Both receiver and source have different positions at each A-scan:
//
//==========================================================================================================================

/* The computational routine - no TGC */
void DnS_1rec_free_pos_alg( double *rf_data,
    						double *source_locations,
	    	                double *receiver_locations,
		    		        double *image_coordinates,
			    			double c, double fsamp,
				    		int Nsrc, int Nt, int Nimg,
		                    double *image)
{
    int srccnt,imgcnt,tcnt;
    int Nsrc2 = 2*Nsrc,     Nimg2 = 2*Nimg;
    double feff = fsamp/c;
    double ximg,yimg,zimg;
	double dist_src_refl,dist_refl_hydr;
    
    for (imgcnt=0; imgcnt<Nimg; imgcnt++){
        ximg = image_coordinates[imgcnt];
        yimg = image_coordinates[imgcnt+Nimg];
        zimg = image_coordinates[imgcnt+Nimg2];
        for (srccnt=0; srccnt<Nsrc; srccnt++){
            dist_refl_hydr = sqrt( squa(ximg - receiver_locations[srccnt]) + 
                                   squa(yimg - receiver_locations[srccnt+Nsrc]) + 
                                   squa(zimg - receiver_locations[srccnt+Nsrc2]) );
            dist_src_refl = sqrt( squa(ximg - source_locations[srccnt]) + 
                                  squa(yimg - source_locations[srccnt+Nsrc]) + 
                                  squa(zimg - source_locations[srccnt+Nsrc2]) );
            tcnt = int(round( (dist_refl_hydr + dist_src_refl)*feff ));
            if (tcnt < Nt){
                image[imgcnt] += rf_data[tcnt + srccnt*Nt];
            }
        }
    }
    return;
}