#include <stdafx.h>
#include <math.h>
#include <stdlib.h>
#include "support_scripts.h"



//==========================================================================================================================
//
// a single, fixed receiver position interrogating pulse-echo signals originating from multiple source locations:
//
//==========================================================================================================================

/* The computational routine - no TGC */
void DMnS_1rec_fixed_pos_alg(   double* rf_data,
                                double* source_locations,
                                double* receiver_location,
                                double* image_coordinates,
                                double c, double fsamp,
                                int Nsrc, int Nt, int Nimg,
                                double* image)
{
    int srccnt, imgcnt, tcnt, icnt, jcnt;
    int Nsrc2 = 2 * Nsrc, Nimg2 = 2 * Nimg;
    double feff = fsamp / c;
    double ximg, yimg, zimg;
    double dist_src_refl, dist_refl_hydr;
    double a, b, product;

    double* delayed;
    delayed = (double*)malloc(Nsrc * sizeof(double));

    for (imgcnt = 0; imgcnt < Nimg; imgcnt++) {
        ximg = image_coordinates[imgcnt];
        yimg = image_coordinates[imgcnt + Nimg];
        zimg = image_coordinates[imgcnt + Nimg2];
        dist_refl_hydr = sqrt(squa(ximg - receiver_location[0]) +
            squa(yimg - receiver_location[1]) +
            squa(zimg - receiver_location[2]));
        for (srccnt = 0; srccnt < Nsrc; srccnt++) {
            dist_src_refl = sqrt(squa(ximg - source_locations[srccnt]) +
                squa(yimg - source_locations[srccnt + Nsrc]) +
                squa(zimg - source_locations[srccnt + Nsrc2]));
            tcnt = int(round((dist_refl_hydr + dist_src_refl) * feff));
            if (tcnt < Nt) { delayed[srccnt] = rf_data[tcnt + srccnt * Nt]; }
            else { delayed[srccnt] = 0.0; }
        }
        for (icnt = 0; icnt < Nsrc - 1; icnt++) {
            //product = delayed[icnt];
            a = delayed[icnt];
            for (jcnt = icnt + 1; jcnt < Nsrc; jcnt++) {
                //product *= delayed[jcnt];
                //product /= sqrt( fabs( product ) );
                b = delayed[jcnt];
                product = sqrt(fabs(a * b));
                product = _copysign(product, a * b);
                image[imgcnt] += product;
            }
        }
    }
    free(delayed);
    return;
}




//==========================================================================================================================
//
// The receiver coincides with the source, which is only approximately true in fibre assemblies:
//
//==========================================================================================================================

/* The computational routine - no TGC */
void DMnS_1rec_as_src_alg(  double* rf_data,
                            double* source_locations,
                            double* image_coordinates,
                            double c, double fsamp,
                            int Nsrc, int Nt, int Nimg,
                            double* image)
{
    int srccnt, imgcnt, tcnt, icnt, jcnt;
    int Nsrc2 = 2 * Nsrc, Nimg2 = 2 * Nimg;
    double feff = fsamp / c;
    double ximg, yimg, zimg;
    double dist_src_refl;
    double a, b, product;

    double* delayed;
    delayed = (double*)malloc(Nsrc * sizeof(double));

    for (imgcnt = 0; imgcnt < Nimg; imgcnt++) {
        ximg = image_coordinates[imgcnt];
        yimg = image_coordinates[imgcnt + Nimg];
        zimg = image_coordinates[imgcnt + Nimg2];
        for (srccnt = 0; srccnt < Nsrc; srccnt++) {
            dist_src_refl = sqrt(squa(ximg - source_locations[srccnt]) +
                squa(yimg - source_locations[srccnt + Nsrc]) +
                squa(zimg - source_locations[srccnt + Nsrc2]));
            tcnt = int(round(2.0 * dist_src_refl * feff));
            if (tcnt < Nt) { delayed[srccnt] = rf_data[tcnt + srccnt * Nt]; }
            else { delayed[srccnt] = 0.0; }
        }
        for (icnt = 0; icnt < Nsrc - 1; icnt++) {
            //product = delayed[icnt];
            a = delayed[icnt];
            for (jcnt = icnt + 1; jcnt < Nsrc; jcnt++) {
                //product *= delayed[jcnt];
                //product /= sqrt( fabs( product ) );
                b = delayed[jcnt];
                product = sqrt(fabs(a * b));
                product = _copysign(product, a * b);
                image[imgcnt] += product;
            }
        }
    }
    free(delayed);
    return;
}




//==========================================================================================================================
//
// Both receiver and source have different positions at each A-scan:
//
//==========================================================================================================================

/* The computational routine - no TGC */
void DMnS_1rec_free_pos_alg(double* rf_data,
                            double* source_locations,
                            double* receiver_locations,
                            double* image_coordinates,
                            double c, double fsamp,
                            int Nsrc, int Nt, int Nimg,
                            double* image)
{
    int srccnt, imgcnt, tcnt, icnt, jcnt;
    int Nsrc2 = 2 * Nsrc, Nimg2 = 2 * Nimg;
    double feff = fsamp / c;
    double ximg, yimg, zimg;
    double dist_src_refl, dist_refl_hydr;
    double a, b, product;

    double* delayed;
    delayed = (double*)malloc(Nsrc * sizeof(double));

    for (imgcnt = 0; imgcnt < Nimg; imgcnt++) {
        ximg = image_coordinates[imgcnt];
        yimg = image_coordinates[imgcnt + Nimg];
        zimg = image_coordinates[imgcnt + Nimg2];
        for (srccnt = 0; srccnt < Nsrc; srccnt++) {
            dist_refl_hydr = sqrt(squa(ximg - receiver_locations[srccnt]) +
                squa(yimg - receiver_locations[srccnt + Nsrc]) +
                squa(zimg - receiver_locations[srccnt + Nsrc2]));
            dist_src_refl = sqrt(squa(ximg - source_locations[srccnt]) +
                squa(yimg - source_locations[srccnt + Nsrc]) +
                squa(zimg - source_locations[srccnt + Nsrc2]));
            tcnt = int(round((dist_refl_hydr + dist_src_refl) * feff));
            if (tcnt < Nt) { delayed[srccnt] = rf_data[tcnt + srccnt * Nt]; }
            else { delayed[srccnt] = 0.0; }
        }
        for (icnt = 0; icnt < Nsrc - 1; icnt++) {
            //product = delayed[icnt];
            a = delayed[icnt];
            for (jcnt = icnt + 1; jcnt < Nsrc; jcnt++) {
                //product *= delayed[jcnt];
                //product /= sqrt( fabs( product ) );
                b = delayed[jcnt];
                product = sqrt(fabs(a * b));
                product = _copysign(product, a * b);
                image[imgcnt] += product;
            }
        }
    }
    return;
}