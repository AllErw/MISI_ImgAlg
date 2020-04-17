#include "stdafx.h"
#include <math.h>
#include <stdlib.h>             // for malloc and calloc
#include "support_scripts.h"    // for "squa" which is faster than "pow"




/* SLSC for the case of a single receiver that moves together with
   the source - approximately true for 2-fibre AOUS probes. */
void SLSC_1rec_fixed_pos_alg(double* rf_data,
    double* source_locations,
    double* receiver_location,
    double* image_coordinates,
    double c, double fsamp,
    int Nsrc, int Nt, int Nimg, int m, int w,
    double* R)
{
    int srccnt, imgcnt, tcnt, mcnt, wcnt, index1, index2;
    int Nsrc2 = 2 * Nsrc, Nimg2 = 2 * Nimg;
    double feff = fsamp / c;
    double ximg, yimg, zimg;
    double dist_src_refl, dist_refl_hydr, Rhat, numerator, denominator, denominator1, denominator2;
    double* rf_data2;
    int* arr_times;

    // Initialise and populate rf_data.^2, and set R to zero:
    rf_data2 = (double*)calloc(Nsrc * Nt, sizeof(double));
    for (tcnt = 0; tcnt < Nsrc * Nt; tcnt++) {
        rf_data2[tcnt] = rf_data[tcnt] * rf_data[tcnt];
    }
    for (tcnt = 0; tcnt < Nimg; tcnt++) { R[tcnt] = 0.0; }

    // Initialise arrival times vector and Rhat:
    arr_times = (int*)calloc(Nsrc, sizeof(int));

    // The actual computation:
    for (imgcnt = 0; imgcnt < Nimg; imgcnt++) {// Loop over all image pixels
        // Extract coordinates of current pixel:
        ximg = image_coordinates[imgcnt];
        yimg = image_coordinates[imgcnt + Nimg];
        zimg = image_coordinates[imgcnt + Nimg2];

        dist_refl_hydr = sqrt( squa(ximg - receiver_location[0]) +
                               squa(yimg - receiver_location[1]) +
                               squa(zimg - receiver_location[2]));

        // Compute the arrival times across the aperture:
        for (srccnt = 0; srccnt < Nsrc; srccnt++) {
            dist_src_refl = sqrt(squa(ximg - source_locations[srccnt]) +
                squa(yimg - source_locations[srccnt + Nsrc]) +
                squa(zimg - source_locations[srccnt + Nsrc2]));
            arr_times[srccnt] = int(round( (dist_refl_hydr + dist_src_refl) * feff));
        }

        for (mcnt = 0; mcnt < m; mcnt++) {
            // Re-initialise Rhat to zero:
            Rhat = 0.0;

            for (srccnt = 0; srccnt < Nsrc - mcnt; srccnt++) {
                numerator = 0.0; denominator1 = 0.0; denominator2 = 0.0;

                for (wcnt = 0; wcnt < w; wcnt++) {
                    index1 = wcnt + arr_times[srccnt] + srccnt * Nt;
                    index2 = wcnt + arr_times[srccnt + mcnt] + (srccnt + mcnt) * Nt;
                    if ((wcnt + arr_times[srccnt + mcnt]) < Nt) {
                        numerator += rf_data[index1] * rf_data[index2];
                        denominator1 += rf_data2[index1];
                        denominator2 += rf_data2[index2];
                    }
                }
                denominator = sqrt(denominator1 * denominator2);

                if (denominator == 0.0) { denominator = 0.000000001; } // avoid divide by 0
                Rhat += numerator / denominator;
            }
            Rhat /= 1.0 * (Nsrc - mcnt);
            R[imgcnt] += Rhat;
        }
    }
    free(rf_data2);
    free(arr_times);
    return;
}

/* SLSC for the case of a single stationary receiver combined with
   an arbitrary number of sources - suitable for the freespace and
   freehand imaging setups. */
void SLSC_1rec_as_src_alg(double* rf_data, 
    double* source_locations,
    double* image_coordinates,
    double c, double fsamp,
    int Nsrc, int Nt, int Nimg, int m, int w,
    double* R)
{
    int srccnt, imgcnt, tcnt, mcnt, wcnt, index1, index2;
    int Nsrc2 = 2 * Nsrc, Nimg2 = 2 * Nimg;
    double feff = fsamp / c;
    double ximg, yimg, zimg;
    double dist_src_refl, Rhat, numerator, denominator, denominator1, denominator2;
    double* rf_data2;
    int* arr_times;
    
    // Initialise and populate rf_data.^2, and set R to zero:
    rf_data2 = (double*)calloc(Nsrc*Nt , sizeof(double));
    for (tcnt = 0; tcnt < Nsrc * Nt; tcnt++) {
        rf_data2[tcnt] = rf_data[tcnt] * rf_data[tcnt];
    }
    for (tcnt = 0; tcnt < Nimg; tcnt++) { R[tcnt] = 0.0; }

    // Initialise arrival times vector and Rhat:
    arr_times = (int*)calloc(Nsrc, sizeof(int));

    // The actual computation:
    for (imgcnt = 0; imgcnt < Nimg; imgcnt++) {// Loop over all image pixels
        // Extract coordinates of current pixel:
        ximg = image_coordinates[imgcnt];
        yimg = image_coordinates[imgcnt + Nimg];
        zimg = image_coordinates[imgcnt + Nimg2];

        // Compute the arrival times across the aperture:
        for (srccnt = 0; srccnt < Nsrc; srccnt++) {
            dist_src_refl = sqrt(squa(ximg - source_locations[srccnt]) +
                                 squa(yimg - source_locations[srccnt + Nsrc]) +
                                 squa(zimg - source_locations[srccnt + Nsrc2]));
            arr_times[srccnt] = int(round(2.0 * dist_src_refl * feff));
        }

        for (mcnt = 0; mcnt < m; mcnt++) {
            // Re-initialise Rhat to zero:
            Rhat = 0.0;

            for (srccnt = 0; srccnt < Nsrc-mcnt; srccnt++) {
                numerator = 0.0; denominator1 = 0.0; denominator2 = 0.0;

                for (wcnt = 0; wcnt < w; wcnt++) {
                    index1 = wcnt + arr_times[srccnt] + srccnt * Nt;
                    index2 = wcnt + arr_times[srccnt+mcnt] + (srccnt+mcnt) * Nt;
                    if ( (wcnt + arr_times[srccnt + mcnt]) < Nt) {
                        numerator += rf_data[index1] * rf_data[index2];
                        denominator1 += rf_data2[index1];
                        denominator2 += rf_data2[index2];
                    }
                }
                denominator = sqrt(denominator1 * denominator2);

                if (denominator == 0.0) { denominator = 0.000000001; } // avoid divide by 0
                Rhat += numerator / denominator;
            }
            Rhat /= 1.0*(Nsrc-mcnt);
            R[imgcnt] += Rhat;
        }
    }
    free(rf_data2);
    free(arr_times);
    return;
}

/* SLSC for the case where for each A-scan the source and detector
   are in different positions. This is the most general, and slowest,
   case that can in principle be used for virtually any imaging system. */
void SLSC_1rec_free_pos_alg(double* rf_data,
    double* source_locations,
    double* receiver_locations,
    double* image_coordinates,
    double c, double fsamp,
    int Nsrc, int Nt, int Nimg, int m, int w,
    double* R)
{
    int srccnt, imgcnt, tcnt, mcnt, wcnt, index1, index2;
    int Nsrc2 = 2 * Nsrc, Nimg2 = 2 * Nimg;
    double feff = fsamp / c;
    double ximg, yimg, zimg;
    double dist_src_refl, dist_refl_hydr, Rhat, numerator, denominator, denominator1, denominator2;
    double* rf_data2;
    int* arr_times;

    // Initialise and populate rf_data.^2, and set R to zero:
    rf_data2 = (double*)calloc(Nsrc * Nt, sizeof(double));
    for (tcnt = 0; tcnt < Nsrc * Nt; tcnt++) {
        rf_data2[tcnt] = rf_data[tcnt] * rf_data[tcnt];
    }
    for (tcnt = 0; tcnt < Nimg; tcnt++) { R[tcnt] = 0.0; }

    // Initialise arrival times vector and Rhat:
    arr_times = (int*)calloc(Nsrc, sizeof(int));

    // The actual computation:
    for (imgcnt = 0; imgcnt < Nimg; imgcnt++) {// Loop over all image pixels
        // Extract coordinates of current pixel:
        ximg = image_coordinates[imgcnt];
        yimg = image_coordinates[imgcnt + Nimg];
        zimg = image_coordinates[imgcnt + Nimg2];

        // Compute the arrival times across the aperture:
        for (srccnt = 0; srccnt < Nsrc; srccnt++) {
            dist_refl_hydr = sqrt(squa(ximg - receiver_locations[srccnt]) +
                squa(yimg - receiver_locations[srccnt + Nsrc]) +
                squa(zimg - receiver_locations[srccnt + Nsrc2]));
            dist_src_refl = sqrt(squa(ximg - source_locations[srccnt]) +
                squa(yimg - source_locations[srccnt + Nsrc]) +
                squa(zimg - source_locations[srccnt + Nsrc2]));
            arr_times[srccnt] = int(round((dist_refl_hydr + dist_src_refl) * feff));
        }

        for (mcnt = 0; mcnt < m; mcnt++) {
            // Re-initialise Rhat to zero:
            Rhat = 0.0;

            for (srccnt = 0; srccnt < Nsrc - mcnt; srccnt++) {
                numerator = 0.0; denominator1 = 0.0; denominator2 = 0.0;

                for (wcnt = 0; wcnt < w; wcnt++) {
                    index1 = wcnt + arr_times[srccnt] + srccnt * Nt;
                    index2 = wcnt + arr_times[srccnt + mcnt] + (srccnt + mcnt) * Nt;
                    if ((wcnt + arr_times[srccnt + mcnt]) < Nt) {
                        numerator += rf_data[index1] * rf_data[index2];
                        denominator1 += rf_data2[index1];
                        denominator2 += rf_data2[index2];
                    }
                }
                denominator = sqrt(denominator1 * denominator2);

                if (denominator == 0.0) { denominator = 0.000000001; } // avoid divide by 0
                Rhat += numerator / denominator;
            }
            Rhat /= 1.0 * (Nsrc - mcnt);
            R[imgcnt] += Rhat;
        }
    }
    free(rf_data2);
    free(arr_times);
    return;
}