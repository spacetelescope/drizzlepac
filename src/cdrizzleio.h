#ifndef CDRIZZLEIO_H
#define CDRIZZLEIO_H

#include "cdrizzleutil.h"
#include "fitsio.h"

/**
Get the geometrical distortion information, either from a file or the
header.

@param[in] coeff_source A string describing where to load the
coefficients from. (\todo Improve description)

@param[in] input_data_file A \a cfitsio file, already opened for
reading.  Only required if \a coeff_source is "header" or "wcs",
otherwise may be NULL.

@param[out] lambda Wavelength in nm.

@param[out] coty The type of coefficients being stored (\todo Add
description)

@param[out] num_coeffs The number of coefficients being stored.

@param[out] x_coeffs An array of x coefficients.  This array should be
            allocated to size MAX_COEFFS, but only the first
            num_coeffs will be valid after returning from this
            function.

@param[out] y_coeffs An array of y coefficients.  This array should be
            allocated to size MAX_COEFFS, but only the first
            num_coeffs will be valid after returning from this
            function.

@param[out] error

@return Non-zero if an error occurred.

was: GETGEO
*/
int
get_geometric_distortion(const char* coeff_source,
                         fitsfile* input_data_file,
                         /* Output arguments */
                         double* lambda,
                         integer_t* coty,
                         integer_t* num_coeffs,
                         double* x_coeffs,
                         double* y_coeffs,
                         struct driz_error_t* error);

#endif /* CDRIZZLEIO_H */
