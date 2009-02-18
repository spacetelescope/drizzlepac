#ifndef CDRIZZLEWCS_H
#define CDRIZZLEWCS_H

#include "cdrizzlemap.h"
#include "cdrizzleutil.h"

/**
Derive best linear transformation coefficients to map
pixel positions from one world coordinate system to another.

was: WCSLIN
*/
int
wcs_derive_linear(struct mapping_param_t* m,
                  /* Output parameters */
                  double* xc, double* yc, double* xs, double* ys,
                  double* xt, double* yt,
                  struct driz_error_t* error);

/**
Update the WCS to include the drizzling transformations.

This is done by applying the transform to a unit square at the centre
pixel in the input whilst retaining the same reference point in the
sky.
*/
int update_wcs(struct driz_param_t* p,
               struct mapping_param_t* m,
               const double wcsin[8],
               /* Output parameters */
               double wcsout[8],
               struct driz_error_t* error);

/**
Update the WCS to include the drizzling transformations

This version is for BLOT and works the other way around to
\a update_wcs.

This is done by applying the transform to a unit square at
the reference pixel in the input whilst retaining the same
reference point on the sky.
*/
int
blot_update_wcs(struct driz_param_t* p,
                struct mapping_param_t* m,
                const double wcsin[8],
                /* Output parameters */
                double wcsout[8],
                struct driz_error_t* error);

#endif /* CDRIZZLEWCS_H */
