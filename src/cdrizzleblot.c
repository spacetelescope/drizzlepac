#define NO_IMPORT_ARRAY
#define NO_IMPORT_ASTROPY_WCS_API

#include "driz_portability.h"
#include "cdrizzlemap.h"
#include "cdrizzleblot.h"

#include <assert.h>
#define _USE_MATH_DEFINES       /* needed for MS Windows to define M_PI */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/**
Procedure to evaluate the bicubic polynomial interpolant.  The array
coeff contains the coefficients of the 2D interpolant.  The procedure
assumes that 1 <= x <= p->dnx and 1 <= y <= p->dny and that
coeff[1+first_point] = datain[1,1]. The interpolant is evaluated using
Everett's central difference formula.

@param[in] coeff Array of shape \a [len_coeff][len_coeff] contains the
coefficients of the 2D interpolant.

@param[in] len_coeff The dimension (each side of the square) of the
coefficient array.

@param[in] firstt Offset of the first data point.  (In practice, this
is always zero.  \todo Remove this parameter?)

@param[in] npts The number of points to calculate.

@param[in] x An array of length \a npts of x values.

@param[in] y An array of length \a npts of y values.

@param[out] zfit An array of length \a npts of interpolated values.

was: IIBIP3
*/
static inline_macro void
ii_bipoly3(const float* coeff /* [len_coeff][len_coeff] */,
           const integer_t len_coeff, const integer_t firstt,
           const integer_t npts,
           const float* x /* [npts] */, const float* y /* [npts] */,
           /* Output parameters */
           float* zfit /* [npts] */) {
  float sx, tx, sx2m1, tx2m1, sy, ty;
  float cd20[4] = {0.0, 0.0, 0.0, 0.0};
  float cd21[4] = {0.0, 0.0, 0.0, 0.0};
  float ztemp[4];
  float cd20y, cd21y;
  integer_t nxold, nyold;
  integer_t nx, ny;
  integer_t firstw, index;
  integer_t i, j;

  assert(coeff);
  assert(x);
  assert(y);
  assert(zfit);

  nxold = nyold = -1;
  for (i = 0; i < npts; ++i) {
    nx = (integer_t)x[i];
    assert(nx >= 0);

    sx = x[i] - (float)nx;
    tx = 1.0f - sx;
    sx2m1 = sx*sx - 1.0f;
    tx2m1 = tx*tx - 1.0f;

    ny = (integer_t)y[i];
    assert(ny >= 0);

    sy = y[i] - (float)ny;
    ty = 1.0f - sy;

    /* Calculate pointer to data[nx, ny-1] */
    firstw = firstt + (ny - 2) * len_coeff + nx - 1;

    /* loop over the 4 surrounding rows of data calculate the central
       differences at each value of y

       If new data point calculate the central differnences in x for
       each y */
    if (nx != nxold || ny != nyold) {
      for (j = 0, index = firstw; j < 4; ++j, index += len_coeff) {
        assert(index > 0 && index < (len_coeff*len_coeff) - 2);

        cd20[j] = 1.0f/6.0f * (coeff[index+1] -
                               2.0f * coeff[index] +
                               coeff[index-1]);
        cd21[j] = 1.0f/6.0f * (coeff[index+2] -
                               2.0f * coeff[index+1] +
                               coeff[index]);
      }
    }

    /* Interpolate in x at each value of y */
    for (j = 0, index = firstw; j < 4; ++j, index += len_coeff) {
      assert(index >= 0 && index < (len_coeff*len_coeff) - 1);

      ztemp[j] = sx * (coeff[index+1] + sx2m1 * cd21[j]) +
                 tx * (coeff[index] + tx2m1 * cd20[j]);
    }

    /* Calculate y central differences */
    cd20y = 1.0f/6.0f * (ztemp[2] - 2.0f * ztemp[1] + ztemp[0]);
    cd21y = 1.0f/6.0f * (ztemp[3] - 2.0f * ztemp[2] + ztemp[1]);

    /* Interpolate in y */
    zfit[i] = sy * (ztemp[2] + (sy * sy - 1.0f) * cd21y) +
              ty * (ztemp[1] + (ty * ty - 1.0f) * cd20y);

    nxold = nx;
    nyold = ny;
  }
}

/**
Procedure to evaluate a biquintic polynomial.  The array coeff
contains the coefficents of the 2D interpolant.  The routine assumes
that 0 <= x < p->dnx and 0 <= y < p->dny. The interpolant is evaluated
using Everett's central difference formula.

@param[in] coeff Array of shape \a [len_coeff][len_coeff] contains the
coefficients of the 2D interpolant.

@param[in] len_coeff The dimension (each side of the square) of the
coefficient array.

@param[in] firstt Offset of the first data point.  (In practice, this
is always zero.  \todo Remove this parameter?)

@param[in] npts The number of points to calculate.

@param[in] x An array of length \a npts of x values.

@param[in] y An array of length \a npts of y values.

@param[out] zfit An array of length \a npts of interpolated values.

was: IIBIP3
*/
static inline_macro void
ii_bipoly5(const float* coeff /* [len_coeff][len_coeff] */,
           const integer_t len_coeff, const integer_t firstt,
           const integer_t npts,
           const float* x /* [npts] */, const float* y /* [npts] */,
           /* Output parameters */
           float* zfit /* [npts] */) {
  integer_t nxold, nyold;
  integer_t nx, ny;
  float sx, sx2, tx, tx2, sy, sy2, ty, ty2;
  float sx2m1, sx2m4, tx2m1, tx2m4;
  float cd20[6], cd21[6], cd40[6], cd41[6];
  float cd20y, cd21y, cd40y, cd41y;
  float ztemp[6];
  integer_t firstw, index;
  integer_t i, j;

  assert(coeff);
  assert(len_coeff > 0);
  assert(npts > 0);
  assert(x);
  assert(y);
  assert(zfit);

  nxold = nyold = -1;
  for (i = 0; i < npts; ++i) {
    nx = (integer_t)x[i];
    ny = (integer_t)y[i];
    assert(nx >= 0);
    assert(ny >= 0);

    sx = x[i] - (float)nx;
    sx2 = sx * sx;
    sx2m1 = sx2 - 1.0f;
    sx2m4 = sx2 - 4.0f;
    tx = 1.0f - sx;
    tx2 = tx * tx;
    tx2m1 = tx2 - 1.0f;
    tx2m4 = tx2 - 4.0f;

    sy = y[i] - (float)ny;
    sy2 = sy * sy;
    ty = 1.0f - sy;
    ty2 = ty * ty;

    /* Calculate value of pointer to data[nx,ny-2] */
    firstw = firstt + (ny - 3)*len_coeff + nx - 1;

    /* Calculate the central differences in x at each value of y */
    if (nx != nxold || ny != nyold) {
      for (j = 0, index = firstw; j < 6; ++j, index += len_coeff) {
        assert(index >= 1 && index < len_coeff*len_coeff - 3);

        cd20[j] = 1.0f/6.0f * (coeff[index+1] -
                               2.0f * coeff[index] +
                               coeff[index-1]);
        cd21[j] = 1.0f/6.0f * (coeff[index+2] -
                               2.0f * coeff[index+1] +
                               coeff[index]);
        cd40[j] = 1.0f/120.0f * (coeff[index-2] -
                                 4.0f * coeff[index-1] +
                                 6.0f * coeff[index] -
                                 4.0f * coeff[index+1] +
                                 coeff[index+2]);
        cd41[j] = 1.0f/120.0f * (coeff[index-1] -
                                 4.0f * coeff[index] +
                                 6.0f * coeff[index+1] -
                                 4.0f * coeff[index+2] +
                                 coeff[index+3]);
      }
    }

    /* Interpolate in x at each value of y */
    for (j = 0, index = firstw; j < 6; ++j, index += len_coeff) {
      assert(index >= 0 && index < len_coeff*len_coeff - 1);

      ztemp[j] = sx * (coeff[index+1] + sx2m1 * (cd21[j] + sx2m4 * cd41[j])) +
        tx * (coeff[index]   + tx2m1 * (cd20[j] + tx2m4 * cd40[j]));
    }

    /* Central differences in y */
    cd20y = 1.0f/6.0f * (ztemp[3] - 2.0f * ztemp[2] + ztemp[1]);
    cd21y = 1.0f/6.0f * (ztemp[4] - 2.0f * ztemp[3] + ztemp[2]);
    cd40y = 1.0f/120.0f * (ztemp[0] -
                           4.0f * ztemp[1] +
                           6.0f * ztemp[2] -
                           4.0f * ztemp[3] +
                           ztemp[4]);
    cd41y = 1.0f/120.0f * (ztemp[1] -
                           4.0f * ztemp[2] +
                           6.0f * ztemp[3] -
                           4.0f * ztemp[4] +
                           ztemp[5]);

    /* Interpolate in y */
    zfit[i] = sy * (ztemp[3] + (sy2 - 1.0f) * (cd21y + (sy2 - 4.0f) * cd41y)) +
      ty * (ztemp[2] + (ty2 - 1.0f) * (cd20y + (ty2 - 4.0f) * cd40y));

    nxold = nx;
    nyold = ny;
  }
}

/** Helper function to look up pixels in a 2D array */
static inline_macro float
data_value(const float* data, const integer_t dnx, const integer_t dny UNUSED_PARAM,
           const integer_t x, const integer_t y) {
  assert(data);
  assert(x >= 0 && x < dnx);
  assert(y >= 0 && y < dny);

  return data[y*dnx + x];
}

/** Macro to pass pointer and dimensions to the data_value function */
#define DATA_VALUE(x, y) (data_value(data, dnx, dny, x, y))

/**
Signature for functions that perform blotting interpolation.
 */
typedef int (interp_function)(const void*,
                              const float*,
                              const integer_t, const integer_t,
                              const float, const float,
                              /* Output parameters */
                              float*,
                              struct driz_error_t*);

/**
A standard set of asserts for all of the interpolation functions
*/
#define INTERPOLATION_ASSERTS \
  assert(data); \
  assert(dnx > 0); \
  assert(dny > 0); \
  assert(x >= 0.0f && x < (float)dnx);      \
  assert(y >= 0.0f && y < (float)dny);      \
  assert(value); \
  assert(error); \

/**
Perform nearest neighbor interpolation.

@param[in] state A pointer to any constant values specific to this
interpolation type.  (For \a interpolate_nearest_neighbor, it should be
NULL).

@param[in] data A 2D data array of shape [dny][dnx]

@param[in] dnx The x dimension of data

@param[in] dny The y dimension of data

@param[in] x The fractional x coordinate

@param[in] y The fractional y coordinate

@param[out] value The resulting value at x, y after interpolating the data

@param[out] error

@return Non-zero if an error occurred
 */
static int
interpolate_nearest_neighbor(const void* state UNUSED_PARAM,
                             const float* data,
                             const integer_t dnx, const integer_t dny,
                             const float x, const float y,
                             /* Output parameters */
                             float* value,
                             struct driz_error_t* error UNUSED_PARAM) {
  assert(state == NULL);
  INTERPOLATION_ASSERTS;

  *value = DATA_VALUE((integer_t)(x + 0.5), (integer_t)(y + 0.5));
  return 0;
}

/**
Perform basic bilinear interpolation.

@param[in] state A pointer to any constant values specific to this
interpolation type.  (For \a interpolate_bilinear, it should be
NULL).

@param[in] data A 2D data array of shape [dny][dnx]

@param[in] dnx The x dimension of data

@param[in] dny The y dimension of data

@param[in] x The fractional x coordinate

@param[in] y The fractional y coordinate

@param[out] value The resulting value at x, y after interpolating the data

@param[out] error

@return Non-zero if an error occurred
 */
static int
interpolate_bilinear(const void* state UNUSED_PARAM,
                     const float* data,
                     const integer_t dnx, const integer_t dny,
                     const float x, const float y,
                     /* Output parameters */
                     float* value,
                     struct driz_error_t* error UNUSED_PARAM) {
  integer_t nx, ny;
  float sx, tx, sy, ty, f00;

  assert(state == NULL);
  INTERPOLATION_ASSERTS;

  nx = (integer_t) x;
  ny = (integer_t) y;

  if (nx < 0 || ny < 0 || nx >= dnx || ny >= dny) {
      driz_error_set_message(error,
          "Bilinear interpolation: point outside of the image.");
      return 1;
  }

  f00 = DATA_VALUE(nx, ny);

  if (nx == (dnx - 1)) {
    if (ny == (dny - 1)) {
      /* This is the last pixel (in x and y). Assign constant value of this pixel. */
      *value = f00;
      return 0;
    }
    /* Interpolate along Y-direction only */
    sy = y - (float)ny;
    *value = (1.0f - sy) * f00 + sy * DATA_VALUE(nx, ny + 1);
  } else if (ny == (dny - 1)) {
    /* Interpolate along X-direction only */
    sx = x - (float)nx;
    *value = (1.0f - sx) * f00 + sx * DATA_VALUE(nx + 1, ny);
  } else {
    /* Bilinear - interpolation */
    sx = x - (float)nx;
    tx = 1.0f - sx;
    sy = y - (float)ny;
    ty = 1.0f - sy;

    *value = tx * ty * f00 +
             sx * ty * DATA_VALUE(nx + 1, ny) +
             sy * tx * DATA_VALUE(nx, ny + 1) +
             sx * sy * DATA_VALUE(nx + 1, ny + 1);
  }

  return 0;
}

/**
Perform cubic polynomial interpolation.

@param[in] state A pointer to any constant values specific to this
interpolation type.  (For \a interpolate_poly3, it should be
NULL).

@param[in] data A 2D data array of shape [dny][dnx]

@param[in] dnx The x dimension of data

@param[in] dny The y dimension of data

@param[in] x The fractional x coordinate

@param[in] y The fractional y coordinate

@param[out] value The resulting value at x, y after interpolating the data

@param[out] error

@return Non-zero if an error occurred
 */
static int
interpolate_poly3(const void* state UNUSED_PARAM,
                  const float* data,
                  const integer_t dnx, const integer_t dny,
                  const float x, const float y,
                  /* Output parameters */
                  float* value,
                  struct driz_error_t* error UNUSED_PARAM) {
  integer_t nx, ny;
  const integer_t rowleh = 4;
  const integer_t nterms = 4;
  float coeff[4][4];
  integer_t i, j;
  integer_t firstw, lastrw;
  float xval, yval;
  float* ci;

  assert(state == NULL);
  INTERPOLATION_ASSERTS;;

  nx = (integer_t)x;
  ny = (integer_t)y;

  ci = &coeff[0][0];
  for (j = ny - 1; j <= ny + 2; ++j) {
    if (j >= 0 && j < dny) {
      for (i = nx - 1; i <= nx + 2; ++i, ++ci) {
        if (i < 0) {
          *ci = 2.0f * DATA_VALUE(0, j) - DATA_VALUE(-i, j);
        } else if (i >= dnx) {
          *ci = 2.0f * DATA_VALUE(dnx - 1, j) - DATA_VALUE(2*dnx - 2 - i, j);
        } else {
          *ci = DATA_VALUE(i, j);
        }
      }
    } else if (j == ny + 2) {
      for (i = nx - 1; i <= nx + 2; ++i, ++ci) {
        if (i < 0) {
          *ci = 2.0f * DATA_VALUE(0, dny - 3) - DATA_VALUE(-i, dny - 3);
        } else if (i >= dnx) {
          *ci = 2.0f * DATA_VALUE(dnx - 1, dny - 3) - DATA_VALUE(2*dnx - 2 - i, dny - 3);
        } else {
          *ci = DATA_VALUE(i, dny - 3);
        }
      }
    } else {
      ci += 4;
    }
  }

  firstw = MAX(0, 1 - ny);
  if (firstw > 0) {
    assert(firstw < nterms);

    for (j = 0; j < firstw; ++j) {
      assert(2*firstw-j >= 0 && 2*firstw-j < nterms);

      weighted_sum_vectors(nterms,
                           &coeff[firstw][0], 2.0,
                           &coeff[2*firstw-j][0], -1.0,
                           &coeff[j][0]);
    }
  }

  lastrw = MIN(nterms - 1, dny - ny);
  if (lastrw < nterms - 1) {
    assert(lastrw >= 0 && lastrw < nterms);

    for (j = lastrw + 1; j <= nterms - 1; ++j) {
      assert(2*lastrw-j >= 0 && 2*lastrw-j < nterms);
      assert(j >= 0 && j < 4);

      weighted_sum_vectors(nterms,
                           &coeff[lastrw][0], 2.0,
                           &coeff[2*lastrw-j][0], -1.0,
                           &coeff[j][0]);
    }
  } else if (lastrw == 1) {
    assert(lastrw >= 0 && lastrw < nterms);

    weighted_sum_vectors(nterms,
                         &coeff[lastrw][0], 2.0,
                         &coeff[3][0], -1.0,
                         &coeff[3][0]);
  } else {
    assert(lastrw >= 0 && lastrw < nterms);
    assert(2*lastrw-3 >= 0 && 2*lastrw-3 < nterms);

    weighted_sum_vectors(nterms,
                         &coeff[lastrw][0], 2.0,
                         &coeff[2*lastrw-3][0], -1.0,
                         &coeff[3][0]);
  }

  xval = 2.0f + (x - (float)nx);
  yval = 2.0f + (y - (float)ny);

  ii_bipoly3(&coeff[0][0], rowleh, 0, 1, &xval, &yval, value);

  return 0;
}

/**
Perform quintic polynomial interpolation.

@param[in] state A pointer to any constant values specific to this
interpolation type.  (For \a interpolate_poly5, it should be
NULL).

@param[in] data A 2D data array of shape [dny][dnx]

@param[in] dnx The x dimension of data

@param[in] dny The y dimension of data

@param[in] x The fractional x coordinate

@param[in] y The fractional y coordinate

@param[out] value The resulting value at x, y after interpolating the data

@param[out] error

@return Non-zero if an error occurred
 */
static int
interpolate_poly5(const void* state UNUSED_PARAM,
                  const float* data,
                  const integer_t dnx, const integer_t dny,
                  const float x, const float y,
                  /* Output parameters */
                  float* value,
                  struct driz_error_t* error UNUSED_PARAM) {
  integer_t nx, ny;
  const integer_t rowleh = 6;
  const integer_t nterms = 6;
  float coeff[6][6];
  integer_t i, j;
  integer_t firstw, lastrw;
  float xval, yval;
  float* ci;

  assert(state == NULL);
  INTERPOLATION_ASSERTS;

  nx = (integer_t)x;
  ny = (integer_t)y;

  ci = &coeff[0][0];
  for (j = ny - 2; j <= ny + 3; ++j) {
    if (j >= 0 && j < dny) {
      for (i = nx - 2; i <= nx + 3; ++i, ++ci) {
        if (i < 0) {
          *ci = 2.0f * DATA_VALUE(0, j) - DATA_VALUE(-i, j);
        } else if (i >= dnx) {
          *ci = 2.0f * DATA_VALUE(dnx - 1, j) - DATA_VALUE(2*dnx - 2 - i, j);
        } else {
          *ci = DATA_VALUE(i, j);
        }
      }
    } else if (j == (ny + 3)) {
      for (i = nx - 2; i <= nx + 3; ++i, ++ci) {
        if (i < 0) {
          *ci = 2.0f * DATA_VALUE(0, dny - 4) - DATA_VALUE(-i, dny - 4);
        } else if (i >= dnx) {
          *ci = 2.0f * DATA_VALUE(dnx - 1, dny - 4) - DATA_VALUE(2 * dnx - 2 - i, dny - 4);
        } else {
          *ci = DATA_VALUE(i, dny - 4);
        }
      }
    } else {
      ci += 6;
    }
  }

  firstw = MAX(0, 2 - ny);
  assert(firstw >= 0 && firstw < nterms);

  if (firstw > 0) {
    for (j = 0; j <= firstw; ++j) {
      assert(2*firstw-j >= 0 && 2*firstw-j < nterms);

      weighted_sum_vectors(nterms,
                           &coeff[firstw][0], 2.0,
                           &coeff[2*firstw-j][0], -1.0,
                           &coeff[j][0]);
    }
  }

  lastrw = MIN(nterms - 1, dny - ny + 1);
  assert(lastrw < nterms);

  if (lastrw < nterms - 1) {
    for (j = lastrw + 1; j <= nterms - 2; ++j) {
      assert(2*lastrw-j >= 0 && 2*lastrw-j < nterms);

      weighted_sum_vectors(nterms,
                           &coeff[lastrw][0], 2.0,
                           &coeff[2*lastrw-j][0], -1.0,
                           &coeff[j][0]);
    }
  } else if (lastrw == 2) {
    weighted_sum_vectors(nterms,
                         &coeff[2][0], 2.0,
                         &coeff[5][0], -1.0,
                         &coeff[5][0]);
  } else {
    assert(2*lastrw - 5 >= 0 && 2*lastrw-5 < nterms);

    weighted_sum_vectors(nterms,
                         &coeff[lastrw][0], 2.0,
                         &coeff[2*lastrw-5][0], -1.0,
                         &coeff[5][0]);
  }

  xval = 3.0f + (x - (float)nx);
  yval = 3.0f + (y - (float)ny);

  ii_bipoly5(&coeff[0][0], rowleh, 0, 1, &xval, &yval, value);

  return 0;
}

/**
A structure to hold parameters for sinc interpolation.
*/
struct sinc_param_t {
  /** The scaling factor for sinc interpolation */
  float sinscl;
};

/**
was: iinisc
*/
#define INTERPOLATE_SINC_NCONV 15

static inline_macro int
interpolate_sinc_(const float* data, const integer_t firstt,
                  const integer_t npts,
                  const float* x /*[npts]*/, const float* y /*[npts]*/,
                  const integer_t len_coeff,
                  const integer_t lenary, const float mindx,
                  const float mindy, const float sinscl,
                  /* Output parameters */
                  float* value,
                  struct driz_error_t* error UNUSED_PARAM) {
  const integer_t nconv = INTERPOLATE_SINC_NCONV;
  const integer_t nsinc = (nconv - 1) / 2;
  /* TODO: This is to match Fortan, but is probably technically less precise */
  const float halfpi = 1.5707963267948966192f; /* M_PI / 2.0; */
  const float sconst = powf((halfpi / (float)nsinc), 2.0f);
  const float a2 = -0.49670f;
  const float a4 = 0.03705f;
  float taper[INTERPOLATE_SINC_NCONV];
  float ac[INTERPOLATE_SINC_NCONV], ar[INTERPOLATE_SINC_NCONV];
  float sdx, dx, dy, dxn, dyn, dx2;
  float ax, ay, px, py;
  float sum, sumx, sumy;
  float tmp;
  integer_t minj, maxj, offj;
  integer_t mink, maxk, offk;
  integer_t nx, ny;
  integer_t i, j, k, m, index;
  integer_t indices[3][4];

  assert(x);
  assert(y);
  assert(value);
  assert(error);

  if ((nsinc % 2) == 0) {
    sdx = 1.0;
    for (j = -nsinc; j <= nsinc; ++j) {
      assert(j + nsinc >= 0 && j + nsinc < INTERPOLATE_SINC_NCONV);

      taper[j + nsinc] = 1.0;
    }
  } else {
    sdx = -1.0;
    errno = 0;
    for (j = -nsinc; j <= nsinc; ++j) {
      assert(j + nsinc >= 0 && j + nsinc < INTERPOLATE_SINC_NCONV);

      dx2 = sconst * (float)j * (float)j;
      tmp = powf(1.0f + a2*dx2 + a4*dx2*dx2, 2.0);
      if (errno != 0) {
        driz_error_set_message(error, "pow failed");
        return 1;
      }
      taper[j + nsinc] = sdx * tmp;

      sdx = -sdx;
    }
  }

  for (i = 0; i < npts; ++i) {
    nx = fortran_round(x[i]);
    ny = fortran_round(y[i]);
    if (nx < 0 || nx >= len_coeff || ny < 0 || ny >= lenary) {
      value[i] = 0.0;
      continue;
    }

    dx = (x[i] - (float)nx) * sinscl;
    dy = (y[i] - (float)ny) * sinscl;

    if (fabsf(dx) < mindx && fabsf(dy) < mindy) {
      index = firstt + (ny - 1) * len_coeff + nx - 1; /* TODO: Base check */
      value[i] = data[index];
      continue;
    }

    dxn = 1.0f + (float)nsinc + dx;
    dyn = 1.0f + (float)nsinc + dy;
    sumx = 0.0f;
    sumy = 0.0f;
    for (j = 0; j < nconv; ++j) {
      /* TODO: These out of range indices also seem to be in Fortran... */
      ax = dxn - (float)j - 1; /* TODO: Base check */
      ay = dyn - (float)j - 1; /* TODO: Base check */
      assert(ax != 0.0);
      assert(ay != 0.0);

      if (ax == 0.0) {
        px = 1.0;
      } else if (dx == 0.0) {
        px = 0.0;
      } else {
        px = taper[j - 1] / ax;
      }

      if (ay == 0.0) {
        py = 1.0;
      } else if (dy == 0.0) {
        py = 0.0;
      } else {
        py = taper[j - 1] / ay;
      }

      /* TODO: These out of range indices also seem to be in Fortran... */
      ac[j - 1] = px;
      ar[j - 1] = py;
      sumx += px;
      sumy += py;
    }

    /* Compute the limits of the convolution */
    minj = MAX(0, ny - nsinc - 1); /* TODO: Bases check */
    maxj = MIN(lenary, ny + nsinc); /* TODO: Bases check */
    offj = nsinc - ny; /* TODO: Bases check */

    mink = MAX(0, nx - nsinc - 1); /* TODO: Bases check */
    maxk = MIN(len_coeff, nx + nsinc); /* TODO: Bases check */
    offk = nsinc - nx; /* TODO: Bases check */

    value[i] = 0.0;

    indices[0][0] = ny - nsinc;
    indices[0][1] = minj - 1;
    indices[0][2] = firstt;
    indices[0][3] = 0;

    indices[1][0] = minj;
    indices[1][1] = maxj;
    indices[1][2] = firstt;
    indices[1][3] = len_coeff;

    indices[2][0] = maxj + 1;
    indices[2][1] = ny + nsinc;
    indices[2][2] = firstt + (lenary - 1) * len_coeff;
    indices[2][3] = 0;

    for (m = 0; m < 3; ++m) {
      for (j = indices[m][0]; j <= indices[m][1]; ++j) {
        sum = 0.0;
        index = indices[m][2] + j * indices[m][3];
        assert(index >= 0 && index < len_coeff*lenary - 1);
        assert(index+len_coeff >= 0 && index+len_coeff < len_coeff*lenary);

        for (k = nx - nsinc; k < mink - 1; ++k) { /* TODO: Bases check */
          assert(k+offk >= 0 && k+offk < INTERPOLATE_SINC_NCONV);

          sum += ac[k+offk] * data[index+1];
        }

        for (k = mink; k <= maxk; ++k) { /* TODO: Bases check */
          assert(k+offk >= 0 && k+offk < INTERPOLATE_SINC_NCONV);
          assert(index+k >= 0 && index+k < len_coeff*lenary);

          sum += ac[k+offk] * data[index+k];
        }

        for (k = maxk + 1; k <= nx + nsinc; ++k) {
          assert(k+offk >= 0 && k+offk < INTERPOLATE_SINC_NCONV);

          sum += ac[k+offk] * data[index+len_coeff];
        }

        assert(j + offj >= 0 && j + offj < INTERPOLATE_SINC_NCONV);

        value[i] += ar[j + offj] * sum;
      }
    }

    assert(sumx != 0.0);
    assert(sumy != 0.0);

    value[i] = value[i] / sumx / sumy;
  }

  return 0;
}

/**
Perform sinc interpolation.

@param[in] state A pointer to any constant values specific to this
interpolation type.  (For \a interpolate_sinc, it must be a pointer
to a \a sinc_param_t object).

@param[in] data A 2D data array of shape [dny][dnx]

@param[in] dnx The x dimension of data

@param[in] dny The y dimension of data

@param[in] x The fractional x coordinate

@param[in] y The fractional y coordinate

@param[out] value The resulting value at x, y after interpolating the data

@param[out] error

@return Non-zero if an error occurred
*/
static int
interpolate_sinc(const void* state,
                 const float* data,
                 const integer_t dnx, const integer_t dny,
                 const float x, const float y,
                 /* Output parameters */
                 float* value,
                 struct driz_error_t* error) {
  const struct sinc_param_t* param = (const struct sinc_param_t*)state;

  assert(state);
  INTERPOLATION_ASSERTS;

  return interpolate_sinc_(data, 0, 1, &x, &y, dnx, dny,
                           0.001f, 0.001f, param->sinscl, value, error);
}

/**
Perform Lanczos interpolation.

@param[in] state A pointer to any constant values specific to this
interpolation type.  (For \a interpolate_lanczos, it must be a pointer
to a \a lanczos_param_t object, already fully filled-in).

@param[in] data A 2D data array of shape [dny][dnx]

@param[in] dnx The x dimension of data

@param[in] dny The y dimension of data

@param[in] x The fractional x coordinate

@param[in] y The fractional y coordinate

@param[out] value The resulting value at x, y after interpolating the data

@param[out] error

@return Non-zero if an error occurred
*/
static int
interpolate_lanczos(const void* state,
                    const float* data,
                    const integer_t dnx, const integer_t dny,
                    const float x, const float y,
                    /* Output parameters */
                    float* value,
                    struct driz_error_t* error UNUSED_PARAM) {
  integer_t ixs, iys, ixe, iye;
  integer_t xoff, yoff;
  float luty, sum;
  integer_t nbox;
  integer_t i, j;
  const struct lanczos_param_t* p = (const struct lanczos_param_t*)state;

  assert(state);
  INTERPOLATION_ASSERTS;

  nbox = p->nbox;

  /* First check for being close to the edge and, if so, return the
     missing value */
  ixs = (integer_t)(x) - nbox;
  ixe = (integer_t)(x) + nbox;
  iys = (integer_t)(y) - nbox;
  iye = (integer_t)(y) + nbox;
  if (ixs < 0 || ixe >= dnx ||
      iys < 0 || iye >= dny) {
    *value = p->misval;
    return 0;
  }

  /* Don't divide-by-zero errors */
  assert(p->space != 0.0);

  /* Loop over the box, which is assumed to be scaled appropriately */
  sum = 0.0;
  for (j = iys; j <= iye; ++j) {
    yoff = (integer_t)(fabs((y - (float)j) / p->space));
    assert(yoff >= 0 && yoff < p->nlut);

    luty = p->lut[yoff];
    for (i = ixs; i <= ixe; ++i) {
      xoff = (integer_t)(fabs((x - (float)i) / p->space));
      assert(xoff >= 0 && xoff < p->nlut);

      sum += DATA_VALUE(i, j) * p->lut[xoff] * luty;
    }
  }

  *value = sum;
  return 0;
}

/**
A mapping from e_interp_t enumeration values to function pointers that actually
perform the interpolation.  NULL elements will raise an "unimplemented" error.
*/
interp_function* interp_function_map[interp_LAST] = {
  &interpolate_nearest_neighbor,
  &interpolate_bilinear,
  &interpolate_poly3,
  &interpolate_poly5,
  NULL,
  &interpolate_sinc,
  &interpolate_sinc,
  &interpolate_lanczos,
  &interpolate_lanczos
};

/* See header file for documentation */
int
doblot(struct driz_param_t* p,
       struct driz_error_t* error) {
  const size_t nlut = 2048;
  const float space = 0.01;
  double *xin = NULL;
  double *xtmp = NULL;
  double *xout = NULL;
  double *yin = NULL;
  double *ytmp = NULL;
  double *yout = NULL;
  integer_t nmiss;
  double dx, dy;
  double yv;
  float xo, yo, v;
  /*float nx, ny;*/
  integer_t i, j;
  interp_function* interpolate;
  struct sinc_param_t sinc;
  void* state = NULL;

  assert(p);
  assert(error);
  assert(space != 0.0);

  /* Some initial settings */
  nmiss = 0;

  /* Select interpolation function */
  assert(p->interpolation >= 0 && p->interpolation < interp_LAST);
  interpolate = interp_function_map[p->interpolation];
  if (interpolate == NULL) {
    driz_error_set_message(error, "Requested interpolation type not implemented.");
    goto doblot_exit_;
  }

  /* Some interpolation functions need some pre-calculated state */
  if (p->interpolation == interp_lanczos3 || p->interpolation == interp_lanczos5) {
    assert(p->kscale != 0.0);
    assert(p->lanczos.lut == NULL);
    if ((p->lanczos.lut = (float*)malloc(nlut * sizeof(float))) == NULL) {
      driz_error_set_message(error, "Out of memory");
      goto doblot_exit_;
    }
    create_lanczos_lut(p->interpolation == interp_lanczos3 ? 3 : 5,
                       nlut, space, p->lanczos.lut);
    p->lanczos.nbox = (integer_t)(3.0 / p->kscale);
    p->kscale2 = 1.0f / (p->kscale * p->kscale);
    p->lanczos.nlut = nlut;
    p->lanczos.space = space;
    p->lanczos.misval = p->misval;
    state = &(p->lanczos);
  } else if (p->interpolation == interp_sinc || p->interpolation == interp_lsinc) {
    sinc.sinscl = p->sinscl;
    state = &sinc;
  } /* Otherwise state is NULL */

  /* Allocate some memory */
  assert(p->onx >= 0);
  assert(p->ony >= 0);

  xin = malloc((size_t)p->onx * sizeof(double));
  if (xin == NULL) {
      driz_error_set_message(error, "Out of memory");
      goto doblot_exit_;
  }

  xtmp = malloc((size_t)p->onx * sizeof(double));
  if (xtmp == NULL) {
      driz_error_set_message(error, "Out of memory");
      goto doblot_exit_;
  }

  xout = malloc((size_t)p->onx * sizeof(double));
  if (xout == NULL) {
      driz_error_set_message(error, "Out of memory");
      goto doblot_exit_;
  }

  yin = malloc((size_t)p->onx * sizeof(double));
  if (yin == NULL) {
      driz_error_set_message(error, "Out of memory");
      goto doblot_exit_;
  }

  ytmp = malloc((size_t)p->onx * sizeof(double));
  if (ytmp == NULL) {
      driz_error_set_message(error, "Out of memory");
      goto doblot_exit_;
  }

  yout = malloc((size_t)p->onx * sizeof(double));
  if (yout == NULL) {
      driz_error_set_message(error, "Out of memory");
      goto doblot_exit_;
  }

  /* In the WCS case, we can't use the scale to calculate the Jacobian,
     so we need to do it.

     Note that we use the center of the image, rather than the reference pixel
     as the reference here.

     This is taken from dobox, except for the inversion of the image order.

     This section applies in WBLOT mode and now contains the addition
     correction to separate the distortion-induced scale change.
  */

  /* Image subset size
  nx = (float)(p->xmax - p->xmin + 1);
  ny = (float)(p->ymax - p->ymin + 1);
  */

  /* Offsets */
  dx = (double)(p->xmin);
  dy = (double)(p->ymin);

  /* Recalculate the area scaling factor */
  assert(p->scale != 0.0);
  p->scale2 = p->scale*p->scale;

  /* Set the X and Y start positions -- most of these don't change
     between iterations. */
  xin[0] = 1.0;
  xin[1] = 0.0;
  yin[1] = 0.0;
  v = 1.0;

  /* Outer look over output image pixels (X, Y) */
  for (j = 0; j < p->ony; ++j) {
    yv = (double)j+1;

    yin[0] = yv;

    /* Transform this vector */
    if (map_value(p, TRUE, p->onx,
                  xin, yin, xtmp, ytmp, xout, yout, error)) {
      goto doblot_exit_;
    }

    /* Loop through the output positions and do the interpolation */
    for (i = 0; i < p->onx; ++i) {
      xo = (float)(xout[i] - dx);
      yo = (float)(yout[i] - dy);

      /* Check it is on the input image */
      if (xo >= 0.0 && xo <= p->dnx &&
          yo >= 0.0 && yo <= p->dny) {

        /* Check for look-up-table interpolation */
        if (interpolate(state, p->data, p->dnx, p->dny, xo, yo, &v, error)) {
          goto doblot_exit_;
        }

        /* TODO: This float cast makes it match Fortran, but technically
           loses more precision */
        *output_data_ptr(p, i, j) = v * p->ef / (float)p->scale2;
      } else {
        /* If there is nothing for us then set the output to missing C
           value flag */
        *output_data_ptr(p, i, j) = p->misval;

        nmiss++;
      }
    }
  }

  /* if (!p->use_wcs) { */
  /*   if (blot_update_wcs(p, m, error)) */
  /*     return 1; */
  /* } */

 doblot_exit_:
  free(p->lanczos.lut); p->lanczos.lut = NULL;
  free(xin); xin = NULL;
  free(xtmp); xtmp = NULL;
  free(xout); xout = NULL;
  free(yin); yin = NULL;
  free(ytmp); ytmp = NULL;
  free(yout); yout = NULL;

  return driz_error_is_set(error);
}
