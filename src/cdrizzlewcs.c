#define NO_IMPORT_ARRAY
#define NO_IMPORT_ASTROPY_WCS_API
#include "driz_portability.h"
#include "cdrizzlewcs.h"
#include "cdrizzlemap.h"

static int
invert_matrix_2x2(const double matin[4],
                  /* Output parameters */
                  double matout[4],
                  struct driz_error_t* error) {
  const double det = matin[0]*matin[1] - matin[2]*matin[3];

  /* Check the determinant for singularity */
  if (det == 0.0) {
    driz_error_set_message(error, "Matrix is singular, cannot update WCS");
    return 1;
  }

  matout[0] = matin[3] / det;
  matout[1] = -matin[1] / det;
  matout[2] = -matin[2] / det;
  matout[3] = matin[0] / det;

  return 0;
}

/**
Convert a pixel position to an equatorial (RA,DEC) position assuming a
tangent projection and WCS with the standard 8 elements.

@param wcs An array of length 8 of WCS parameters

@param x The input x coordinate

@param y The input y coordinate

@param[out] ra The output ra coordinate

@param[out] dec The output dec coordinate

@param[out] error
*/
static void
xy2rd(const double wcs[8], const double x, const double y,
      /* Output parameters */
      double* ra, double* dec) {
  const double pi_by = M_PI / 180.0;
  double xi, eta, ra0, dec0;

  assert(wcs);
  assert(ra);
  assert(dec);

  /* First convert pixel coordinates to tangent plane in radians */
  xi = (wcs[4]*(x - wcs[0]) + wcs[6]*(y - wcs[2])) * pi_by;
  eta = (wcs[5]*(x - wcs[0]) + wcs[7]*(y - wcs[2])) * pi_by;

  /* And convert the reference point on the sky to radians */
  ra0 = wcs[1] * pi_by;
  dec0 = wcs[3] * pi_by;

  /* Now go to equatorial from tangent plane */
  *ra = (atan2(xi, cos(dec0) - eta*sin(dec0)) + ra0) / pi_by;
  *dec = (atan2(eta*cos(dec0) + sin(dec0),
                sqrt(pow(cos(dec0) - eta*sin(dec0), 2) +
                     xi*xi))) / pi_by;

  /* Convert back to degrees and check the range */
  if (*ra < 0.0) {
    *ra += 360.0;
  }
}

/**
Convert an equatorial (ra, dec) position to a pixel position assuming
a tangent projection and WCS with the standard 8 elements.

@param wcs An array of length 8 of WCS parameters

@param r The input ra coordinate

@param d The input dec coordinate

@param[out] x The output x coordinate

@param[out] y The output y coordinate

@param[out] error

@return Non-zero if an error occurred, particularly if the wcs matrix
is non-invertible.
*/
static int
rd2xy(const double wcs[8],
      const double r, const double d,
      /* Output parameters */
      double* x, double* y,
      struct driz_error_t* error) {
  const double pi_by = M_PI / 180.0; /* TODO: This is more accurate than the Fortran */
  double ra, dec, xi, eta, ra0, dec0;
  double det, cdinv[2][2], bottom;

  assert(wcs);
  assert(x);
  assert(y);
  assert(error);

  /* First, invert the CD matrix */
  det = wcs[4] * wcs[7] - wcs[6] * wcs[5];

  if (det == 0.0) {
    driz_error_set_message(error, "Non-invertible matrix");
    return 1;
  }

  cdinv[0][0] = wcs[7] / det;
  cdinv[0][1] = -wcs[6] / det;
  cdinv[1][0] = -wcs[5] / det;
  cdinv[1][1] = -wcs[4] / det;

  /* Translate from RA, DEC to X, Y */
  ra0 = wcs[1] * pi_by;
  dec0 = wcs[3] * pi_by;

  ra = r * pi_by;
  dec = d * pi_by;

  bottom = sin(dec)*sin(dec0) + cos(dec)*cos(dec0)*cos(ra-ra0);
  if (bottom == 0.0) {
    driz_error_set_message(error, "Non-invertible matrix");
    return 1;
  }

  /* Calculate tangent plane position and convert to degrees */
  xi = cos(dec)*sin(ra-ra0) / bottom / pi_by;
  eta = (sin(dec)*cos(dec0) - cos(dec)*sin(dec0)*cos(ra-ra0)) / bottom / pi_by;

  /* Convert back to pixels using the inverse of the CD matrix */
  *x = cdinv[0][0]*xi + cdinv[0][1]*eta + wcs[0];
  *y = cdinv[1][0]*xi + cdinv[1][1]*eta + wcs[2];

  return 0;
}

#define MAX_MATRIX_ORDER 10

static inline_macro double*
mat_ptr(double* mat, const integer_t norder,
        const integer_t i, const integer_t j) {
  assert(mat);
  assert(i >= 0 && i < norder);
  assert(j >= 0 && j < norder);
  return (mat + (j * norder + i));
}

#define MAT_PTR(i, j) (mat_ptr(mat, norder, i, j))

static void
invert_matrix_accumulate(double* mat,
                         const integer_t norder,
                         const integer_t k,
                         const double amax,
                         /* Input/Output parameters */
                         double* det) {
  integer_t i, j;

  assert(mat);
  assert(det);
  assert(norder <= MAX_MATRIX_ORDER);
  assert(k < norder);
  assert(amax != 0.0);

  /* Accumulate elements of inverse matrix */
  for (i = 0; i < norder; ++i) {
    if (i - k != 0) {
      *MAT_PTR(i, k) /= -amax;
    }
  }

  for (i = 0; i < norder; ++i) {
    for (j = 0; j < norder; ++j) {
      /* TODO: Why is this not i != k && j != k ? */
      if (i - k != 0 && j - k != 0) {
        *MAT_PTR(i, j) = *MAT_PTR(i, j) + *MAT_PTR(i, k) * *MAT_PTR(k, j);
      }
    }
  }

  for (j = 0; j < norder; ++j) {
    if (j - k != 0) {
      *MAT_PTR(k, j) /= amax;
    }
  }

  *MAT_PTR(k, k) = 1.0 / amax;
  *det *= amax;
}

/**
Invert a symmetric matrix and calculate its determinant.

was: MATINV
*/
static int
invert_matrix(double* mat /* Input/output */,
              const integer_t norder,
              /* Output parameters */
              double* det,
              struct driz_error_t* error UNUSED_PARAM) {
  integer_t i, j, k, l;
  double amax, save;
  integer_t ik[MAX_MATRIX_ORDER], jk[MAX_MATRIX_ORDER];

  assert(mat);
  assert(det);
  assert(error);
  assert(norder <= MAX_MATRIX_ORDER);

  for (i = 0; i < norder; ++i) {
    ik[i] = 0;
    jk[i] = 0;
  }

  *det = 1.0;

  amax = 0.0;
  /* Find largest element array in rest of matrix */
  for (k = 0; k < norder; ++k) {
    for (i = k; i < norder; ++i) {
      for (j = k; j < norder; ++j) {
        if (fabs(amax) - fabs(*MAT_PTR(i, j)) <= 0.0) {
          amax = *MAT_PTR(i, j);
          ik[k] = i;
          jk[k] = j;
        }
      }
    }

    /* Interchange rows and columns to put amax in array[k][k] */
    if (amax < 0.0) {
      *det = 0.0;
    } else {
      i = ik[k];
      if (i - k < 0) {
        continue;
      } else if (i - k == 0) {
        j = jk[k];
        if (j - k < 0) {
          continue;
        } else if (j - k == 0) {
          invert_matrix_accumulate(mat, norder, k, amax, det);
        } else { /* j - k > 0 */
          for (i = 0; i < norder; ++i) {
            /* Save and invert */
            save = *MAT_PTR(i, k);
            *MAT_PTR(i, k) = *MAT_PTR(i, j);
            *MAT_PTR(i, j) = -save;
          }
          invert_matrix_accumulate(mat, norder, k, amax, det);
        }
      } else { /* i - k > 0 */
        for (j = 0; j < norder; ++j) {
          /* Swap and invert */
          save = *MAT_PTR(k, j);
          *MAT_PTR(k, j) = *MAT_PTR(i, j);
          *MAT_PTR(i, j) = -save;
        }

        invert_matrix_accumulate(mat, norder, k, amax, det);
      }
    }
    amax = 0.0;
  }

  /* Restore ordering of matrix */
  for (l = 0; l < norder; ++l) {
    k = norder - l - 1; /* TODO: Base check */
    assert(k >= 0 && k < norder);

    j = ik[k];
    if (j - k > 0) {
      for (i = 0; i < norder; ++i) {
        save = *MAT_PTR(i, k);
        *MAT_PTR(i, k) = -*MAT_PTR(i, j);
        *MAT_PTR(i, j) = save;
      }
    }

    i = jk[k];
    if (i - k > 0) {
      for (j = 0; j < norder; ++j) {
        save = *MAT_PTR(k, j);
        *MAT_PTR(k, j) = -*MAT_PTR(i, j);
        *MAT_PTR(i, j) = save;
      }
    }
  }

  return 0;
}

/**
Fit a linear transformation (with six parameters, equivalent
to the linear WCS) to two sets of points.

This uses the standard least-squares linear equations method.
*/
static int
fit_linear_transformation(const double* xo /* [n] */, const double* yo /* [n] */,
                          const double* x /* [n] */, const double* y /* [n] */,
                          const integer_t n,
                          /* Output parameters */
                          double* x0, double* y0,
                          double* a, double* b, double* c, double* d,
                          struct driz_error_t* error) {
  integer_t i;
  double mat[9];
  const integer_t norder = 3;
  double xorg, yorg, xoorg, yoorg;
  double sigxox = 0.0, sigxoy = 0.0, sigxo = 0.0;
  double sigyox = 0.0, sigyoy = 0.0, sigyo = 0.0;
  double dx, dy, dox, doy;
  double det;

  assert(xo);
  assert(yo);
  assert(x);
  assert(y);
  assert(x0);
  assert(y0);
  assert(a);
  assert(b);
  assert(c);
  assert(d);
  assert(error);

  /* Initialize the matrix (3x3) */
  for (i = 0; i < 9; ++i) {
    mat[i] = 0.0;
  }

  /* Take an offset */
  xorg = x[0];
  yorg = y[0];
  xoorg = xo[0];
  yoorg = yo[0];

  /* Setup the normal equations */
  for (i = 0; i < n; ++i) {
    dx = x[i] - xorg;
    dy = y[i] - yorg;
    dox = xo[i] - xoorg;
    doy = yo[i] - yoorg;

    *MAT_PTR(0, 0) += dx*dx;
    *MAT_PTR(0, 1) += dx*dy;
    *MAT_PTR(0, 2) += dx;

    *MAT_PTR(1, 1) += dy*dy;
    *MAT_PTR(1, 2) += dy;

    sigxox += dox*dx;
    sigxoy += dox*dy;
    sigxo += dx;

    sigyox += doy*dx;
    sigyoy += doy*dy;
    sigyo += dy;
  }

  /* Use symmetry (the matrix is diagonal) */
  *MAT_PTR(2, 2) = (double)n;
  *MAT_PTR(1, 0) = *MAT_PTR(0, 1);
  *MAT_PTR(2, 0) = *MAT_PTR(0, 2);
  *MAT_PTR(2, 1) = *MAT_PTR(1, 2);

  /* Invert the matrix (we check it isn't singular) */
  if (invert_matrix(mat, 3, &det, error))
    return 1;

  if (det == 0.0) {
    driz_error_set_message(error, "Linear transformation matrix is singular");
    return 1;
  }

  /* Multiply the inverse by the vector */
  *a = sigxox * *MAT_PTR(0, 0) + sigxoy * *MAT_PTR(0, 1) + sigxo * *MAT_PTR(0, 2);
  *b = sigxox * *MAT_PTR(1, 0) + sigxoy * *MAT_PTR(1, 1) + sigxo * *MAT_PTR(1, 2);
  *x0 = sigxox * *MAT_PTR(2, 0) + sigxoy * *MAT_PTR(2, 1) + sigxo * *MAT_PTR(2, 2);

  *c = sigyox * *MAT_PTR(0, 0) + sigyoy * *MAT_PTR(0, 1) + sigyo * *MAT_PTR(0, 2);
  *d = sigyox * *MAT_PTR(1, 0) + sigyoy * *MAT_PTR(1, 1) + sigyo * *MAT_PTR(1, 2);
  *y0 = sigyox * *MAT_PTR(2, 0) + sigyoy * *MAT_PTR(2, 1) + sigyo * *MAT_PTR(2, 2);

  /* Note that x0 and y0 haven't been corrected for the offsets
     Normally they are not used
  */
  return 0;
}

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
                  struct driz_error_t* error) {
  bool_t new_reference;
  double xdref, ydref, xdoff, ydoff, xcen, ycen;
  double xin[4], yin[4], xout[4], yout[4], x[4], y[4];
  double ra, dec;
  double x0, y0, a, b, c, d;
  integer_t coeff_type;
  integer_t i;

  assert(m);
  assert(xc);
  assert(yc);
  assert(xs);
  assert(ys);
  assert(xt);
  assert(yt);
  assert(error);

  xcen = m->xcen;
  ycen = m->ycen;

  /* First check for the presence of "refpix" additional information
     in the coefficients. If it is, set a flag and offset again */
  if (m->coeff_type > COEFF_OFFSET) {
    new_reference = 1;
    m->coeff_type -= COEFF_OFFSET;
    xdref = m->x_coeffs[m->num_coeffs-1];
    ydref = m->y_coeffs[m->num_coeffs-1];
    m->num_coeffs--;
  } else {
    new_reference = 0;
    xdref = 0;
    ydref = 0;
  }

  /* Set up a square at the reference pixel of the input
     image (WCSIN(1)=CRPIX1 and WCSIN(3)=CRPIX2) */
  xin[0] = m->wcs[0];
  xin[1] = m->wcs[0];
  xin[2] = m->wcs[0] + 1.0;
  xin[3] = m->wcs[0] + 1.0;
  yin[0] = m->wcs[2];
  yin[1] = m->wcs[2] + 1.0;
  yin[2] = m->wcs[2] + 1.0;
  yin[3] = m->wcs[2];

  /* Transform these points onto the sky and then back out again
     using the target WCS - all double precision */
  for (i = 0; i < 4; ++i) {
    xy2rd(m->wcs, xin[i], yin[i], &ra, &dec);
    if (rd2xy(m->wcsout, ra, dec, &xout[i], &yout[i], error))
      return 1;
  }

  /* Check for different reference pixel */
  if (new_reference) {
    xdoff = xcen - xdref;
    ydoff = ycen - ydref;
  } else {
    xdoff = 0.0;
    ydoff = 0.0;
  }

  /* Now we apply the geometric distortion to the input points so that
     the linear transformation which we derive is appropriate after
     the distortion is corrected.

     Just use LINEAR terms */
  coeff_type = m->coeff_type;

  /* Why limit the evaluation to only linear terms???  This fit is
     only a linear fit, therefore, only linear terms need to be
     considered.  When all terms are used, it introduces slight
     non-orthogonality of the CD matrix after correction, as well as
     additional offsets between the chips.  AK, WJH 1-Aug-2006
  */
  if (coeff_type > 1) coeff_type = 1;

  if (coeff_type == 3) {
    for (i = 0; i < 4; ++i) {
      x[i] = eval3(xin[i] - xcen + xdoff, yin[i] - ycen + ydoff, m->x_coeffs) - xdoff;
      y[i] = eval3(xin[i] - xcen + xdoff, yin[i] - ycen + ydoff, m->y_coeffs) - ydoff;
    }
  } else if (coeff_type == 4) {
    for (i = 0; i < 4; ++i) {
      x[i] = eval4(xin[i] - xcen + xdoff, yin[i] - ycen + ydoff, m->x_coeffs) - xdoff;
      y[i] = eval4(xin[i] - xcen + xdoff, yin[i] - ycen + ydoff, m->y_coeffs) - ydoff;
    }
  } else if (coeff_type == 5) {
    for (i = 0; i < 4; ++i) {
      x[i] = eval5(xin[i] - xcen + xdoff, yin[i] - ycen + ydoff, m->x_coeffs) - xdoff;
      y[i] = eval5(xin[i] - xcen + xdoff, yin[i] - ycen + ydoff, m->y_coeffs) - ydoff;
    }
  } else if (coeff_type >= 6 || coeff_type == 1 || coeff_type == 2) {
    for (i = 0; i < 4; ++i) {
      if (evaln(xin[i] - xcen + xdoff, yin[i] - ycen + ydoff, m->x_coeffs, coeff_type, &x[i], error) ||
          evaln(xin[i] - xcen + xdoff, yin[i] - ycen + ydoff, m->y_coeffs, coeff_type, &y[i], error))
        return 1;
      x[i] -= xdoff;
      y[i] -= ydoff;
    }
  } else if (coeff_type == -3) {
    for (i = 0; i < 4; ++i) {
      rad3(xin[i] - xcen + xdoff, yin[i] - ycen + ydoff, m->x_coeffs, &x[i], &y[i]);
      x[i] -= xdoff;
      y[i] -= ydoff;
    }
  }

  /* Now we have the inputs and outputs and can derive the linear
     transform between them.  This is now done in a general way using
     least squares */
  if (fit_linear_transformation(xout, yout, x, y, 4, &x0, &y0, &a, &b, &c, &d, error))
    return 1;

  /* We change the sign here to fit in with the convention later */
  b *= -1.0;

  /* And now the linear offset */
  *xt = xout[0] - a*x[0] + b*y[0];
  *yt = yout[0] - c*x[0] + b*y[0];

  *xc = a;
  *ys = b;
  *xs = c;
  *yc = d;

  /* Before returning, reset the offsets, if there are any */
  if (new_reference) {
    m->coeff_type += COEFF_OFFSET;
    m->num_coeffs++;
  }

  return 0;
}

int
update_wcs(struct driz_param_t* p,
           struct mapping_param_t* m,
           const double wcsin[8],
           /* Output parameters */
           double wcsout[8],
           struct driz_error_t* error) {
  double xin[3], yin[3], xtmp[3], ytmp[3], xout[3], yout[3];
  bool_t old_use_distortion_image;
  integer_t old_coeff_type;
  double mat[4];
  double tmp[4];
  integer_t i;

  assert(p);
  assert(m);
  assert(wcsin);
  assert(wcsout);
  assert(error);

  /* If we have the WCS already, just return */
  if (m->use_wcs)
    return 0;

  /* Set up a single point at the reference pixel to map the reference
     point */
  xin[0] = wcsin[0];
  yin[0] = wcsin[2];

  /* Verify that reference pixel falls within distortion correction
     image.  If not, turn off use of distortion correction image */
  old_use_distortion_image = m->use_distortion_image;
  if (xin[0] < 0.0 ||
      xin[0] >= m->dnx ||
      yin[0] < 0.0 ||
      yin[0] >= m->dny) {
    m->use_distortion_image = FALSE;
  }

  /* Transform */
  if (map_value(p, FALSE, 1,
                xin, yin, xtmp, ytmp, xout, yout, error))
    return 1;

  m->use_distortion_image = old_use_distortion_image;

  /* We can immediately set the reference point on the sky */
  wcsout[0] = xout[0];
  wcsout[2] = yout[0];
  wcsout[1] = wcsin[1];
  wcsout[3] = wcsin[3];

  /* Set up a 1x1 box at the centre pixel (three sides only) to allow us
     to update the WCS */
  xin[1] = xin[0] + 1.0;
  yin[1] = yin[0];
  xin[2] = xin[0];
  yin[2] = yin[0] + 1.0;

  /* Transform.  Use only LINEAR terms and ignore distortion images */
  old_coeff_type = m->coeff_type;
  m->coeff_type = 1;
  old_use_distortion_image = m->use_distortion_image;
  m->use_distortion_image = FALSE;
  if (map_value(p, FALSE, 3,
                xin, yin, xtmp, ytmp, xout, yout, error)) {
    return 1;
  }

  /* Restore order */
  m->coeff_type = old_coeff_type;
  m->use_distortion_image = old_use_distortion_image;

  /* Now work out the effective CD matrix of the transformation */
  mat[0] = (xout[1] - xout[0]);
  mat[1] = (xout[2] - xout[0]);
  mat[2] = (yout[1] - yout[0]);
  mat[3] = (yout[2] - yout[0]);

  /* Invert the matrix */
  if (invert_matrix_2x2(mat, mat, error)) {
    return 1;
  }

  /* Correct the CD matrix */
  tmp[0] = mat[0]*wcsin[4] + mat[2]*wcsin[6];
  tmp[1] = mat[0]*wcsin[5] + mat[2]*wcsin[7];
  tmp[2] = mat[1]*wcsin[4] + mat[3]*wcsin[6];
  tmp[3] = mat[1]*wcsin[5] + mat[2]*wcsin[7];

  for (i = 0; i < 4; ++i) {
    wcsout[i+4] = tmp[i];
  }

  return 0;
}

int
blot_update_wcs(struct driz_param_t* p,
                struct mapping_param_t* m,
                const double wcsin[8],
                /* Output parameters */
                double wcsout[8],
                struct driz_error_t* error) {
  const double off = 0.1;
  double xin[3], yin[3], xtmp[3], ytmp[3], xout[3], yout[3];
  double r[3], d[3];
  integer_t i;

  assert(p);
  assert(m);
  assert(wcsin);
  assert(wcsout);
  assert(error);

  /* If we have the WCS already, just return */
  if (m->use_wcs)
    return 0;

  /* Set up a 1x1 box at the centre of the output image (three sides only) */
  xin[0] = (double)(m->dnx * 0.5);
  yin[0] = (double)(m->dny * 0.5);
  xin[1] = (double)(m->dnx * 0.5) + off;
  yin[1] = (double)(m->dny * 0.5);
  xin[2] = (double)(m->dnx * 0.5);
  xin[2] = (double)(m->dny * 0.5) + off;

  /* Transform */
  if (map_value(p, FALSE, 3,
                xin, yin, xtmp, ytmp, xout, yout, error)) {
    return 1;
  }

  /* Convert the output pixel position to a sky position using the WCS */
  for (i = 0; i < 3; ++i) {
    xy2rd(wcsin, xout[i], yout[i], &r[i], &d[i]);
  }

  /* We can immediately set the reference point on the sky */
  wcsout[0] = xin[0];
  wcsout[2] = yin[0];
  wcsout[1] = r[0];
  wcsout[3] = d[0];

  /* Now work out the effective CD matrix of the transformation */
  wcsout[4] = cos(d[0] * M_PI / 180.0) * (r[1] - r[0]) / off;
  wcsout[5] = (d[1] - d[0]) / off;
  wcsout[6] = cos(d[0] * M_PI / 180.0) * (r[2] - r[0]) / off;
  wcsout[7] = (d[2] - d[0]) / off;

  return 0;
}
