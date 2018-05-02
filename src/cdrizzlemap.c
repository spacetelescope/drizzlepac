#define NO_IMPORT_ARRAY
#define NO_IMPORT_ASTROPY_WCS_API
#include "driz_portability.h"
#include "astropy_wcs_api.h"

#define _USE_MATH_DEFINES       /* needed for MS Windows to define M_PI */
#include <math.h>
#include <string.h>
#include <time.h>

#include "cdrizzlemap.h"
#include "cdrizzlewcs.h"


static inline_macro int
drizzle_polynomial(void* state,
                   const double xd, const double yd,
                   const integer_t n,
                   const double* xin /*[n]*/, const double* yin /*[n]*/,
                   /* Output parameters */
                   double* xout, double* yout,
                   struct driz_error_t* error) {
  struct mapping_param_t* m = (struct mapping_param_t*)state;
  double xdoff, ydoff;
  integer_t i;
  bool_t new_reference;

  /* Check for the presence of "refpix" additional information in the
     coefficients.  If it is, set a flag and offset again */
  if (m->coeff_type > COEFF_OFFSET / 2) {
    new_reference = TRUE;
    m->coeff_type -= COEFF_OFFSET;
    xdoff = m->xcen - m->x_coeffs[m->num_coeffs - 1] + 1.0;
    ydoff = m->ycen - m->y_coeffs[m->num_coeffs - 1] + 1.0;
    m->num_coeffs--;
  } else {
    new_reference = FALSE;
    xdoff = 2.0;
    ydoff = 2.0;
  }

  if (m->coeff_type == 3) {
    for (i = 0; i < n; ++i) {
      xout[i] = eval3(xin[i] + xdoff, yin[i] + ydoff, m->x_coeffs) - xdoff;
      yout[i] = eval3(xin[i] + xdoff, yin[i] + ydoff, m->y_coeffs) - ydoff;
    }
  } else if (m->coeff_type == 4) {
    for (i = 0; i < n; ++i) {
      xout[i] = eval4(xin[i] + xdoff, yin[i] + ydoff, m->x_coeffs) - xdoff;
      yout[i] = eval4(xin[i] + xdoff, yin[i] + ydoff, m->y_coeffs) - ydoff;
    }
  } else if (m->coeff_type == 5) {
    for (i = 0; i < n; ++i) {
      xout[i] = eval5(xin[i] + xdoff, yin[i] + ydoff, m->x_coeffs) - xdoff;
      yout[i] = eval5(xin[i] + xdoff, yin[i] + ydoff, m->y_coeffs) - ydoff;
    }
  } else if (m->coeff_type >= 6 || m->coeff_type == 1 || m->coeff_type == 2) {
    for (i = 0; i < n; ++i) {
      if (evaln(xin[i] + xdoff, yin[i] + ydoff, m->x_coeffs, m->coeff_type, &xout[i], error) ||
          evaln(xin[i] + xdoff, yin[i] + ydoff, m->y_coeffs, m->coeff_type, &yout[i], error))
        return 1;
      xout[i] -= xdoff;
      yout[i] -= ydoff;
    }
  } else if (m->coeff_type == -3) {
    for (i = 0; i < n; ++i) {
      rad3(xin[i] + xdoff, yin[i] + ydoff, m->x_coeffs, &xout[i], &yout[i]);
      xout[i] -= xdoff;
      yout[i] -= ydoff;
    }
  } else {
    driz_error_format_message(error, "Invalid coefficient type %d", m->coeff_type);
    return 1;
  }

  if (new_reference) {
    m->coeff_type += COEFF_OFFSET;
    m->num_coeffs++;
  }

  return 0;
}

static inline_macro int
drizzle_distortion_image(void* state,
                         const double xd, const double yd,
                         const integer_t n,
                         const double* xin /*[n]*/, const double* yin /*[n]*/,
                         /* Output parameters */
                         double* xout, double* yout,
                         struct driz_error_t* error) {
  struct mapping_param_t* m = (struct mapping_param_t*)state;
  integer_t ix, iy;
  integer_t i;

  if (m->use_distortion_image) {
    for (i = 0; i < n; ++i) {
      ix = (integer_t)(xin[i] + 1.0 - xd + m->xcen) - 1;
      iy = (integer_t)(yin[i] + 1.0 - yd + m->ycen) - 1;
      xout[i] += (double)(*x_distortion_ptr(m, ix, iy));
      yout[i] += (double)(*y_distortion_ptr(m, ix, iy));
    }
  }

  return 0;
}

static inline_macro int
drizzle_alpha_beta(void* state,
                   const double xd, const double yd,
                   const integer_t n,
                   const double* xin /*[n]*/, const double* yin /*[n]*/,
                   /* Output parameters */
                   double* xout, double* yout,
                   struct driz_error_t* error) {
  /* TODO: There is definitely something fishy in the original here
     since the x in in the first line will have an effect in the
     second.  We use xout and yout temporaries here to eliminate that. */
  struct mapping_param_t* m = (struct mapping_param_t*)state;
  double xt, yt, xa, ya;
  integer_t i;

  if (m->alpha != 0.0 || m->beta != 0.0) {
    for (i = 0; i < n; ++i) {
      xa = (xin[i] + 1.0 - m->xcen) / m->xcen;
      ya = (yin[i] + 1.0 - m->ycen) / m->ycen;

      xt = m->beta * xa + m->alpha * ya - 1.0;
      yt = m->beta * ya + m->alpha * xa - 1.0;
      xout[i] += xt;
      yout[i] += yt;
    }
  }

  return 0;
}

static inline_macro int
drizzle_linear_or_secondary(void* state,
                            const double xd, const double yd,
                            const integer_t n,
                            const double* xin /*[n]*/, const double* yin /*[n]*/,
                            /* Output parameters */
                            double* xout, double* yout,
                            struct driz_error_t* error) {
  struct mapping_param_t* m = (struct mapping_param_t*)state;
  double xs, ys, xc, yc, xt, yt;
  double xs2, ys2, xc2, yc2, xp2, yp2, xt2, yt2;
  double xtmp, ytmp;
  double sinth, costh, xf, yf, xoff, yoff;
  integer_t i;

  sinth = sin(m->rotation);
  costh = cos(m->rotation);
  xoff = m->x_shift / m->scale;
  yoff = m->y_shift / m->scale;
  xf = 1.0 / m->scale;
  yf = 1.0 / m->scale;

  xs = xf * sinth;
  xc = xf * costh;
  ys = yf * sinth;
  yc = yf * costh;
  xt = xoff + m->xp;
  yt = yoff + m->yp;

  if (m->use_wcs) {
    if (wcs_derive_linear(m, &xc, &yc, &xs, &ys, &xt, &yt, error))
      return 1;
    m->do_shift_first = shift_input;
  }

  if (m->has_secondary_parameters) {
    xs2 = ys2 = sin(m->rotation2);
    xc2 = yc2 = cos(m->rotation2);
    xp2 = m->xp;
    yp2 = m->yp;
    xt2 = xp2 + m->x_shift2;
    yt2 = yp2 + m->y_shift2;

    for (i = 0; i < n; ++i) {
      /* Apply the linear transform
         There are two ways this can be done - shift then
         rotate or rotate then shift */
      if (m->do_shift_first == shift_input) {
        xtmp = xc * (xin[i] + 1.0) - ys * (yin[i] + 1.0) + xt - m->xp;
        ytmp = xs * (xin[i] + 1.0) + yc * (yin[i] + 1.0) + yt - m->yp;
      } else {
        xtmp = xc * (xin[i] + 1.0 + m->x_shift) - ys * (yin[i] + 1.0 + m->y_shift);
        ytmp = xs * (xin[i] + 1.0 + m->x_shift) + yc * (yin[i] + 1.0 + m->y_shift);
      }

      /* Apply the secondary transform */
      if (m->do_shift2_first == shift_output) {
        xout[i] = xc2 * xtmp - ys2 * ytmp + xt2 - 1.0;
        yout[i] = xs2 * xtmp + yc2 * ytmp + yt2 - 1.0;
      } else {
        xout[i] = xc2 * (xtmp + m->x_shift2) -
          ys2 * (ytmp + m->y_shift2) + xp2 - 1.0;
        yout[i] = xs2 * (xtmp + m->x_shift2) +
          yc2 * (ytmp + m->y_shift2) + yp2 - 1.0;
      }
    }
  } else {
    if (m->do_shift_first == shift_output) {
      for (i = 0; i < n; ++i) {
        xout[i] = xc * (xin[i] + 1.0) - ys * (yin[i] + 1.0) + xt - 1.0;
        yout[i] = xs * (xin[i] + 1.0) + yc * (yin[i] + 1.0) + yt - 1.0;
      }
    } else {
      for (i = 0; i < n; ++i) {
        xout[i] = xc * (xin[i] + 1.0 + m->x_shift) - ys * (yin[i] + 1.0 + m->y_shift) + m->xp - 1.0;
        yout[i] = xs * (xin[i] + 1.0 + m->x_shift) + yc * (yin[i] + 1.0 + m->y_shift) + m->yp - 1.0;
      }
    }
  }

  return 0;
}

int
map_value(struct driz_param_t* p,
          const bool_t regular,
          const integer_t n,
          const double* xin /*[n]*/, const double* yin /*[n]*/,
          /* Output parameters */
          double* xtmp /*[n]*/, double* ytmp /*[n]*/,
          double* xout /*[n]*/, double* yout /*[n]*/,
          struct driz_error_t* error) {
  double x, y, xd, yd;
  integer_t i;

  assert(p);
  assert(p->mapping_callback);
  assert(xin);
  assert(yin);
  assert(xtmp);
  assert(ytmp);
  assert(xout);
  assert(yout);
  assert(xin != xout);
  assert(yin != yout);
  assert(xin != xtmp);
  assert(yin != ytmp);
  assert(xtmp != xout);
  assert(ytmp != yout);
  assert(error);

  if (regular) {
    /* x = xin[0] - p->x_scale + 1.0; */
    /* y = yin[0] + yin[1] + 2.0; */
    x = xin[0];
    y = yin[0];
    xd = xin[0];
    yd = yin[1];

    for (i = 0; i < n; ++i) {
      xtmp[i] = x;
      ytmp[i] = y;
      x += p->x_scale;
    }
  } else {
    xd = yd = 0.0;

    if (xtmp != xin) {
        memcpy(xtmp, xin, sizeof(double) * n);
    }
    if (ytmp != yin) {
        memcpy(ytmp, yin, sizeof(double) * n);
    }
  }

  if (p->mapping_callback(p->mapping_callback_state, xd, yd, n,
                          xtmp, ytmp, xout, yout, error))
    return 1;

  return 0;
}

static int
default_wcsmap_direct(struct wcsmap_param_t* m,
                      const double xd, const double yd,
                      const integer_t n,
                      double* xin /*[n]*/, double* yin /*[n]*/,
                      /* Output parameters */
                      double* xout, double* yout,
                      struct driz_error_t* error) {

  integer_t  i;
  int        status;
  double    *memory = NULL;
  double    *ptr    = NULL;
  double    *xyin   = NULL;
  double    *skyout = NULL;
  double    *xyout  = NULL;
  double    *imgcrd = NULL;
  double    *phi    = NULL;
  double    *theta  = NULL;
  int       *stat   = NULL;

  /* Allocate memory for new 2-D array */
  ptr = memory = (double *) malloc(n * 10 * sizeof(double));
  if (memory == NULL) return 1;

  xyin = ptr;
  ptr += n * 2;
  xyout = ptr;
  ptr += n * 2;
  skyout = ptr;
  ptr += n * 2;
  imgcrd = ptr;
  ptr += n * 2;
  phi = ptr;
  ptr += n;
  theta = ptr;

  stat = (int *)malloc(n * sizeof(int));
  if (stat == NULL) {
      free(memory);
      return 1;
  }

  /* The input arrays need to be converted to 2-D arrays for input
     to the PyWCS (and related) functions. */

  /* Populate new 2-D array with values from x and y input arrays */
  for (i = 0; i < n; ++i) {
    xyin[2*i] = xin[i];
    xyin[2*i+1] = yin[i];
  }

  /*
    Apply pix2sky() transformation from PyWCS
  */

  wcsprm_python2c(m->input_wcs->wcs);
  status = pipeline_all_pixel2world(m->input_wcs, n, 2, xyin, skyout);
  wcsprm_c2python(m->input_wcs->wcs);
  if (status) {
    free(memory);
    free(stat);
    return 1;
  }

  /*
    Finally, call wcs_sky2pix() for the output object.
  */
  wcsprm_python2c(m->output_wcs->wcs);
  status = wcss2p(m->output_wcs->wcs, n, 2,
                  skyout, phi, theta, imgcrd, xyout, stat);
  wcsprm_c2python(m->output_wcs->wcs);
  if (status) {
    free(memory);
    free(stat);
    return 1;
  }

  /*
    Transform results back to 2 1-D arrays, like the input.
  */
  for (i = 0; i < n; ++i){
    xout[i] = xyout[2*i];
    yout[i] = xyout[2*i+1];
  }

  /* Free memory allocated to internal 2-D arrays */
  free(memory);
  free(stat);
  return 0;
}

static int
default_wcsmap_interpolate(struct wcsmap_param_t* m,
                           const double xd, const double yd,
                           const integer_t n,
                           double* xin /*[n]*/, double* yin /*[n]*/,
                           /* Output parameters */
                           double* xout, double* yout,
                           struct driz_error_t* error) {

  int     i;
  double *xiptr;
  double *yiptr;
  double *xoptr;
  double *yoptr;
  double *table;
  double  x, y;
  int     xi, yi;
  double  xf, yf, ixf, iyf;
  double  tabx00, tabx01, tabx10, tabx11;

  /* do the bilinear interpolation */
  xiptr = xin;
  yiptr = yin;
  xoptr = xout;
  yoptr = yout;
  table = m->table;

#define TABLE_X(x, y) (table[((y)*m->snx + (x))*2])
#define TABLE_Y(x, y) (table[((y)*m->snx + (x))*2 + 1])

  for (i = 0; i < n; ++i) {
    x = *xiptr++ / m->factor;
    y = *yiptr++ / m->factor;
    xi = (int)floor(x);
    yi = (int)floor(y);
    xf = x - (double)xi;
    yf = y - (double)yi;
    ixf = 1.0 - xf;
    iyf = 1.0 - yf;

    tabx00 = TABLE_X(xi, yi);
    tabx10 = TABLE_X(xi+1, yi);
    tabx01 = TABLE_X(xi, yi+1);
    tabx11 = TABLE_X(xi+1, yi+1);

    /* Account for interpolating across 360-0 boundary */
    if ((tabx00 - tabx10) > 359) {
      tabx00 -= 360.0;
      tabx01 -= 360.0;
    } else if ((tabx00 - tabx10) < -359) {
      tabx10 -= 360.0;
      tabx11 -= 360.0;
    }

    *xoptr++ =
      tabx00 * ixf * iyf +
      tabx10 * xf * iyf +
      tabx01 * ixf * yf +
      tabx11 * xf * yf;

    *yoptr++ =
      TABLE_Y(xi, yi)     * ixf * iyf +
      TABLE_Y(xi+1, yi)   * xf * iyf +
      TABLE_Y(xi, yi+1)   * ixf * yf +
      TABLE_Y(xi+1, yi+1) * xf * yf;
  }

#undef TABLE_X
#undef TABLE_Y

  return 0;
}



/*

Default WCS mapping code

*/
int
default_wcsmap(void* state,
               const double xd, const double yd,
               const integer_t n,
               double* xin /*[n]*/, double* yin /*[n]*/,
               /* Output parameters */
               double* xout, double* yout,
               struct driz_error_t* error) {

  struct wcsmap_param_t* m = (struct wcsmap_param_t*)state;

  if (m->factor == 0) {
    return default_wcsmap_direct(m, xd, yd, n, xin, yin, xout, yout, error);
  } else {
    return default_wcsmap_interpolate(m, xd, yd, n, xin, yin, xout, yout, error);
  }
}

int
default_wcsmap_init(struct wcsmap_param_t* m,
                    pipeline_t* input,
                    pipeline_t* output,
                    int nx, int ny,
                    double factor,
                    struct driz_error_t* error) {
  int     n;
  int     table_size;
  double *pixcrd = NULL;
  double *ptr    = NULL;
  double *tmp    = NULL;
  double *phi    = NULL;
  double *theta  = NULL;
  double *imgcrd = NULL;
  int    *stat   = NULL;
  int     snx = nx + 2;
  int     sny = ny + 2;
  int     i;
  int     j;
  int     istat;

  assert(m);
  assert(input);
  assert(output);
  assert(m->input_wcs == NULL);
  assert(m->output_wcs == NULL);
  assert(m->table == NULL);

  if (factor > 0) {
    snx = (int)((double)nx / factor) + 2;
    sny = (int)((double)ny / factor) + 2;

    n = (snx) * (sny);
    table_size = n << 1;

    pixcrd = malloc(table_size * sizeof(double));
    if (pixcrd == NULL) {
      driz_error_set_message(error, "Out of memory");
      goto exit;
    }

    m->table = malloc(table_size * sizeof(double));
    if (m->table == NULL) {
      driz_error_set_message(error, "Out of memory");
      goto exit;
    }

    tmp = malloc(table_size * sizeof(double));
    if (tmp == NULL) {
      driz_error_set_message(error, "Out of memory");
      goto exit;
    }

    phi = malloc(n * sizeof(double));
    if (phi == NULL) {
      driz_error_set_message(error, "Out of memory");
      goto exit;
    }

    theta = malloc(n * sizeof(double));
    if (theta == NULL) {
      driz_error_set_message(error, "Out of memory");
      goto exit;
    }

    imgcrd = malloc(table_size * sizeof(double));
    if (imgcrd == NULL) {
      driz_error_set_message(error, "Out of memory");
      goto exit;
    }

    stat = malloc(n * sizeof(int));
    if (stat == NULL) {
      driz_error_set_message(error, "Out of memory");
      goto exit;
    }

    ptr = pixcrd;
    for (j = 0; j < sny; ++j) {
      for (i = 0; i < snx; ++i) {
        *ptr++ = (double)i * factor;
        *ptr++ = (double)j * factor;
      }
    }

    wcsprm_python2c(input->wcs);
    istat = pipeline_all_pixel2world(input, n, 2, pixcrd, tmp);
    wcsprm_c2python(input->wcs);

    if (istat) {
      free(m->table);
      m->table = NULL;
      driz_error_set_message(error, wcslib_get_error_message(istat));
      goto exit;
    }

    wcsprm_python2c(output->wcs);
    istat = wcss2p(output->wcs, n, 2, tmp, phi, theta, imgcrd, m->table, stat);
    wcsprm_c2python(output->wcs);

    if (istat) {
      free(m->table);
      m->table = NULL;
      driz_error_set_message(error, wcslib_get_error_message(istat));
      goto exit;
    }
  } /* End if_then for factor > 0 */

  m->input_wcs = input;
  m->output_wcs = output;

  m->nx = nx;
  m->ny = ny;
  m->snx = snx;
  m->sny = sny;
  m->factor = factor;

 exit:

  free(pixcrd);
  free(tmp);
  free(phi);
  free(theta);
  free(imgcrd);
  free(stat);

  return 0;
}

void
wcsmap_param_dump(struct wcsmap_param_t* m) {
  assert(m);

  printf("WCS MAPPING PARAMETERS:\n"
         "input WCS:            \n");
  wcsprt(m->input_wcs->wcs);
  printf("output WCS:           \n");
  wcsprt(m->input_wcs->wcs);
}

void
wcsmap_param_free(struct wcsmap_param_t* m) {
  free(m->table);
  wcsmap_param_init(m);
}

void
wcsmap_param_init(struct wcsmap_param_t* m) {
  assert(m);

  /* Pointers to the PyWCS objects */
  m->input_wcs = NULL;
  m->output_wcs = NULL;
  m->table = NULL;
}

/*

Default pixel-based mapping code:DefaultMapping

*/

int
default_mapping(void* state,
                const double xd, const double yd,
                const integer_t n,
                double* xin /*[n]*/, double* yin /*[n]*/,
                /* Output parameters */
                double* xout, double* yout,
                struct driz_error_t* error) {
  struct mapping_param_t* m = (struct mapping_param_t*)state;
  integer_t i;

  /* The built-in "default" mapping needs some pre-processing on its
     input values. */
  for (i = 0; i < n; ++i) {
    xin[i] = m->x_scale * (xin[i] + 1.0 - m->xcen) - 1.0;
    yin[i] = m->y_scale * (yin[i] + 1.0 - m->ycen) - 1.0;
  }

  if (drizzle_polynomial(state, xd, yd, n, xin, yin, xout, yout, error))
    return 1;

  if (drizzle_distortion_image(state, xd, yd, n, xin, yin, xout, yout, error))
    return 1;

  /* These next two build on the values already in xout/yout, so they
     are passed as inputs.  This is not a mistake. */
  if (drizzle_alpha_beta(state, xd, yd, n, xout, yout, xout, yout, error))
    return 1;

  if (drizzle_linear_or_secondary(state, xd, yd, n, xout, yout, xout, yout, error))
    return 1;

  return 0;
}
