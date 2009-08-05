#define NO_IMPORT_ARRAY
#define NO_IMPORT_PYWCS_API
#include "pywcs_api.h"

#include <math.h>
#include <string.h>
#include <time.h>

#include "cdrizzleio.h"
#include "cdrizzlemap.h"
#include "cdrizzlewcs.h"


static inline int
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

static inline int
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

static inline int
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

static inline int
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

/*#define WCSMAP_ORIGINAL_SLOW */

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

  integer_t  i;
  int        status;
  int        result = 1;
  double    *memory = NULL;
  double    *ptr    = NULL;
  double    *xyin   = NULL;
  double    *skyout = NULL;
  double    *xyout  = NULL;
  double    *imgcrd = NULL;
  double    *phi    = NULL;
  double    *theta  = NULL;
  int       *stat   = NULL;
  double    *iptr;
  double    *optr;
  double    *table;
  double     x, y;
  int        xi, yi;
  double     xf, yf, ixf, iyf;

  /* Call PyWCS methods here to perform the transformation... */
  /* The input arrays need to be converted to 2-D arrays for input
      to the PyWCS (and related) functions. */

  /* Allocate memory for new 2-D array */
  ptr = memory = (double *) malloc(n * 10 * sizeof(double));
  if (memory == NULL) goto exit;
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
  if (stat == NULL) goto exit;

  /* Populate new 2-D array with values from x and y input arrays */
  for (i = 0; i < n; ++i) {
      xyin[2*i] = xin[i];
      xyin[2*i+1] = yin[i];
  }

  /* Start by checking to see whether DET2IM correction needs to
  be applied and applying it as appropriate. */

#ifdef WCSMAP_ORIGINAL_SLOW
  /*
  Apply pix2sky() transformation from PyWCS
  */

  wcsprm_python2c(m->input_wcs->wcs);
  status = pipeline_all_pixel2world(m->input_wcs, n, 2, xyin, skyout);
  if (status)
    return 1;
  wcsprm_c2python(m->input_wcs->wcs);

#else
  /* TODO: Apply only distortion and sip in this way, and the rest using the lookup table */
  /* status = pipeline_pix2foc(m->input_wcs, n, 2, xyin, xyin); */
  /* if (status) */
  /*   return 1; */

  /* do the bilinear interpolation */
  iptr = xyin;
  optr = skyout;
  table = m->table;

#define TABLE_X(x, y) (table[((y)*m->snx + (x))*2])
#define TABLE_Y(x, y) (table[((y)*m->snx + (x))*2 + 1])

  for (i = 0; i < n; ++i) {
    x = *iptr++ / m->factor;
    y = *iptr++ / m->factor;
    xi = (int)floor(x);
    yi = (int)floor(y);
    xf = x - (double)xi;
    yf = y - (double)yi;
    ixf = 1.0 - xf;
    iyf = 1.0 - yf;

    *optr++ =
      TABLE_X(xi, yi)     * ixf * iyf +
      TABLE_X(xi+1, yi)   * xf * iyf +
      TABLE_X(xi, yi+1)   * ixf * yf +
      TABLE_X(xi+1, yi+1) * xf * yf;

    *optr++ =
      TABLE_Y(xi, yi)     * ixf * iyf +
      TABLE_Y(xi+1, yi)   * xf * iyf +
      TABLE_Y(xi, yi+1)   * ixf * yf +
      TABLE_Y(xi+1, yi+1) * xf * yf;
  }

#undef TABLE_X
#undef TABLE_Y

#endif

  /*
  Finally, call wcs_sky2pix() for the output object.
  */
  wcsprm_python2c(m->output_wcs->wcs);
  status = wcss2p(m->output_wcs->wcs, n, 2,
                  skyout, phi, theta, imgcrd, xyout, stat);
  wcsprm_c2python(m->output_wcs->wcs);

  /*
  Transform results back to 2 1-D arrays, like the input.
  */
  for (i = 0; i < n; ++i){
      xout[i] = xyout[2*i];
      yout[i] = xyout[2*i+1];
  }

  result = 0;

 exit:
  /*
  Free memory allocated to internal 2-D arrays
  */
  free(memory);
  free(stat);

  return result;
}

int
default_wcsmap_init(struct wcsmap_param_t* m,
                    pipeline_t* input,
                    pipeline_t* output,
                    int nx, int ny,
                    double factor,
                    struct driz_error_t* error) {
  int     table_size;
  double *pixcrd = NULL;
  double *ptr  = NULL;
  int     snx;
  int     sny;
  int     i;
  int     j;
  int     status = 1;
  int     istat;

  assert(m);
  assert(input);
  assert(output);
  assert(m->input_wcs == NULL);
  assert(m->output_wcs == NULL);
  assert(m->table == NULL);

#ifndef WCSMAP_ORIGINAL_SLOW
  snx = (int)((double)nx / factor) + 2;
  sny = (int)((double)ny / factor) + 2;

  table_size = (snx + 1) * (sny + 1) * 2;

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

  ptr = pixcrd;
  for (j = 0; j < sny; ++j) {
    for (i = 0; i < snx; ++i) {
      *ptr++ = (double)i * factor;
      *ptr++ = (double)j * factor;
    }
  }

  wcsprm_python2c(input->wcs);
  istat = pipeline_all_pixel2world(input, table_size / 2, 2, pixcrd,
                                   m->table);
  wcsprm_c2python(input->wcs);

  if (istat) {
    free(m->table);
    driz_error_set_message(error, wcslib_get_error_message(istat));
    goto exit;
  }
#endif

  m->input_wcs = input;
  m->output_wcs = output;

  m->nx = nx;
  m->ny = ny;
  m->snx = snx;
  m->sny = sny;
  m->factor = factor;

  status = 0;

 exit:

  free(pixcrd);

  return status;
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

static int
default_mapping_init_wcs(struct mapping_param_t* m,
                         struct driz_error_t* error) {
  struct driz_param_t p;
  double xin[4], yin[4], xtmp[4], ytmp[4], xout[4], yout[4];
  double scall, scdis;

  driz_param_init(&p);

  xin[0] = m->xcen;
  xin[1] = m->xcen;
  xin[2] = m->xcen + 1.0;
  xin[3] = m->xcen + 1.0;

  yin[0] = m->ycen;
  yin[1] = m->ycen + 1.0;
  yin[2] = m->ycen + 1.0;
  yin[3] = m->ycen;

  p.x_scale = m->x_scale;
  p.y_scale = m->y_scale;

  if (map_value(&p, FALSE, 4, xin, yin, xtmp, ytmp, xout, yout, error)) {
    return 1;
  }

  scall = sqrt(1.0 / abs(0.5 * (xout[1] - xout[3]) * (yout[0] - yout[2]) -
                         (xout[0] - xout[2]) * (yout[1] - yout[3])));

  /* Now calculate how much of this if from the geometric distortion */
  m->x_shift = 0.0;
  m->y_shift = 0.0;
  m->rotation = 0.0;
  m->scale = 1.0;
  m->has_secondary_parameters = FALSE;
  m->use_wcs = FALSE;

  if (map_value(&p, FALSE, 4, xin, yin, xtmp, ytmp, xout, yout, error)) {
    return 1;
  }

  scdis = sqrt(1.0 / fabs(0.5 * (xout[1] - xout[3]) * (yout[0] - yout[2]) -
                          (xout[0] - xout[2]) * (yout[1] - yout[3])));
  m->use_wcs = TRUE;
  m->scale = scall / scdis;

  return 0;
}


int
default_mapping_init(struct mapping_param_t* m,
                     const integer_t dnx,
                     const integer_t dny,
                     const integer_t onx,
                     const integer_t ony,
                     const double xsh,
                     const double ysh,
                     const enum e_shift_t shftfr,
                     const enum e_shift_t shftun,
                     const double drot,
                     const double scale,
                     const double xsh2,
                     const double ysh2,
                     const double xscale,
                     const double yscale,
                     const double rot2,
                     enum e_shift_t shfr2,
                     const float* pxg,
                     const float* pyg,
                     const integer_t xgdim,
                     const integer_t ygdim,
                     const enum e_align_t align,
                     const char* coeffs,
                     const double wcs[8],
                     const double alpha,
                     const double beta,
                     struct driz_error_t* error) {
  integer_t i;
  float align_offset;

  mapping_param_init(m);
  /* Get the geometric distortion coefficient information */
  /* Get geometric distortion coefficients */
  if (get_geometric_distortion(coeffs,
                               NULL,
                               &m->lambda,
                               &m->coeff_type,
                               &m->num_coeffs,
                               m->x_coeffs,
                               m->y_coeffs,
                               error))
    return 1;

  m->x_shift = xsh;
  m->y_shift = ysh;
  m->do_shift_first = shftfr;
  m->shift_units = shftun;
  m->scale = scale;
  m->x_shift2 = xsh2;
  m->y_shift2 = ysh2;
  m->x_scale = xscale;
  m->y_scale = yscale;
  m->do_shift2_first = shfr2;
  m->x_distortion = pxg;
  m->y_distortion = pyg;
  m->x_dist_dim = xgdim;
  m->y_dist_dim = ygdim;
  m->align = align;
  m->dnx = dnx;
  m->dny = dny;
  if (wcs != NULL) {
    for (i = 0; i < 8; ++i) {
      m->wcs[i] = wcs[i];
    }
  }

  /* Set disim logical based on whether distortion images have
     been passed in for use or not. */
  m->use_distortion_image = !(m->x_dist_dim == 2 && m->y_dist_dim == 2);

  /* Convert the shift units if necessary */
  if (m->shift_units == shift_output) {
    m->x_shift *= m->scale;
    m->y_shift *= m->scale;
  }

  /* Convert the rotation to radians */
  /* TODO: Use double precision here */
  m->rotation = drot * (float)M_PI / 180.0f;

  /* Convert the secondary parameters into suitable units */
  /* TODO: Use double precision here */
  m->rotation2 = rot2 * (float)M_PI / 180.0f;

  /* Give a warning message if secondary parameters will have an effect */
  m->has_secondary_parameters =
    (m->x_scale != 1.0 || m->y_scale != 1.0 ||
     m->x_shift2 != 0.0 || m->y_shift2 != 0.0 ||
     m->rotation2 != 0.0);

  /* Note that the xcen and xp values are the opposite of what they
     are in dobox */
  align_offset = (m->align == align_corner) ? 0.5f : 1.0f;
  m->xcen = (double)(dnx / 2) + align_offset;
  m->ycen = (double)(dny / 2) + align_offset;
  m->xp = (double)(onx / 2) + align_offset;
  m->yp = (double)(ony / 2) + align_offset;

  m->use_wcs = FALSE;

  if (m->use_wcs) {
    if (default_mapping_init_wcs(m, error)) {
      return 1;
    }
  }

  return 0;
}

void
mapping_param_dump(struct mapping_param_t* m) {
  assert(m);

  printf("MAPPING PARAMETERS:\n"
         "num_coeffs:           %d\n"
         "coeff_type:           %d\n"
         "lambda:               %f\n"
         "x_distortion:         %p\n"
         "x_dist_dim:           %d\n"
         "y_distortion:         %p\n"
         "y_dist_dim:           %d\n"
         "use_distortion_image: %s\n"
         "scale:                %f\n"
         "rotation:             %f\n"
         "x_shift:              %f\n"
         "y_shift:              %f\n"
         "do_shift_first:       %s\n"
         "shift_units:          %s\n"
         "align:                %s\n"
         "x_scale:              %f\n"
         "y_scale:              %f\n"
         "x_shift2:             %f\n"
         "y_shift2:             %f\n"
         "rotation2:            %f\n"
         "do_shift_first2:      %s\n"
         "alpha:                %f\n"
         "beta:                 %f\n\n"
         "xcen:                 %f\n"
         "ycen:                 %f\n"
         "xp:                   %f\n"
         "yp:                   %f\n"
         "dnx:                  %d\n"
         "dny:                  %d\n",

         m->num_coeffs,
         m->coeff_type,
         m->lambda,
         m->x_distortion,
         m->x_dist_dim,
         m->y_distortion,
         m->y_dist_dim,
         bool2str(m->use_distortion_image),
         m->scale,
         m->rotation,
         m->x_shift,
         m->y_shift,
         shift_enum2str(m->do_shift_first),
         shift_enum2str(m->shift_units),
         align_enum2str(m->align),
         m->x_scale,
         m->y_scale,
         m->x_shift2,
         m->y_shift2,
         m->rotation2,
         shift_enum2str(m->do_shift2_first),
         m->alpha,
         m->beta,
         m->xcen,
         m->ycen,
         m->xp,
         m->yp,
         m->dnx,
         m->dny
         );
}

void
mapping_param_init(struct mapping_param_t* m) {
  assert(m);

  /* Geometric distortion coefficients */
  m->num_coeffs = 0;
  m->coeff_type = 0;
  m->lambda = 0.0;

  /* Distortion image arrays */
  m->x_distortion = NULL;
  m->x_dist_dim = 0;
  m->y_distortion = NULL;
  m->y_dist_dim = 0;
  m->use_distortion_image = 0;

  /* Primary linear transformation parameters */
  m->scale = 1.0;
  m->rotation = 0.0;
  m->x_shift = 0.0;
  m->y_shift = 0.0;
  m->do_shift_first = shift_input;
  m->shift_units = shift_input;
  m->align = align_center;

  /* Secondary transformation parameters */
  m->x_scale = 1.0;
  m->y_scale = 1.0;
  m->x_shift2 = 0.0;
  m->y_shift2 = 0.0;
  m->rotation2 = 0.0;
  m->do_shift2_first = shift_input;

  /* WCS mode parameters */
  m->use_wcs = FALSE;
  m->output_scale = 1.0;
  m->ra_ref = 0.0;
  m->dec_ref = 0.0;
  m->x_ref_pixel = 0.0;
  m->y_ref_pixel = 0.0;
  m->orientation = 0.0;

  m->alpha = 0.0;
  m->beta = 0.0;

  m->xcen = 0.0;
  m->ycen = 0.0;
  m->xp = 0.0;
  m->yp = 0.0;
  m->dnx = 0.0;
  m->dny = 0.0;
}



