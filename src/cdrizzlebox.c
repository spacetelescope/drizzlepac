#define NO_IMPORT_ARRAY
#define NO_IMPORT_ASTROPY_WCS_API

#include "driz_portability.h"
#include "cdrizzlemap.h"
#include "cdrizzlebox.h"
#include "cdrizzleutil.h"

#include <assert.h>
#define _USE_MATH_DEFINES       /* needed for MS Windows to define M_PI */
#include <math.h>
#include <stdlib.h>

/**
Check how much of a line will overlap an output image, if any,
after applying the standard drizzle transformation.

This is intended to allow the number of points which are needlessly
drizzled outside the output image to be minimized.

was: CHOVER
*/
#define CHECK_OVER_NPOINT 21

inline_macro static void
update_data(struct driz_param_t* p, const integer_t ii, const integer_t jj,
            const float d, const float vc, const float dow) {

  const double vc_plus_dow = vc + dow;

  /* Just a simple calculation without logical tests */
  if (vc == 0.0) {
    set_pixel(p->output_data, ii, jj, d);

  } else if (vc_plus_dow != 0.0) {
    double value;
    value = (get_pixel(p->output_data, ii, jj) * vc + dow * d) / (vc_plus_dow);
    set_pixel(p->output_data, ii, jj, value);
  }

  set_pixel(p->output_counts, ii, jj, vc_plus_dow);
}

/**
To calculate area under a line segment within unit square at origin.
This is used by BOXER.

NOTE: This is the single most frequently called function.  Ripe
for optimization.
*/
static inline_macro double
sgarea(const double x1, const double y1, const double x2, const double y2) {
  double m, c, dx, dy, xlo, xhi, ylo, yhi, xtop;
  int negdx;

  dy = y2 - y1;

  dx = x2 - x1;
  /* Trap vertical line */
  if (dx == 0.0)
    return 0.0;

  negdx = (int)(dx < 0.0);
  if (negdx) {
    xlo = x2;
    xhi = x1;
  } else {
    xlo = x1;
    xhi = x2;
  }

  /* And determine the bounds ignoring y for now */
  if (xlo >= 1.0 || xhi <= 0.0)
    return 0.0;

  xlo = MAX(xlo, 0.0);
  xhi = MIN(xhi, 1.0);

  /* Now look at y */
  m = dy / dx;
  assert(m != 0.0);
  c = y1 - m * x1;
  ylo = m * xlo + c;
  yhi = m * xhi + c;

  /* Trap segment entirely below axis */
  if (ylo <= 0.0 && yhi <= 0.0)
    return 0.0;

  /* Adjust bounds if segment crosses axis (to exclude anything below
     axis) */
  if (ylo < 0.0) {
    ylo = 0.0;
    xlo = -c / m;
  }

  if (yhi < 0.0) {
    yhi = 0.0;
    xhi = -c / m;
  }

  /* There are four possibilities: both y below 1, both y above 1 and
     one of each. */
  if (ylo >= 1.0 && yhi >= 1.0) {
    /* Line segment is entirely above square */
    if (negdx) {
      return xlo - xhi;
    } else {
      return xhi - xlo;
    }
  }

  if (ylo <= 1.0) {
    if (yhi <= 1.0) {
      /* Segment is entirely within square */
      if (negdx) {
        return 0.5 * (xlo - xhi) * (yhi + ylo);
      } else {
        return 0.5 * (xhi - xlo) * (yhi + ylo);
      }
    }

    /* Otherwise, it must cross the top of the square */
    xtop = (1.0 - c) / m;

    if (negdx) {
      return -(0.5 * (xtop - xlo) * (1.0 + ylo) + xhi - xtop);
    } else {
      return 0.5 * (xtop - xlo) * (1.0 + ylo) + xhi - xtop;
    }
  }

  xtop = (1.0 - c) / m;

  if (negdx) {
    return -(0.5 * (xhi - xtop) * (1.0 + yhi) + xtop - xlo);
  } else {
    return 0.5 * (xhi - xtop) * (1.0 + yhi) + xtop - xlo;
  }

  /* Shouldn't ever get here */
  assert(FALSE);
  return 0.0;
}

/**
 compute area of box overlap

 Calculate the area common to input clockwise polygon x(n), y(n) with
 square (is, js) to (is+1, js+1).
 This version is for a quadrilateral.

 Used by do_square_kernel.
*/
static inline_macro double
boxer(double is, double js,
      const double x[4], const double y[4]) {
  integer_t i;
  double sum;
  double px[4], py[4];

  assert(x);
  assert(y);

  is -= 0.5;
  js -= 0.5;
  /* Set up coords relative to unit square at origin Note that the
     +0.5s were added when this code was included in DRIZZLE */

  for (i = 0; i < 4; ++i) {
    px[i] = x[i] - is;
    py[i] = y[i] - js;
  }

  /* For each line in the polygon (or at this stage, input
     quadrilateral) calculate the area common to the unit square
     (allow negative area for subsequent `vector' addition of
     subareas). */
  sum = 0.0;
  for (i = 0; i < 4; ++i) {
    sum += sgarea(px[i], py[i], px[(i+1) & 0x3], py[(i+1) & 0x3]);
  }

  return sum;
}

/**
Calculate overlap between an arbitrary rectangle, aligned with the
axes, and a pixel.

This is a simplified version of the BOXER code.

Used by do_kernel_turbo
*/
static inline_macro double
over(const integer_t i, const integer_t j,
     const double xmin, const double xmax,
     const double ymin, const double ymax) {
  double dx, dy;

  assert(xmin <= xmax);
  assert(ymin <= ymax);

  dx = MIN(xmax, (double)(i) + 0.5) - MAX(xmin, (double)(i) - 0.5);
  dy = MIN(ymax, (double)(j) + 0.5) - MAX(ymin, (double)(j) - 0.5);

  if (dx > 0.0 && dy > 0.0)
    return dx*dy;

  return 0.0;
}

/***************************************************************************
 KERNEL HANDLERS
*/

static int
do_kernel_point(struct driz_param_t* p, const integer_t j,
                const integer_t x1, const integer_t x2) {
  integer_t i, ii, jj;
  float vc, d, dow;
  double dx, dy;

  dx = 1.0;
  dy = 1.0;

  /* Offset within the subset */
  for (i = x1; i < x2; ++i) {
    ii = fortran_round(get_pixmap(p->pixmap, i, j)[0]);
    jj = fortran_round(get_pixmap(p->pixmap, i, j)[1]);

    /* Check it is on the output image */
    if (ii >= p->xmin && ii < p->xmax &&
        jj >= p->ymin && jj < p->ymax) {
      vc = get_pixel(p->output_counts, ii, jj);

      /* Allow for stretching because of scale change */
      d = get_pixel(p->data, i, j) * (float)p->scale2;

      /* Scale the weighting mask by the scale factor.  Note that we
         DON'T scale by the Jacobian as it hasn't been calculated */
      if (p->weights) {
        dow = get_pixel(p->weights, i, j) * p->weight_scale;
      } else {
        dow = 1.0;
      }

      /* If we are creating of modifying the context image,
         we do so here. */
      if (p->output_context && dow > 0.0) {
        set_bit(p->output_context, ii, jj, p->bv);
      }

      update_data(p, ii, jj, d, vc, dow);
    } else {

      ++ p->nmiss;
    }
  }

  return 0;
}

static int
do_kernel_tophat(struct driz_param_t* p, const integer_t j,
                 const integer_t x1, const integer_t x2) {
  integer_t i, ii, jj, nhit, nxi, nxa, nyi, nya;
  float vc, d, dow;
  double xx, yy, xxi, xxa, yyi, yya, dx, dy, ddx, ddy, r2;

  dx = 1.0;
  dy = 1.0;

  for (i = x1; i < x2; ++i) {
    /* Offset within the subset */
    xx = get_pixmap(p->pixmap, i, j)[0];
    yy = get_pixmap(p->pixmap, i, j)[1];

    xxi = xx - p->pfo;
    xxa = xx + p->pfo;
    yyi = yy - p->pfo;
    yya = yy + p->pfo;

    nxi = MAX(fortran_round(xxi), p->xmin);
    nxa = MIN(fortran_round(xxa), p->xmax-1);
    nyi = MAX(fortran_round(yyi), p->ymin);
    nya = MIN(fortran_round(yya), p->ymax-1);

    nhit = 0;

    /* Allow for stretching because of scale change */
    d = get_pixel(p->data, i, j) * (float)p->scale2;

    /* Scale the weighting mask by the scale factor and inversely by
       the Jacobian to ensure conservation of weight in the output */
    if (p->weights) {
      dow = get_pixel(p->weights, i, j) * p->weight_scale;
    } else {
      dow = 1.0;
    }

    /* Loop over output pixels which could be affected */
    for (jj = nyi; jj <= nya; ++jj) {
      ddy = yy - (double)jj;

      /* Check it is on the output image */
      for (ii = nxi; ii <= nxa; ++ii) {
        ddx = xx - (double)ii;

        /* Radial distance */
        r2 = ddx*ddx + ddy*ddy;

        /* Weight is one within the specified radius and zero outside.
           Note: weight isn't conserved in this case */
        if (r2 <= p->pfo2) {
          /* Count the hits */
          nhit++;
          vc = get_pixel(p->output_counts, ii, jj);

          /* If we are create or modifying the context image,
             we do so here. */
          if (p->output_context && dow > 0.0) {
            set_bit(p->output_context, ii, jj, p->bv);
          }

          update_data(p, ii, jj, d, vc, dow);
        }
      }
    }

    /* Count cases where the pixel is off the output image */
    if (nhit == 0) ++ p->nmiss;
  }

  return 0;
}

static int
do_kernel_gaussian(struct driz_param_t* p, const integer_t j,
                   const integer_t x1, const integer_t x2) {
  integer_t i, ii, jj, nxi, nxa, nyi, nya, nhit;
  float vc, d, dow;
  double xx, yy, xxi, xxa, yyi, yya, w, dx, dy, ddx, ddy, r2, dover;

  dx = 1.0;
  dy = 1.0;

  for (i = x1; i < x2; ++i) {
    xx = get_pixmap(p->pixmap, i, j)[0];
    yy = get_pixmap(p->pixmap, i, j)[1];

    xxi = xx - p->pfo;
    xxa = xx + p->pfo;
    yyi = yy - p->pfo;
    yya = yy + p->pfo;

    nxi = MAX(fortran_round(xxi), p->xmin);
    nxa = MIN(fortran_round(xxa), p->xmax-1);
    nyi = MAX(fortran_round(yyi), p->ymin);
    nya = MIN(fortran_round(yya), p->ymax-1);

    nhit = 0;

    /* Allow for stretching because of scale change */
    d = get_pixel(p->data, i, j) * (float)p->scale2;

    /* Scale the weighting mask by the scale factor and inversely by
       the Jacobian to ensure conservation of weight in the output */
    if (p->weights) {
      w = get_pixel(p->weights, i, j) * p->weight_scale;
    } else {
      w = 1.0;
    }

    /* Loop over output pixels which could be affected */
    for (jj = nyi; jj <= nya; ++jj) {
      ddy = yy - (double)jj;
      for (ii = nxi; ii <= nxa; ++ii) {
        ddx = xx - (double)ii;
        /* Radial distance */
        r2 = ddx*ddx + ddy*ddy;

        /* Weight is a scaled Gaussian function of radial
           distance */
        dover = p->gaussian.es * exp(-r2 * p->gaussian.efac);

        /* Count the hits */
        ++nhit;

        vc = get_pixel(p->output_counts, ii, jj);
        dow = (float)dover * w;

        /* If we are create or modifying the context image, we do so
           here. */
        if (p->output_context && dow > 0.0) {
          set_bit(p->output_context, ii, jj, p->bv);
        }

        update_data(p, ii, jj, d, vc, dow);
      }
    }

    /* Count cases where the pixel is off the output image */
    if (nhit == 0) ++ p->nmiss;
  }

  return 0;
}

static int
do_kernel_lanczos(struct driz_param_t* p, const integer_t j,
                  const integer_t x1, const integer_t x2) {
  integer_t i, ii, jj, nxi, nxa, nyi, nya, nhit, ix, iy;
  float vc, d, dow;
  double xx, yy, xxi, xxa, yyi, yya, w, dx, dy, dover;

  dx = 1.0;
  dy = 1.0;

  for (i = x1; i < x2; ++i) {
    xx = get_pixmap(p->pixmap, i, j)[0];
    yy = get_pixmap(p->pixmap, i, j)[1];

    xxi = xx - dx - p->pfo;
    xxa = xx - dx + p->pfo;
    yyi = yy - dy - p->pfo;
    yya = yy - dy + p->pfo;

    nxi = MAX(fortran_round(xxi), p->xmin);
    nxa = MIN(fortran_round(xxa), p->xmax-1);
    nyi = MAX(fortran_round(yyi), p->ymin);
    nya = MIN(fortran_round(yya), p->ymax-1);

    nhit = 0;

    /* Allow for stretching because of scale change */
    d = get_pixel(p->data, i, j) * (float)p->scale2;

    /* Scale the weighting mask by the scale factor and inversely by
       the Jacobian to ensure conservation of weight in the output */
    if (p->weights) {
      w = get_pixel(p->weights, i, j) * p->weight_scale;
    } else {
      w = 1.0;
    }

    /* Loop over output pixels which could be affected */
    for (jj = nyi; jj <= nya; ++jj) {
      for (ii = nxi; ii <= nxa; ++ii) {
        /* X and Y offsets */
        ix = fortran_round(fabs(xx - (double)ii) * p->lanczos.sdp) + 1;
        iy = fortran_round(fabs(yy - (double)jj) * p->lanczos.sdp) + 1;

        /* Weight is product of Lanczos function values in X and Y */
        dover = p->lanczos.lut[ix] * p->lanczos.lut[iy];

        /* Count the hits */
        ++nhit;

        /* VALGRIND REPORTS: Address is 1 bytes after a block of size
           435 */
        vc = get_pixel(p->output_counts, ii, jj);
        dow = (float)(dover * w);

        /* If we are create or modifying the context image, we do so
           here. */
        if (p->output_context && dow > 0.0) {
          set_bit(p->output_context, ii, jj, p->bv);
        }

        update_data(p, ii, jj, d, vc, dow);
      }
    }

    /* Count cases where the pixel is off the output image */
    if (nhit == 0) ++ p->nmiss;
  }

  return 0;
}

static int
do_kernel_turbo(struct driz_param_t* p, const integer_t j,
                const integer_t x1, const integer_t x2) {
  integer_t i, ii, jj, nxi, nxa, nyi, nya, nhit, iis, iie, jjs, jje;
  float vc, d, dow;
  double xxi, xxa, yyi, yya, w, dover, xoi, yoi;

  nhit = 0;

  for (i = x1; i < x2; ++i) {
    /* Offset within the subset */
    xoi = get_pixmap(p->pixmap, i, j)[0];
    yoi = get_pixmap(p->pixmap, i, j)[1];
    xxi = xoi - p->pfo;
    xxa = xoi + p->pfo;
    yyi = yoi - p->pfo;
    yya = yoi + p->pfo;

    nxi = fortran_round(xxi);
    nxa = fortran_round(xxa);
    nyi = fortran_round(yyi);
    nya = fortran_round(yya);
    iis = MAX(nxi, p->xmin);  /* Needed to be set to 0 to avoid edge effects */
    iie = MIN(nxa, p->xmax-1);
    jjs = MAX(nyi, p->ymin);  /* Needed to be set to 0 to avoid edge effects */
    jje = MIN(nya, p->ymax-1);

    nhit = 0;

    /* Allow for stretching because of scale change */
    d = get_pixel(p->data, i, j) * (float)p->scale2;

    /* Scale the weighting mask by the scale factor and inversely by
       the Jacobian to ensure conservation of weight in the output. */
    if (p->weights) {
      w = get_pixel(p->weights, i, j) * p->weight_scale;
    } else {
      w = 1.0;
    }

    /* Loop over the output pixels which could be affected */
    for (jj = jjs; jj <= jje; ++jj) {
      for (ii = iis; ii <= iie; ++ii) {
        /* Calculate the overlap using the simpler "aligned" box
           routine */
        dover = over(ii, jj, xxi, xxa, yyi, yya);

        if (dover > 0.0) {
          /* Correct for the pixfrac area factor */
          dover *= p->scale2 * p->ac;

          /* Count the hits */
          ++nhit;

          vc = get_pixel(p->output_counts, ii, jj);
          dow = (float)(dover * w);

          /* If we are create or modifying the context image,
             we do so here. */
          if (p->output_context && dow > 0.0) {
            set_bit(p->output_context, ii, jj, p->bv);
          }

          update_data(p, ii, jj, d, vc, dow);
        }
      }
    }

    /* Count cases where the pixel is off the output image */
    if (nhit == 0) ++ p->nmiss;
  }

  return 0;
}

int
do_kernel_square(struct driz_param_t* p, const integer_t j, 
                 const integer_t x1, const integer_t x2) {

  integer_t i, nhit, ii, jj, min_ii, max_ii, min_jj, max_jj, n;
  float vc, d, dow;
  double dh, jaco, tem, dover, w;
  double xyin[4][2], xyout[2], xout[4], yout[4];

  /* TODO: These are constant across calls -- perhaps cache??? */
  dh = 0.5 * p->pixel_fraction;
  n = x2 - x1;

  /* Next the "classic" drizzle square kernel...  this is different
     because we have to transform all four corners of the shrunken
     pixel */

  /* Set the input corner positions */

  xyin[0][1] = (double) j + dh;
  xyin[1][1] = (double) j + dh;
  xyin[2][1] = (double) j - dh;
  xyin[3][1] = (double) j - dh;

  for (i = x1; i < x2; ++i) {
    xyin[0][0] = (double) i - dh;
    xyin[1][0] = (double) i + dh;
    xyin[2][0] = (double) i + dh;
    xyin[3][0] = (double) i - dh;

    for (ii = 0; ii < 4; ++ii) {
      map_point(p->pixmap, xyin[ii], xyout);
      xout[ii] = xyout[0];
      yout[ii] = xyout[1];
    }

    /* Work out the area of the quadrilateral on the output grid.
       Note that this expression expects the points to be in clockwise
       order */
    
    jaco = 0.5f * ((xout[1] - xout[3]) * (yout[0] - yout[2]) -
                   (xout[0] - xout[2]) * (yout[1] - yout[3]));

    if (jaco < 0.0) {
      jaco *= -1.0;
      /* Swap */
      tem = xout[1]; xout[1] = xout[3]; xout[3] = tem;
      tem = yout[1]; yout[1] = yout[3]; yout[3] = tem;
    }

    nhit = 0;

    /* Allow for stretching because of scale change */
    d = get_pixel(p->data, i, j) * (float)p->scale2;

    /* Scale the weighting mask by the scale factor and inversely by
       the Jacobian to ensure conservation of weight in the output */
    if (p->weights) {
      w = get_pixel(p->weights, i, j) * p->weight_scale;
    } else {
      w = 1.0;
    }

    /* Loop over output pixels which could be affected */
    min_jj = MAX(fortran_round(min_doubles(yout, 4)), p->ymin);
    max_jj = MIN(fortran_round(max_doubles(yout, 4)), p->ymax-1);
    min_ii = MAX(fortran_round(min_doubles(xout, 4)), p->xmin);
    max_ii = MIN(fortran_round(max_doubles(xout, 4)), p->xmax-1);

    for (jj = min_jj; jj <= max_jj; ++jj) {
      for (ii = min_ii; ii <= max_ii; ++ii) {
        /* Call boxer to calculate overlap */
        dover = boxer((double)ii, (double)jj, xout, yout);

        if (dover > 0.0) {
          /* Re-normalise the area overlap using the Jacobian */
          dover /= jaco;

          /* Count the hits */
          ++nhit;

          vc = get_pixel(p->output_counts, ii, jj);
          dow = (float)(dover * w);

          /* If we are creating or modifying the context image we do
             so here */
          if (p->output_context && dow > 0.0) {
            set_bit(p->output_context, ii, jj, p->bv);
          }

          update_data(p, ii, jj, d, vc, dow);
        }
      }
    }

    /* Count cases where the pixel is off the output image */
    if (nhit == 0) ++ p->nmiss;
  }

  return 0;
}

static kernel_handler_t
kernel_handler_map[] = {
  do_kernel_square,
  do_kernel_gaussian,
  do_kernel_point,
  do_kernel_tophat,
  do_kernel_turbo,
  do_kernel_lanczos,
  do_kernel_lanczos
};

/**
This module does the actual mapping of input flux to output images
using "boxer", a code written by Bill Sparks for FOC geometric
distortion correction, rather than the "drizzling" approximation.

This works by calculating the positions of the four corners of a
quadrilateral on the output grid corresponding to the corners of the
input pixel and then working out exactly how much of each pixel in the
output is covered, or not.

In V1.6 this was simplified to use the DRIVAL routine and also to
include some limited multi-kernel support.
*/
int
dobox(struct driz_param_t* p) {

  const double nsig = 2.5;
  const size_t nlut = 512;
  const float del = 0.01;
  kernel_handler_t kernel_handler = NULL;
  integer_t j, np;
  float inv_exposure_time;
  int kernel_order;
  int margin;
  size_t bit_no;
  integer_t xbounds[2], ybounds[2];

  assert(p);

  /* The bitmask, trimmed to the appropriate range */
  np = (p->uuid - 1) / 32 + 1;
  bit_no = (size_t)(p->uuid - 1 - (32 * (np - 1)));
  assert(bit_no < 32);
  p->bv = (integer_t)(1 << bit_no);

  assert(p->pixel_fraction != 0.0);
  p->ac = 1.0 / (p->pixel_fraction * p->pixel_fraction);

  /* Recalculate the area scaling factor */
  p->scale2 = p->scale * p->scale;

  /* Half pixfrac on output */
  assert(p->scale != 0.0);
  p->pfo = p->pixel_fraction / p->scale / 2.0;

  switch (p->kernel) {
  /* Some Gaussian related numbers */
  case kernel_gaussian:
    p->gaussian.efac = (2.3548*2.3548) * p->scale2 * p->ac / 2.0;
    p->gaussian.es = p->gaussian.efac / M_PI;
    p->pfo = nsig * p->pixel_fraction / 2.3548 / p->scale;
    /* Added in V2.9 - make sure this doesn't get less than 1.2
       divided by the scale so that there are never holes in the
       output */
    p->pfo = CLAMP_ABOVE(p->pfo, 1.2 / p->scale);
    break;
  case kernel_lanczos2:
  case kernel_lanczos3:
    kernel_order = (p->kernel == kernel_lanczos2) ? 2 : 3;
    p->lanczos.nlut = nlut;
    assert(p->lanczos.lut == NULL);
    if ((p->lanczos.lut = malloc(nlut * sizeof(float))) == NULL) {
      driz_error_set_message(p->error, "Out of memory");
      goto dobox_exit_;
    }
    /* Set up a look-up-table for Lanczos-style interpolation
       kernels */
    create_lanczos_lut(kernel_order, nlut, del, p->lanczos.lut);
    p->pfo = (double)kernel_order * p->pixel_fraction / p->scale;
    p->lanczos.sdp = p->scale / del / p->pixel_fraction;
    break;

  default:
    break;
  }

  p->pfo2 = p->pfo*p->pfo;

  /* Set up a function pointer to handle the appropriate kernel */
  if (p->kernel >= kernel_LAST) {
    driz_error_set_message(p->error, "Invalid kernel type");
    goto dobox_exit_;
  }
  
  kernel_handler = kernel_handler_map[p->kernel];
  if (kernel_handler == NULL) {
    driz_error_set_message(p->error, "Invalid kernel type");
    goto dobox_exit_;
  }


  /* If the input image is not in CPS we need to divide by the
     exposure */
  if (p->in_units != unit_cps) {
    if (p->exposure_time == 0.0) {
      driz_error_set_message(p->error, "Invalid exposure time");
      goto dobox_exit_;
    }

    assert(p->exposure_time != 0.0);
    inv_exposure_time = 1.0f / p->exposure_time;
    scale_image(p->data, inv_exposure_time);
  }

  DRIZLOG("-Drizzling using kernel = %s\n",kernel_enum2str(p->kernel));

  /* This is the outer loop over all the lines in the input image */

  margin = 5;
  check_image_overlap(p, margin, ybounds);

  p->nskip = (p->ymax - p->ymin) - (ybounds[1] - ybounds[0]);
  p->nmiss = p->nskip * (p->ymax - p->ymin);
  
  for (j = ybounds[0]; j < ybounds[1]; ++j) {
    /* Check the overlap with the output */
    check_line_overlap(p, margin, j, xbounds);
    
    /* We know there may be some misses */
    p->nmiss += (p->xmax - p->xmin) - (xbounds[1] - xbounds[0]);
    if (xbounds[0] == xbounds[1]) ++ p->nskip;
    
    if (kernel_handler(p, j, xbounds[0], xbounds[1])) {
      goto dobox_exit_;
    }
  }

 dobox_exit_:
  free(p->lanczos.lut); p->lanczos.lut = NULL;

  return driz_error_is_set(p->error);
}
