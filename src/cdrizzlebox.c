#define NO_IMPORT_ARRAY
#define NO_IMPORT_ASTROPY_WCS_API

#include "driz_portability.h"
#include "cdrizzlemap.h"
#include "cdrizzlebox.h"
#include "cdrizzlewcs.h"
#include "cdrizzleutil.h"

#include <assert.h>
#define _USE_MATH_DEFINES       /* needed for MS Windows to define M_PI */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static inline_macro double*
mapping_4_ptr(struct driz_param_t* p, double* arr, integer_t i0, integer_t i1) {
  assert(p);
  assert(arr);
  assert(i0 >= 0 && i0 <= p->dnx);
  assert(i1 >= 0 && i1 < 4);
  return (arr + ((i1 * p->dnx) + i0));
}

static inline_macro double*
mapping_ptr(struct driz_param_t* p, double* arr, integer_t i0) {
  assert(p);
  assert(arr);
  assert(i0 >= 0 && i0 <= p->dnx);
  return (arr + i0);
}

/**
Check how much of a line will overlap an output image, if any,
after applying the standard drizzle transformation.

This is intended to allow the number of points which are needlessly
drizzled outside the output image to be minimized.

was: CHOVER
*/
#define CHECK_OVER_NPOINT 21

static int
check_over(struct driz_param_t* p, const integer_t y, const integer_t margin,
           /* Output parameters */
           double* ofrac, integer_t* x1, integer_t* x2,
           struct driz_error_t* error) {
  const integer_t npoint = CHECK_OVER_NPOINT;
  double xval[CHECK_OVER_NPOINT], yval[CHECK_OVER_NPOINT];
  double xtmp[CHECK_OVER_NPOINT], ytmp[CHECK_OVER_NPOINT];
  double xout[CHECK_OVER_NPOINT], yout[CHECK_OVER_NPOINT];
  integer_t logo[CHECK_OVER_NPOINT];
  integer_t step, first, last;
  integer_t nhit, nmiss;
  integer_t i, np;

  assert(p);
  assert(ofrac);
  assert(x1);
  assert(x2);
  assert(error);

  if (p->dnx < npoint)
    step = 1;
  else
    step = p->dnx / (npoint / 2);

  for (i = 1, np = 0; i <= p->dnx; i += step, ++np) {
    assert(np < npoint);
    xval[np] = (double)i;
    yval[np] = (double)y;
  }
  assert(np < npoint);

  /* Check end point */
  if (xval[np - 1] < (double)p->dnx) {
    xval[np] = (double)p->dnx;
    yval[np] = (double)y;
    ++np;
  }

  /* Transform */
  if (map_value(p, FALSE, np, xval, yval, xtmp, ytmp, xout, yout, error))
    return 1;

  /* Check where the overlap starts and ends */
  for (i = 0; i < np; ++i) {
    logo[i] = 0;
  }

  for (i = 0; i < np - 1; ++i) {
    if (MAX(xout[i], xout[i+1]) >= 1.0 - (double)margin &&
        MIN(xout[i], xout[i+1]) < (double)(p->onx + margin) &&
        MAX(yout[i], yout[i+1]) >= 1.0 - (double)margin &&
        MIN(yout[i], yout[i+1]) < (double)(p->ony + margin)) {
      logo[i] = 1;
      logo[i+1] = 1;
    }
  }

  nhit = 0;
  for (i = 0; i < np; ++i) {
    if (logo[i]) {
      ++nhit;
    }
  }

  nmiss = np - nhit;

  if (nhit == 0) {
    *ofrac = 0.0;
    *x1 = 0;
    *x2 = 0;
    return 0;
  }

  first = 0;
  last = 0;

  for (i = 0; i < np; ++i) {
    if (logo[i]) {
      first = i;
      break;
    }
  }

  for (i = np - 1; i >= 0; --i) {
    if (logo[i]) {
      last = i;
      break;
    }
  }

  if (nhit == 0 && nmiss == 0) {
    *ofrac = 0.0;
  } else {
    *ofrac = (double)nhit / (double)(nhit + nmiss);
  }
  *x1 = (first > 0) ? (integer_t)xval[first - 1] : (integer_t)xval[0];
  *x2 = (last < np - 1) ? (integer_t)xval[last + 1] : (integer_t)xval[np - 1];

  assert(*x1 > 0 && *x1 <= p->dnx);
  assert(*x2 > 0 && *x2 <= p->dnx);

  return 0;
}

/**
Match up a context against a context table Note that after drizzle
V0.40 (EIS) the context table itself has an extra second column for
the counter.
*/
static bool_t
match(const integer_t* a /*[n]*/, const integer_t* b /*[n]*/,
      const integer_t n UNUSED_PARAM) {
  integer_t i, an, bn;

  assert(a);
  assert(b);

  an = a[0];
  bn = b[0];

  if (an != bn) {
    return 0;
  } else {
    assert(an <= n - 2);
    for (i = 1; i <= an; ++i) {
      if (a[i] != b[i+1])
        return 0;
    }
  }

  return 1;
}

/* A comparision predicate for use with qsort, used by update_context */
static int
_compare_integer_t(const void* a, const void* b) {
  integer_t ad = *(integer_t*)a;
  integer_t bd = *(integer_t*)b;
  if (ad < bd)
    return -1;
  if (ad == bd)
    return 0;
  return 1;
}

/**
Update the context image.

This routine is called in the heart of Drizzle when the context image
needs to be updated.

was: UPCON
*/
static int
update_context_image(struct driz_param_t* p, const integer_t ii, const integer_t jj,
                     /* Input/output parameters */
                     integer_t* oldcon,
                     /* Output parameters */
                     integer_t* newcon, struct driz_error_t* error) {
  integer_t icon, nn, k;
  integer_t newma[100];

  assert(p);
  assert(oldcon);
  assert(newcon);
  assert(error);
  assert(ii >= 0 && ii < p->nsx);
  assert(jj >= 0 && jj < p->nsy);

  /* Look up the current context value */
  icon = *output_context_ptr(p, ii, jj);

  /* If it is the same as the last one, we don't need to go further */
  if (icon == *oldcon) {
    *output_context_ptr(p, ii, jj) = *newcon;
    goto upcon_continue;
  }

  /* Combine with the new one */
  if (icon == 0) {
    nn = 0;
    newma[0] = 0;
    newma[1] = p->uuid;
  } else {
    nn = *intab_ptr(p, 0, icon);

    /* Check for whether this image is already on this pixel */
    for (k = 2; k <= nn+1; ++k)
      if (p->uuid == *intab_ptr(p, k, icon))
        goto upcon_continue;

    /* If not create the new context by adding the current image */
    newma[0] = nn + 1;
    for (k = 1; k <= nn; ++k)
      newma[k] = *intab_ptr(p, k+1, icon);

    /* Check for too many images at a given context */
    if (nn > MAXIM - 3) {
      driz_error_set_message(error, "Too many images.  Context table overloaded.");
      return 1;
    }

    newma[nn+2] = p->uuid; /* TODO: Base check */
  }

  /* Before matching sort the context array */
  if (nn > 0) {
    qsort(newma, (size_t)nn, sizeof(integer_t), &_compare_integer_t);
  }

  /* See whether we have had this context before */
  for (k = p->nen; k >= 0; --k) {
    if (match(newma, intab_ptr(p, 0, k), MAXIM)) {
      *output_context_ptr(p, ii, jj) = k;
      goto upcon_continue;
    }
  }

  /* No context match found: make a new one */
  ++(p->nen);

  /* Check for full table */
  if (p->nen == MAXEN) {
    driz_error_set_message(error, "Context table full");
    return 1;
  }

  *output_context_ptr(p, ii, jj) = p->nen;
  *intab_ptr(p, 0, p->nen) = newma[0];
  *intab_ptr(p, 1, p->nen) = 0;
  for (k = 2; k < nn+4; ++k) { /* TODO: Base check */
    *intab_ptr(p, k, p->nen) = newma[k - 1];
  }

  /* Tell the user about this new context */
  /* TODO */

 upcon_continue:
  /* Save the old values for quick comparison */
  *oldcon = icon;
  *newcon = *output_context_ptr(p, ii, jj);

  /* Lastly, we update the counter */
  if (*oldcon != *newcon) {
    if (*oldcon > 0) {
      *intab_ptr(p, 1, *oldcon) -= 1;
    }
    *intab_ptr(p, 1, *newcon) += 1;
  }

  *output_done_ptr(p, ii, jj) = 1;

  return 0;
}

inline_macro static int
update_context(struct driz_param_t* p, const integer_t ii, const integer_t jj,
               const double dow,
               /* Input/output parameters */
               integer_t* oldcon,
               /* Output parameters */
               integer_t* newcon, struct driz_error_t* error) {
  if (p->output_context && dow > 0.0) {
    if (p->output_done == NULL) {
      *output_context_ptr(p, ii, jj) |= p->bv;
    } else {
      if (*output_done_ptr(p, ii, jj) == 0) {
        /* TODO: The error case seems to be ignored here in original */
        if (update_context_image(p, ii, jj, oldcon, newcon, error))
          return 1;
      }
    }
  }

  return 0;
}

inline_macro static void
update_data(struct driz_param_t* p, const integer_t ii, const integer_t jj,
            const float d, const float vc, const float dow) {
  const double vc_plus_dow = vc + dow;

  /* Just a simple calculation without logical tests */
  if (vc == 0.0) {
    *output_data_ptr(p, ii, jj) = d;
  } else {
    if (vc_plus_dow != 0.0) {
      *output_data_ptr(p, ii, jj) =
        (*output_data_ptr(p, ii, jj) * vc + dow * d) / (vc_plus_dow);
    }
  }

  *output_counts_ptr(p, ii, jj) = vc_plus_dow;
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

typedef int (*kernel_handler_t)(struct driz_param_t*, const integer_t,
                                const integer_t, const integer_t,
                                double*, double*,
                                integer_t*, integer_t*, integer_t*,
                                struct driz_error_t*);

static int
do_kernel_point(struct driz_param_t* p, const integer_t j,
                const integer_t x1, const integer_t x2,
                double* xo, double* yo,
                /* Input/output parameters */
                integer_t* oldcon, integer_t* newcon, integer_t* nmiss,
                /* Output parameters */
                struct driz_error_t* error) {
  integer_t i, ii, jj;
  float vc, d, dow;
  double dx, dy;
  integer_t xarr,yarr;

  dx = (double)(p->xmin);
  dy = (double)(p->ymin);

  /* Offset within the subset */
  for (i = x1; i <= x2; ++i) {
    ii = fortran_round(*mapping_ptr(p, xo, i) - dx);
    jj = fortran_round(*mapping_ptr(p, yo, i) - dy);

    /* Check it is on the output image */
    if (ii >= 0 && ii < p->nsx &&
        jj >= 0 && jj < p->nsy) {
      vc = *output_counts_ptr(p, ii, jj);
    /* Convert i,j 1-based pixel positions into 0-based
       indices for accessing data array. */
    xarr = i-1;
    yarr = j-1;

      /* Allow for stretching because of scale change */
      d = *data_ptr(p, xarr, yarr) * (float)p->scale2;

      /* Scale the weighting mask by the scale factor.  Note that we
         DON'T scale by the Jacobian as it hasn't been calculated */
      if (p->weights) {
        dow = *weights_ptr(p, xarr, yarr) * p->weight_scale;
      } else {
        dow = 1.0;
      }

      /* If we are creating of modifying the context image,
         we do so here. */
      if (update_context(p, ii, jj, dow, oldcon, newcon, error)) {
        return 1;
      }

      update_data(p, ii, jj, d, vc, dow);
    } else {

      ++(*nmiss);
    }
  }

  return 0;
}

static int
do_kernel_tophat(struct driz_param_t* p, const integer_t j,
                 const integer_t x1, const integer_t x2,
                 double* xo, double* yo,
                 /* Input/output parameters */
                 integer_t* oldcon, integer_t* newcon, integer_t* nmiss,
                 struct driz_error_t* error) {
  integer_t i, ii, jj, nhit, nxi, nxa, nyi, nya;
  float vc, d, dow;
  double xx, yy, xxi, xxa, yyi, yya, dx, dy, ddx, ddy, r2;
  integer_t xarr,yarr;

  dx = (double)(p->xmin);
  dy = (double)(p->ymin);

  for (i = x1; i <= x2; ++i) {
    /* Offset within the subset */
    xx = *mapping_ptr(p, xo, i) - dx;
    yy = *mapping_ptr(p, yo, i) - dy;

    xxi = xx - p->pfo;
    xxa = xx + p->pfo;
    yyi = yy - p->pfo;
    yya = yy + p->pfo;

    nxi = MAX(fortran_round(xxi), 0);
    nxa = MIN(fortran_round(xxa), p->nsx - 1);
    nyi = MAX(fortran_round(yyi), 0);
    nya = MIN(fortran_round(yya), p->nsy - 1);

    nhit = 0;
    /* Convert i,j 1-based pixel positions into 0-based
       indices for accessing data array. */
    xarr = i-1;
    yarr = j-1;

    /* Allow for stretching because of scale change */
    d = *data_ptr(p, xarr, yarr) * (float)p->scale2;

    /* Scale the weighting mask by the scale factor and inversely by
       the Jacobian to ensure conservation of weight in the output */
    if (p->weights) {
      dow = *weights_ptr(p, xarr, yarr) * p->weight_scale;
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
          vc = *output_counts_ptr(p, ii, jj);

          /* If we are create or modifying the context image,
             we do so here. */
          if (update_context(p, ii, jj, dow, oldcon, newcon, error)) {
            return 1;
          }

          update_data(p, ii, jj, d, vc, dow);
        }
      }
    }

    /* Count cases where the pixel is off the output image */
    if (nhit == 0) ++(*nmiss);
  }

  return 0;
}

static int
do_kernel_gaussian(struct driz_param_t* p, const integer_t j,
                   const integer_t x1, const integer_t x2,
                   double* xo, double* yo,
                   /* Input/output parameters */
                   integer_t* oldcon, integer_t* newcon, integer_t* nmiss,
                   struct driz_error_t* error) {
  integer_t i, ii, jj, nxi, nxa, nyi, nya, nhit;
  float vc, d, dow;
  double xx, yy, xxi, xxa, yyi, yya, w, dx, dy, ddx, ddy, r2, dover;
  integer_t xarr,yarr;

  dx = (double)(p->xmin);
  dy = (double)(p->ymin);

  for (i = x1; i <= x2; ++i) {
    xx = *mapping_ptr(p, xo, i) - dx;
    yy = *mapping_ptr(p, yo, i) - dy;

    xxi = xx - p->pfo;
    xxa = xx + p->pfo;
    yyi = yy - p->pfo;
    yya = yy + p->pfo;

    nxi = MAX(fortran_round(xxi), 0);
    nxa = MIN(fortran_round(xxa), p->nsx - 1);
    nyi = MAX(fortran_round(yyi), 0);
    nya = MIN(fortran_round(yya), p->nsy - 1);

    nhit = 0;
    /* Convert i,j 1-based pixel positions into 0-based
       indices for accessing data array. */
    xarr = i-1;
    yarr = j-1;


    /* Allow for stretching because of scale change */
    d = *data_ptr(p, xarr, yarr) * (float)p->scale2;

    /* Scale the weighting mask by the scale factor and inversely by
       the Jacobian to ensure conservation of weight in the output */
    if (p->weights) {
      w = *weights_ptr(p, xarr, yarr) * p->weight_scale;
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

        vc = *output_counts_ptr(p, ii, jj);
        dow = (float)dover * w;

        /* If we are create or modifying the context image, we do so
           here. */
        if (update_context(p, ii, jj, dow, oldcon, newcon, error)) {
          return 1;
        }

        update_data(p, ii, jj, d, vc, dow);
      }
    }

    /* Count cases where the pixel is off the output image */
    if (nhit == 0) ++(*nmiss);
  }

  return 0;
}

static int
do_kernel_lanczos(struct driz_param_t* p, const integer_t j,
                  const integer_t x1, const integer_t x2,
                  double* xo, double *yo,
                  /* Input/output parameters */
                  integer_t* oldcon, integer_t* newcon, integer_t* nmiss,
                  struct driz_error_t* error) {
  integer_t i, ii, jj, nxi, nxa, nyi, nya, nhit, ix, iy;
  float vc, d, dow;
  double xx, yy, xxi, xxa, yyi, yya, w, dx, dy, dover;
  integer_t xarr,yarr;

  dx = (double)(p->xmin);
  dy = (double)(p->ymin);

  for (i = x1; i <= x2; ++i) {
    xx = *mapping_ptr(p, xo, i) - dx;
    yy = *mapping_ptr(p, yo, i) - dy;

    xxi = xx - p->pfo;
    xxa = xx + p->pfo;
    yyi = yy - p->pfo;
    yya = yy + p->pfo;

    nxi = MAX(fortran_round(xxi), 0);
    nxa = MIN(fortran_round(xxa), p->nsx - 1);
    nyi = MAX(fortran_round(yyi), 0);
    nya = MIN(fortran_round(yya), p->nsy - 1);

    nhit = 0;
    /* Convert i,j 1-based pixel positions into 0-based
       indices for accessing data array. */
    xarr = i-1;
    yarr = j-1;


    /* Allow for stretching because of scale change */
    d = *data_ptr(p, xarr, yarr) * (float)p->scale2;

    /* Scale the weighting mask by the scale factor and inversely by
       the Jacobian to ensure conservation of weight in the output */
    if (p->weights) {
      w = *weights_ptr(p, xarr, yarr) * p->weight_scale;
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
        vc = *output_counts_ptr(p, ii, jj);
        dow = (float)(dover * w);

        /* If we are create or modifying the context image, we do so
           here. */
        if (update_context(p, ii, jj, dow, oldcon, newcon, error)) {
          return 1;
        }

        update_data(p, ii, jj, d, vc, dow);
      }
    }

    /* Count cases where the pixel is off the output image */
    if (nhit == 0) ++(*nmiss);
  }

  return 0;
}

static int
do_kernel_turbo(struct driz_param_t* p, const integer_t j,
                const integer_t x1, const integer_t x2,
                double* xo, double *yo,
                /* Input/output parameters */
                integer_t* oldcon, integer_t* newcon, integer_t* nmiss,
                struct driz_error_t* error) {
  integer_t i, ii, jj, nxi, nxa, nyi, nya, nhit, iis, iie, jjs, jje;
  float vc, d, dow;
  double xxi, xxa, yyi, yya, w, dx, dy, dover,xoi,yoi;
  integer_t xarr,yarr;

  dx = (double)(p->xmin);
  dy = (double)(p->ymin);

  nhit = 0;

  for (i = x1; i <= x2; ++i) {
    /* Offset within the subset */
    xoi = *mapping_ptr(p, xo, i);
    yoi = *mapping_ptr(p, yo, i);
    xxi = xoi - dx - p->pfo;
    xxa = xoi - dx + p->pfo;
    yyi = yoi - dy - p->pfo;
    yya = yoi - dy + p->pfo;

    nxi = fortran_round(xxi);
    nxa = fortran_round(xxa);
    nyi = fortran_round(yyi);
    nya = fortran_round(yya);
    iis = MAX(nxi, 0);  /* Needed to be set to 0 to avoid edge effects */
    iie = MIN(nxa, p->nsx - 1);
    jjs = MAX(nyi, 0);  /* Needed to be set to 0 to avoid edge effects */
    jje = MIN(nya, p->nsy - 1);

    nhit = 0;

    /* Convert i,j 1-based pixel positions into 0-based
       indices for accessing data array. */
    xarr = i-1;
    yarr = j-1;

    /* Allow for stretching because of scale change */
    d = *data_ptr(p, xarr, yarr) * (float)p->scale2;

    /* Scale the weighting mask by the scale factor and inversely by
       the Jacobian to ensure conservation of weight in the output. */
    if (p->weights) {
      w = *weights_ptr(p, xarr, yarr) * p->weight_scale;
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

          vc = *output_counts_ptr(p, ii, jj);
          dow = (float)(dover * w);

          /* If we are create or modifying the context image,
             we do so here. */
          if (update_context(p, ii, jj, dow, oldcon, newcon, error)) {
            return 1;
          }

          update_data(p, ii, jj, d, vc, dow);
        }
      }
    }

    /* Count cases where the pixel is off the output image */
    if (nhit == 0) ++(*nmiss);
  }

  return 0;
}

static int
do_kernel_square(struct driz_param_t* p,
                 const integer_t j, double y,
                 const integer_t x1, const integer_t x2,
                 const integer_t last_x1, const integer_t last_x2,
                 /* Input/output parameters */
                 double* xi, double* yi,
                 double* xtmp, double* ytmp,
                 double* xo, double* yo,
                 integer_t* oldcon, integer_t* newcon, integer_t* nmiss,
                 struct driz_error_t* error) {
  integer_t i, nhit, ii, jj, min_ii, max_ii, min_jj, max_jj, n;
  float vc, d, dow;
  double dh, jaco, tem, dover, dx, dy, w;
  double xout[4], yout[4];

  /* TODO: These are constant across calls -- perhaps cache??? */
  dh = 0.5 * p->pixel_fraction;
  dx = (double)(p->xmin) - 1;
  dy = (double)(p->ymin) - 1;
  n = x2 - x1 + 1;

  /* Next the "classic" drizzle square kernel...  this is different
     because we have to transform all four corners of the shrunken
     pixel */
  /* Set the start corner positions */

  *mapping_4_ptr(p, xi, x1, 0) = (double)x1 - dh;
  *mapping_4_ptr(p, xi, x1, 1) = (double)x1 + dh;
  *mapping_4_ptr(p, xi, x1, 2) = (double)x1 + dh;
  *mapping_4_ptr(p, xi, x1, 3) = (double)x1 - dh;

  *mapping_4_ptr(p, yi, x1, 0) = y + dh;
  *mapping_4_ptr(p, yi, x1, 1) = y + dh;
  *mapping_4_ptr(p, yi, x1, 2) = y - dh;
  *mapping_4_ptr(p, yi, x1, 3) = y - dh;

  *mapping_4_ptr(p, yi, x1+1, 0) = dh;
  *mapping_4_ptr(p, yi, x1+1, 1) = dh;
  *mapping_4_ptr(p, yi, x1+1, 2) = -dh;
  *mapping_4_ptr(p, yi, x1+1, 3) = -dh;

  /* Transform onto the output grid */
  for (i = 0; i < 4; ++i) {
    if (map_value(p, TRUE, n,
                  mapping_4_ptr(p, xi, x1, i), mapping_4_ptr(p, yi, x1, i),
                  xtmp, ytmp,
                  mapping_4_ptr(p, xo, x1, i), mapping_4_ptr(p, yo, x1, i),
                  error)) {
      return 1;
    }
  }

  for (i = x1; i <= x2; ++i) {
    /* Offset within the subset */
    for (ii = 0; ii < 4; ++ii) {
      /* The offset by 1 here is needed to match the alignment in the
         output frame generated by the other kernels (such as turbo).
      */
      xout[ii] = *mapping_4_ptr(p, xo, i, ii) - dx - 1;
      yout[ii] = *mapping_4_ptr(p, yo, i, ii) - dy - 1;
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
    d = *data_ptr(p, i-1, j) * (float)p->scale2;

    /* Scale the weighting mask by the scale factor and inversely by
       the Jacobian to ensure conservation of weight in the output */
    if (p->weights) {
      w = *weights_ptr(p, i-1, j) * p->weight_scale;
    } else {
      w = 1.0;
    }

    /* Loop over output pixels which could be affected */
    min_jj = MAX(fortran_round(min_doubles(yout, 4)), 0);
    max_jj = MIN(fortran_round(max_doubles(yout, 4)), p->nsy - 1);
    min_ii = MAX(fortran_round(min_doubles(xout, 4)), 0);
    max_ii = MIN(fortran_round(max_doubles(xout, 4)), p->nsx - 1);

    for (jj = min_jj; jj <= max_jj; ++jj) {
      for (ii = min_ii; ii <= max_ii; ++ii) {
        /* Call boxer to calculate overlap */
        dover = boxer((double)ii, (double)jj, xout, yout);

        if (dover > 0.0) {
          /* Re-normalise the area overlap using the Jacobian */
          dover /= jaco;

          /* Count the hits */
          ++nhit;

          vc = *output_counts_ptr(p, ii, jj);
          dow = (float)(dover * w);

          /* If we are creating or modifying the context image we do
             so here */
          if (update_context(p, ii, jj, dow, oldcon, newcon, error)) {
            return 1;
          }

          update_data(p, ii, jj, d, vc, dow);
        }
      }
    }
    /* Count cases where the pixel is off the output image */
    if (nhit == 0) ++(*nmiss);
  }

  return 0;
}

static kernel_handler_t
kernel_handler_map[] = {
  NULL,
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
dobox(struct driz_param_t* p, const integer_t ystart,
      /* Output parameters */
      integer_t* nmiss, integer_t* nskip, struct driz_error_t* error) {
  const double nsig = 2.5;
  const size_t nlut = 512;
  const float del = 0.01;
  integer_t j, x1, x2, last_x1, last_x2;
  double y, dh, ofrac;
  kernel_handler_t kernel_handler = NULL;
  integer_t oldcon, newcon;
  integer_t np;
  double* xi = NULL;
  double* yi = NULL;
  double* xtmp = NULL;
  double* ytmp = NULL;
  double* xo = NULL;
  double* yo = NULL;
  float inv_exposure_time;
  float* data_begin, *data_end;
  int kernel_order;
  size_t new_buffer_size;
  size_t bit_no;

  assert(p);
  assert(nmiss);
  assert(nskip);
  assert(error);

  /* We skip all this if there is no overlap */
  if (p->no_over) {
    /* If there is no overlap at all, set appropriate values */
    *nskip = p->dny;
    *nmiss = p->dnx * p->dny;
    return 0;
  }

  /* Some initial settings - note that the reference pixel position is
     determined by the value of ALIGN */
  oldcon = -1;

  /* The bitmask, trimmed to the appropriate range */
  np = (p->uuid - 1) / 32 + 1;
  bit_no = (size_t)(p->uuid - 1 - (32 * (np - 1)));
  assert(bit_no < 32);
  p->bv = (integer_t)(1 << bit_no);

  /* Before we start we can fill the X arrays as they don't change
     with Y */
  new_buffer_size = (size_t)((p->kernel == kernel_square) ? p->dnx*4 : p->dnx);

  /* Image subset size */
  p->nsx = p->xmax - p->xmin + 1;
  p->nsy = p->ymax - p->ymin + 1;
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
      driz_error_set_message(error, "Out of memory");
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

  /* assert(p->output_done == NULL); */
  /* p->output_done = malloc(p->nsx * p->nsy * sizeof(integer_t)); */
  /* if (p->output_done == NULL) { */
  /*   driz_error_set_message(error, "Out of memory"); */
  /*   goto dobox_exit_; */
  /* } */
  /* for (i = 0; i < p->nsx * p->nsy; ++i) { */
  /*   p->output_done[i] = 0; */
  /* } */

  xi = malloc(new_buffer_size * sizeof(double));
  if (xi == NULL) {
    driz_error_set_message(error, "Out of memory");
    goto dobox_exit_;
  }

  yi = malloc(new_buffer_size * sizeof(double));
  if (yi == NULL) {
    driz_error_set_message(error, "Out of memory");
    goto dobox_exit_;
  }

  xtmp = malloc(new_buffer_size * sizeof(double));
  if (xtmp == NULL) {
    driz_error_set_message(error, "Out of memory");
    goto dobox_exit_;
  }

  ytmp = malloc(new_buffer_size * sizeof(double));
  if (ytmp == NULL) {
    driz_error_set_message(error, "Out of memory");
    goto dobox_exit_;
  }

  xo = malloc((new_buffer_size + 1) * sizeof(double));
  if (xo == NULL) {
    driz_error_set_message(error, "Out of memory");
    goto dobox_exit_;
  }

  yo = malloc((new_buffer_size + 1) * sizeof(double));
  if (yo == NULL) {
    driz_error_set_message(error, "Out of memory");
    goto dobox_exit_;
  }

  if (p->kernel == kernel_square) {
    dh = 0.5 * p->pixel_fraction;
    *mapping_4_ptr(p, xi, 1, 0) = 1.0 - dh;
    *mapping_4_ptr(p, xi, 1, 1) = 1.0 + dh;
    *mapping_4_ptr(p, xi, 1, 2) = 1.0 + dh;
    *mapping_4_ptr(p, xi, 1, 3) = 1.0 - dh;
  } else {
    *mapping_ptr(p, xi, 0) = 1.0;

    /* Set up a function pointer to handle the appropriate kernel */
    if (p->kernel >= kernel_LAST) {
      driz_error_set_message(error, "Invalid kernel type");
      goto dobox_exit_;
    }
    kernel_handler = kernel_handler_map[p->kernel];
    if (kernel_handler == NULL) {
      driz_error_set_message(error, "Invalid kernel type");
      goto dobox_exit_;
    }
  }

  /* If the input image is not in CPS we need to divide by the
     exposure */
  if (p->in_units != unit_cps) {
    if (p->exposure_time == 0.0) {
      driz_error_set_message(error, "Invalid exposure time");
      goto dobox_exit_;
    }
    assert(p->exposure_time != 0.0);
    inv_exposure_time = 1.0f / p->exposure_time;
    /* TODO: Removing this printf causes the results to be
       less accurate.  Frustrating Heisenbug */
    /*printf("%f\n", inv_exposure_time); */
    data_begin = p->data;
    data_end = data_begin + (p->ny * p->dnx);
    for (; data_begin != data_end; ++data_begin) {
      *data_begin *= inv_exposure_time;
    }
  }

  DRIZLOG("-Drizzling using kernel = %s\n",kernel_enum2str(p->kernel));

  /* This is the outer loop over all the lines in the input image */
  last_x1 = p->dnx;
  last_x2 = 0;
  y = (double)ystart;
  for (j = 0; j < p->ny; ++j) {
    y += 1.0;
    /* Check the overlap with the output */
    if (check_over(p, (integer_t)y, 5, &ofrac, &x1, &x2, error)) {
      goto dobox_exit_;
    }

    /* If the line falls completely off the output, then skip it */
    if (ofrac != 0.0) {
      assert(x1 > 0 && x1 <= p->dnx);
      assert(x2 > 0 && x2 <= p->dnx);

      /* We know there may be some misses */
      *nmiss += p->dnx - (x2 - x1 + 1);

      /* Don't read past the edge of the image
      if (x2 == p->dnx) {
          x2 -= 1;
      }
      */
      /* At this point we can handle the different kernels separately.
         First the cases where we just transform a single point rather
         than four - every case except the "classic" square-pixel
         kernel */
      if (p->kernel != kernel_square) {
        *mapping_ptr(p, xi, x1) = (double)x1;

        *mapping_ptr(p, yi, x1) = y;
        *mapping_ptr(p, yi, x1+1) = 0.0;


        if (map_value(p, TRUE, x2 - x1 + 1,
                      mapping_ptr(p, xi, x1), mapping_ptr(p, yi, x1),
                      xtmp, ytmp,
                      mapping_ptr(p, xo, x1), mapping_ptr(p, yo, x1), error)) {
          goto dobox_exit_;
        }

        if (kernel_handler(p, y, x1, x2, xo, yo,
                           &oldcon, &newcon, nmiss, error)) {
          goto dobox_exit_;
        }
      } else {
        if (do_kernel_square(p, j, y, x1, x2, last_x1, last_x2,
                             xi, yi, xtmp, ytmp, xo, yo,
                             &oldcon, &newcon, nmiss, error)) {
          goto dobox_exit_;
        }
      }
      last_x1 = x1;
      last_x2 = x2;
    } else {
      /* If we are skipping a line, count it */
      ++(*nskip);
      *nmiss += p->dnx;
      last_x1 = p->dnx;
      last_x2 = 0;
    }
  }

 dobox_exit_:
  free(p->lanczos.lut); p->lanczos.lut = NULL;
  free(p->output_done); p->output_done = NULL;
  free(xi); xi = NULL;
  free(yi); yi = NULL;
  free(xo); xo = NULL;
  free(yo); yo = NULL;
  free(xtmp); xtmp = NULL;
  free(ytmp); ytmp = NULL;

  return driz_error_is_set(error);
}
