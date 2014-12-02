#include <assert.h>
#include <math.h>
#include <string.h>
#include <Python.h>
#include <numpy/arrayobject.h>

#include "driz_portability.h"
#include "cdrizzlemap.h"
#include "cdrizzleutil.h"

/* Line segment structure, used for computing overlap
 * The first index on line is the endpoint
 * The second index is {x, y) coordinate of the point
 * The valid flag is non-zero if it does not intersect the image
 */

struct segment {
    double  point[2][2];
    int     side;
};

static inline_macro void
initialize_segment(struct segment *self, double x1, double y1, double x2, double y2) {
  self->point[0][0] = x1;
  self->point[0][1] = y1;
  self->point[1][0] = x2;
  self->point[1][1] = y2;
  self->side = 0;

  return;
}

void
map_point(PyArrayObject *pixmap,
          const double xyin[2],
          double xyout[2]
         ) {

  int        k;
  integer_t  pix[2];
  double     frac[2];
  integer_t  xydim[2];

  /* Divide the position into an integral pixel position
   * plus a fractional offset within the pixel */

  get_dimensions(pixmap, xydim);

  for (k = 0; k < 2; ++k) {
    frac[k] = xyin[k];
    pix[k]  = frac[k];
    pix[k]  = CLAMP(pix[k], 0, xydim[k] - 2);
    frac[k] = frac[k] - pix[k];
  }
  
  if (frac[0] == 0.0 && frac[1] == 0.0) {
    /* No interpolation needed if input position has no fractional part */
    for (k = 0; k < 2; ++k) {
      xyout[k] = get_pixmap(pixmap, pix[0], pix[1])[k];
    }

  } else {
    /* Bilinear interpolation btw pixels, see Wikipedia for explanation */
    for (k = 0; k < 2; ++k) {
      xyout[k] = (1.0 - frac[0]) * (1.0 - frac[1]) * get_pixmap(pixmap, pix[0], pix[1])[k] +
                 frac[0] * (1.0 - frac[1]) * get_pixmap(pixmap, pix[0]+1, pix[1])[k] +
                 (1.0 - frac[0]) * frac[1] * get_pixmap(pixmap, pix[0], pix[1]+1)[k] +
                 frac[0] * frac[1] * get_pixmap(pixmap, pix[0]+1, pix[1]+1)[k];
    }
  }

  return;
}


/* Clip a line segment to the limits of an image */

void
clip_bounds(PyArrayObject *pixmap, int jdim, struct segment *xylimit, struct segment *xybounds) {
    int ipoint, idim;
    
    xybounds->side = 1; /* Track if bounds are both outside the image */
    
    for (ipoint = 0; ipoint < 2; ++ipoint) {
        int m = 21;         /* maximum iterations */
        int side = 0;       /* flag indicating which side moved last */

        double xyin[2], xyout[2];
        double a, b, c, fa, fb, fc;
        
        /* starting values at endpoints of interval */
    
        for (idim = 0; idim < 2; ++idim) {
          xyin[idim] = xybounds->point[0][idim];
        }
        
        map_point(pixmap, xyin, xyout);
        a = xybounds->point[0][jdim];
        fa = xyout[jdim] - xylimit->point[ipoint][jdim];
        
        for (idim = 0; idim < 2; ++idim) {
          xyin[idim] = xybounds->point[1][idim];
        }

        map_point(pixmap, xyin, xyout);
        c = xybounds->point[1][jdim];
        fc = xyout[jdim] - xylimit->point[ipoint][jdim];

        /* Solution via the method of false position (regula falsi) */

        if (fa * fc < 0.0) {
            int n; /* for loop limit is just for safety's sake */
            
            for (n = 0; n < m; n++) {
                b = (fa * c - fc * a) / (fa - fc);
                
                /* Solution is exact if within a pixel because linear interpolation */
                if (floor(a) == floor(c)) break;
                
                xyin[jdim] = b;
                map_point(pixmap, xyin, xyout);
                fb = xyout[jdim] - xylimit->point[ipoint][jdim];

                /* Maintain the bound by copying b to the variable
                 * with the same sign as b
                 */
                
                if (fb * fc > 0.0) {
                    c = b;
                    fc = fb;
                    if (side == -1) {
                        fa /= 2.0;
                    }
                    side = -1;
        
                } else if (fa * fb > 0.0) {
                    a = b;
                    fa = fb;
                    if (side == +1) {
                        fc /= 2.0;
                    }
                    side = +1;
                    
                } else {
                    /* if the product is zero, we have converged */
                    break;
                } 
            }

            assert(n < m);
            xybounds->side = 0;
            xybounds->point[ipoint][jdim] = b;

        } else {
            /* No bracket, so track which side the bound lies on */
            xybounds->side *= fa > 0.0 ? +1 : -1;
        }
    }

    if (xybounds->side > 0) {
        /* Positive means both bounds are outside the image */
        xybounds->point[1][jdim] = xybounds->point[0][jdim]; 
    } else {
        /* Negative means both bounds are inside, which is not a problem */
        xybounds->side = 0;
    }
    
    return;
}

void
check_line_overlap(PyArrayObject *pixmap, PyArrayObject *output_data,
                   const int margin, const integer_t iy, integer_t *xbounds) {

  struct segment xylimit, xybounds;
  integer_t isize[2], osize[2];
  int idim;
    
  get_dimensions(pixmap, isize);
  get_dimensions(output_data, osize);

  initialize_segment(&xylimit, - margin, - margin, osize[0] + margin, osize[1] + margin);
  initialize_segment(&xybounds, 0.0, (double) iy, isize[0], (double) iy);

  for (idim = 0; idim < 2; ++idim) {
    clip_bounds(pixmap, idim, &xylimit, &xybounds);
    if (xybounds.side) break;
  }
 
  xbounds[0] = floor(xybounds.point[0][0]);
  xbounds[1] = xybounds.side ? floor(xybounds.point[1][0]) : ceil(xybounds.point[1][0]);

  return;
}

