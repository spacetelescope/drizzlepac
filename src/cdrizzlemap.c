#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <Python.h>
#include <numpy/arrayobject.h>

#include "driz_portability.h"
#include "cdrizzlemap.h"
#include "cdrizzleutil.h"

/*----------------------------------------------------------------------------------------------------
 * Initialize a segment structure to contain the two points (x1, y1) and (x2, y2)
 * the valid flag is initially set to 0 (valid)
 */

void
initialize_segment(struct segment *self, integer_t x1, integer_t y1, integer_t x2, integer_t y2) {
  self->point[0][0] = x1;
  self->point[0][1] = y1;
  self->point[1][0] = x2;
  self->point[1][1] = y2;
  self->valid = 0;

  return;
}

/*----------------------------------------------------------------------------------------------------
 * Generate a string representation of a segment for debugging
 *
 * self: the segment
 * str:  the string representation, at least 64 chars (output)
 */

void
show_segment(struct segment *self, char *str) {
  sprintf(str, "(%10f,%10f) - (%10f,%10f) [%2d]",
          self->point[0][0], self->point[0][1],
          self->point[1][0], self->point[1][1],
          self->valid);

  return;
}

/*----------------------------------------------------------------------------------------------------
 * Take the the union of several line segments along a dimension.
 * That is, the result is the combined range of all the segments along a dimension
 *
 * npoint:   number of line segments to combine
 * jdim:     dimension to take union along, x (0) or y (1)
 * sybounds: the array of line segments
 * bounds:   union of the segment range (output)
 */

void
union_of_segments(int npoint, int jdim, struct segment xybounds[], integer_t bounds[2]) {
  int ipoint;
  int none = 1;
  
  for (ipoint = 0; ipoint < npoint; ++ipoint) {
    if (xybounds[ipoint].valid == 0) {
      integer_t lo = floor(xybounds[ipoint].point[0][jdim]);
      integer_t hi = ceil(xybounds[ipoint].point[1][jdim]);

      if (none == 0) {
        if (lo < bounds[0]) bounds[0] = lo;
        if (hi > bounds[1]) bounds[1] = hi;

      } else {
        none = 0;
        bounds[0] = lo;
        bounds[1] = hi;
      }
    }
  }
  
  if (none) {
    bounds[1] = bounds[0];
  }

  return;
}

/*----------------------------------------------------------------------------------------------------
 * Map a point on the input image to the output image using
 * a mapping of the pixel centers between the two by interpolating
 * between the centers in the mapping
 *
 * pixmap: The mapping of the pixel centers from input to output image
 * xyin:   An (x,y) point on the input image
 * xyout:  The same (x, y) point on the output image (output)
 */

void
map_point(PyArrayObject *pixmap,
          const double xyin[2],
          double xyout[2]
         ) {

  int        idim;
  integer_t  pix[2];
  double     frac[2];
  integer_t  xydim[2];

  /* Divide the position into an integral pixel position
   * plus a fractional offset within the pixel */

  get_dimensions(pixmap, xydim);

  for (idim = 0; idim < 2; ++idim) {
    frac[idim] = xyin[idim];
    pix[idim]  = frac[idim];
    pix[idim]  = CLAMP(pix[idim], 0, xydim[idim] - 2);
    frac[idim] = frac[idim] - pix[idim];
  }
  
  if (frac[0] == 0.0 && frac[1] == 0.0) {
    /* No interpolation needed if input position has no fractional part */
    for (idim = 0; idim < 2; ++idim) {
      xyout[idim] = get_pixmap(pixmap, pix[0], pix[1])[idim];
    }

  } else {
    /* Bilinear interpolation btw pixels, see Wikipedia for explanation */
    for (idim = 0; idim < 2; ++idim) {
      xyout[idim] = (1.0 - frac[0]) * (1.0 - frac[1]) * get_pixmap(pixmap, pix[0], pix[1])[idim] +
                    frac[0] * (1.0 - frac[1]) * get_pixmap(pixmap, pix[0]+1, pix[1])[idim] +
                    (1.0 - frac[0]) * frac[1] * get_pixmap(pixmap, pix[0], pix[1]+1)[idim] +
                    frac[0] * frac[1] * get_pixmap(pixmap, pix[0]+1, pix[1]+1)[idim];
    }
  }

  return;
}


/*----------------------------------------------------------------------------------------------------
 * Clip a line segment from an input image to the limits of an output image along one dimension
 *
 * pixmap:   the mapping between input and output images
 * jdim:     the dimension to do the clipping to, x (0) or y (1)
 * xylimit:  the limits  of the output image
 * xybounds: the clipped line segment (output)
 * 
 */

void
clip_bounds(PyArrayObject *pixmap, int jdim, struct segment *xylimit, struct segment *xybounds) {
  int ipoint, idim;
  
  xybounds->valid = 1; /* Track if bounds are both outside the image */
  
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
      xybounds->valid = 0;
      xybounds->point[ipoint][jdim] = b;

    } else {
      /* No bracket, so track which side the bound lies on */
      xybounds->valid *= fa > 0.0 ? +1 : -1;
    }
  }

  if (xybounds->valid > 0) {
      /* Positive means both bounds are outside the image */
      xybounds->point[1][jdim] = xybounds->point[0][jdim]; 
  } else {
      /* Negative means both bounds are inside, which is not a problem */
      xybounds->valid = 0;
  }
  
  return;
}

/*----------------------------------------------------------------------------------------------------
 * Determine the range of pixels in a specified line of an input image
 * which are inside the output image. Range is one-sided, that is, the second
 * value returned is one greater than the last pixel that is on the image.
 *
 * p:       the stucture containing the image pointers
 * margin:  a margin in pixels added to the limits
 * j:       the index of the line in the input image whose range is computed
 * xbounds: the input pixels bounding the overlap (output)
 */

void
check_line_overlap(struct driz_param_t* p, int margin, integer_t j, integer_t *xbounds) {

  struct segment xylimit, xybounds;
  integer_t isize[2];
  int idim;
    
  get_dimensions(p->pixmap, isize);
  initialize_segment(&xybounds, 0, j, isize[0], j);
  initialize_segment(&xylimit, p->xmin - margin, p->ymin - margin,
                               p->xmax + margin, p->ymax + margin);

  for (idim = 0; idim < 2; ++idim) {
    clip_bounds(p->pixmap, idim, &xylimit, &xybounds);
    if (xybounds.valid) break;
  }
 
  xbounds[0] = floor(xybounds.point[0][0]);
  xbounds[1] = ceil(xybounds.point[1][0]);
  if (xybounds.valid) xbounds[1] = xbounds[0];

  return;
}

/*----------------------------------------------------------------------------------------------------
 * Determine the range of lines in the input image that overlap the output image
 * Range is one-sided, that is, the second value returned is one greater than the
 * last line that is on the image.
 * 
 * p:       the stucture containing the image pointers
 * margin:  a margin in pixels added to the limits
 * ybounds: the input lines bounding the overlap (output)
 */

void
check_image_overlap(struct driz_param_t* p, const int margin, integer_t *ybounds) {

  struct segment xylimit, xybounds[2];
  integer_t isize[2], osize[2];
  int ipoint, idim;
  
  get_dimensions(p->pixmap, isize);

  ybounds[0] = 0;
  ybounds[1] = isize[0];
  
  initialize_segment(&xylimit, p->xmin - margin, p->ymin - margin,
                               p->xmax + margin, p->ymax + margin);

  for (ipoint = 0; ipoint < 2; ++ipoint) {
    initialize_segment(&xybounds[ipoint], ybounds[ipoint], 0, ybounds[ipoint], isize[1]);
    
    for (idim = 0; idim < 2; ++idim) {
      clip_bounds(p->pixmap, idim, &xylimit, &xybounds[ipoint]);
      if (xybounds[ipoint].valid) break;
    }
  }

  union_of_segments(2, 1, xybounds, ybounds);
  return;
}

