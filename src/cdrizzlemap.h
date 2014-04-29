#ifndef CDRIZZLEDRIZ_H
#define CDRIZZLEDRIZ_H

#include "driz_portability.h"
#include "astropy_wcs_api.h"
#include "cdrizzleutil.h"
/**

Declarations for supporting the DefaultWCSMapping (WCS-based)
transformations.

*/
struct wcsmap_param_t {
  /* Pointers to PyWCS objects for input and output WCS */
  pipeline_t* input_wcs;
  pipeline_t* output_wcs;
  double*     table;
  int         nx, ny;
  int         snx, sny;
  double      factor;
};

/**
Initialize all of the members of the mapping_param_t to sane default
values, mostly zeroes.  Note, these are not *meaningful* values, just
ones that help with memory management etc.  It is up to users of the
struct, e.g. cdrizzle_, to fill the struct with valid parameters.
*/
void
wcsmap_param_init(struct wcsmap_param_t* m);

void
wcsmap_param_dump(struct wcsmap_param_t* m);

void
wcsmap_param_free(struct wcsmap_param_t* m);

int
default_wcsmap(void* state,
                const double xd, const double yd,
                const integer_t n,
                double* xin /*[n]*/, double* yin /*[n]*/,
                /* Output parameters */
                double* xout, double* yout,
                struct driz_error_t* error);
int
default_wcsmap_init(struct wcsmap_param_t* m,
                    pipeline_t* input,
                    pipeline_t* output,
                    int nx, int ny, double factor,
                    /* Output parameters */
                    struct driz_error_t* error);

/**

Declarations for supporting the DefaultMapping (pixel-based)
transformations.

*/
struct mapping_param_t {
  /* Geometric distortion coefficients */
  double x_coeffs[MAX_COEFFS]; /* was: XCO */
  double y_coeffs[MAX_COEFFS]; /* was: YCO */
  integer_t num_coeffs; /* was: CONUM */
  integer_t coeff_type; /* TODO: Define this variable */
  double lambda; /* was: LAM */

  /* Distortion image arrays */
  const float* x_distortion; /* [y_dist_dim][x_dist_dim] was: PXG */
  integer_t x_dist_dim; /* was: XGDIM */
  const float* y_distortion; /* [y_dist_dim][x_dist_dim] was: PYG */
  integer_t y_dist_dim; /* was: YGDIM */
  bool_t use_distortion_image; /* was: DISIM */

  /* Primary linear transformation parameters */
  double scale;
  double rotation; /* was: ROT */
  double x_shift; /* was: XSH */
  double y_shift; /* was: YSH */
  enum e_shift_t do_shift_first; /* was: SHFTFR */
  enum e_shift_t shift_units; /* was: SHFTUN */
  enum e_align_t align;

  /* Secondary transformation parameters */
  double x_scale;
  double y_scale;
  double x_shift2; /* was: XSH2 */
  double y_shift2; /* was: YSH2 */
  double rotation2; /* was: ROT2 */
  enum e_shift_t do_shift2_first; /* was SHFR2 */

  bool_t rotate_first2; /* was: ROTF2 */
  bool_t has_secondary_parameters; /* was: SECPAR */

  /* Never seem to get set, but are here for completeness */
  double alpha;
  double beta;

  /* WCS mode parameters of the output (only when use_wcs == TRUE) */
  bool_t use_wcs;
  double output_scale; /* was: OUTSCL */
  double ra_ref; /* was RACEN */
  double dec_ref; /* was DECCEN */
  double x_ref_pixel; /* was XREFP */
  double y_ref_pixel; /* was YREFP */
  double orientation; /* was: ORIENT */

  /* WCS coefficients, based on WCS mode parameters in public
     section. */
  double wcsout[8];
  double wcs[8];
  bool_t got_wcs;
  char ctype1[16];
  char ctype2[16];

  double xcen;
  double ycen;
  double xp;
  double yp;
  integer_t dnx;
  integer_t dny;
};


static inline_macro const float*
x_distortion_ptr(struct mapping_param_t* m, integer_t i0, integer_t i1) {
  assert(m);
  assert(m->x_distortion);
  assert(i0 >= 0 && i0 < m->x_dist_dim);
  assert(i1 >= 0 && i1 < m->y_dist_dim);
  return (m->x_distortion + (i1 * m->x_dist_dim) + i0);
}

static inline_macro const float*
y_distortion_ptr(struct mapping_param_t* m, integer_t i0, integer_t i1) {
  assert(m);
  assert(m->y_distortion);
  assert(i0 >= 0 && i0 < m->x_dist_dim);
  assert(i1 >= 0 && i1 < m->y_dist_dim);
  return (m->y_distortion + (i1 * m->x_dist_dim) + i0);
}

/**

This function will be used by both the pixel-based and WCS-based
mapping functions; namely, DefaultMapping and DefaultWCSMapping.

Apply the standard Drizzle transformation from input to output pixel
coordinates.  This may optionally include a polynomial distortion
solution or be specified using input and output WCS.
*/
int
map_value(struct driz_param_t* p,
          const bool_t regular,
          const integer_t n,
          const double* xin /*[n]*/, const double* yin /*[n]*/,
          /* Output parameters */
          double* xtmp /*[n]*/, double* ytmp /*[n]*/,
          double* xout /*[n]*/, double* yout /*[n]*/,
          struct driz_error_t* error);


#endif /* CDRIZZLEDRIZ_H */
