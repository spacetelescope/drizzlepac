#ifndef CDRIZZLEUTIL_H
#define CDRIZZLEUTIL_H
#include "driz_portability.h"

#include <assert.h>
#include <errno.h>
#define _USE_MATH_DEFINES       /* needed for MS Windows to define M_PI */ 
#include <math.h>
#if __STDC_VERSION__ >= 199901L
#include <stdint.h>
#endif
#include <stdlib.h>

/*****************************************************************
 ERROR HANDLING
*/
#define MAX_DRIZ_ERROR_LEN 512

struct driz_error_t {
  char last_message[MAX_DRIZ_ERROR_LEN];
};

void driz_error_init(struct driz_error_t* error);
void driz_error_set_message(struct driz_error_t* error, const char* message);
void driz_error_format_message(struct driz_error_t* error, const char* format, ...);
const char* driz_error_get_message(struct driz_error_t* error);
int driz_error_is_set(struct driz_error_t* error);
void driz_error_unset(struct driz_error_t* error);

/*****************************************************************
 LOGGING
*/

typedef void (*driz_log_func_t)(const char *, ...);
extern driz_log_func_t driz_log_func;
void driz_default_log_func(const char *, ...);

#define DRIZLOG(format, ...) ((*driz_log_func)(format, __VA_ARGS__))

/*****************************************************************
 CONVENIENCE MACROS
*/
#if !defined(MIN)
  #define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#endif

#if !defined(MAX)
  #define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#endif

#define CLAMP(x, low, high)  (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))
#define CLAMP_ABOVE(x, low)  (((x) < low) ? (low) : (x))
#define CLAMP_BELOW(x, high)  (((x) > high) ? (high) : (x))

#ifdef __GNUC__
#define UNUSED_PARAM __attribute__((unused))
#else
#define UNUSED_PARAM
#endif

#define MAX_DOUBLE 1.7976931348623158e+308
#define MIN_DOUBLE 2.2250738585072014e-308

#define MAX_COEFFS 128
#define COEFF_OFFSET 100

/* This parameter controls the maximum number of entries in
   the context table. This value is rather arbitrary */
#define MAXEN 1000

/* The following parameter controls the number of images which
   may be combined when contexts are in use. */
#define MAXIM 128

#undef TRUE
#define TRUE 1

#undef FALSE
#define FALSE 0

/*****************************************************************
 DATA TYPES
*/
typedef unsigned int length_t;
typedef int integer_t;
#if __STDC_VERSION__ >= 199901L
typedef int_fast8_t bool_t;
#else
typedef unsigned char bool_t;
#endif

enum e_shift_t {
  shift_input,
  shift_output
};

enum e_align_t {
  align_center,
  align_corner
};

enum e_kernel_t {
  kernel_square,
  kernel_gaussian,
  kernel_point,
  kernel_tophat,
  kernel_turbo,
  kernel_lanczos2,
  kernel_lanczos3,
  kernel_LAST
};

enum e_unit_t {
  unit_counts,
  unit_cps
};

enum e_interp_t {
  interp_nearest,
  interp_bilinear,
  interp_poly3,
  interp_poly5,
  interp_spline3,
  interp_sinc,
  interp_lsinc,
  interp_lanczos3,
  interp_lanczos5,
  interp_LAST
};

/* Lanczos values */
struct lanczos_param_t {
  size_t nlut;
  float* lut;
  double sdp;
  integer_t nbox;
  float space;
  float misval;
};

typedef int (*mapping_callback_t) \
  (void* state,
   const double, const double,
   const integer_t,
   double*, double*,
   /* Output parameters */
   double* /*[n]*/, double* /*[n]*/,
   struct driz_error_t*);

struct driz_param_t {
  /* Drizzle callback to perform the actual drizzling */
  mapping_callback_t mapping_callback;
  void* mapping_callback_state;

  /* Kernel shape and size */
  enum e_kernel_t kernel;
  double pixel_fraction; /* was: PIXFRAC */

  /* Exposure time */
  float exposure_time; /* was: EXPIN */

  /* Weight scale */
  float weight_scale; /* was: WTSCL */

  /* Filling */
  float fill_value; /* was: FILVAL */
  bool_t do_fill; /* was: FILL */

  /* CPS / counts */
  enum e_unit_t in_units; /* was: INCPS, either counts or CPS */
  enum e_unit_t out_units; /* was: INCPS, either counts or CPS */

  /* Input data */
  integer_t dny;
  integer_t dnx;
  integer_t ny;
  float* data; /* [dny][dnx] */
  float* weights; /* [dny][dnx] */

  /* Output data */
  integer_t onx;
  integer_t ony;
  float* output_data; /* [ony][onx] */
  float* output_counts; /* [ony][onx] was: COU */
  integer_t* output_context; /* [ony][onx] was: CONTIM */

  /* Blotting-specific parameters */
  enum e_interp_t interpolation; /* was INTERP */
  float ef; /* TODO: Rename these variables */
  float misval;
  float sinscl;
  float kscale;
  float kscale2;

  double ox;
  double oy;

  integer_t uuid; /* was: UNIQID */

  integer_t xmin;
  integer_t xmax;
  integer_t ymin;
  integer_t ymax;

  bool_t sub;
  bool_t no_over;

  integer_t nsx;
  integer_t nsy;

  integer_t intab[MAXEN*MAXIM]; /* [maxen][maxim] */
  integer_t nen; /* TODO: Rename me */

  integer_t bv;
  double ac;
  double pfo;
  double pfo2;

  integer_t* output_done; /* [nsy][nsx] */

  /* Stuff specific to certain kernel types */
  /* Gaussian values */
  struct {
    double efac;
    double es;
  } gaussian;
  struct lanczos_param_t lanczos;

  /* Scaling */
  enum e_align_t align;
  double scale;
  double scale2;
  double x_scale;
  double y_scale;
};

/**
Initialize all of the members of the drizzle_param_t to sane default
values, mostly zeroes.  Note, these are not *meaningful* values, just
ones that help with memory management etc.  It is up to users of the
struct, e.g. cdrizzle_, to fill the struct with valid parameters.
*/
void
driz_param_init(struct driz_param_t* p);

void
driz_param_dump(struct driz_param_t* p);

/****************************************************************************/
/* ARRAY ACCESSORS */
static inline_macro float*
data_ptr(struct driz_param_t* p, integer_t x, integer_t y) {
  assert(p);
  assert(p->data);
  assert(x >= 0 && x < p->dnx);
  assert(y >= 0 && y < p->dny);
  return (p->data + (y * p->dnx) + x);
}

static inline_macro const float*
weights_ptr(struct driz_param_t* p, integer_t x, integer_t y) {
  assert(p);
  assert(p->weights);
  assert(x >= 0 && x < p->dnx);
  assert(y >= 0 && y < p->dny);
  return (p->weights + (y * p->dnx) + x);
}

static inline_macro float*
output_data_ptr(struct driz_param_t* p, integer_t x, integer_t y) {
  assert(p);
  assert(p->output_data);
  assert(x >= 0 && x < p->onx);
  assert(y >= 0 && y < p->ony);
  return (p->output_data + (y * p->onx) + x);
}

static inline_macro float*
output_counts_ptr(struct driz_param_t* p, integer_t x, integer_t y) {
  assert(p);
  assert(p->output_counts);
  assert(x >= 0 && x < p->onx);
  assert(y >= 0 && y < p->ony);
  return (p->output_counts + (y * p->onx) + x);
}

static inline_macro integer_t*
output_context_ptr(struct driz_param_t* p, integer_t x, integer_t y) {
  assert(p);
  assert(p->output_context);
  assert(x >= 0 && x < p->onx);
  assert(y >= 0 && y < p->ony);
  return (p->output_context + (y * p->onx) + x);
}

static inline_macro integer_t*
output_done_ptr(struct driz_param_t* p, integer_t x, integer_t y) {
  assert(p);
  assert(p->output_done);
  assert(x >= 0 && x < p->onx);
  assert(y >= 0 && y < p->ony);
  return (p->output_done + (y * p->onx) + x);
}

static inline_macro integer_t*
intab_ptr(struct driz_param_t* p, integer_t i0, integer_t i1) {
  assert(p);
  assert(p->intab);
  assert(i0 >= 0 && i0 < MAXIM);
  assert(i1 >= 0 && i1 < MAXEN);
  return (p->intab + (i1 * MAXIM) + i0);
}

/*****************************************************************
 STRING TO ENUMERATION CONVERSIONS
*/
int
shift_str2enum(const char* s, enum e_shift_t* result, struct driz_error_t* error);

int
align_str2enum(const char* s, enum e_align_t* result, struct driz_error_t* error);

int
kernel_str2enum(const char* s, enum e_kernel_t* result, struct driz_error_t* error);

int
unit_str2enum(const char* s, enum e_unit_t* result, struct driz_error_t* error);

int
interp_str2enum(const char* s, enum e_interp_t* result, struct driz_error_t* error);

const char*
shift_enum2str(enum e_shift_t value);

const char*
align_enum2str(enum e_align_t value);

const char*
kernel_enum2str(enum e_kernel_t value);

const char*
unit_enum2str(enum e_unit_t value);

const char*
interp_enum2str(enum e_interp_t value);

const char*
bool2str(bool_t value);

/*****************************************************************
 NUMERICAL UTILITIES
*/
/**
Fill up a look-up-table of Lanczos interpolation kernel values for
rapid weighting determination for kernel == kernel_lanczos.

@param kernel_order the order of the kernel.
@param npix the size of the lookup table
@param del the spacings of the sampling of the function
@param lanczos_lut 1d array of lookup values.  This is a single-sided Lanczos
   function with lanczos_lut[0] being the central value.

Note that no checking is done to see whether the values are sensible.

was: FILALU
*/
void
create_lanczos_lut(const int kernel_order, const size_t npix,
                   const float del, float* lanczos_lut);

void
put_fill(struct driz_param_t* p, const float fill_value);

/**
 Calculate the refractive index of MgF2 for a given C wavelength (in
 nm) using the formula given by Trauger (1995)
*/
double
mgf2(double lambda);

/**
Weighted sum of 2 real vectors.

was: WSUMR
*/
static inline_macro void
weighted_sum_vectors(const integer_t npix,
                     const float* a /*[npix]*/, const float w1,
                     const float* b /*[npix]*/, const float w2,
                     /* Output arguments */
                     float* c /*[npix]*/) {
  float* c_end = c + npix;

  assert(a);
  assert(b);
  assert(c);

  while(c != c_end)
    *(c++) = *(a++) * w1 + *(b++) * w2;
}

/**
 Round to nearest integer in a way that mimics fortrans NINT
*/
static inline_macro integer_t
fortran_round(const double x) {
  return (x >= 0) ? (integer_t)floor(x + .5) : (integer_t)-floor(.5 - x);
}

static inline_macro double
min_doubles(const double* a, const integer_t size) {
  const double* end = a + size;
  double value = MAX_DOUBLE;
  for ( ; a != end; ++a)
    if (*a < value)
      value = *a;
  return value;
}

static inline_macro double
max_doubles(const double* a, const integer_t size) {
  const double* end = a + size;
  double value = MIN_DOUBLE;
  for ( ; a != end; ++a)
    if (*a > value)
      value = *a;
  return value;
}

/**
Evaluate a cubic geometric distortion in 2D

@param x The x coordinate

@param y The y coordinate

@param co An array of length 10 of coefficients

@return The distorted value
*/
static inline_macro double
eval3(const double x, const double y, const double* co) {

  register double x2 = x * x;
  register double y2 = y * y;

  assert(co);

  return
    co[0] +
    co[1] * x +
    co[2] * y +
    co[3] * x2 +
    co[4] * x * y +
    co[5] * y2 +
    co[6] * x2 * x +
    co[7] * x2 * y +
    co[8] * x * y2 +
    co[9] * y * y2;
}

/**
Evaluate a 4th order (quartic) geometric distortion in 2d

@param x The x coordinate

@param y The y coordinate

@param co An array of length 15 of coefficients

@return The distorted value
*/
static inline_macro double
eval4(const double x, const double y, const double* co) {

  register double x2 = x * x;
  register double y2 = y * y;
  register double x3 = x2 * x;
  register double y3 = y2 * y;

  assert(co);

  return
    co[0] +
    co[1] * x +
    co[2] * y +
    co[3] * x2 +
    co[4] * x * y +
    co[5] * y2 +
    co[6] * x3 +
    co[7] * x2 * y +
    co[8] * x * y2 +
    co[9] * y3 +
    co[10] * x2 * x2 +
    co[11] * x3 * y +
    co[12] * x2 * y2 +
    co[13] * x * y3 +
    co[14] * y2 * y2;
}

/**
Evaluate a 5th order geometric distortion in 2d

@param x The x coordinate

@param y The y coordinate

@param co An array of length 21 of coefficients

@return The distorted value
*/
static inline_macro double
eval5(const double x, const double y, const double* co) {

  register double x2 = x * x;
  register double y2 = y * y;
  register double x3 = x2 * x;
  register double y3 = y2 * y;
  register double xy = x * y;

  assert(co);

  return
    co[0] +
    co[1] * x +
    co[2] * y +
    co[3] * x2 +
    co[4] * xy +
    co[5] * y2 +
    co[6] * x3 +
    co[7] * x2 * y +
    co[8] * x * y2 +
    co[9] * y3 +
    co[10] * x2 * x2 +
    co[11] * x3 * y +
    co[12] * x2 * y2 +
    co[13] * x * y3 +
    co[14] * y2 * y2 +
    co[15] * x2 * x3 +
    co[16] * x3 * xy +
    co[17] * x3 * y2 +
    co[18] * x2 * y3 +
    co[19] * xy * y3 +
    co[20] * y2 * y3;
}

/**
Evaluate the value of a general 2d polynomial in X and Y
complete with "half" cross-terms.

For orders lower than 7 it is slightly more efficient to
use the specific (eval4 etc) routines instead.

@param x The x coordinate

@param y The y coordinate

@param co An array of length 10 of coefficients

@param order The order of the polynomial

@param result The distorted value

@return 1 if error occurred
*/
static inline_macro int
evaln(const double x, const double y, const double* co, const integer_t order,
      /* Output parameters */
      double* result, struct driz_error_t* error) {
  integer_t n, m;
  double t;
  const double* c;
  double xco, yco;

  assert(co);

  t = 0.0;
  c = co;
  errno = 0;
  for (n = 1; n <= order + 1; ++n) {
    for (m = 1; m <= n; ++m) {
      xco = pow(x, (double)(n - m));
      if (errno != 0) {
        driz_error_set_message(error, "pow failed");
        return 1;
      }
      yco = pow(y, (double)(m - 1));
      if (errno != 0) {
        driz_error_set_message(error, "pow failed");
        return 1;
      }
      t += *(c++) * xco * yco;
    }
  }

  *result = t;

  return 0;
}

/**
Evaluate a 3rd order radial geometric distortion in 2d
X version. Note that there is no zero order coefficient
as this is physically meaningless.

@param x The x coordinate

@param y The y coordinate

@param co An array of length 4 of coefficients

@param[out] xo The distorted x coordinate

@param[out] yo The distorted y coordinate
*/
static inline_macro void
rad3(const double x, const double y, const double* co,
     /* Output parameters */
     double* xo, double* yo) {
  double r, f;

  assert(co);
  assert(xo);
  assert(yo);

  r = sqrt(x*x + y*y);

  f = 1.0 + co[0] + co[1]*r + co[2]*r*r;
  *xo = f*x;
  *yo = f*y;
}

#endif /* CDRIZZLEUTIL_H */
