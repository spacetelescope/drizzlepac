
#include <Python.h>

#define _USE_MATH_DEFINES       /* needed for MS Windows to define M_PI */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "astropy_wcs_api.h"
#include "astropy_wcs.h"

#include "cdrizzleblot.h"
#include "cdrizzlebox.h"
#include "cdrizzlemap.h"
#include "cdrizzleutil.h"
#include "cdrizzlewcs.h"

static PyObject *gl_Error;


/*
 A mapping callback that delegates to a Python-based drizzle
 callback.

 The Python mapping callback must have the following signature:

    def mapping(xin, yin):
        return xout, yout

 xin, yin are the input coordinates, and xout and yout are the output
 coordinates.  All are 1-dimensional Numpy DOUBLE arrays of the same
 length.
*/
static int
py_mapping_callback(void* state,
                    const double xd, const double yd,
                    const integer_t n,
                    double* xin /*[n]*/, double* yin /*[n]*/,
                    /* Output parameters */
                    double* xout, double* yout,
                    struct driz_error_t* error) {
  PyObject* callback = (PyObject*)state;
  npy_intp dims = n;
  PyArrayObject* py_xin = NULL;
  PyArrayObject* py_yin = NULL;
  PyObject* py_xout_obj = NULL;
  PyObject* py_yout_obj = NULL;
  PyArrayObject* py_xout = NULL;
  PyArrayObject* py_yout = NULL;
  PyObject* callback_result = NULL;
  PyObject* callback_tuple = NULL;
  int result = TRUE;

  py_xin = (PyArrayObject*)PyArray_SimpleNewFromData(1, &dims, NPY_FLOAT64, (double*)xin);
  if (py_xin == NULL)
    goto _py_mapping_callback_exit;

  py_yin = (PyArrayObject*)PyArray_SimpleNewFromData(1, &dims, NPY_FLOAT64, (double*)yin);
  if (py_yin == NULL)
    goto _py_mapping_callback_exit;

  callback_result = PyObject_CallFunctionObjArgs(callback, py_xin, py_yin, NULL);

  if (callback_result == NULL)
    goto _py_mapping_callback_exit;

  callback_tuple = PySequence_Tuple(callback_result);
  if (callback_tuple == NULL)
    goto _py_mapping_callback_exit;

  if (!PyArg_UnpackTuple(callback_tuple, "result", 2, 2, &py_xout_obj, &py_yout_obj))
    goto _py_mapping_callback_exit;

  py_xout = (PyArrayObject*)PyArray_ContiguousFromAny(py_xout_obj, NPY_FLOAT64, 1, 1);
  if (py_xout == NULL)
    goto _py_mapping_callback_exit;

  py_yout = (PyArrayObject*)PyArray_ContiguousFromAny(py_yout_obj, NPY_FLOAT64, 1, 1);
  if (py_yout == NULL)
    goto _py_mapping_callback_exit;

  if (PyArray_DIM(py_xout, 0) != n) {
    PyErr_Format(PyExc_ValueError,
                 "Returned arrays must be same dimension as passed-in arrays.  Expected '%d', got '%d'",
                 (int)n, (int)PyArray_DIM(py_xout, 0));
    goto _py_mapping_callback_exit;
  }

  if (PyArray_DIM(py_yout, 0) != n) {
    PyErr_Format(PyExc_ValueError,
                 "Returned arrays must be same dimension as passed-in arrays.  Expected '%d', got '%d'",
                 (int)n, (int)PyArray_DIM(py_yout, 0));
    goto _py_mapping_callback_exit;
  }

  memcpy(xout, PyArray_DATA(py_xout), n*sizeof(double));
  memcpy(yout, PyArray_DATA(py_yout), n*sizeof(double));

  result = FALSE;

 _py_mapping_callback_exit:
  Py_XDECREF(py_xin);
  Py_XDECREF(py_yin);
  Py_XDECREF(callback_result);
  Py_XDECREF(callback_tuple);
  Py_XDECREF(py_xout);
  Py_XDECREF(py_yout);

  if (result)
    driz_error_set_message(error, "<PYTHON>");

  return result;
}

/**

Code to implement the WCS-based C interface for the mapping.

It uses the same py_mapping_callback as the DefaultMapping
(pixel-based default transformation) uses, since we want to
retain the support for providing Python-based transformations
in addition to the C-based transformations.

*/
typedef struct {
  PyObject_HEAD
  struct wcsmap_param_t m;
  PyObject* py_input;
  PyObject* py_output;
} PyWCSMap;

static void
PyWCSMap_dealloc(PyWCSMap* self)
{
  /* Deal with our reference-counted members */
  Py_XDECREF(self->py_input);  self->py_input = NULL;
  Py_XDECREF(self->py_output); self->py_output = NULL;
  wcsmap_param_free(&self->m);

  Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject *
PyWCSMap_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  PyWCSMap *self;

  self = (PyWCSMap *)type->tp_alloc(type, 0);
  self->py_input = NULL;
  self->py_output = NULL;
  if (self != NULL) {
    wcsmap_param_init(&self->m);
  }

  return (PyObject *)self;
}

static int
PyWCSMap_init(PyWCSMap *self, PyObject *args, PyObject *kwds)
{
  /* Arguments in the order they appear */
  PyObject *input_obj = NULL;
  PyObject *output_obj = NULL;
  int nx, ny;
  double factor;
  int status = -1;

  /* Other miscellaneous local variables */
  struct driz_error_t error;
  int istat = 1;

  driz_error_init(&error);

  /* TODO: Make factor a kwarg */
  if (! PyArg_ParseTuple(args, "OOiid:DefaultWCSMapping.__init__",
                         &input_obj, &output_obj, &nx, &ny, &factor)){
    goto exit;
  }

  /* Create the C struct from all of these mapping parameters */
  istat = default_wcsmap_init(
      &self->m,
      &((Wcs*)input_obj)->x, &((Wcs*)output_obj)->x,
      nx, ny, factor,
      &error);

  if (istat || driz_error_is_set(&error)) {
    if (strcmp(driz_error_get_message(&error), "<PYTHON>") != 0)
      PyErr_SetString(PyExc_Exception, driz_error_get_message(&error));
    goto exit;
  }

  Py_INCREF(input_obj);
  Py_INCREF(output_obj);
  self->py_input = input_obj;
  self->py_output = output_obj;

  status = 0;

 exit:

  return status;
}

static PyObject*
PyWCSMap_call(PyWCSMap* self, PyObject* args, PyObject* kwargs)
{
  PyObject*           py_xin_obj = NULL;
  PyObject*           py_yin_obj = NULL;
  PyArrayObject*      py_xin     = NULL;
  PyArrayObject*      py_yin     = NULL;
  PyArrayObject*      py_xout    = NULL;
  PyArrayObject*      py_yout    = NULL;
  npy_intp            dims;
  PyObject*           result     = NULL;

  struct driz_error_t error;

  driz_error_init(&error);

  if (!PyArg_ParseTuple(args, "OO", &py_xin_obj, &py_yin_obj)) {
    return NULL;
  }

  py_xin = (PyArrayObject*)PyArray_ContiguousFromAny(py_xin_obj, NPY_FLOAT64, 1, 1);
  if (py_xin == NULL) {
    goto _py_wcsmap_call_exit;
  }

  py_yin = (PyArrayObject*)PyArray_ContiguousFromAny(py_yin_obj, NPY_FLOAT64, 1, 1);
  if (py_yin == NULL) {
    goto _py_wcsmap_call_exit;
  }

  if (PyArray_DIM(py_xin, 0) != PyArray_DIM(py_yin, 0)) {
    PyErr_Format(PyExc_ValueError,
                 "Passed in arrays must have the same dimensions.  Got '%d' and '%d'",
                 (int)PyArray_DIM(py_xin, 0), (int)PyArray_DIM(py_yin, 0));
    goto _py_wcsmap_call_exit;
  }

  dims = PyArray_DIM(py_xin, 0);

  /* Generate output arrays */
  py_xout = (PyArrayObject*)PyArray_SimpleNew(1, &dims, NPY_FLOAT64);
  if (py_xout == NULL) {
    goto _py_wcsmap_call_exit;
  }

  py_yout = (PyArrayObject*)PyArray_SimpleNew(1, &dims, NPY_FLOAT64);
  if (py_yout == NULL) {
    goto _py_wcsmap_call_exit;
  }

  if (default_wcsmap(&self->m, 0, 0, (integer_t)dims,
                      PyArray_DATA(py_xin), PyArray_DATA(py_yin),
                      PyArray_DATA(py_xout), PyArray_DATA(py_yout),
                      &error)) {
    goto _py_wcsmap_call_exit;
  }

  result = Py_BuildValue("OO", py_xout, py_yout);
  if (result == NULL) {
    goto _py_wcsmap_call_exit;
  }

 _py_wcsmap_call_exit:
  Py_XDECREF(py_xin);
  Py_XDECREF(py_yin);
  Py_XDECREF(py_xout);
  Py_XDECREF(py_yout);

  if (driz_error_is_set(&error)) {
    PyErr_SetString(PyExc_Exception, driz_error_get_message(&error));
  }

  return result;
}

static PyTypeObject WCSMapType = {
  PyVarObject_HEAD_INIT(NULL, 0)
  (char *) "cdriz.DefaultWCSMapping",              /*tp_name*/
  sizeof(PyWCSMap),                                /*tp_basicsize*/
  0,                                               /*tp_itemsize*/
  (destructor) PyWCSMap_dealloc,                   /*tp_dealloc*/
  0,                                               /*tp_print*/
  0,                                               /*tp_getattr*/
  0,                                               /*tp_setattr*/
  0,                                               /*tp_compare*/
  0,                                               /*tp_repr*/
  0,                                               /*tp_as_number*/
  0,                                               /*tp_as_sequence*/
  0,                                               /*tp_as_mapping*/
  0,                                               /*tp_hash */
  (ternaryfunc) PyWCSMap_call,                     /*tp_call*/
  0,                                               /*tp_str*/
  0,                                               /*tp_getattro*/
  0,                                               /*tp_setattro*/
  0,                                               /*tp_as_buffer*/
  (long) Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
  (char *) "DefaultWCSMapping(input,output)",      /* tp_doc */
  0,                                               /* tp_traverse */
  0,                                               /* tp_clear */
  0,                                               /* tp_richcompare */
  0,                                               /* tp_weaklistoffset */
  0,                                               /* tp_iter */
  0,                                               /* tp_iternext */
  0,                                               /* tp_methods */
  0,                                               /* tp_members */
  0,                                               /* tp_getset */
  0,                                               /* tp_base */
  0,                                               /* tp_dict */
  0,                                               /* tp_descr_get */
  0,                                               /* tp_descr_set */
  0,                                               /* tp_dictoffset */
  (initproc)PyWCSMap_init,                         /* tp_init */
  0,                                               /* tp_alloc */
  PyWCSMap_new,                                    /* tp_new */
};

static PyObject *
tdriz(PyObject *obj UNUSED_PARAM, PyObject *args)
{
  /* Arguments in the order they appear */
  PyObject *oimg, *owei, *oout, *owht, *ocon;
  long uniqid, ystart, xmin, ymin, dny;
  double scale, xscale, yscale;
  char *align_str;
  double pfract;
  char *kernel_str, *inun_str;
  float expin, wtscl;
  char *fillstr;
  integer_t nmiss, nskip, vflag;
  PyObject *callback_obj;

  /* Derived values */
  PyArrayObject *img = NULL, *wei = NULL, *out = NULL, *wht = NULL, *con = NULL;
  enum e_align_t align;
  enum e_kernel_t kernel;
  enum e_unit_t inun;
  integer_t nx, ny, onx, ony;
  char *fillstr_end;
  bool_t do_fill;
  float fill_value;
  mapping_callback_t callback = NULL;
  void* callback_state = NULL;
  int istat = 0;
  struct driz_error_t error;
  struct driz_param_t p;

  /* struct wcsmap_param_t* m = NULL; */
  /* clock_t start_t, end_t; */
  /* double delta_time; */

  driz_error_init(&error);

  if (!PyArg_ParseTuple(args,"OOOOOllllldddsdssffsiiiO:tdriz",
                        &oimg, &owei, &oout, &owht, &ocon, &uniqid, &ystart,
                        &xmin, &ymin, &dny, &scale, &xscale, &yscale,
                        &align_str, &pfract, &kernel_str, &inun_str,
                        &expin, &wtscl, &fillstr, &nmiss,&nskip, &vflag,
                        &callback_obj)) {
    return PyErr_Format(gl_Error, "cdriz.tdriz: Invalid Parameters.");
  }

  /* Check for invalid scale */
  if (scale == 0.0) {
    driz_error_format_message(&error, "Invalid scale %f (must be non-zero)", scale);
    goto _exit;
  }

  if (pfract < 0.0) {
    driz_error_format_message(&error, "Invalid pfract %f (must be greater than or equal to 0.0)", scale);
    goto _exit;
  }

  if (expin <= 0.0) {
    driz_error_format_message(&error, "Invalid expin %f (must be greater than 0.0)", scale);
    goto _exit;
  }

  if (PyObject_TypeCheck(callback_obj, &WCSMapType)) {
    /* If we're using the default mapping, we can set things up to avoid
       the Python/C bridge */
    callback = default_wcsmap;
    callback_state = (void *)&(((PyWCSMap *)callback_obj)->m);
    /*scale = ((PyWCSMap *)callback_obj)->m.scale; */
  } else {
    callback = py_mapping_callback;
    callback_state = (void *)callback_obj;
  }

  /* Get raw C-array data */
  img = (PyArrayObject *)PyArray_ContiguousFromAny(oimg, NPY_FLOAT32, 2, 2);
  if (!img) {
    driz_error_set_message(&error, "Invalid input array");
    goto _exit;
  }

  wei = (PyArrayObject *)PyArray_ContiguousFromAny(owei, NPY_FLOAT32, 2, 2);
  if (!wei) {
    driz_error_set_message(&error, "Invalid weights array");
    goto _exit;
  }

  out = (PyArrayObject *)PyArray_ContiguousFromAny(oout, NPY_FLOAT32, 2, 2);
  if (!out) {
    driz_error_set_message(&error, "Invalid output array");
    goto _exit;
  }

  wht = (PyArrayObject *)PyArray_ContiguousFromAny(owht, NPY_FLOAT32, 2, 2);
  if (!wht) {
    driz_error_set_message(&error, "Invalid array");
    goto _exit;
  }

  con = (PyArrayObject *)PyArray_ContiguousFromAny(ocon, NPY_INT32, 2, 2);
  if (!con) {
    driz_error_set_message(&error, "Invalid context array");
    goto _exit;
  }

  /* Convert strings to enumerations */
  if (align_str2enum(align_str, &align, &error) ||
      kernel_str2enum(kernel_str, &kernel, &error) ||
      unit_str2enum(inun_str, &inun, &error)) {
    goto _exit;
  }
  if (pfract <= 0.001){
    printf("kernel reset to POINT due to pfract being set to 0.0...\n");
    kernel_str2enum("point", &kernel, &error);
  }

  /* Convert the fill value string */
  if (fillstr == NULL ||
      *fillstr == 0 ||
      strncmp(fillstr, "INDEF", 6) == 0 ||
      strncmp(fillstr, "indef", 6) == 0)
  {
    do_fill = 0;
    fill_value = 0.0;
  } else {
    do_fill = 1;
#ifdef _WIN32
    fill_value = atof(fillstr);
#else
    fill_value = strtof(fillstr, &fillstr_end);
    if (fillstr == fillstr_end || *fillstr_end != '\0') {
      driz_error_format_message(&error, "Could not convert fill value '%s'",
                                fillstr);
      goto _exit;
    }
#endif
  }

  nx = PyArray_DIMS(img)[1];
  ny = PyArray_DIMS(img)[0];
  onx = PyArray_DIMS(out)[1];
  ony = PyArray_DIMS(out)[0];

  nmiss = 0;
  nskip = 0;

  driz_param_init(&p);

  p.data = PyArray_DATA(img);
  p.weights = PyArray_DATA(wei);
  p.output_data = PyArray_DATA(out);
  p.output_counts = PyArray_DATA(wht);
  p.output_context = PyArray_DATA(con);
  p.uuid = uniqid;
  p.xmin = xmin;
  p.ymin = ymin;
  p.dnx = nx;
  p.dny = ny;
  p.ny = dny;
  p.onx = p.xmax = onx;
  p.ony = p.ymax = ony;
  p.scale = scale;
  p.x_scale = xscale;
  p.y_scale = yscale;
  p.align = align;
  p.pixel_fraction = pfract;
  p.kernel = kernel;
  p.in_units = inun;
  p.exposure_time = expin;
  p.weight_scale = wtscl;
  p.mapping_callback = callback;
  p.mapping_callback_state = callback_state;

  /* Setup reasonable defaults for drizzling */
  p.no_over = FALSE;

  /*
  start_t = clock();
  */
  /* Do the drizzling */
  if (dobox(&p, ystart, &nmiss, &nskip, &error)) {
    goto _exit;
  }
  /*
  end_t = clock();
  delta_time = difftime(end_t, start_t)/1e+6;
  printf("==> Finished dobox() in %0.3f seconds\n",delta_time);

  start_t = clock();
  */
  /* Put in the fill values (if defined) */
  if (do_fill) {
    put_fill(&p, fill_value);
  }
  /*
  if (callback == default_wcsmap){
    m = (struct wcsmap_param_t *)p.mapping_callback_state;
    printf("==> Coordinate transformation times: total=%0.3f sec, d2im= %0.3f, sip= %0.3f.\n",m->dtime_coord,m->dtime_d2im,m->dtime_dgeosip);
  }
  */

 _exit:
  Py_XDECREF(con);
  Py_XDECREF(img);
  Py_XDECREF(wei);
  Py_XDECREF(out);
  Py_XDECREF(wht);

  if (istat || driz_error_is_set(&error)) {
    if (strcmp(driz_error_get_message(&error), "<PYTHON>") != 0)
      PyErr_SetString(PyExc_Exception, driz_error_get_message(&error));
    return NULL;
  } else {
    return Py_BuildValue("sii", "Callable C-based DRIZZLE Version 0.8 (20th May 2009)", nmiss, nskip);
  }
}

/*
static PyObject *
twdriz(PyObject *obj, PyObject *args)
{

    PyObject *oimg, *owei, *oout, *owht, *owcsin, *owcsout;
    PyArrayObject *img, *wei, *out, *wht, *wcsin, *wcsout;
    PyObject *opxg, *opyg;
    PyArrayObject *pxg, *pyg;
    double pfract;
    long xmin, ymin, ystart;
    long nmiss, nskip;
    char *kernel;
    char *coeffs;
    char vers[80];
    long vflag;
    float istat;
    long nx,ny,onx,ony,dny;
    long xgdim, ygdim;
    char *fillstr;
    long kernel_len, coeffs_len;
    long vers_len, fillstr_len;

        extern doublereal twdriz_(real *data, real *wei, real *ndat, real *ncou,
    integer *ystart, integer *nx, integer *ny, integer *dny, integer *onx, integer *ony,
        doublereal *wcs, doublereal *wcsout, real *pxg, real *pyg, integer *xgdim,
    integer *ygdim, doublereal *pfract, char *kernel,
        char *coeffs, char *filstr, integer *vflag, integer *clen, integer *nmiss,
        integer *nskip, char *vers, ftnlen kernel_len,
    ftnlen coeffs_len, ftnlen filstr_len, ftnlen vers_len);

    if (!PyArg_ParseTuple(args,"OOOOllllOOOOdssslll",
            &oimg,&owei,&oout,&owht, &ystart,&xmin,&ymin,&dny, &owcsin, &owcsout,
            &opxg, &opyg, &pfract, &kernel,&coeffs, &fillstr,
            &nmiss, &nskip, &vflag)){
         return PyErr_Format(gl_Error, "cdriz.twdriz: Invalid Parameters.");
    }

    img = (PyArrayObject *)NA_InputArray(oimg, tFloat32, C_ARRAY);
    wei = (PyArrayObject *)NA_InputArray(owei, tFloat32, C_ARRAY);
    out = (PyArrayObject *)NA_IoArray(oout, tFloat32, 0);
    wht = (PyArrayObject *)NA_IoArray(owht, tFloat32, 0);

    wcsin   = (PyArrayObject *)NA_IoArray(owcsin, tFloat64, 0);
    wcsout  = (PyArrayObject *)NA_IoArray(owcsin, tFloat64, 0);
    pxg = (PyArrayObject *)NA_InputArray(opxg, tFloat32, C_ARRAY);
    pyg = (PyArrayObject *)NA_InputArray(opyg, tFloat32, C_ARRAY);

    nx = img->dimensions[1];
    ny = img->dimensions[0];
    onx = out->dimensions[1];
    ony = out->dimensions[0];
    xgdim = pxg->dimensions[1];
    ygdim = pxg->dimensions[0];

    kernel_len = 8;
    coeffs_len = strlen(coeffs)+1;
    vers_len = 44;
    fillstr_len = strlen(fillstr) + 1;

    istat = twdriz_(NA_OFFSETDATA(img), NA_OFFSETDATA(wei),
                    NA_OFFSETDATA(out),NA_OFFSETDATA(wht),
                    &ystart,&nx,&ny, &dny, &onx,&ony,
                    NA_OFFSETDATA(wcsin), NA_OFFSETDATA(wcsout),
                    NA_OFFSETDATA(pxg),NA_OFFSETDATA(pyg),&xgdim, &ygdim,
                    &pfract, kernel, coeffs, fillstr,
                    &vflag, &coeffs_len,
                    &nmiss, &nskip, vers,
                    kernel_len, coeffs_len,
                    fillstr_len, vers_len);

    Py_DECREF(img);
    Py_DECREF(wei);
    Py_DECREF(out);
    Py_DECREF(wht);

    Py_DECREF(wcsout);
    Py_DECREF(wcsin);
    Py_DECREF(pxg);
    Py_DECREF(pyg);

    return Py_BuildValue("s#ii",vers,vers_len,nmiss,nskip);
}*/


static PyObject *
tblot(PyObject *obj, PyObject *args)
{
  /* Arguments in the order they appear */
  PyObject *oimg, *oout;
  long xmin, xmax, ymin, ymax;
  double scale;
  float kscale;
  double xscale, yscale;
  char *align_str, *interp_str;
  float ef, misval, sinscl;
  long vflag;
  PyObject *callback_obj = NULL;

  PyArrayObject *img = NULL, *out = NULL;
  enum e_align_t align;
  enum e_interp_t interp;
  mapping_callback_t callback = NULL;
  void *callback_state = NULL;
  long nx,ny,onx,ony;
  int istat = 0;
  struct driz_error_t error;
  struct driz_param_t p;

  driz_error_init(&error);

  if (!PyArg_ParseTuple(args,"OOlllldfddssffflO:tblot", &oimg, &oout, &xmin,
                        &xmax, &ymin, &ymax, &scale, &kscale, &xscale,
                        &yscale, &align_str, &interp_str, &ef, &misval,
                        &sinscl, &vflag, &callback_obj)){
    return PyErr_Format(gl_Error, "cdriz.tblot: Invalid Parameters.");
  }

  /* Check for invalid scale */
  if (scale == 0.0) {
    driz_error_format_message(&error, "Invalid scale %f (must be non-zero)", scale);
    goto _exit;
  }

  if (kscale == 0.0) {
    driz_error_format_message(&error, "Invalid kscale %f (must be non-zero)", scale);
    goto _exit;
  }

  callback = py_mapping_callback;
  callback_state = (void *)callback_obj;

  img = (PyArrayObject *)PyArray_ContiguousFromAny(oimg, NPY_FLOAT32, 2, 2);
  if (!img) {
    driz_error_set_message(&error, "Invalid input array");
    goto _exit;
  }
  out = (PyArrayObject *)PyArray_ContiguousFromAny(oout, NPY_FLOAT32, 2, 2);
  if (!out) {
    driz_error_set_message(&error, "Invalid output array");
    goto _exit;
  }

  if (align_str2enum(align_str, &align, &error) ||
      interp_str2enum(interp_str, &interp, &error)) {
    goto _exit;
  }

  nx = PyArray_DIMS(img)[1];
  ny = PyArray_DIMS(img)[0];
  onx = PyArray_DIMS(out)[1];
  ony = PyArray_DIMS(out)[0];

  driz_param_init(&p);

  p.data = PyArray_DATA(img);
  p.output_data = PyArray_DATA(out);
  p.xmin = xmin;
  p.xmax = xmax;
  p.ymin = ymin;
  p.ymax = ymax;
  p.dnx = nx;
  p.dny = ny;
  p.onx = onx;
  p.ony = ony;
  p.scale = scale;
  p.kscale = kscale;
  p.x_scale = xscale;
  p.y_scale = yscale;
  p.in_units = unit_cps;
  p.align = align;
  p.interpolation = interp;
  p.ef = ef;
  p.misval = misval;
  p.sinscl = sinscl;
  p.mapping_callback = callback;
  p.mapping_callback_state = callback_state;

  istat = doblot(&p, &error);

 _exit:
  Py_DECREF(img);
  Py_DECREF(out);

  if (istat || driz_error_is_set(&error)) {
    if (strcmp(driz_error_get_message(&error), "<PYTHON>") != 0)
      PyErr_SetString(PyExc_Exception, driz_error_get_message(&error));
    return NULL;
  } else {
    return Py_BuildValue("i",istat);
  }
}


/* To replace the default prinf log; instead log to a pythonic log */
void cdriz_log_func(const char *format, ...) {
  static PyObject *logging = NULL;
  va_list args;
  PyObject *logger;
  PyObject *string;
  char msg[256];
  int n;

  va_start(args, format);

  if (logging == NULL) {
    logging = PyImport_ImportModuleNoBlock("logging");
    if (logging == NULL) return;
  }

  n = PyOS_vsnprintf(msg, sizeof(msg), format, args);

  va_end(args);

  if (n < 0) {
    /* XXX: An error occurred in string formatting; just ignore for now */
    return;
  }

  /* XXX: Provide a way to specify the log level to use */
  string = Py_BuildValue("s", msg);
  if (string == NULL) return;

  logger = PyObject_CallMethod(logging, "getLogger", "s",
                               "drizzlepac.cdriz");
  if (logger == NULL) {
      Py_XDECREF(string);
      return;
  }

  PyObject_CallMethod(logger, "info", "O", string);

  Py_XDECREF(logger);
  Py_XDECREF(string);
  return;
}


static PyObject *
arrmoments(PyObject *obj, PyObject *args)
{
  /* Arguments in the order they appear */
  PyObject *oimg;
  long p,q ;

  /* Derived values */
  PyArrayObject *img = NULL;
  long x,y;
  double moment = 0.0;
  integer_t i,j;

  if (!PyArg_ParseTuple(args,"Oll:arrmoments", &oimg, &p, &q)){
    return PyErr_Format(gl_Error, "cdriz.arrmoments: Invalid Parameters.");
  }

  img = (PyArrayObject *)PyArray_ContiguousFromAny(oimg, NPY_FLOAT32, 2, 2);
  if (!img) {
    goto _exit;
  }

  x = PyArray_DIMS(img)[1];
  y = PyArray_DIMS(img)[0];

  /* Perform computation */
  for (i = 0; i < y; i++)
    for (j = 0; j < x; j++)
      moment += pow(i,p) * pow(j,q)* (*(float *)((char*)PyArray_DATA(img) +
                j*PyArray_STRIDES(img)[1] + i*PyArray_STRIDES(img)[0]));

 _exit:
  Py_DECREF(img);

  return Py_BuildValue("d", moment);
}

static PyObject *
arrxyround(PyObject *obj, PyObject *args)
{
  /* Arguments (mostly) in the order they appear */
  PyObject *oimg, *oker2d;
  double x0,y0, skymode ;
  double xsigsq,ysigsq, datamin,datamax;

  /* Derived values */
  PyArrayObject *img = NULL;
  PyArrayObject *ker2d = NULL;
  long nyk,nxk;

  double xc, yc, round;
  long j, k;
  double xhalf, yhalf;
  long xmiddle, ymiddle;
  double sg, sd, wt;
  double sumgd, sumgsq, sumg, sumd, sumdx;
  double sdgdx, sdgdxsq, sddgdx, sgdgdx;
  double p;
  float pixval, ker2dval;
  long px,py;
  long n;
  double dxk, dgdx, dyj;
  double hx, hy, hx1, hy1, dy, dx, dy1, skylvl;

  if (!PyArg_ParseTuple(args,"OdddOdddd:arrxyround", &oimg, &x0, &y0, &skymode,
                        &oker2d, &xsigsq, &ysigsq, &datamin, &datamax)){
    return PyErr_Format(gl_Error, "cdriz.arrxyround: Invalid Parameters.");
  }

  img = (PyArrayObject *)PyArray_ContiguousFromAny(oimg, NPY_FLOAT32, 2, 2);
  if (!img) return Py_BuildValue("");

  ker2d = (PyArrayObject *)PyArray_ContiguousFromAny(oker2d, NPY_FLOAT64, 2, 2);
  if (!ker2d) {
    Py_DECREF(img);
    return Py_BuildValue("");
  }

  nxk = PyArray_DIMS(ker2d)[1];
  nyk = PyArray_DIMS(ker2d)[0];

  /* Perform computation */
  /*These are all 0-based indices */
  xhalf = (nxk / 2.0) - 0.5;
  yhalf = (nyk / 2.0) - 0.5;
  xmiddle = (int)floor(nxk / 2);
  ymiddle = (int)floor(nyk / 2);

  /* Initialize the x fit. */
  sumgd = 0.0;
  sumgsq = 0.0;
  sumg = 0.0;
  sumd = 0.0;
  sumdx = 0.0;
  sdgdx = 0.0;
  sdgdxsq = 0.0;
  sddgdx = 0.0;
  sgdgdx = 0.0;
  p = 0.0;
  n = 0;

  /* Compute the sums required for the x fit. */
  for (k=0; k < nxk; k++){
      sg = 0.0;
      sd = 0.0;
      for (j=0; j < nyk; j++){
          wt = (float)(ymiddle+1 - labs (j - ymiddle));
          px = x0-xmiddle+k;
          py = y0-ymiddle+j;
          pixval = *(float *)((char*)PyArray_DATA(img) + px*PyArray_STRIDES(img)[1] + py*PyArray_STRIDES(img)[0]);
          /* pixval = data[y0-ymiddle+j,x0-xmiddle+k]; */
          if ((pixval < datamin) || (pixval > datamax)){
              sg=DBL_MIN;
              break;
          }

          sd += (pixval - skymode) * wt;
          ker2dval = *(double *)((char*)PyArray_DATA(ker2d) + k*PyArray_STRIDES(ker2d)[1] + j*PyArray_STRIDES(ker2d)[0]);
          sg += ker2dval * wt;
      }
      if (sg == DBL_MIN){
          break;
      }
      dxk = xmiddle-k;
      wt = (float)(xmiddle+1L - labs(xmiddle-k));
      sumgd += wt * sg * sd;
      sumgsq += wt * pow(sg, 2);
      sumg += wt * sg;
      sumd += wt * sd;
      sumdx += wt * sd * dxk;
      p += wt;
      n += 1;
      dgdx = sg * dxk;
      sdgdxsq += wt * pow(dgdx, 2);
      sdgdx += wt * dgdx;
      sddgdx += wt * sd * dgdx;
      sgdgdx += wt * sg * dgdx;
  }
  /*
  Need at least three points to estimate the x height, position
  and local sky brightness of the star.
  */
  if ( (sg == DBL_MIN) || ((n <= 2) || (p <= 0.0))){
    Py_DECREF(img);
    Py_DECREF(ker2d);
    return Py_BuildValue("");
  }

  /*
  Solve for the height of the best-fitting gaussian to the
  xmarginal. Reject the star if the height is non-positive.
  */
  hx1 = sumgsq - (pow(sumg,2)) / p;
  if (hx1 <= 0.0){
    Py_DECREF(img);
    Py_DECREF(ker2d);
    return Py_BuildValue("");
  }

  hx = (sumgd - sumg * sumd / p) / hx1;
  if (hx <= 0.0){
    Py_DECREF(img);
    Py_DECREF(ker2d);
    return Py_BuildValue("");
  }

  /* Solve for the new x centroid. */
  skylvl = (sumd - hx * sumg) / p;
  dx = (sgdgdx - (sddgdx - sdgdx * (hx * sumg + skylvl * p))) / (hx * sdgdxsq / xsigsq);

  if (fabs(dx) > xhalf){
      if (sumd == 0.0){
          dx = 0.0;
      } else{
          dx = sumdx / sumd;
      }
      if (fabs(dx) > xhalf){
          dx = 0.0;
      }
  }
  xc = (int)floor(x0) + dx;

  /* Initialize y fit. */
  sumgd = 0.0;
  sumgsq = 0.0;
  sumg = 0.0;
  sumd = 0.0;
  sumdx = 0.0;
  sdgdx = 0.0;
  sdgdxsq = 0.0;
  sddgdx = 0.0;
  sgdgdx = 0.0;
  p = 0.0;
  n = 0;

  for (j=0; j < nyk; j++){
      sg = 0.0;
      sd = 0.0;
      for (k=0; k<nxk; k++){
          wt = (float)(xmiddle+1L - labs(k - xmiddle));
          px = x0-xmiddle+k;
          py = y0-ymiddle+j;
          pixval = *(float *)((char*)PyArray_DATA(img) + px*PyArray_STRIDES(img)[1] + py*PyArray_STRIDES(img)[0]);
          /* pixval = data[y0-ymiddle+j,x0-xmiddle+k]; */
          if ((pixval < datamin) || (pixval > datamax)){
              sg = DBL_MIN;
              break;
          }
          sd += (pixval - skymode) * wt;
          ker2dval = *(double *)((char*)PyArray_DATA(ker2d) + k*PyArray_STRIDES(ker2d)[1] + j*PyArray_STRIDES(ker2d)[0]);
          sg += ker2dval * wt;
      }
      if (sg == DBL_MIN){
          break;
      }
      dyj = ymiddle-j;
      wt = (float)(ymiddle+1L - labs(j - ymiddle));

      sumgd += wt * sg * sd;
      sumgsq += wt * pow(sg, 2);
      sumg += wt * sg;
      sumd += wt * sd;
      sumdx += wt * sd * dyj;
      p = p + wt;
      n = n + 1;
      dgdx = sg * dyj;
      sdgdx += wt * dgdx;
      sdgdxsq += wt * pow(dgdx, 2);
      sddgdx += wt * sd * dgdx;
      sgdgdx += wt * sg * dgdx;
  }
  /*
   Need at least three points to estimate the y height, position
   and local sky brightness of the star.
  */

  if ((sg == DBL_MIN) || ((n <= 2) || (p <= 0.0))){
    Py_DECREF(img);
    Py_DECREF(ker2d);
    return Py_BuildValue("");
  }
  /*
  Solve for the height of the best-fitting gaussian to the
  y marginal. Reject the star if the height is non-positive.
  */

  hy1 = sumgsq - (pow(sumg, 2) / p);
  if (hy1 <= 0.0){
    Py_DECREF(img);
    Py_DECREF(ker2d);
    return Py_BuildValue("");
  }

  hy = (sumgd - ((sumg * sumd) / p)) / hy1;
  if (hy <= 0.0){
    Py_DECREF(img);
    Py_DECREF(ker2d);
    return Py_BuildValue("");
  }

  /* Solve for the new x centroid. */
  skylvl = (sumd - hy * sumg) / p;
  dy1 = sgdgdx - (sddgdx - (sdgdx * ((hy * sumg) + (skylvl * p))));
  dy = dy1 / (hy * sdgdxsq / ysigsq);

  if (fabs(dy) > yhalf) {
      if (sumd == 0.0)
        dy = 0.0;
      else
        dy = sumdx / sumd;

      if (fabs(dy) > yhalf)
        dy = 0.0;
  }
  yc = (int)floor(y0) + dy;

  round = 2.0 * (hx - hy) / (hx + hy);

  Py_DECREF(img);
  Py_DECREF(ker2d);
  return Py_BuildValue("ddd", xc, yc, round);
}

/* ==== Allocate a double *vector (vec of pointers) ======================
    Memory is Allocated!  See void free_Carray(double ** )                  */
double **ptrvector(long n)  {
    double **v;
    v=(double **)malloc((size_t) (n*sizeof(double)));
    if (!v)   {
        printf("In **ptrvector. Allocation of memory for double array failed.");
        exit(0);
    }
    return v;
}


/* ==== Create Carray from PyArray ======================
    Assumes PyArray is contiguous in memory.
    Memory is allocated!                                    */
double **pymatrix_to_Carrayptrs(PyArrayObject *arrayin)  {
    double **c, *a;
    integer_t i,n,m;

    n=PyArray_DIMS(arrayin)[0];
    m=PyArray_DIMS(arrayin)[1];

    c=(double **)ptrvector(n);
    a=(double *) PyArray_DATA(arrayin);  /* pointer to arrayin data as double */
    for ( i=0; i<n; i++)
      c[i]=a+i*m;
    return c;
}

/* ==== Free a double *vector (vec of pointers) ========================== */
void free_Carrayptrs(double **v)  {
    free((char*) v);
}

static PyObject *
arrxyzero(PyObject *obj, PyObject *args)
{
  /* Arguments (mostly) in the order they appear */
  PyObject *oimgxy, *orefxy;
  double searchrad;

  /* Derived values */
  PyArrayObject *imgxy = NULL;
  PyArrayObject *refxy = NULL;
  PyArrayObject *ozpmat = NULL;
  double **zpmat = NULL;

  integer_t imgnum, refnum;
  npy_intp dimensions[2];
  integer_t xind, yind;
  double dx, dy;
  integer_t j, k;

  if (!PyArg_ParseTuple(args,"OOd:arrxyzero", &oimgxy, &orefxy, &searchrad)){
    return PyErr_Format(gl_Error, "cdriz.arrxyzero: Invalid Parameters.");
  }

  imgxy = (PyArrayObject *)PyArray_ContiguousFromAny(oimgxy, NPY_FLOAT32, 2, 2);
  if (!imgxy)
    goto _exit;

  refxy = (PyArrayObject *)PyArray_ContiguousFromAny(orefxy, NPY_FLOAT32, 2, 2);
  if (!refxy)
    goto _exit;

  dimensions[0] = (integer_t)(searchrad*2) + 1;
  dimensions[1] = (integer_t)(searchrad*2) + 1;
  ozpmat = (PyArrayObject *) PyArray_SimpleNew(2, dimensions, NPY_DOUBLE);
  if (!ozpmat)
    goto _exit;
  PyArray_FILLWBYTE(ozpmat, 0);

  /* Allocate memory for return matrix */
  zpmat = pymatrix_to_Carrayptrs(ozpmat);

  imgnum = PyArray_DIMS(imgxy)[0];
  refnum = PyArray_DIMS(refxy)[0];

  /* For each entry in the input image...*/
  for (j=0; j< imgnum; j++)
    /* compute the delta relative to each source in ref image */
    for (k = 0; k < refnum; k++){
        dx = *(float *)((char*)PyArray_DATA(imgxy) + j*PyArray_STRIDES(imgxy)[0]) -
             *(float *)((char*)PyArray_DATA(refxy) + k*PyArray_STRIDES(refxy)[0]);
        dy = *((float *)((char*)PyArray_DATA(imgxy) +
               j*PyArray_STRIDES(imgxy)[0]+ PyArray_STRIDES(imgxy)[1]))
           - *((float *)((char*)PyArray_DATA(refxy) +
               k*PyArray_STRIDES(refxy)[0]+ PyArray_STRIDES(refxy)[1]));
        if ((fabs(dx) < searchrad) && (fabs(dy) < searchrad)) {
            xind = (integer_t)(dx+searchrad);
            yind = (integer_t)(dy+searchrad);
            zpmat[yind][xind] += 1;
        }
    }

 _exit:
  Py_DECREF(imgxy);
  Py_DECREF(refxy);
  free_Carrayptrs(zpmat);

  return PyArray_Return(ozpmat);
}

static PyMethodDef cdriz_methods[] =
  {
    {"tdriz",  tdriz, METH_VARARGS, "tdriz(image, weight, output, outweight, context, uniqid, ystart, xmin, ymin, dny, scale, xscale, yscale, align, pfrace, kernel, inun, expin, wtscl, fill, nmiss, nskip, vflag, callback)"},
    /*{"twdriz",  tdriz, METH_VARARGS, "triz(image, weight, output, outweight, ystart, xmin, ymin, dny, wcsin, wcsout,pxg,pyg,pfract, kernel, coeffs, fillstr,nmiss,nskip,vflag)"},*/
    {"tblot",  tblot, METH_VARARGS, "tblot(image, output, xmin, xmax, ymin, ymax, scale, kscale, xscale, yscale, align, interp, ef, misval, sinscl, vflag, callback)"},
    {"arrmoments", arrmoments, METH_VARARGS, "arrmoments(image, p, q)"},
    {"arrxyround", arrxyround, METH_VARARGS, "arrxyround(data,x0,y0,skymode,ker2d,xsigsq,ysigsq,datamin,datamax)"},
    {"arrxyzero", arrxyzero, METH_VARARGS, "arrxyzero(imgxy,refxy,searchrad,zpmat)"},
    {0, 0, 0, 0}                             /* sentinel */
  };

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "cdriz",
    NULL,
    -1,
    cdriz_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyObject *PyInit_cdriz(void)
{
  PyObject* m;
  driz_log_func = &cdriz_log_func;

  if (PyType_Ready(&WCSMapType) < 0) {
    return NULL;
  }
  m = PyModule_Create(&moduledef);
  if (m == NULL) {
    return NULL;
  }

  import_array();
  import_astropy_wcs();

  Py_INCREF(&WCSMapType);
  PyModule_AddObject(m, "DefaultWCSMapping", (PyObject *)&WCSMapType);

  return m;
}
