
#include <Python.h>

#define _USE_MATH_DEFINES       /* needed for MS Windows to define M_PI */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <stdio.h>

#include <numpy/arrayobject.h>

#include "astropy_wcs_api.h"
#include "astropy_wcs.h"

#include "cdrizzleblot.h"
#include "cdrizzlebox.h"
#include "cdrizzlemap.h"
#include "cdrizzleutil.h"
#include "tests/drizzletest.h"

static PyObject *gl_Error;

/**
 
 The Python mapping callback must have the following signature:

    def mapping(xin, yin):
        return xout, yout

 xin, yin are the input coordinates, and xout and yout are the output
 coordinates.  All are 1-dimensional Numpy DOUBLE arrays of the same
 length.
*/

static PyObject *
tdriz(PyObject *obj UNUSED_PARAM, PyObject *args)
{
  /* Arguments in the order they appear */
  PyObject *oimg, *owei, *oout, *owht, *ocon;
  long uniqid, ystart, xmin, xmax, ymin, ymax;
  double scale, xscale, yscale;
  char *align_str;
  double pfract;
  char *kernel_str, *inun_str;
  float expin, wtscl;
  char *fillstr;
  integer_t nmiss, nskip, vflag;
  PyObject *pixmap;

  /* Derived values */
  PyArrayObject *img = NULL, *wei = NULL, *out = NULL, *wht = NULL, *con = NULL, *map = NULL;
  enum e_align_t align;
  enum e_kernel_t kernel;
  enum e_unit_t inun;
  char *fillstr_end;
  bool_t do_fill;
  float fill_value;
  int istat = 0;
  struct driz_error_t error;
  struct driz_param_t p;

  /* clock_t start_t, end_t; */
  /* double delta_time; */

  driz_error_init(&error);

  if (!PyArg_ParseTuple(args,"OOOOOllllldddsdssffsiiiO:tdriz",
                        &oimg, &owei, &oout, &owht, &ocon, &uniqid, 
                        &xmin, &xmax, &ymin, &ymax, &scale, &xscale, &yscale,
                        &align_str, &pfract, &kernel_str, &inun_str,
                        &expin, &wtscl, &fillstr, &nmiss, &nskip, &vflag,
                        &pixmap)) {
    return PyErr_Format(gl_Error, "cdriz.tdriz: Invalid Parameters.");
  }

  /* Get raw C-array data */
  img = (PyArrayObject *)PyArray_ContiguousFromAny(oimg, PyArray_FLOAT, 2, 2);
  if (!img) {
    driz_error_set_message(&error, "Invalid input array");
    goto _exit;
  }

  wei = (PyArrayObject *)PyArray_ContiguousFromAny(owei, PyArray_FLOAT, 2, 2);
  if (!wei) {
    driz_error_set_message(&error, "Invalid weights array");
    goto _exit;
  }

  out = (PyArrayObject *)PyArray_ContiguousFromAny(oout, PyArray_FLOAT, 2, 2);
  if (!out) {
    driz_error_set_message(&error, "Invalid output array");
    goto _exit;
  }

  wht = (PyArrayObject *)PyArray_ContiguousFromAny(owht, PyArray_FLOAT, 2, 2);
  if (!wht) {
    driz_error_set_message(&error, "Invalid array");
    goto _exit;
  }

  con = (PyArrayObject *)PyArray_ContiguousFromAny(ocon, PyArray_INT32, 2, 2);
  if (!con) {
    driz_error_set_message(&error, "Invalid context array");
    goto _exit;
  }

  map = (PyArrayObject *)PyArray_ContiguousFromAny(pixmap, PyArray_DOUBLE, 3, 3);

  if (!map) {
    driz_error_set_message(&error, "Invalid pixmap array");
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

  nmiss = 0;
  nskip = 0;

  /* Setup reasonable defaults for drizzling */
  driz_param_init(&p);

  p.data = img;
  p.weights = wei;
  p.output_data = out;
  p.output_counts = wht;
  p.output_context = con;
  p.uuid = uniqid;
  p.xmin = xmin;
  p.ymin = ymin;
  p.xmax = xmax;
  p.ymax = ymax;
  p.scale = scale;
  p.x_scale = xscale;
  p.y_scale = yscale;
  p.align = align;
  p.pixel_fraction = pfract;
  p.kernel = kernel;
  p.in_units = inun;
  p.exposure_time = expin;
  p.weight_scale = wtscl;
  p.pixmap = map;
 
  /*
  start_t = clock();
  */
  /* Do the drizzling */
  if (dobox(&p, &nmiss, &nskip, &error)) {
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

 _exit:
  Py_XDECREF(con);
  Py_XDECREF(img);
  Py_XDECREF(wei);
  Py_XDECREF(out);
  Py_XDECREF(wht);
  Py_XDECREF(map);

  if (istat || driz_error_is_set(&error)) {
    if (strcmp(driz_error_get_message(&error), "<PYTHON>") != 0)
      PyErr_SetString(PyExc_Exception, driz_error_get_message(&error));
    return NULL;
  } else {
    return Py_BuildValue("sii", "Callable C-based DRIZZLE Version 0.8 (20th May 2009)", nmiss, nskip);
  }
}

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
  PyObject *pixmap;

  PyArrayObject *img = NULL, *out = NULL, *map = NULL;
  enum e_align_t align;
  enum e_interp_t interp;
  int istat = 0;
  struct driz_error_t error;
  struct driz_param_t p;
  double maxdiff = 0.0;

  driz_error_init(&error);

  if (!PyArg_ParseTuple(args,"OOlllldfddssffflO:tblot", &oimg, &oout, &xmin,
                        &xmax, &ymin, &ymax, &scale, &kscale, &xscale,
                        &yscale, &align_str, &interp_str, &ef, &misval,
                        &sinscl, &vflag, &pixmap)){
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

  img = (PyArrayObject *)PyArray_ContiguousFromAny(oimg, PyArray_FLOAT, 2, 2);
  if (!img) {
    driz_error_set_message(&error, "Invalid input array");
    goto _exit;
  }
  
  out = (PyArrayObject *)PyArray_ContiguousFromAny(oout, PyArray_FLOAT, 2, 2);
  if (!out) {
    driz_error_set_message(&error, "Invalid output array");
    goto _exit;
  }

  map = (PyArrayObject *)PyArray_ContiguousFromAny(pixmap, PyArray_DOUBLE, 3, 3);
  if (!map) {
    driz_error_set_message(&error, "Invalid pixmap array");
    goto _exit;
  }
  
  if (align_str2enum(align_str, &align, &error) ||
      interp_str2enum(interp_str, &interp, &error)) {
    goto _exit;
  }

  driz_param_init(&p);

  p.data = img;
  p.output_data = out;
  p.xmin = xmin;
  p.xmax = xmax;
  p.ymin = ymin;
  p.ymax = ymax;
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
  p.pixmap = map;

  istat = doblot(&p, &error);

 _exit:
  Py_DECREF(img);
  Py_DECREF(out);
  Py_DECREF(map);

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
    if (logging == NULL) {
      return;
    }
  }

  n = PyOS_vsnprintf(msg, sizeof(msg), format, args);

  va_end(args);

  if (n < 0) {
    /* XXX: An error occurred in string formatting; just ignore for now */
    return;
  }

  /* XXX: Provide a way to specify the log level to use */
  string = Py_BuildValue("s", msg);
  if (string == NULL) {
    goto cleanup;
  }

  logger = PyObject_CallMethod(logging, "getLogger", "s",
                               "drizzlepac.cdriz");
  if (logger == NULL) {
      goto cleanup;
  }

  PyObject_CallMethod(logger, "info", "O", string);

cleanup:
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
  double moment;
  integer_t i,j;
  double val;

  if (!PyArg_ParseTuple(args,"Oll:arrmoments", &oimg, &p, &q)){
    return PyErr_Format(gl_Error, "cdriz.arrmoments: Invalid Parameters.");
  }

  img = (PyArrayObject *)PyArray_ContiguousFromAny(oimg, PyArray_FLOAT, 2, 2);
  if (!img) {
    goto _exit;
  }

  x = img->dimensions[1];
  y = img->dimensions[0];

  moment = 0.0;
  /* Perform computation */
  for (i = 0; i < y; i++) {
    for (j = 0; j < x; j++) {
      val = *(float *)(img->data + j*img->strides[1] + i*img->strides[0]);
      moment += pow(i,p)*pow(j,q)*val;
    }
  }

 _exit:
  Py_DECREF(img);

  return Py_BuildValue("d",moment);
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

  double xc,yc,round;
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
  integer_t return_val = 0;


  if (!PyArg_ParseTuple(args,"OdddOdddd:arrxyround", &oimg, &x0, &y0, &skymode,
                        &oker2d, &xsigsq, &ysigsq, &datamin, &datamax)){
    return PyErr_Format(gl_Error, "cdriz.arrxyround: Invalid Parameters.");
  }

  img = (PyArrayObject *)PyArray_ContiguousFromAny(oimg, PyArray_FLOAT, 2, 2);
  if (!img) {
    goto _exit;
  }

  ker2d = (PyArrayObject *)PyArray_ContiguousFromAny(oker2d, PyArray_DOUBLE, 2, 2);
  if (!ker2d) {
    goto _exit;
  }

  nxk = ker2d->dimensions[1];
  nyk = ker2d->dimensions[0];

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
          wt = (float)(ymiddle+1 - abs (j - ymiddle));
          px = x0-xmiddle+k;
          py = y0-ymiddle+j;
          pixval = *(float *)(img->data + px*img->strides[1] + py*img->strides[0]);
          /* pixval = data[y0-ymiddle+j,x0-xmiddle+k]; */
          if ((pixval < datamin) || (pixval > datamax)){
              sg=DBL_MIN;
              break;
          }

          sd += (pixval - skymode) * wt;
          ker2dval = *(double *)(ker2d->data + k*ker2d->strides[1] + j*ker2d->strides[0]);
          sg += ker2dval * wt;
      }
      if (sg == DBL_MIN){
          break;
      }
      dxk = xmiddle-k;
      wt = (float)(xmiddle+1 - abs(dxk));
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
      return_val = -1;
      goto _exit;
  }

  /*
  Solve for the height of the best-fitting gaussian to the
  xmarginal. Reject the star if the height is non-positive.
  */
  hx1 = sumgsq - (pow(sumg,2)) / p;
  if (hx1 <= 0.0){
      return_val = -1;
      goto _exit;
  }
  hx = (sumgd - sumg * sumd / p) / hx1;
  if (hx <= 0.0){
      return_val = -1;
      goto _exit;
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
          wt = (float)(xmiddle+1 - abs(k - xmiddle));
          px = x0-xmiddle+k;
          py = y0-ymiddle+j;
          pixval = *(float *)(img->data + px*img->strides[1] + py*img->strides[0]);
          /* pixval = data[y0-ymiddle+j,x0-xmiddle+k]; */
          if ((pixval < datamin) || (pixval > datamax)){
              sg = DBL_MIN;
              break;
          }
          sd += (pixval - skymode) * wt;
          ker2dval = *(double *)(ker2d->data + k*ker2d->strides[1] + j*ker2d->strides[0]);
          sg += ker2dval * wt;
      }
      if (sg == DBL_MIN){
          break;
      }
      dyj = ymiddle-j;
      wt = (float)(ymiddle+1 - abs(j - ymiddle));

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
      return_val = -1;
      goto _exit;
  }
  /*
  Solve for the height of the best-fitting gaussian to the
  y marginal. Reject the star if the height is non-positive.
  */

  hy1 = sumgsq - (pow(sumg,2) / p);
  if (hy1 <= 0.0){
      return_val = -1;
      goto _exit;
  }
  hy = (sumgd - ((sumg * sumd) / p)) / hy1;
  if (hy <= 0.0){
      return_val = -1;
      goto _exit;
  }

  /* Solve for the new x centroid. */
  skylvl = (sumd - hy * sumg) / p;
  dy1 = sgdgdx - (sddgdx - (sdgdx * ((hy * sumg) + (skylvl * p))));
  dy = dy1 / (hy * sdgdxsq / ysigsq);

  if (fabs(dy) > yhalf){
      if (sumd == 0.0){
          dy = 0.0;
      } else {
          dy = sumdx / sumd;
      }
      if (fabs(dy) > yhalf){
          dy = 0.0;
      }
  }
  yc = (int)floor(y0) + dy;

  round = 2.0 * (hx - hy) / (hx + hy);

 _exit:
  Py_DECREF(img);
  Py_DECREF(ker2d);

  if (return_val < 0){
      return Py_BuildValue("");
  } else {
      return Py_BuildValue("ddd",xc,yc,round);
  }
}

/* ==== Allocate a double *vector (vec of pointers) ======================
    Memory is Allocated!  See void free_Carray(double ** )                  */
double **ptrvector(long n)  {
    double **v;
    v=(double **)malloc((size_t) (n*sizeof(double)));
    if (!v)   {
        printf("In **ptrvector. Allocation of memory for double array failed.");
        exit(0);  }
    return v;
}


/* ==== Create Carray from PyArray ======================
    Assumes PyArray is contiguous in memory.
    Memory is allocated!                                    */
double **pymatrix_to_Carrayptrs(PyArrayObject *arrayin)  {
    double **c, *a;
    long i,n,m;

    n=arrayin->dimensions[0];
    m=arrayin->dimensions[1];
    c=(double **)ptrvector(n);
    a=(double *) arrayin->data;  /* pointer to arrayin data as double */
    for ( i=0; i<n; i++)  {
        c[i]=a+i*m;  }
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
  double **zpmat;
  long *a;

  long imgnum, refnum;
  integer_t dimensions[2];
  integer_t xind, yind;
  double dx, dy;
  long j, k;
  long nsource = 0;

  if (!PyArg_ParseTuple(args,"OOd:arrxyzero", &oimgxy, &orefxy, &searchrad)){
    return PyErr_Format(gl_Error, "cdriz.arrxyzero: Invalid Parameters.");
  }

  imgxy = (PyArrayObject *)PyArray_ContiguousFromAny(oimgxy, PyArray_FLOAT, 2, 2);
  if (!imgxy) {
    goto _exit;
  }

  refxy = (PyArrayObject *)PyArray_ContiguousFromAny(orefxy, PyArray_FLOAT, 2, 2);
  if (!refxy) {
    goto _exit;
  }

  dimensions[0] = (integer_t)(searchrad*2) + 1;
  dimensions[1] = (integer_t)(searchrad*2) + 1;
  ozpmat = (PyArrayObject *)PyArray_FromDims(2, dimensions, NPY_DOUBLE);
  if (!ozpmat) {
    goto _exit;
  }
  /* Allocate memory for return matrix */
  zpmat=pymatrix_to_Carrayptrs(ozpmat);

  imgnum = imgxy->dimensions[0];
  refnum = refxy->dimensions[0];

  /* For each entry in the input image...*/
  for (j=0; j< imgnum; j++){
    /* compute the delta relative to each source in ref image */
    for (k = 0; k < refnum; k++){
        dx = *(float *)(imgxy->data + j*imgxy->strides[0]) - *(float *)(refxy->data + k*refxy->strides[0]);
        dy = *(float *)(imgxy->data + j*imgxy->strides[0]+ imgxy->strides[1]) -
             *(float *)(refxy->data + k*refxy->strides[0]+ refxy->strides[1]);
        if ((fabs(dx) < searchrad) && (fabs(dy) < searchrad)) {
            xind = (integer_t)(dx+searchrad);
            yind = (integer_t)(dy+searchrad);
            zpmat[yind][xind] += 1;
        }
    }
  }

 _exit:
  Py_DECREF(imgxy);
  Py_DECREF(refxy);
  free_Carrayptrs(zpmat);

  return PyArray_Return(ozpmat);
}

static PyObject *
test_cdrizzlepac(PyObject *self, PyObject *args)
{
  PyObject *data, *weights, *pixmap, *output_data, *output_counts, *output_context;
  PyArrayObject *dat, *wei, *map, *odat, *ocnt, *ocon;

  int argc = 1;
  char *argv[] = {"utest_cdrizzlepac", NULL};
  
  if (!PyArg_ParseTuple(args,"OOOOOO:test_cdrizzlepac", &data, &weights, &pixmap,
                                          &output_data, &output_counts, &output_context)) {
    return PyErr_Format(gl_Error, "cdriz.test_cdrizzlepac: Invalid Parameters.");
  }

  dat = (PyArrayObject *)PyArray_ContiguousFromAny(data, PyArray_FLOAT, 2, 2);
  if (! dat) {
    return PyErr_Format(gl_Error, "cdriz.test_cdrizzlepac: Invalid Data Array.");
  }

  wei = (PyArrayObject *)PyArray_ContiguousFromAny(weights, PyArray_FLOAT, 2, 2);
  if (! wei) {
    return PyErr_Format(gl_Error, "cdriz.test_cdrizzlepac: Invalid Weghts Array.");
  }

  map = (PyArrayObject *)PyArray_ContiguousFromAny(pixmap, PyArray_DOUBLE, 3, 3);
  if (! map) {
    return PyErr_Format(gl_Error, "cdriz.test_cdrizzlepac: Invalid Pixmap.");
  }
  
  odat = (PyArrayObject *)PyArray_ContiguousFromAny(output_data, PyArray_FLOAT, 2, 2);
  if (! odat) {
    return PyErr_Format(gl_Error, "cdriz.test_cdrizzlepac: Invalid Output Data Array.");
  }

  ocnt = (PyArrayObject *)PyArray_ContiguousFromAny(output_counts, PyArray_FLOAT, 2, 2);
  if (! ocnt) {
    return PyErr_Format(gl_Error, "cdriz.test_cdrizzlepac: Invalid Output Counts Array.");
  }

  ocon = (PyArrayObject *)PyArray_ContiguousFromAny(output_context, PyArray_INT32, 2, 2);
  if (! ocon) {
    return PyErr_Format(gl_Error, "Invalid context array");
  }

  set_test_arrays(dat, wei, map, odat, ocnt, ocon);
  utest_cdrizzlepac(argc, argv);
  
  return Py_BuildValue("");

}

static PyMethodDef cdriz_methods[] =
  {
    {"tdriz",  tdriz, METH_VARARGS, "tdriz(image, weight, output, outweight, context, uniqid, ystart, xmin, ymin, scale, xscale, yscale, align, pfrace, kernel, inun, expin, wtscl, fill, nmiss, nskip, vflag, pixmap)"},
    {"tblot",  tblot, METH_VARARGS, "tblot(image, output, xmin, xmax, ymin, ymax, scale, kscale, xscale, yscale, align, interp, ef, misval, sinscl, vflag, pixmap)"},
    {"arrmoments", arrmoments, METH_VARARGS, "arrmoments(image, p, q)"},
    {"arrxyround", arrxyround, METH_VARARGS, "arrxyround(data,x0,y0,skymode,ker2d,xsigsq,ysigsq,datamin,datamax)"},
    {"arrxyzero", arrxyzero, METH_VARARGS, "arrxyzero(imgxy,refxy,searchrad,zpmat)"},
    {"test_cdrizzlepac", test_cdrizzlepac, METH_VARARGS, "test_cdrizzlepac(data, weights, pixmap, output_data, output_counts)"},
    {0, 0, 0, 0}                             /* sentinel */
  };

void initcdriz(void)
{
  PyObject* m;

  driz_log_func = &cdriz_log_func;

  m = Py_InitModule("cdriz", cdriz_methods);
  if (m == NULL)
    return;

  import_array();
  import_astropy_wcs();

}
