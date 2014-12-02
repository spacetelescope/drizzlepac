/* Declarations for test functions */

#include <Python.h>
#include <numpy/arrayobject.h>

int do_kernel_square(struct driz_param_t* p, const integer_t j, const integer_t x1, const integer_t x2,
                     integer_t* nmiss, struct driz_error_t* error);

void set_test_arrays(PyArrayObject *dat, PyArrayObject *wei, PyArrayObject *map,
                     PyArrayObject *odat, PyArrayObject *ocnt, PyArrayObject *ocon);

int utest_cdrizzlepac(int argc, char* argv[]);
