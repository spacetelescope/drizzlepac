#include <Python.h>
#include <assert.h>
#define _USE_MATH_DEFINES       /* needed for MS Windows to define M_PI */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "driz_portability.h"
#include "cdrizzleutil.h"

#ifdef _WIN32
#define ssize_t SSIZE_T
#endif

#ifdef WITH_CFITSIO
/**
Read Spitzer-format WCS distortion coefficients from an image header
and populate the distortion arrays.

@param[in] input_data_file A \a cfitsio data file pointer, already opened

@param[out] coeff_type The type of coefficients being stored (TODO: Document
me in a central place)

@param[out] num_coeffs The number of coefficients being stored

@param[out] x_coeffs An array of x coefficients.  This array should be
            allocated to size MAX_COEFFS, but only the first
            num_coeffs will be valid after returning from this
            function.

@param[out] y_coeffs An array of y coefficients.  This array should be
            allocated to size MAX_COEFFS, but only the first
            num_coeffs will be valid after returning from this
            function.

@param[out] error

@return Non-zero if an error occurred.

was: GTSPCO
*/
static int
get_spitzer_coefficients(fitsfile* input_data_file,
                         /* Output arguments */
                         integer_t* coeff_type,
                         integer_t* num_coeffs,
                         double* x_coeffs /* [num_coeffs] */,
                         double* y_coeffs /* [num_coeffs] */,
                         struct driz_error_t* error) {
  char* ctype1;
  char* ctype2;
  double crpix1, crpix2;
  int x_order, y_order, order;
  int nc, n, m;
  char key[16];
  int status = 0;

  assert(input_data_file);
  assert(coeff_type);
  assert(num_coeffs);
  assert(x_coeffs);
  assert(y_coeffs);
  assert(error);

  /* First try to get the CTYPE keywords and confirm they are as
     expected */
  if (fits_read_key(input_data_file, TSTRING, "CTYPE1", &ctype1, NULL, &status)) {
    driz_error_set_message(error, "Could not read CTYPE1 record");
    return 1;
  }

  if (strncmp(ctype1, "RA---TAN-SIP", 16) != 0) {
    driz_error_format_message(error,
                              "SIP (Spitzer) style CTYPE1 not found (got '%s')",
                              ctype1);
    return 1;
  }

  if (fits_read_key(input_data_file, TSTRING, "CTYPE2", &ctype2, NULL, &status)) {
    driz_error_set_message(error, "Could not read CTYPE2 record");
    return 1;
  }

  if (strncmp(ctype2, "DEC--TAN-SIP", 16) != 0) {
    /* TODO: Improve exception */
    driz_error_format_message(error,
                              "SIP (Spitzer) style CTYPE2 not found (got '%s')",
                              ctype2);
    return 1;
  }

  /* Get the reference pixel */
  if (fits_read_key(input_data_file, TDOUBLE, "CRPIX1", &crpix1, NULL, &status)) {
    driz_error_set_message(error, "Could not read CRPIX1 record");
    return 1;
  }

  if (fits_read_key(input_data_file, TDOUBLE, "CRPIX2", &crpix2, NULL, &status)) {
    driz_error_set_message(error, "Could not read CRPIX2 record");
    return 1;
  }

  /* Now look for the order in X and Y */
  if (fits_read_key(input_data_file, TINT, "A_ORDER", &x_order, NULL, &status)) {
    driz_error_set_message(error, "Could not read A_ORDER record");
    return 1;
  }

  if (fits_read_key(input_data_file, TINT, "B_ORDER", &y_order, NULL, &status)) {
    driz_error_set_message(error, "Could not read B_ORDER record");
    return 1;
  }

  /* Set the order to the higher of the two */
  order = MAX(x_order, y_order);

  /* Offset to show we have an explicit refpix */
  *coeff_type = MAX_COEFFS + order;
  *num_coeffs = (order + 1) * (order + 2) / 2 + 1;
  if (*num_coeffs > MAX_COEFFS) {
    driz_error_set_message(error, "Too many coefficients");
    return 1;
  }

  /* And set the refpix to the CRPIXs */
  x_coeffs[*num_coeffs-1] = crpix1;
  y_coeffs[*num_coeffs-1] = crpix2;

  /* Read all the coefficients and set any ommitted ones to zero */
  /* X first */
  nc = 0;
  for (n = 1; n <= order + 1; ++n) {
    for (m = 1; m <= n; ++m) {
      snprintf(key, 16, "A_%d_%d", n-m, m-1);
      if (fits_read_key(input_data_file, TDOUBLE, key, x_coeffs + nc, NULL, &status)) {
        driz_error_format_message(error, "Could not read record '%s'", key);
        return 1;
      }
      snprintf(key, 16, "B_%d_%d", n-m, m-1);
      if (fits_read_key(input_data_file, TDOUBLE, key, y_coeffs + nc, NULL, &status)) {
        driz_error_format_message(error, "Could not read record '%s'", key);
        return 1;
      }
      ++nc;
    }
  }

  /* The Spitzer coefficients are offsets, so we need to add an extra X and Y */
  x_coeffs[1] += 1.0;
  y_coeffs[2] += 1.0;

  return 0;
}

/**
This routine is invoked when it has been requested that the
coefficients file name for geometric distortion, and the associated
wavelength is to be extracted from the header.

It is currently handles \a WFPC2, \a STIS and \a NICMOS and assumes
coefficients in a directory \a drizzle$coeffs.

@param[in] input_data_file A cfitsio data file pointer, already opened

@param[in] coeff_filename A buffer that will be filled with the
filename of the coefficients file on return.

@param[in] coeff_filename_length The length of the coeff_filename
buffer.  It must be at least 255.

@param[out] lambda The wavelength, in nm.

@param[out] error

@return Non-zero if an error has occurred.

was: GTCOEF
*/
static int
get_geometric_distortion_filename_from_header(fitsfile* input_data_file,
                                              /* Output arguments */
                                              char* coeff_filename,
                                              const size_t coeff_filename_length,
                                              double* lambda,
                                              struct driz_error_t* error) {
  int status = 0;
  int detector;
  double plam;
  char* instrument;
  char* detector_name;

  assert(input_data_file);
  assert(coeff_filename);
  assert(coeff_filename_length >= 255);
  assert(lambda);
  assert(error);

  /* First try to get the instrument name */
  if (fits_read_key(input_data_file, TSTRING, "INSTRUME", &instrument, NULL, &status)) {
    driz_error_set_message(error, "Could not read INSTRUME record");
    return 1;
  }

  /* Check instrument can be handled */
  if (strncmp(instrument, "WFPC2", 8) == 0) {
    if (fits_read_key(input_data_file, TINT, "DETECTOR", &detector, NULL, &status)) {
      driz_error_set_message(error, "Could not read DETECTOR record");
      return 1;
    }

    if (detector < 1 || detector > 4) {
      driz_error_format_message(error,
                                "Invalid detector number in header (must be 1-4, got %d)",
                                detector);
      return 1;
    }

    if (fits_read_key(input_data_file, TDOUBLE, "PHOTPLAM", &plam, NULL, &status)) {
      driz_error_set_message(error, "Could not read PHOTPLAM record");
      return 1;
    }

    /* We have to convert the wavelength from Angstroms to nm */
    *lambda = plam / 1.0;

    /* Now we build the full "Trauger style" name for the coefficients */
    switch (detector) {
    case 1:
      strncpy(coeff_filename, "drizzle$coeffs/pc1-trauger", coeff_filename_length);
      break;
    case 2:
      strncpy(coeff_filename, "drizzle$coeffs/wf2-trauger", coeff_filename_length);
      break;
    case 3:
      strncpy(coeff_filename, "drizzle$coeffs/wf3-trauger", coeff_filename_length);
      break;
    case 4:
      strncpy(coeff_filename, "drizzle$coeffs/wf4-trauger", coeff_filename_length);
      break;
    }
  } else if (strncmp(instrument, "NICMOS", 8) == 0) {
    if (fits_read_key(input_data_file, TINT, "CAMERA", &detector, NULL, &status)) {
      driz_error_set_message(error, "Could not read CAMERA record");
      return 1;
    }

    if (detector < 1 || detector > 3) {
      driz_error_format_message(error,
                                "Invalid camera number in header (must be 1-3, got %d)",
                                detector);
      return 1;
    }

    snprintf(coeff_filename, coeff_filename_length,
             "drizzle$coeffs/nic-%d", detector);
  } else if (strncmp(instrument, "STIS", 8) == 0) {
    if (fits_read_key(input_data_file, TSTRING, "DETECTOR", &detector_name, NULL, &status)) {
      driz_error_set_message(error, "Could not read DETECTOR record");
      return 1;
    }

    if (strncmp(detector_name, "CCD", 8) == 0) {
      strncpy(coeff_filename, "drizzle$coeffs/stis-ccd", coeff_filename_length);
    } else if (strncmp(detector_name, "NUV-MAMA", 8) == 0) {
      strncpy(coeff_filename, "drizzle$coeffs/stis-nuv-mama", coeff_filename_length);
    } else if (strncmp(detector_name, "FUV-MAMA", 8) == 0) {
      strncpy(coeff_filename, "drizzle$coeffs/stis-fuv-mama", coeff_filename_length);
    } else {
      driz_error_format_message(error,
                                "Invalid STIS detector name in header '%s'",
                                detector_name);
      return 1;
    }
  } else {
    driz_error_format_message(error,
                              "Invalid instrument in header.  Currently "
                              "handled: WFPC2, NICMOS, and STIS, got '%s'",
                              instrument);
    return 1;
  }

  return 0;
}

#endif
#ifndef BUFSIZ
#define BUFSIZ (1 << 16)
#endif

/*
This is a work-alike of the getline function in GNU libc borrowed from
X11.  This allows our code to be portable to non-GNU libc platforms.

It has the following original license, which appears to be compatible
with our license.

Copyright (c) 1996  X Consortium

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and sell copies of the Software, and to permit
persons to whom the Software is furnished to do so, subject to the
following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE X CONSORTIUM BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

Except as contained in this notice, the name of the X Consortium shall
not be used in advertising or otherwise to promote the sale, use or
other dealings in this Software without prior written authorization
from the X Consortium.
*/
static int
__getline(char **pbuf, size_t *plen, FILE *f) {
  int c;
  size_t i;

  i = 0;
  while(1) {
    if (i+2 > *plen) {
      if (*plen) {
        *plen *= 2;
      } else {
        *plen = BUFSIZ;
      }
      if (*pbuf) {
        *pbuf = (char *) realloc (*pbuf, *plen + 1);
      } else {
        *pbuf = (char *) malloc (*plen + 1);
      }
      if (! *pbuf) {
        return 0;
      }
    }
    c = getc (f);
    if (c == EOF) {
      break;
    }
    (*pbuf)[i++] = c;
    if (c == '\n') {
      i--;
      break;
    }
  }
  (*pbuf)[i] = '\0';
  return i;
}

/**
A wrapper around getline that automatically ignores empty and
commented lines.

See "man 3 getline" for usage info.
*/
static inline_macro int
get_line(FILE *fd,
         char** line,
         size_t* line_size) {
  ssize_t num_read;
  char first_char;

  assert(fd);
  assert(line);
  assert(line_size);

  while (1) {
    num_read = __getline(line, line_size, fd);
    first_char = **line;
    if (num_read == -1 ||
        !(first_char == 0 ||
          first_char == '\n' ||
          first_char == '\r' ||
          first_char == '#')) {
      return num_read;
    }
  }

  /* Shouldn't ever get here */
  assert(FALSE);
  return -1;
}

/**
Reads \a n doubles separated by whitespace from \a s and writes them into
\a array.

@param[in,out] s A string containing \a n numbers at its beginning.
Upon return, this will point to the character immediately following
the \a n read numbers.

@param[in] n The number of numbers to read.

@param[out] array An array to store the numbers.  Must be at least
length \a n.

@param[out] error

@return Non-zero if an error occurred.
*/
int
read_numbers(char** s, const integer_t n,
             /* Output parameters */
             double* array /* [n] */,
             struct driz_error_t* error) {
  integer_t i;
  char* end;

  assert(s);
  assert(*s);
  assert(**s);
  assert(n > 0);
  assert(array);
  assert(error);

  for (i = 0; i < n; ++i) {
    array[i] = strtod(*s, &end);
    if (end == *s) {
      driz_error_set_message(error, "Error reading numbers");
      return 1;
    }
    *s = end;
  }

  return 0;
}

/**
Reads the entire contents of the file (following the current seek
pointer) into a buffer.  The buffer will be allocated on the heap to
the correct size, and it is the caller's responsibility to free it.

@param[in] fd A C \a FILE pointer that has already been opened and is
\a seeked to the beginning of the data to be read.

@param[out] buffer A pointer to the beginning of the buffer containing
the file data.

@param[out] error

@return Non-zero if an error has occurred.
*/
int
read_all_file(FILE* fd,
              /* Output parameters */
              char** buffer,
              struct driz_error_t* error) {
  long file_pos;
  long file_size;
  long bytes_read;

  assert(fd);
  assert(buffer);
  assert(*buffer == NULL);

  file_pos = ftell(fd);
  if (fseek(fd, 0, SEEK_END)) {
    driz_error_set_message(error, "Unable to determine size of file");
    return 1;
  }
  file_size = ftell(fd) - file_pos;

  if (fseek(fd, file_pos, SEEK_SET)) {
    driz_error_set_message(error, "Unable to determine size of file");
    return 1;
  }

  *buffer = malloc(file_size);
  if (*buffer == NULL) {
    driz_error_set_message(error, "Out of memory");
    return 1;
  }

  bytes_read = fread(*buffer, 1, file_size, fd);
  if (bytes_read != file_size) {
    driz_error_set_message(error, "Could not read all contents of file");
    free(*buffer); *buffer = NULL;
    return 1;
  }

  return 0;
}

/**
Read geometric distortion coefficients in free format from a text file
which has already been opened and return them.

@param[in] fd A C \a FILE pointer, already opened.

@param[out] lambda The wavelength in nm.

@param[out] coeff_type The type of coefficients being stored (TODO: Add description).

@param[out] num_coeffs The number of coefficients being stored.

@param[out] x_coeffs An array of x coefficients.  This array should be
            allocated to size MAX_COEFFS, but only the first
            num_coeffs will be valid after returning from this
            function.

@param[out] y_coeffs An array of y coefficients.  This array should be
            allocated to size MAX_COEFFS, but only the first
            num_coeffs will be valid after returning from this
            function.

@param[out] error

@return Non-zero if an error occurred.

was: GETCO
*/
int
get_coefficients_from_file(FILE* fd,
                           /* Output arguments */
                           double* lambda,
                           integer_t* coeff_type,
                           integer_t* num_coeffs,
                           double* x_coeffs /* [num_coeffs] */,
                           double* y_coeffs /* [num_coeffs] */,
                           struct driz_error_t* error) {
  char* line;
  char* end;
  size_t line_size = BUFSIZ;
  ssize_t num_read;
  double n;
  double a[5];
  double value;
  integer_t j;
  char* buffer = NULL;
  char* cursor = NULL;
  double xdref = 0.0, ydref = 0.0;
  bool_t newref = 0;
  long order;
  int stat = 1;

  assert(fd);
  assert(lambda);
  assert(coeff_type);
  assert(num_coeffs);
  assert(x_coeffs);
  assert(y_coeffs);
  assert(error);

  line = malloc(line_size);
  if (line == NULL) {
    driz_error_set_message(error, "Out of memory");
    return 1;
  }

  while (1) {
    num_read = get_line(fd, &line, &line_size);
    if (num_read == -1) break;

    /* First "Trauger" style coefficients -- 60 of them */
    if (num_read > 7 && memcmp(line, "trauger", 7) == 0) {
      /* Now we need to calculate the coefficients which are a
         function of wavelength */
      n = mgf2(*lambda);

      j = 0;
      while (j < 20) {
        num_read = get_line(fd, &line, &line_size);
        if (num_read == -1) break;

        cursor = line;
        if (read_numbers(&cursor, 3, a, error)) break;

        value = a[0] + a[1] * (n - 1.5) + pow(a[2] * (n - 1.5), 2);

        if (j < 10) {
          x_coeffs[j] = value;
        } else {
          y_coeffs[j-10] = value;
        }
      }

      *coeff_type = 3;
      *num_coeffs = 10;
      stat = 0;
      break;
    } else if (num_read >= 5 &&
               (memcmp(line, "cubic", 5) == 0 ||
                memcmp(line, "poly3", 5) == 0)) {
      if (read_all_file(fd, &buffer, error)) break;

      cursor = buffer;
      if (read_numbers(&cursor, 10, x_coeffs, error)) break;
      if (read_numbers(&cursor, 10, y_coeffs, error)) break;

      *coeff_type = 3;
      *num_coeffs = 10;
      stat = 0;

      break;
    } else if ((num_read >= 5 && memcmp(line, "poly4", 5) == 0) ||
               (num_read >= 7 && memcmp(line, "quartic", 7) == 0)) {
      if (read_all_file(fd, &buffer, error)) break;

      cursor = buffer;
      if (read_numbers(&cursor, 15, x_coeffs, error)) break;
      if (read_numbers(&cursor, 15, y_coeffs, error)) break;

      *coeff_type = 4;
      *num_coeffs = 15;
      stat = 0;

      break;
    } else if ((num_read >= 5 && memcmp(line, "poly5", 5) == 0) ||
               (num_read >= 7 && memcmp(line, "quintic", 7) == 0)) {
      if (read_all_file(fd, &buffer, error)) break;

      cursor = buffer;
      if (read_numbers(&cursor, 21, x_coeffs, error)) break;
      if (read_numbers(&cursor, 21, y_coeffs, error)) break;

      *coeff_type = 5;
      *num_coeffs = 21;
      stat = 0;

      break;
    } else if ((num_read >= 5 && memcmp(line, "poly", 4) == 0)) {
      order = strtol(line + 4, &end, 10);
      if (order <= 0 || order > MAX_COEFFS) {
        driz_error_set_message(error, "Invalid number of coefficients");
        break;
      }

      *coeff_type = (integer_t)order;
      *num_coeffs = (*coeff_type + 1) * (*coeff_type + 2) / 2;

      if (*num_coeffs > MAX_COEFFS) {
        driz_error_set_message(error, "Too many of coefficients");
        break;
      }

      if (read_all_file(fd, &buffer, error)) break;

      cursor = buffer;
      if (read_numbers(&cursor, *num_coeffs, x_coeffs, error)) {
        break;
      }
      if (read_numbers(&cursor, *num_coeffs, y_coeffs, error)) {
        break;
      }

      stat = 0;
      break;
    } else if ((num_read >= 6 && memcmp(line, "radial", 6) == 0) ||
               (num_read >= 4 && memcmp(line, "rad3", 4) == 0)) {
      if (read_all_file(fd, &buffer, error)) break;

      *num_coeffs = 3;
      *coeff_type = -3;

      cursor = buffer;
      if (read_numbers(&cursor, *num_coeffs, x_coeffs, error)) {
        break;
      }

      stat = 0;
      break;
    } else if ((num_read >= 6 && memcmp(line, "refpix", 6) == 0)) {
      cursor = line + 6;
      if (read_numbers(&cursor, 2, a, error)) {
        break;
      }
      xdref = a[0];
      ydref = a[1];
      newref = 1;
    } else {
      driz_error_set_message(error, "Unknown coefficient set type");
      break;
    }
  }

  /* If we have a modified reference pixel we encode this in the
     standard coefficients by adding an offset and incrementing the
     number of coefficients by one */
  if (newref) {
    *coeff_type += COEFF_OFFSET;
    x_coeffs[*num_coeffs] = xdref;
    y_coeffs[*num_coeffs] = ydref;
    (*num_coeffs)++;
  }

  free(line); line = NULL;
  free(buffer); buffer = NULL;
  return stat;
}
