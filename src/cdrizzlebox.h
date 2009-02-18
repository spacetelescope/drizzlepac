#ifndef CDRIZZLEBOX_H
#define CDRIZZLEBOX_H

#include "cdrizzleutil.h"

/**
dobox

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
dobox(struct driz_param_t* p, const integer_t ystart, integer_t* nmiss,
      integer_t* nskip, struct driz_error_t* error);

#endif /* CDRIZZLEBOX_H */
