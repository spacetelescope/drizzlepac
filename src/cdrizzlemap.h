#ifndef CDRIZZLEMAP_H
#define CDRIZZLEMAP_H

#include "driz_portability.h"
#include "cdrizzleutil.h"

void
map_point(PyArrayObject * pixmap,
          const double xyin[2],
          double xyout[2]
         );



#endif /* CDRIZZLEMAP_H */
