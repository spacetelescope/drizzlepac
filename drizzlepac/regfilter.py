"""
`regfilter` provides functions for filtering out regions outside the image space.

:Authors: Mihai Cara

:License: :doc:`LICENSE`

"""
# This is specifically NOT intended to match the package-wide version information.
__version__ = '0.1'
__version_date__ = '17-Nov-2013'
__author__ = 'Mihai Cara'


def fast_filter_outer_regions(reglist, width, height, origin=1):
    # fast_filter_outer_regions filters regions that are outside a rectangle
    # ('image's rectangle') of width 'width' and height 'height' that is has the
    # bottom-left corner at (origin,origin). This function is based on checking
    # for the intersection of the image's rectangle with the bounding box of the
    # regions and therefore it is approximate (some regions that in reality are
    # not within the image's rectangle will not be filtered out if their
    # bounding box still intersects the image rectangle even though the shape
    # itself does not intersect the image's rectangle.
    for k in range(len(reglist)-1,-1,-1):
        reg = reglist[k]
        regname = reg.name.lower()
        if regname[:3] == 'cir' or regname == 'annulus':
            blc, trc = _get_bb_circle(reg)
        elif regname[-3:] == 'box':
            blc, trc = _get_bb_box(reg)
        elif regname == 'ellipse':
            blc, trc = _get_bb_ellipse(reg)
        elif regname == 'polygon':
            blc, trc = _get_bb_polygon(reg)
        elif regname == 'point':
            x = reg.coord_list[0]
            y = reg.coord_list[1]
            if not _is_point_inside(width, height, x, y, origin=origin):
                del reglist[k]
                continue
        elif regname[:4] == 'rect':
            blc, trc = _get_bb_rect(reg)
        elif regname == 'panda':
            blc, trc = _get_bb_circle(reg, True)
        elif regname == 'epanda':
            blc, trc = _get_bb_ellipse(reg, True)
        elif regname == 'bpanda':
            blc, trc = _get_bb_box(reg, True)
        else:
            continue

        if not _is_rect_inside(width, height, blc, trc, origin=origin):
            del reglist[k]
            continue


def _is_rect_inside(w1, h1, blc2, trc2, origin=1):
    pad = 0.5
    o = origin-pad
    return ((o < trc2[0]) and (o + w1 > blc2[0]) \
        and (o < trc2[1]) and (o + h1 > blc2[1]))


def _is_point_inside(w, h, x, y, origin=1):
    pad = 0.5
    o = origin-pad
    return (o < x and (o + w > x) and (o < y) and (o + h > y))


def _get_bb_rect(shape):
    # CIAO rectangles
    return (shape.coord_list[0],shape.coord_list[2]), \
        (shape.coord_list[1],shape.coord_list[3])


def _get_bb_box(shape, bpanda=False):
    from math import sin, cos, radians

    # check if angle is provided:
    rem = len(shape.coord_list) % 2
    # check if bpanda:
    pnd = 1 if bpanda else 0

    xc = shape.coord_list[0]
    yc = shape.coord_list[1]
    w  = shape.coord_list[-2-rem-pnd] / 2.0
    h  = shape.coord_list[-1-rem-pnd] / 2.0
    th = radians(shape.coord_list[-1]) if rem > 0 else 0.0
    cs = cos(th)
    sn = sin(th)
    xm = max(abs(w*cs-h*sn),abs(w*cs+h*sn))
    ym = max(abs(w*sn+h*cs),abs(w*sn-h*cs))
    return (xc-xm,yc-ym),(xc+xm,yc+ym)


def _get_bb_circle(shape, panda=False):
    # check if panda:
    pnd = 1 if panda else 0

    xc = shape.coord_list[0]
    yc = shape.coord_list[1]
    r  = shape.coord_list[-1-pnd]
    return (xc-r,yc-r),(xc+r,yc+r)


def _get_bb_ellipse(shape, epanda=False):
    from math import sin, cos, radians, sqrt

    # check if angle is provided:
    rem = len(shape.coord_list) % 2
    # check if epanda:
    pnd = 1 if epanda else 0

    xc = shape.coord_list[0]
    yc = shape.coord_list[1]
    a  = shape.coord_list[-2-rem-pnd]
    b  = shape.coord_list[-1-rem-pnd]
    th = radians(shape.coord_list[-1]) if rem > 0 else 0.0
    cs = cos(th)
    sn = sin(th)
    xm = sqrt( (a*cs)**2 + (b*sn)**2 )
    ym = sqrt( (a*sn)**2 + (b*cs)**2 )
    return (xc-xm,yc-ym),(xc+xm,yc+ym)


def _get_bb_point(shape):
    xc = shape.coord_list[0]
    yc = shape.coord_list[1]
    return (xc-0.5,yc-0.5),(xc+0.5,yc+0.5)


def _get_bb_polygon(shape):
    xs = shape.coord_list[0::2]
    ys = shape.coord_list[1::2]
    minx = min(xs)
    maxx = max(xs)
    miny = min(ys)
    maxy = max(ys)
    return (minx,miny),(maxx,maxy)
