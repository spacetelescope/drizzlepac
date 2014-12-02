import numpy as np

def calc_pixmap(first_wcs, second_wcs):
    """
    Calculate a mapping between the pixels of two images
    """

    first_naxis1 = first_wcs._naxis1
    first_naxis2 = first_wcs._naxis2

    # We add one to the pixel co-ordinates before the transformation and subtract
    # it afterwards because wcs co-ordinates are one based, while pixel co-ordinates
    # are zero based, The result is the final values in pixmap give a tranformation
    # between the pixel co-ordinates in the first image to pixel co-ordinates in the
    # co-ordinate system of the second.

    one = np.ones(2, dtype='float64')
    idxmap = np.indices((first_naxis1, first_naxis2), dtype='float64')
    idxmap = idxmap.transpose() + one
    
    idxmap = idxmap.reshape(first_naxis2 * first_naxis1, 2)
        
    worldmap = first_wcs.all_pix2world(idxmap, 1)

    if second_wcs.sip is None:
        pixmap = second_wcs.wcs_world2pix(worldmap, 1)
    else:
        pixmap = second_wcs.all_world2pix(worldmap, 1)

    pixmap = pixmap.reshape(first_naxis2, first_naxis1, 2)
    pixmap = pixmap - one

    # Check to see if two images do not overlap,
    # return None if that is the case as a flag
    # to not perform the drizzle

    xmin = pixmap[:,:,0].min()
    xmax = pixmap[:,:,0].max()
    ymin = pixmap[:,:,1].min()
    ymax = pixmap[:,:,1].max() 

    if (xmax < -1 or xmin > second_wcs._naxis1 + 1 or
        ymax < -1 or ymin > second_wcs._naxis2 + 1):
        return None
    
    return pixmap
