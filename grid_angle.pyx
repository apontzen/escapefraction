cimport numpy as np
cimport cython
cimport libc.math as cmath
from libc.math cimport atan, pow, sqrt
from libc.stdlib cimport malloc, free
import numpy as np

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def grid_pixel_angles(double z_cam, double width, int npix) :
    """Work out the solid angle for each pixel of an image with the camera at z_cam
    and the total image width in the z=0 plane being width"""
    
    cdef np.ndarray[double, ndim=2] output = np.empty((npix,npix),dtype='float64')
    cdef double delta2 = (width/npix)**2
    cdef double x2,y2
    cdef int ix, iy
    cdef double z_cam2 = z_cam*z_cam
    cdef double total_dist2
    
    with nogil: 
        for ix in xrange(npix/2+1) :
            for iy in xrange(npix/2+1) :
                R2 = delta2*(npix/2-ix)*(npix/2-ix) + delta2*(npix/2-iy)*(npix/2-iy)
                total_dist2 =  R2+z_cam2
                output[ix,iy] = (delta2*(z_cam/sqrt(total_dist2)))/total_dist2
                output[npix-1-ix,iy] = output[ix,iy]
                output[ix,npix-1-iy] = output[ix,iy]
                output[npix-1-ix,npix-1-iy] = output[ix,iy]
    return output
                

