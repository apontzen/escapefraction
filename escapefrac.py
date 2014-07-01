import pynbody
import numpy as np
import healpy as hp
import pylab as p
import math
import contextlib
import sys
import os

@pynbody.derived_array
def rhoHI(f) :
    return f['rho']*f['HI']



class quiet(object) :
    def __enter__(self) :
        self.stdout = sys.stdout
        self.stderr = sys.stderr
        sys.stdout = open(os.devnull, 'w')
        sys.stderr = sys.stdout

    def __exit__(self, *args) :
        sys.stdout = self.stdout
        sys.stderr = self.stderr

    
    
def escape_fraction_from_position(f, pos, to_distance=100.0,nside=64,plot=True) :
    with translate(f,-pos) :
        im = pynbody.sph.render_spherical_image(f.gas,nside=nside,distance=to_distance,
                                            qty='rhoHI',kernel=pynbody.sph.Kernel2D(),
                                            out_units="m_p cm^-2",denoise=False,threaded=False)
    
    if plot :
        hp.mollview(np.log10(im))

    esc_frac_per_pixel = np.exp(-6.3e-18*im)
    return esc_frac_per_pixel.mean()


def dirty_escape_fraction_from_position(f, pos, to_distance=100.0,nside=64,plot=True,
                                        skysize=0.1) :
    """A quick and dirty escape fraction choosing just a random
    portion of the sky to include.  This will give you very wrong
    results for an individual star, but averaged over many stars it
    comes out to the correct answer."""

    cos_theta_cut = 1-skysize

    with translate(f,-pos) :
        randvec = np.random.normal(0,1,3)
        randvec/=np.linalg.norm(randvec)
        cos_theta = np.dot(randvec,np.asarray(f.gas['pos']).T)/f.gas['r']
        use_particles = cos_theta>cos_theta_cut
        f_sub = f.gas[use_particles]
        im = pynbody.sph.render_spherical_image(f_sub,nside=nside,distance=to_distance,
                                                qty='rhoHI',kernel=pynbody.sph.Kernel2D(),
                                                out_units="m_p cm^-2",denoise=False,threaded=False)

        # find the pixels well within the disc
        use_pixels = hp.query_disc(nside, randvec, math.acos(cos_theta_cut)*0.95)

        if plot :
            hp.mollview(np.log10(im))
            im2 = np.zeros_like(im)
            im2[use_pixels]=im[use_pixels]
            hp.mollview(np.log10(im2))


        im = im[use_pixels]
        esc_frac_per_pixel = np.exp(-6.3e-18*im)
    return esc_frac_per_pixel.mean()
    
def escape_fraction_from_youngstars(h,max_age='50 Myr',min_age='10 Myr',
                                    max_num_stars=100) :
    """Returns the mean escape fraction from a selection of young stars
    chosen according to the criteria min_age<age<max_age, then randomly subsampled
    so there are a maximum of max_num_stars."""
    
    st = h.stars[pynbody.filt.BandPass('age',min_age,max_age)]
    print "There are",len(st),"stars"
    
    myorder = np.random.permutation(np.arange(len(st)))[:max_num_stars]
    print "Truncated to",len(myorder),"stars"

    efrac = np.empty(len(myorder))
    for i,index in enumerate(myorder) :
        print "Calculating for star %d (%d in snapshot)"%(i,index)
        efrac[i]=escape_fraction_from_position(h,st['pos'][index],to_distance=1000.0,plot=max_num_stars<3)
        print " --> estimate is ",efrac[i]

    return efrac.mean()


def dirty_escape_fraction_from_youngstars(h,max_age='50 Myr',min_age='10 Myr',
                                          max_num_stars=100,skysize=0.2) :
    
    
    st = h.stars[pynbody.filt.BandPass('age',min_age,max_age)]
    print "There are",len(st),"stars"
    
    myorder = np.random.permutation(np.arange(len(st)))[:max_num_stars]
    print "Truncated to",len(myorder),"stars"

    efrac = np.empty(len(myorder))
    for i,index in enumerate(myorder) :
        print "Calculating for star %d (number %d in snapshot)"%(i,index)
        with quiet() :
            efrac[i]=dirty_escape_fraction_from_position(h,st['pos'][index],to_distance=1000.0,plot=max_num_stars<3,skysize=skysize)
        print " --> estimate is ",efrac[i]

    return efrac.mean()
