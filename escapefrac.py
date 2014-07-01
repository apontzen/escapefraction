import pynbody
import numpy as np
import healpy as hp
import pylab as p
import math

@pynbody.derived_array
def rhoHI(f) :
    return f['rho']*f['HI']

@pynbody.derived_array
def stellar_age(f) :
    return f['tform'].max()-f['tform']
    
def escape_fraction_from_position(f, pos, to_distance=100.0,nside=64,plot=True) :
    f['pos']-=pos
    im = pynbody.sph.render_spherical_image(f.gas,nside=nside,distance=to_distance,
                                            qty='rhoHI',kernel=pynbody.sph.Kernel2D(),
                                            out_units="m_p cm^-2",denoise=False,threaded=False)
    
    if plot :
        hp.mollview(np.log10(im))

    esc_frac_per_pixel = np.exp(-6.3e-18*im)
    f['pos']+=pos
    return esc_frac_per_pixel.mean()


def dirty_escape_fraction_from_position(f, pos, to_distance=100.0,nside=64,plot=True,
                                        skysize=0.1) :
    """A quick and dirty escape fraction choosing just a random
    portion of the sky to include.  This will give you very wrong
    results for an individual star, but averaged over many stars it
    comes out to the correct answer."""

    cos_theta_cut = 1-skysize
    print "cos_theta_cut=",cos_theta_cut
    
    f['pos']-=pos
    randvec = np.random.normal(0,1,3)
    randvec/=np.linalg.norm(randvec)
    print "randvec=",randvec
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
    f['pos']+=pos
    return esc_frac_per_pixel.mean()
    
def escape_fraction_from_youngstars(h,max_age='50 Myr',min_age='10 Myr',
                                    max_num_stars=100) :
    """Returns the mean escape fraction from a selection of young stars
    chosen according to the criteria min_age<age<max_age, then randomly subsampled
    so there are a maximum of max_num_stars."""
    
    st = h.stars[pynbody.filt.BandPass('age',min_age,max_age)]
    print "There are",len(st),"stars"
    if len(st)>max_num_stars :
        st = st[np.random.permutation(np.arange(len(st)))[:max_num_stars]]
        print "Truncated to",len(st),"stars"

    efrac = np.empty(len(st))
    for i in xrange(len(st)) :
        print "Calculating for star",i
        efrac[i]=escape_fraction_from_position(h,st['pos'][i],to_distance=1000.0,plot=False)

    return efrac.mean()


def dirty_escape_fraction_from_youngstars(h,max_age='50 Myr',min_age='10 Myr',
                                          max_num_stars=100,skysize=0.2) :
    
    st = h.stars[pynbody.filt.BandPass('age',min_age,max_age)]
    print "There are",len(st),"stars"
    if len(st)>max_num_stars :
        st = st[np.random.permutation(np.arange(len(st)))[:max_num_stars]]
        print "Truncated to",len(st),"stars"

    efrac = np.empty(len(st))
    for i in xrange(len(st)) :
        print "Calculating for star",i
        efrac[i]=dirty_escape_fraction_from_position(h,st['pos'][i],to_distance=1000.0,plot=False,skysize=skysize)

    return efrac.mean()

