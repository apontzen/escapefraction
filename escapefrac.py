import pynbody
import numpy as np
import healpy as hp
import pylab as p
import math
import contextlib
import sys
import os
import pyximport
pyximport.install(setup_args={'include_dirs': np.get_include()},reload_support=True)
import grid_angle

from pynbody.transformation import translate, rotate_x, rotate_y, rotate_z

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

def cuboid_escape_fraction(f, pos, to_distance=100.0, resolution=500) :
    pixel_solid_angles = grid_angle.grid_pixel_angles(to_distance,to_distance*2,resolution)
    
    def mkimage(X,star_point) :
        im =  pynbody.plot.sph.image(X,width=to_distance*2,qty='rhoHI', units='m_p cm^-2',noplot=False,z_camera=star_point,resolution=resolution)
        im = np.exp(-6.3e-18*im).view(np.ndarray)
        return im

    def get_sides(f) :
        state = f.gas
        try :
            state = translate(state,[0.0,0.0,to_distance])
            side_a = mkimage(f.gas[pynbody.filt.BandPass('z',0,to_distance)],to_distance)
            
            state = translate(state,[0.0,0.0,-2*to_distance])
            side_b = mkimage(f.gas[pynbody.filt.BandPass('z',-to_distance,0)],-to_distance)

        finally :
            state.revert()
            
        return side_a,side_b
    
    with translate(f,-pos) :
        state = f.gas
        try:
            side_a,side_b = get_sides(f)
            state = rotate_x(state,90)
            side_c,side_d = get_sides(f)
            state = rotate_y(rotate_x(state,-90,defer=True),90)

            with f.gas.rotate_y(90) :
                side_e,side_f = get_sides(f)


    return (side_a*pixel_solid_angles+side_b*pixel_solid_angles+side_c*pixel_solid_angles+
            side_d*pixel_solid_angles+side_e*pixel_solid_angles+side_f*pixel_solid_angles).sum()/(6*pixel_solid_angles.sum())
        

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
