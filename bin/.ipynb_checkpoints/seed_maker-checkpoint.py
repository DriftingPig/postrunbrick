import numpy as np
import sys
sys.path.append("./")
from grid import GlassDistribution
from RegularGrid import RectGrid, HexGrid
from astropy.table import vstack,Table
import astropy.io.fits as fits
from astrometry.util.fits import fits_table
from legacypipe.survey import wcs_for_brick,LegacySurveyData
import multiprocessing as mp
class seed_maker(object):
    def __init__(self, ndraws=None, ramin=None, ramax=None, decmin=None, decmax=None, outdir=None, rotation=None, seed=None, surveybricks = None, bricklist = None):
        self.ndraws = ndraws
        self.ramin = ramin
        self.ramax = ramax
        self.decmin = decmin
        self.decmax = decmax
        self.outdir = outdir
        self.rotation = 0.
        self.seed = seed
        self.surveybricks = surveybricks
        self.bricklist = bricklist
    def __call__(self):
        print("grid_sampler")
        self.grid_sampler('glass')
        print("grid_transform")
        self.grid_transform()
        print("make_input")
        self.make_input()
        print("split_to_bricks")
        self.split_to_bricks()
        print("maskbits")
        self.maskbits()
        print("mask")
        self.mask()
    def grid_sampler(self,grid_type):
        """
        sampling a set of grids ranging within [[0,1],[0,1]]
        """
        if grid_type == 'glass':
            self.distrib = GlassDistribution(npoints = self.ndraws)
            self.distrib()
        elif grid_type == 'rect':
            self.distrib = RectGrid(shape=self.ndraws, rotation=self.rotation)
        elif grid_type == 'hex':
            self.distrib = HexGrid(shape=self.ndraws, rotation=self.rotation)
        else:
            raise ValueError("unknown grid type %s"%grid_type)
        self.x = self.distrib.positions[:,0]
        self.y = self.distrib.positions[:,1]
        
    def grid_transform(self):
        """
        transform grid to ra,dec coordinate
        """
        cmin = np.sin(self.decmin*np.pi/180)
        cmax = np.sin(self.decmax*np.pi/180)
        RA   = self.ramin + self.x*(self.ramax - self.ramin)
        DEC  = 90-np.arccos(cmin + self.y*(cmax - cmin))*180./np.pi
        self.ra = RA
        self.dec = DEC
        
    def make_input(self):
        """
        add seed info to current catalog
        """
        random_state=np.random.RandomState()
        seed = self.seed
        T = fits_table()
        T.set('id',np.arange(self.ndraws))
        T.set('ra',self.ra)
        T.set('dec',self.dec)
        ids = random_state.randint(low=0,high=len(seed),size=self.ndraws)
        T.set('g',seed['g'][ids])
        T.set('r',seed['r'][ids])
        T.set('z',seed['z'][ids])
        T.set('n',seed['n'][ids])
        T.set('rhalf',seed['rhalf'][ids])
        T.set('id_sample', seed['id_sample'][ids])
        T.set('w1',seed['w1'][ids])
        T.set('w2',seed['w2'][ids])
        T.set('redshift', np.ones_like(ids))
        T.set('e1',seed['e1'][ids])
        T.set('e2',seed['e2'][ids])
        self.T = T
        fn = outdir+'randoms.fits'
        T.writeto(fn,overwrite=True)
    def split_to_bricks(self):
        """
        split the current catalog to a per brick level
        """
        p = mp.Pool(30)
        p.map(self._split_to_bricks, self.bricklist)
    def _split_to_bricks(self, brickname):
            surveybrick_i = self.surveybricks[surveybricks['BRICKNAME']==brickname]
            ra1 = surveybrick_i['RA1'][0]
            ra2 = surveybrick_i['RA2'][0]
            dec1 = surveybrick_i['DEC1'][0]
            dec2 = surveybrick_i['DEC2'][0]
            T_sel = (self.T.ra>=ra1)&(self.T.ra<=ra2)&(self.T.dec>=dec1)&(self.T.dec<=dec2)
            T_new = self.T[T_sel]
            fn = self.outdir+'/divided_randoms_2/bricks_%s.fits'%brickname
            T_new.writeto(fn)
        
        
    def maskbits(self):
        """
        add maskbits to catalog, this needs to be done in docker
        """
        p = mp.Pool(30)
        p.map(self._maskbits, self.bricklist)
    def _maskbits(self, brickname):
            fn = self.outdir+'divided_randoms_2/bricks_%s.fits'%brickname
            T = fits_table(fn)
            survey = LegacySurveyData(survey_dir="/global/cfs/cdirs/cosmo/work/legacysurvey/dr9/")
            brick = survey.get_brick_by_name(brickname)
            brickwcs = wcs_for_brick(brick)
            W, H, pixscale = brickwcs.get_width(), brickwcs.get_height(), brickwcs.pixel_scale()
            targetwcs = wcs_for_brick(brick, W=W, H=H, pixscale=pixscale)
            flag, target_x, target_y = targetwcs.radec2pixelxy(T.ra, T.dec)
            T.set('bx',target_x)
            T.set('by', target_y)
            
            maskbits_dr9 = fits.getdata("/global/cfs/cdirs/cosmo/work/legacysurvey/dr9/south/coadd/%s/%s/legacysurvey-%s-maskbits.fits.fz"%(brickname[:3],brickname,brickname))
            nexp_g = fits.getdata("/global/cfs/cdirs/cosmo/work/legacysurvey/dr9/south/coadd/%s/%s/legacysurvey-%s-nexp-g.fits.fz"%(brickname[:3],brickname,brickname))
            nexp_r = fits.getdata("/global/cfs/cdirs/cosmo/work/legacysurvey/dr9/south/coadd/%s/%s/legacysurvey-%s-nexp-r.fits.fz"%(brickname[:3],brickname,brickname))
            nexp_z = fits.getdata("/global/cfs/cdirs/cosmo/work/legacysurvey/dr9/south/coadd/%s/%s/legacysurvey-%s-nexp-z.fits.fz"%(brickname[:3],brickname,brickname))
            bx = (T.bx+0.5).astype(int)
            by = (T.by+0.5).astype(int)
            mask_flag = maskbits_dr9[(by),(bx)]
            mask_flag_g = nexp_g[(by),(bx)]
            mask_flag_r = nexp_r[(by),(bx)]
            mask_flag_z = nexp_z[(by),(bx)]
            T.set('maskbits',mask_flag)
            T.set('nobs_g',mask_flag_g)
            T.set('nobs_r',mask_flag_r)
            T.set('nobs_z',mask_flag_z)
            fn = self.outdir+'/divided_randoms_2/bricks_%s.fits'%brickname
            T.writeto(fn,overwrite=True)
        
    def mask(self):
        """
        mask out un-used bits, to save time in future process
        """
        p = mp.Pool(30)
        p.map(self._mask, self.bricklist)
        
    def _mask(self,brickname):
            fn = self.outdir+'/divided_randoms_2/bricks_%s.fits'%brickname
            randoms = fits_table(fn)
            maskbits = (randoms.maskbits&2**1)|(randoms.maskbits&2**8)|(randoms.maskbits&2**9)|(randoms.maskbits&2**12)|(randoms.maskbits&2**13)
            sel_obs = (randoms.nobs_g>0)&(randoms.nobs_r>0)&(randoms.nobs_z>0)
            randoms = randoms[(maskbits==0)&(sel_obs)]
            randoms.writeto(fn, overwrite = True)
    
if __name__ == "__main__":
    ndraws = 50*2000
    ramin = 146.5
    ramax = 157.5
    decmin = -5.5
    decmax = 5.5
    outdir = "/global/cscratch1/sd/huikong/Obiwan/dr9_LRG/obiwan_out/deep1/"
    seed = fits.getdata(outdir+'/seed.fits')
    surveybricks = fits.getdata("/global/cfs/cdirs/cosmo/work/legacysurvey/dr9/survey-bricks.fits.gz")
    bricklist = np.loadtxt("/global/cscratch1/sd/huikong/Obiwan/dr9_LRG/obiwan_out/deep1/bricklist.txt",dtype = np.str)
    seeds = seed_maker(ndraws=ndraws, ramin=ramin, ramax=ramax, decmin=decmin, decmax=decmax, outdir=outdir, rotation=0., seed=seed, surveybricks = surveybricks, bricklist = bricklist)
    seeds()
    
#shifter --module=mpich-cle6 --image=legacysurvey/legacypipe:DR9.7.1 /bin/bash    