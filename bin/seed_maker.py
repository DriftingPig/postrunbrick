import numpy as np
import sys
sys.path.append("./")
from grid import GlassDistribution
from RegularGrid import RectGrid, HexGrid
from astropy.table import vstack,Table
import astropy.io.fits as fits
try:
    from astrometry.util.fits import fits_table, merge_tables
    from legacypipe.survey import wcs_for_brick,LegacySurveyData
except:
    pass
import multiprocessing as mp
import os
import subprocess
class seed_maker(object):
    def __init__(self, ndraws=None, outdir=None, rotation=None, seed=None, surveybricks = None, bricklist = None, grid_type = None ,seed_num = None, rand_fn='divided_randoms_2'):
        self.bits = [1, 5, 6, 7, 11, 12, 13]# BRIGHT, ALLMASK_G,R,Z, MEDIUM, GALAXY, CLUSTER
        self.ndraws = ndraws
        self.outdir = outdir
        self.rotation = 0.
        self.seed = seed
        self.grid_type = grid_type
        self.surveybricks = surveybricks
        self.bricklist = bricklist
        self.seed_num = seed_num
        self.rand_fn = rand_fn
        fn = self.outdir+'/'+self.rand_fn
        if not os.path.isdir(fn):
            print(fn)
            subprocess.call(["mkdir", "-p", fn])
            
    def __call__(self):
        print("grid_sampler")
        self.grid_sampler(self.grid_type)
        
    def grid_sampler(self, grid_type):
        self.grid_type = grid_type
        p = mp.Pool(30)
        p.map(self._grid_sampler, self.bricklist)
    def _grid_sampler(self,brickname):
        """
        sampling a set of grids ranging within [W,H]
        """
        grid_type = self.grid_type
        survey = LegacySurveyData(survey_dir="/global/cfs/cdirs/cosmo/work/legacysurvey/dr9/")
        brick = survey.get_brick_by_name(brickname)
        brickwcs = wcs_for_brick(brick)
        W, H, pixscale = brickwcs.get_width(), brickwcs.get_height(), brickwcs.pixel_scale()
        side = 180 # 60*0.262 = 15.72 arcsec
        seed_num = self.seed_num + brick.brickcol    
        rng = np.random.RandomState(seed=seed_num)
        offset = rng.uniform(0.,1.,size=2)
        
        if grid_type == 'glass':
            distrib = GlassDistribution(npoints = self.ndraws)
            distrib()
            #assuming W=H
            distrib.positions *= W
            x = distrib.positions[:,0]
            y = distrib.positions[:,1]
        
        elif grid_type == 'rect':
            distrib = RectGrid(spacing = side, shape = (W,H), shift=brick.brickcol % 2)
            distrib.positions += offset + side/2
            distrib.positions = distrib._mask(distrib.positions)
            x = distrib.positions[:,0]
            y = distrib.positions[:,1]
            
        elif grid_type == 'hex':
            # for HexGrid spacing is defined as radius... here we rather want space along x and y, so divide by correct factor
            # we also alternate start of horizontal lines depending on brick column, to allow for better transition between bricks
            distrib = HexGrid(spacing=side/(2.*np.tan(np.pi/6)),shape=(W,H),shift=brick.brickcol % 2)
            distrib.positions += offset + side/2 # we add random pixel fraction offset, then side/2 because grid.positions start at 0
            distrib.positions = distrib._mask(distrib.positions)
            x = distrib.positions[:,0]
            y = distrib.positions[:,1]
    
        elif grid_type == "rand":
            #assuming W=H
            x, y = (W-1) * rng.uniform(size=(2, self.ndraws) )
        
        else:
            raise ValueError("unknown grid type %s"%grid_type)
        
        T = fits_table()
        N = len(x)
        T.set('id', np.arange(N))
        T.set('bx', x)
        T.set('by', y)
        
        #set ra,dec
        ra,dec = self.grid_transform(brickname, x/float(W), y/float(H))
        T.set('ra',ra)
        T.set('dec', dec)
        
        #set other properties from truth file
        truth = self.seed
        ids = rng.randint(low=0,high=len(truth),size=N)
        T.set('g',truth['g'][ids])
        T.set('r',truth['r'][ids])
        T.set('z',truth['z'][ids])
        T.set('n',truth['n'][ids])
        T.set('rhalf',truth['rhalf'][ids])
        T.set('id_sample', truth['id_sample'][ids])
        T.set('w1',truth['w1'][ids])
        T.set('w2',truth['w2'][ids])
        T.set('redshift', np.ones_like(ids))
        T.set('e1',truth['e1'][ids])
        T.set('e2',truth['e2'][ids])
        
        mask_flag, mask_flag_g, mask_flag_r, mask_flag_z = self.maskbits(brickname, x, y)
        T.set('maskbits',mask_flag)
        T.set('nobs_g',mask_flag_g)
        T.set('nobs_r',mask_flag_r)
        T.set('nobs_z',mask_flag_z)
        obs = (T.nobs_g>0)&(T.nobs_r>0)&(T.nobs_z>0)
        masked = self.mask(brickname, T.maskbits, obs)
        T = T[masked]
        fn = self.outdir+'/'+self.rand_fn+'/bricks_%s.fits'%brickname
        T.writeto(fn,overwrite=True)
        print(fn)
        
        
    def grid_transform(self, brickname, x, y):
        """
        transform grid to ra,dec coordinate
        """
        surveybrick_i = self.surveybricks[surveybricks['BRICKNAME']==brickname]
        ra1 = surveybrick_i['RA1'][0]
        ra2 = surveybrick_i['RA2'][0]
        dec1 = surveybrick_i['DEC1'][0]
        dec2 = surveybrick_i['DEC2'][0]
        
        cmin = np.sin(dec1*np.pi/180)
        cmax = np.sin(dec2*np.pi/180)
        RA   = ra1 + x*(ra2-ra1)
        DEC  = 90-np.arccos(cmin + y*(cmax - cmin))*180./np.pi
        return RA, DEC
        
        
    def maskbits(self, brickname, bx, by):    
            maskbits_dr9 = fits.getdata("/global/cfs/cdirs/cosmo/work/legacysurvey/dr9/south/coadd/%s/%s/legacysurvey-%s-maskbits.fits.fz"%(brickname[:3],brickname,brickname))
            nexp_g = fits.getdata("/global/cfs/cdirs/cosmo/work/legacysurvey/dr9/south/coadd/%s/%s/legacysurvey-%s-nexp-g.fits.fz"%(brickname[:3],brickname,brickname))
            nexp_r = fits.getdata("/global/cfs/cdirs/cosmo/work/legacysurvey/dr9/south/coadd/%s/%s/legacysurvey-%s-nexp-r.fits.fz"%(brickname[:3],brickname,brickname))
            nexp_z = fits.getdata("/global/cfs/cdirs/cosmo/work/legacysurvey/dr9/south/coadd/%s/%s/legacysurvey-%s-nexp-z.fits.fz"%(brickname[:3],brickname,brickname))
            bx = (bx+0.5).astype(int)
            by = (by+0.5).astype(int)
            mask_flag = maskbits_dr9[(by),(bx)]
            mask_flag_g = nexp_g[(by),(bx)]
            mask_flag_r = nexp_r[(by),(bx)]
            mask_flag_z = nexp_z[(by),(bx)]
            return mask_flag, mask_flag_g, mask_flag_r, mask_flag_z
        
    def mask(self,brickname, maskbits, obs):
            mb = np.ones_like(maskbits, dtype='?')
            for bit in self.bits:
                mb &= ((maskbits & 2**bit)==0)
            mb &= obs
            return mb
    def stack(self,ofn):
        TT = []
        for brickname in self.bricklist:
            TT.append(fits_table(self.outdir+'/'+self.rand_fn+'/bricks_%s.fits'%brickname))
        T = merge_tables(TT)
        T.writeto(ofn, overwrite = True)
        print("written %s"%ofn)
        
def corr_test(fn1,fn2):
    import treecorr
    dat1 = fits.getdata(fn1)
    dat2 = fits.getdata(fn2)
    cat1 = treecorr.Catalog(ra = dat1['ra'], dec = dat1['dec'],ra_units='degrees', dec_units='degrees')
    cat2 = treecorr.Catalog(ra = dat2['ra'], dec = dat2['dec'], ra_units='degrees', dec_units='degrees')
    thetamin=0.1
    thetamax=2
    nthetabins=20
    bin_slop=0
    bin_type="Log"
    nn = treecorr.NNCorrelation(min_sep=thetamin, max_sep=thetamax, nbins = nthetabins, bin_slop=bin_slop, bin_type=bin_type,sep_units='degrees')
    dr = treecorr.NNCorrelation(min_sep=thetamin, max_sep=thetamax, nbins = nthetabins, bin_slop=bin_slop, bin_type=bin_type,sep_units='degrees')
    rr = treecorr.NNCorrelation(min_sep=thetamin, max_sep=thetamax, nbins = nthetabins, bin_slop=bin_slop, bin_type=bin_type,sep_units='degrees')   
    rd = treecorr.NNCorrelation(min_sep=thetamin, max_sep=thetamax, nbins = nthetabins, bin_slop=bin_slop, bin_type=bin_type,sep_units='degrees')  
    nn.process(cat1)
    dr.process(cat1,cat2)
    rr.process(cat2)
    rd.process(cat2,cat1)
    xi,varxi = nn.calculateXi(rr,dr,rd)
    x = nn.meanr
    y = xi
    np.savetxt('./corr.txt',np.array([x,y,varxi]).transpose())
    
    
if __name__ == "__main__":
    ndraws = 50
    #ramin = 146.5
    #ramax = 157.5
    #decmin = -5.5
    #decmax = 5.5
    grid_type = 'rand' #hex, rect, glass
    outdir = "/global/cscratch1/sd/huikong/Obiwan/dr9_LRG/obiwan_out/deep1/"
    seed = fits.getdata(outdir+'/seed.fits')
    surveybricks = fits.getdata("/global/cfs/cdirs/cosmo/work/legacysurvey/dr9/survey-bricks.fits.gz")
    bricklist = np.loadtxt("/global/cscratch1/sd/huikong/Obiwan/dr9_LRG/obiwan_out/deep1/bricklist.txt",dtype = np.str)
    seeds = seed_maker(ndraws=ndraws,outdir=outdir, rotation=0., seed=seed, surveybricks = surveybricks, bricklist = bricklist, grid_type=grid_type, seed_num = 42)
    seeds()
    ##made for correlation function testing
    #seeds.stack(outdir+'randoms_2.fits')
    #fn1 = outdir+'randoms_2.fits'
    #fn2 = outdir+'randoms.fits'
    #corr_test(fn1,fn2)
    
#shifter --module=mpich-cle6 --image=legacysurvey/legacypipe:DR9.7.1 /bin/bash    