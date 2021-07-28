import numpy as np
import astropy.io.fits as fits

#also selecting dec>-30 for decals south only
region = 'south'
topdir = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/'+region+'/sweep/9.0/'
fns = glob.glob(topdir+'sweep*')
sweep_keys= ['BRICKNAME','RA','DEC','TYPE','OBJID','EBV','FLUX_G','FLUX_R','FLUX_Z','FLUX_W1','FLUX_W2','FLUX_IVAR_G','FLUX_IVAR_R','FLUX_IVAR_Z','FLUX_IVAR_W1','FLUX_IVAR_W2','MW_TRANSMISSION_G','MW_TRANSMISSION_R','MW_TRANSMISSION_Z','MW_TRANSMISSION_W1','MW_TRANSMISSION_W2','NOBS_G','NOBS_R','NOBS_Z','NOBS_W1','NOBS_W2','SHAPE_R','SHAPE_E1','SHAPE_E2','FIBERFLUX_G','FIBERFLUX_R','FIBERFLUX_Z','MASKBITS','SERSIC','DCHISQ','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z','PSFDEPTH_G','PSFDEPTH_R','PSFDEPTH_Z','GALDEPTH_G','GALDEPTH_R','GALDEPTH_Z','WISEMASK_W1','WISEMASK_W2','ANYMASK_G','ANYMASK_R','ANYMASK_Z','BX','BY','GAIA_PHOT_G_MEAN_MAG']


class isLRG_like(object):
    def __init__(self,method):
        exec("self.%s()"%method)
    def LRG_basline_like(self):#1.4% LRG lost rate, south
        self.cut1 = 'zmag - w1mag > 0.8 * (rmag-zmag) - 0.6'
        self.cut2 = 'zfibermag < 21.5+0.13'
        self.cut3 = 'rmag - w1mag > 1.1-0.04'
        self.cut4 = 'rmag - w1mag > (w1mag - 17.22) * 1.8 -0.42'
        self.cut5 = 'rmag - w1mag > w1mag - 16.37-0.18'
        self.structure = '(%s)&(%s)&(%s)&(%s)&(%s)'%(self.cut1, self.cut2, self.cut3, self.cut4, self.cut5)
    def LRG_sv3_like(self):#south 
        #adding 0.01 to 0s
        self.cut1 = 'zmag - w1mag > 0.8 * (rmag - zmag) - 0.6-0.010000'
        self.cut2 = 'zfibermag < 21.7+0.130000'
        self.cut3_1 = 'gmag - rmag > 1.3-0.970000'
        self.cut3_2 = '(gmag - rmag) > -1.55 * (rmag - w1mag)+3.13-0.010000'
        self.cut3_3 = 'rmag - w1mag > 1.8-0.470000'
        self.cut4_1 = 'rmag - w1mag > (w1mag - 17.26) * 1.8-0.940000'
        self.cut4_2 = 'rmag - w1mag > (w1mag - 16.36) * 1.-0.070000'
        self.cut4_3 = 'rmag - w1mag > 3.29-1.960000'
        self.structure = '(%s)&(%s)&((%s)&(%s)|(%s))&((%s)&(%s)|(%s))'%(self.cut1, self.cut2, self.cut3_1, self.cut3_2, self.cut3_3, self.cut4_1, self.cut4_2,  self.cut4_3)
    def isLRGlike_color(self,data):
        #assume the data are sweep files, things are in upper case
        rflux = data['FLUX_R']/data['MW_TRANSMISSION_R']
        zflux = data['FLUX_Z']/data['MW_TRANSMISSION_Z']
        w1flux = data['FLUX_W1']/data['MW_TRANSMISSION_W1']
        zfiberflux=data['FIBERFLUX_Z']/data['MW_TRANSMISSION_Z']
        gnobs = data['NOBS_G']
        rnobs = data['NOBS_R']
        znobs = data['NOBS_Z']
        rfluxivar=data['FLUX_IVAR_R']
        zfluxivar=data['FLUX_IVAR_Z']
        w1fluxivar=data['FLUX_IVAR_W1']
        gaiagmag=data['GAIA_PHOT_G_MEAN_MAG']
        maskbits = data['MASKBITS']
        zfibertotflux = data['FIBERTOTFLUX_Z']/data['MW_TRANSMISSION_Z']
        wisemask = data['WISEMASK_W1']
        R = locals()
        R.pop('self');R.pop('data')
        lrg = self.notinLRG_mask(primary=None, **R)
        assert(hasattr(self,'structure'))
        lrg &= eval(self.structure)
        return lrg
    @staticmethod
    def notinLRG_mask(primary=None, rflux=None, zflux=None, w1flux=None,
                  zfiberflux=None, gnobs=None, rnobs=None, znobs=None,
                  rfluxivar=None, zfluxivar=None, w1fluxivar=None,
                  gaiagmag=None, maskbits=None, zfibertotflux=None,wisemask=None):
        """See :func:`~desitarget.cuts.isLRG` for details.
        Returns
        -------
        :class:`array_like`
            ``True`` if and only if the object is NOT masked for poor quality.
        """
        if primary is None:
            primary = np.ones_like(rflux, dtype='?')
        lrg = primary.copy()

        # ADM to maintain backwards-compatibility with mocks.
        if zfiberflux is None:
            log.warning('Setting zfiberflux to zflux!!!')
            zfiberflux = zflux.copy()

        lrg &= (rfluxivar > 0) & (rflux > 0)   # ADM quality in r.
        lrg &= (zfluxivar > 0) & (zflux > 0) & (zfiberflux > 0)   # ADM quality in z.
        lrg &= (w1fluxivar > 0) & (w1flux > 0)  # ADM quality in W1.

        lrg &= (gaiagmag == 0) | (gaiagmag > 18)  # remove bright GAIA sources

        # ADM remove stars with zfibertot < 17.5 that are missing from GAIA.
        lrg &= zfibertotflux < 10**(-0.4*(17.5-22.5))

        # ADM observed in every band.
        lrg &= (gnobs > 0) & (rnobs > 0) & (znobs > 0)

        # ADM default mask bits from the Legacy Surveys not set.
        lrg &= (maskbits==0)
    
        lrg &= (wisemask==0)

        return lrg


def process_one_set(X):
    (fns_i,set_num,rank,savedir)=X
    cut_method = isLRG_like('LRG_sv3_like')
    idx=0
    for fn in fns_i:
        data = Table.read(fn)
        sel = cut_method.isLRGlike_color(data)
        source = Table()
        for key in sweep_keys:
            source[key]=dat_i[source_sel][key]
        source.write(savedir+'/desi_source_rank%d_set%d_%d.fits'%(rank,set_num,idx),overwrite=True)
        idx+=1
        
    

def main():
    from mpi4py import MPI
    node_tot = 10
    rank = MPI.COMM_WORLD.Get_rank()
    node_sub_fns = np.array_split(fns,node_tot)
    this_sub_fn = node_sub_fns[rank]
    N=16
    sub_fns = np.array_split(this_sub_fn,N)
    commands=[]
    obj = 'lrg_sv3'
    savedir = '/global/cscratch1/sd/huikong/Obiwan/dr9_LRG/obiwan_out/truth/sweep_%s/'%obj
    import subprocess
    subprocess.call(["mkdir",'-p',savedir])
    for i in range(N):
        print(len(sub_fns[i]))
        commands.append((sub_fns[i],i,rank, savedir))
    with mp.Pool(N) as p:
        p.map(process_one_set,commands)
        
main()