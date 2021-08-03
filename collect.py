import os
import glob
import numpy as n
import numpy as np
from astropy.io import fits
from math import *
from astropy.table import Column, Table
from astropy.coordinates import SkyCoord
from astropy import units as u
import subprocess
from astropy.table import vstack,Table
from filesystem import LegacySimData
import subprocess
from SurveySource import BaseSource
import os

class Collect(LegacySimData):
    def __init__(self, survey_dir=None, outdir=None, subsection=None, brick=None, bricklist=None, **kwargs):
        super(Collect,self).__init__(survey_dir=survey_dir, outdir=outdir, subsection=subsection, brick=brick)
    
    def brick_match(self, threads = None, bricklist = None, mp=None, subsection=None, startid=None, nobj=None, angle=1.5/3600, mode='sim',tracer=None):
        inputs = []
        self.subsection = subsection
        bricklist_split = np.array_split(bricklist,threads)
        for i in range(threads):
            inputs.append((bricklist_split[i],startid, nobj, angle, mode, True,tracer))
        results = list(mp.map(self._brick_match_core, inputs))
        final_tab = None
        if len(results)>0:
            for tab in results:
                if tab is None:
                    continue
                if final_tab is not None:
                    final_tab = vstack((final_tab,tab))
                else:
                    final_tab = tab
            topdir = os.path.dirname(os.path.abspath(self.outdir))+'/subset'
            subprocess.call(["mkdir","-p",topdir])
            final_tab.write(topdir+'/subset_%s.fits'%(self.subsection),overwrite=True)
            print('written %s/subset_%s.fits'%(topdir,self.subsection))
            
    def brick_match_dr9(self, bricklist, tracer = None,region='south',south=True, threads = None, mp=None):
        surveydir = "/global/cfs/cdirs/cosmo/work/legacysurvey/dr9/%s/"%region
        outdir = "/global/cfs/cdirs/cosmo/work/legacysurvey/dr9/%s/"%region
        final_tab = None
        bricklist_split = np.array_split(bricklist,threads)
        inputs = []
        for i in range(threads):
            inputs.append((bricklist_split[i],outdir,tracer,south))
        results = mp.map(self._brick_match_dr9, inputs)
        final_tab = None
        for result in results:
            if final_tab is None:
                final_tab = result
            else:
                final_tab = vstack((final_tab,result))
        topdir = os.path.dirname(os.path.abspath(self.outdir))+'/subset'
        final_tab.write(topdir+'/subset_dr9_%s.fits'%tracer,overwrite=True)    
        print('written %s/subset_dr9_%s.fits'%(topdir,tracer))    
    
    def _brick_match_dr9(self, X):
        (bricklist, outdir,tracer,south) = X
        final_tab = None
        for brickname in bricklist:
            catalog = BaseSource(filetype='tractor', survey_dir=outdir, outdir=outdir,subsection=None, brick=brickname)
            tracer_sel = catalog.target_selection(target=tracer,south=south)
            tractor_fn = catalog.find_file('tractor')
            T = Table.read(tractor_fn)
            T_new = T[tracer_sel]
            if final_tab is not None:
                final_tab = vstack((final_tab,T_new))
            else:
                final_tab = T_new
        return final_tab
    
    def brick_match_random(self, bricklist, threads = None, mp = None):
        outdir = os.path.dirname(os.path.abspath(self.outdir))+'/randoms/'
        final_tab = None
        
        bricklist_split = np.array_split(bricklist,threads)
        
        
        inputs = []
        for i in range(threads):
            inputs.append((bricklist_split[i],outdir))
        results = mp.map(self._brick_match_random, inputs)
        final_tab = None
        for result in results:
            if final_tab is None:
                final_tab = result
            else:
                final_tab = vstack((final_tab,result))
        topdir = os.path.dirname(os.path.abspath(self.outdir))+'/subset'
        final_tab.write(topdir+'/subset_random.fits',overwrite=True)    
        print('written %s/subset_random.fits'%(topdir)) 
        
    def _brick_match_random(self, X):
        (bricklist, outdir) = X
        final_tab = None
        for brickname in bricklist:
            tractor_fn = outdir+'random-%s.fits'%brickname
            T = Table.read(tractor_fn)
            if final_tab is not None:
                final_tab = vstack((final_tab,T))
            else:
                final_tab = T
        return final_tab
        
    def _brick_match_core(self, X):
        '''
        mode:
        -- sim: match outputs to sim, output shares same dimension as sim
        -- tractor: match outputs to tractor, output shares same dimension as tractor
        '''
        (bricklist, startid, nobj, angle, mode, MS_star, tracer) = X
        final_tab = None
        for brickname in bricklist:
            
            catalog = BaseSource(filetype='tractor', survey_dir=self.survey_dir, outdir=self.outdir,subsection=self.subsection, brick=brickname)
            tracer_sel = catalog.target_selection(tracer,south=True)
            
            
            sim = Table.read(self.find_file('simcat',brick=brickname))
            tractor = Table.read(self.find_file('tractor',brick=brickname))
            tractor[tracer] = tracer_sel
            #MS stars
            try:
                original_sim = Table.read(self.find_file('simorigin',brick=brickname))[startid:startid+nobj] 
            except:
                original_sim = Table.read(self.find_file('simorigin',brick=brickname))
            
            c1 = SkyCoord(ra=sim['ra']*u.degree, dec=sim['dec']*u.degree)
            c2 = SkyCoord(ra=np.array(tractor['ra'])*u.degree, dec=np.array(tractor['dec'])*u.degree)
            c3 = SkyCoord(ra=original_sim['ra']*u.degree, dec=original_sim['dec']*u.degree)
            if mode == 'sim':
                idx1, d2d, d3d = c1.match_to_catalog_sky(c2)
                idx2, d2d2, d3d2 = c1.match_to_catalog_sky(c3)
                matched = d2d.value <= angle
                distance = d2d.value
                tc = tractor[idx1]
                ors = original_sim[idx2]
            
            elif mode == 'tractor':
                idx1, d2d, d3d = c2.match_to_catalog_sky(c1)
                idx2, d2d, d3d = c2.match_to_catalog_sky(c3)
                matched = d2d.value <= angle
                distance = d2d.value
                tc = tractor
                sim = sim[idx1]
                ors = original_sim[idx2]
    
            tc.add_column(sim['ra'],name = 'sim_ra')
            tc.add_column(sim['dec'],name = 'sim_dec')
            tc.add_column(sim['gflux'],name = 'sim_gflux')
            tc.add_column(sim['rflux'],name='sim_rflux')
            tc.add_column(sim['zflux'],name='sim_zflux')
            tc.add_column(ors['w1'],name='sim_w1')
            tc.add_column(ors['id_sample'],name='id_sample')
            tc.add_column(ors['redshift'],name='sim_redshift')
            tc.add_column(sim['rhalf'],name='sim_rhalf')
            tc.add_column(sim['e1'],name='sim_e1')
            tc.add_column(sim['e2'],name='sim_e2')
            tc.add_column(sim['x'],name='sim_bx')
            tc.add_column(sim['y'],name='sim_by')
            tc['angle'] = np.array(d2d.value*3600.,dtype=np.float)
            tc['matched'] = np.array(matched,dtype=np.bool)
            tc.add_column(sim['n'],name='sim_sersic_n')
            if MS_star:
                #match closest MS star to the output
                metric = fits.getdata(self.find_file('ref-sources',brick=brickname))
                c1 = SkyCoord(ra=tc['ra'],dec=tc['dec'])
                c2 = SkyCoord(ra=metric['ra']*u.degree, dec=metric['dec']*u.degree)
                idx1, d2d, d3d = c1.match_to_catalog_sky(c2)
                mt = metric[idx1]
                distance = d2d.value
                tc['star_distance'] = np.array(distance,dtype=np.float)
                tc['star_radius'] = np.array(mt['radius'], dtype=np.float)
                tc['MS_delta_ra'] = np.array(mt['ra']-tc['ra'],dtype=np.float)
                tc['MS_delta_dec']= np.array(mt['dec']-tc['dec'],dtype=np.float)
            if final_tab is None:
                final_tab = tc
            else:
                final_tab = vstack((final_tab,tc))
        return final_tab

    
    
def add_more_TS(prefix, survey_dir, outdir, startid, set_num, south=True):
    #TS on cards like 'flux_g%s'%prefix
    catalog = BaseSource(filetype='processed_one_rsp', survey_dir=survey_dir, outdir=outdir,subsection='rs%d_cosmos%d'%(startid,set_num),brick=None, force_construct=True)
    LRG_rsp = catalog.target_selection('LRG_sv3',south=south, prefix='_rsp')
    fn = catalog.find_file('processed_one_rsp')
    T = Table.read(fn)
    T['lrg_sv3%s'%prefix] = LRG_rsp
    T.write(fn, overwrite=True)
    print("written %s"%fn)
if __name__ == "__main__":
    survey_dir = "/global/project/projectdirs/cosmo/data/legacysurvey/dr9/" 
    outdir = "/global/cscratch1/sd/huikong/Obiwan/dr9_LRG/obiwan_out/cosmos_resampled_seed/output/"
    startid = 9999
    set_num = 80
    
    add_more_TS(prefix = "_rsp",survey_dir = survey_dir, outdir = outdir, startid = startid, set_num = set_num)