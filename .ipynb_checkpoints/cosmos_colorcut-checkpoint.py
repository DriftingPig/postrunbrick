import subprocess
from glob import glob
import astropy.io.fits as fits
from astropy.table import Table
from SurveySource import BaseSource
import os
import numpy as np
'''
currently considering baseline LR selection
https://github.com/desihub/desitarget/blob/master/py/desitarget/sv1/sv1_cuts.py#L341
'''
def flux2mag(flux,mwtransmission):
        mag= 22.5 - 2.5 * np.log10(flux / mwtransmission)
        return mag

def get_scatter_ratio(equation, target_type='lrgir'):
        #scatter info for cosmos 
        topdir = '/global/cscratch1/sd/huikong/Obiwan/dr9_LRG/obiwan_out/cosmos_subsets/cosmos_all_stacked/'
        for n in range(10):
            dat = fits.getdata(topdir+'cosmos_set%d.fits'%n)
            #forgot to apply wisemask previously, do it here...
            sel = (dat['WISEMASK_W1']==0)
            dat = dat[sel]
            cut1_list = []
            ratio = []
            for i in range(10):
                if i==n:
                    continue
                sel_target = dat['set_%d_%s'%(n,target_type)]&dat['set_%d_matched'%i]
                zmag = flux2mag(dat[sel_target]['set_%d_flux_z'%i],dat[sel_target]['MW_TRANSMISSION_Z'])
                w1mag = flux2mag(dat[sel_target]['set_%d_flux_w1'%i],dat[sel_target]['MW_TRANSMISSION_W1'])
                rmag = flux2mag(dat[sel_target]['set_%d_flux_r'%i],dat[sel_target]['MW_TRANSMISSION_R'])
                gmag = flux2mag(dat[sel_target]['set_%d_flux_g'%i],dat[sel_target]['MW_TRANSMISSION_G'])
                zfibermag = flux2mag(dat[sel_target]['set_%d_fiberflux_z'%i],dat[sel_target]['MW_TRANSMISSION_Z'])
                cut1 = ~eval(equation)
                ratio.append(cut1.sum()/len(cut1))
            ave_ratio = np.mean(ratio)
            return ave_ratio
        

def get_cotamin_ratio(equation, equation_origin):
        #scatter info for cosmos 
        topdir = '/global/cscratch1/sd/huikong/Obiwan/dr9_LRG/obiwan_out/cosmos_subsets/cosmos_all_stacked/'
        for n in range(10):
            dat = fits.getdata(topdir+'cosmos_set%d.fits'%n)
            #forgot to apply wisemask previously, do it here...
            sel = (dat['WISEMASK_W1']==0)
            dat = dat[sel]
            cut1_list = []
            ratio = []
            for i in range(10):
                if i==n:
                    continue
                sel_target = dat['set_%d_matched'%i]
                zmag = flux2mag(dat[sel_target]['set_%d_flux_z'%i],dat[sel_target]['MW_TRANSMISSION_Z'])
                w1mag = flux2mag(dat[sel_target]['set_%d_flux_w1'%i],dat[sel_target]['MW_TRANSMISSION_W1'])
                rmag = flux2mag(dat[sel_target]['set_%d_flux_r'%i],dat[sel_target]['MW_TRANSMISSION_R'])
                gmag = flux2mag(dat[sel_target]['set_%d_flux_g'%i],dat[sel_target]['MW_TRANSMISSION_G'])
                zfibermag = flux2mag(dat[sel_target]['set_%d_fiberflux_z'%i],dat[sel_target]['MW_TRANSMISSION_Z'])
                cut1 = eval(equation)
                cut2 = eval(equation_origin)
                ratio_i = (cut1&~cut2).sum()/cut2.sum()
                ratio.append(ratio_i)
            ave_ratio = np.mean(ratio)
            return ave_ratio
        

def cosmos_colorcut(step=0.01, change=0.002):
        #a module dedicated to checking LRG color cuts using cosmos repeats data
        cut1 = 'zmag - w1mag > 0.8 * (rmag-zmag) - 0.6'
        cut2 = 'zfibermag < 21.5'
        cut3 = 'rmag - w1mag > 1.1'
        cut4 = 'rmag - w1mag > (w1mag - 17.22) * 1.8 '
        cut5 = 'rmag - w1mag > w1mag - 16.37'
        cut_all = '(%s)&(%s)&(%s)&(%s)&(%s)'%(cut1, cut2, cut3, cut4, cut5)
        #the original color cut
        ratio1_old = get_scatter_ratio(cut1)
        ratio2_old = get_scatter_ratio(cut2)
        ratio3_old = get_scatter_ratio(cut3)
        ratio4_old = get_scatter_ratio(cut4)
        ratio5_old = get_scatter_ratio(cut5)
        ratio_all_old = get_scatter_ratio(cut_all)
        ratio1_new = ratio1_old
        ratio2_new = ratio2_old
        ratio3_new = ratio3_old
        ratio4_new = ratio4_old
        ratio5_new = ratio5_old
        
        cut1_old = cut1
        cut2_old = cut2
        cut3_old = cut3
        cut4_old = cut4
        cut5_old = cut5
        cutall_old = cut_all
        
        cutlist1 = [ratio1_old]
        cutlist2 = [ratio2_old]
        cutlist3 = [ratio3_old]
        cutlist4 = [ratio4_old]
        cutlist5 = [ratio5_old]
        cutlist6 = [ratio_all_old]
        #cotamination ratio
        cutlist7 = [0]
        
        cut1_N=0;cut2_N=0;cut3_N=0;cut4_N=0;cut5_N=0
        l_cut1 = [cut1_N];l_cut2 = [cut2_N];l_cut3 = [cut3_N];l_cut4 = [cut4_N];l_cut5 = [cut5_N]
        while(ratio_all_old>0.01):
            #import pdb;pdb.set_trace()
            print(ratio_all_old)
            max_ratio = max(ratio1_old, ratio2_old, ratio3_old, ratio4_old, ratio5_old)
            if max_ratio == ratio1_old:
                while(ratio1_old-ratio1_new<change):
                    cut1_N+=1
                    cut1_new = '%s - %f'%(cut1, step*cut1_N)
                    ratio1_new = get_scatter_ratio(cut1_new)    
                cut1_old = cut1_new
                ratio1_old = ratio1_new
            if max_ratio == ratio2_old:
                while(ratio2_old - ratio2_new < change):
                    cut2_N += 1
                    cut2_new = '%s + %f'%(cut2, step*cut2_N)
                    ratio2_new = get_scatter_ratio(cut2_new)
                ratio2_old = ratio2_new
                cut2_old = cut2_new
            if max_ratio == ratio3_old:
                while(ratio3_old-ratio3_new<change):
                    cut3_N+=1
                    cut3_new = '%s - %f'%(cut3, step*cut3_N)
                    ratio3_new = get_scatter_ratio(cut3_new)    
                cut3_old = cut3_new
                ratio3_old = ratio3_new
            if max_ratio == ratio4_old:
                while(ratio4_old-ratio4_new<change):
                    cut4_N+=1
                    cut4_new = '%s - %f'%(cut4, step*cut4_N)
                    ratio4_new = get_scatter_ratio(cut4_new)    
                cut4_old = cut4_new
                ratio4_old = ratio4_new
            if max_ratio == ratio5_old:
                while(ratio5_old-ratio5_new<change):
                    cut5_N+=1
                    cut5_new = '%s - %f'%(cut5, step*cut5_N)
                    ratio5_new = get_scatter_ratio(cut5_new)    
                cut5_old = cut5_new
                ratio5_old = ratio5_new
            cut_all = '(%s)&(%s)&(%s)&(%s)&(%s)'%(cut1_old, cut2_old, cut3_old, cut4_old, cut5_old)
            ratio_all_old = get_scatter_ratio(cut_all)
            ratio_cotamin = get_cotamin_ratio(cut_all,cutall_old)
            print(ratio1_old,ratio2_old,ratio3_old,ratio4_old,ratio5_old,ratio_all_old,ratio_cotamin)
            cutlist1.append(ratio1_old)
            cutlist2.append(ratio2_old)
            cutlist3.append(ratio3_old)
            cutlist4.append(ratio4_old)
            cutlist5.append(ratio5_old)
            cutlist6.append(ratio_all_old)
            cutlist7.append(ratio_cotamin)
            l_cut1.append(cut1_N)
            l_cut2.append(cut2_N)
            l_cut3.append(cut3_N)
            l_cut4.append(cut4_N)
            l_cut5.append(cut5_N)
        t = Table()
        t['cut1']=np.array(cutlist1)
        t['cut2']=np.array(cutlist2)
        t['cut3']=np.array(cutlist3)
        t['cut4']=np.array(cutlist4)
        t['cut5']=np.array(cutlist5)
        t['cut6']=np.array(cutlist6)
        t['cut7']=np.array(cutlist7)
        t['cut1_n']=np.array(l_cut1)
        t['cut2_n']=np.array(l_cut2)
        t['cut3_n']=np.array(l_cut3)
        t['cut4_n']=np.array(l_cut4)
        t['cut5_n']=np.array(l_cut5)
        t.write('/global/cscratch1/sd/huikong/Obiwan/dr9_LRG/obiwan_out/cosmos_subsets/cuts_info.fits',overwrite=True)
        
cosmos_colorcut()