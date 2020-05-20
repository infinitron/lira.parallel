#!/usr/bin/env python

import yaml
from rebin_img import rebin_img
from rebin_psf import rebin_psf
from ciao_contrib.runtool import *
from sherpa.astro.ui import *
from sherpa_contrib.chart import *
from sim_image import sim_image
import os

def prepare_the_params(params):
    #parse the params and return the 
    param_def = get_param_definition()
    
    #go through all the params and edit the if a user supplies it
    for key,defn in param_def.items():
        if(defn['type']!=type(params[key]).__name__): raise("Incorrect value for {0}".format(key))
        if(defn['required']) and not key in params: raise("{0} is required".format(key))
        if key in params: param_def[key]['value'] = params[key]
    
    return param_def
                
def get_param_definition():
    param_list = ['evt_file','binsize','core_reg','bkg_reg','n_psf_sims','n_null_sims','inp_size','psf_size','nH']
    param_types = ['str','float','str','str','int','int','int','int','float']
    param_req = [True,False,True,False,True,False,False,True,False]
    default_params = [None, 0.5,None,None,50,50,64,32,None]

    param_definition = {}
    for i in range(0,len(param_list)):
        param_definition[param_list[i]] = {'type':param_types[i],'required':param_req[i],'value':default_params[i]} 
        
    return param_definition
 
 
def create_the_inputs(params):
    
    #create the image
    image_file,core_cen = create_image(params['evt_file']['value'],
                 params['core_reg']['value'],
                 params['binsize']['value'],
                 params['inp_size']['value'])

    #simulate the psf             
    simulated_psf_file = sim_core_psf(params['evt_file']['value'],
                 params['n_psf_sims']['value'],
                 params['core_reg']['value'],
                 params['bkg_reg']['value'],
                 params['binsize']['value'],
                 min(params['binsize']['value'],params['psf_size']['value']),
                 core_cen[0],core_cen[1],
                 params['nH']['value'])
#
    #psf_file,_cen = create_image('{0}[energy=500:7000]'.format(simulated_psf_file),
    #             params['core_reg']['value'],
    #             params['binsize']['value'],
    #             params['psf_size']['value'],outfile='core_psf.fits')

    #create_image(simulated_psf_file
    #    ,params['core_reg']['value']
    #    ,params['binsize']['value']
    #    ,params['psf_size']['value']
    #    ,outfile='core_psf.fits')

    #rebin the PSF to match the image grid
    rebin_psf('{0}[energy=500:7000]'.format(simulated_psf_file),
            outfile='core_psf.fits',
            binsize=params['binsize']['value'],
            nsize=params['psf_size']['value'],
            xcen=core_cen[0],
            ycen=core_cen[1])
    

    #create the baseline image and simulate images from it
    simulate_null_images(image_file,'core_psf.fits',params['n_null_sims']['value'])
        
        
def get_centroid_physical(evt_file,reg_file):
    dmstat.punlearn()
    dmstat("{0}[sky=region({1})][cols sky]".format(evt_file,reg_file))
    
    vals = [float(x) for x in dmstat.out_mean.split(',')]
    xval=vals[0]
    yval=vals[1]
    
    return (xval,yval)

def get_centroid_ra_dec(evt_file,reg_file):
    dmcoords.punlearn()
    cntrd_phys = get_centroid_physical(evt_file,reg_file)
    dmcoords(evt_file,option='sky',x=cntrd_phys[0],y=cntrd_phys[1],celfmt='deg')
    return (dmcoords.ra,dmcoords.dec)
 
def create_image(evt_file,reg_file,binsize,nsize,outfile=None):
    dmstat.punlearn()
    #get the centroid of the image
    dmstat("{0}[sky=region({1})][cols sky]".format(evt_file,reg_file))
    
    vals = [float(x) for x in dmstat.out_mean.split(',')]
    xval=vals[0]
    yval=vals[1]

    if(outfile is None):outfile = "img_{0}x{1}_{2}.fits".format(nsize,nsize,binsize)
    
    print("Creating the input image file")
    centroid = rebin_img(infile="{0}[energy=500:7000]".format(evt_file),outfile=outfile,
              binsize=binsize, nsize=nsize, xcen=xval,ycen=yval)
    
    return (outfile,vals)
    
def create_psf_image(evt_file,reg_file,binsize,nsize,outfile=None):
    dmstat.punlearn()
    #get the centroid of the image
    dmstat("{0}[sky=region({1})][cols sky]".format(evt_file,reg_file))
    
    vals = [float(x) for x in dmstat.out_mean.split(',')]
    xval=vals[0]
    yval=vals[1]

    if(outfile is None):outfile = "img_{0}x{1}_{2}.fits".format(nsize,nsize,binsize)
    
    print("Creating the input image file")
    centroid = rebin_psf(infile="{0}[energy=500:7000]".format(evt_file),outfile=outfile,
              binsize=binsize, nsize=nsize, xcen=xval,ycen=yval)
    
    return (outfile,centroid)

def sim_core_psf(evt_file,npsf_sims,core_reg,bkg_reg,binsize,psf_size,xcen,ycen,nH=None):
    #extract the spectrum
    specextract.punlearn()
    print('Extracting the spectrum')
    specextract(infile="{0}[sky=region({1})]".format(evt_file,core_reg),outroot="core_spectrum",bkgfile="{0}[sky=region({1})]".format(evt_file,bkg_reg),clobber=True)
    
    #fit the spectrum
    print('Fitting the spectrum')
    load_pha("core_spectrum.pi")
    set_source(xsphabs.abs1*powlaw1d.srcp1)
    if (nH is not None) and (nH > 0.0):
        abs1.nH=nH
        freeze(abs1.nH)
    group_counts(1,10)
    notice(0.5,7)
    fit()
    covar()
    save_chart_spectrum("core_flux_chart.dat", elow=0.5, ehigh=7.0)
    clean()
    
    #get the ra dec from sky coords
    ra_dec = get_centroid_ra_dec(evt_file,core_reg)
    #dmcoords.punlearn()
    #dmcoords(evt_file,option='sky',x=xcen,y=ycen,celfmt='deg')
    
    #simulate the psf
    print('Simulating the psf')
    simulate_psf.punlearn()
    simulate_psf(infile=evt_file,outroot="core_psf_sim",spectrum='core_flux_chart.dat',numiter=npsf_sims,ra=ra_dec[0],dec=ra_dec[1],binsize=binsize
    #,minsize=psf_size
    )
    return 'core_psf_sim_projrays.fits'

            
def simulate_null_images(infile,psffile,num_sims,mcmciter=5000):
    print('Creating the null file')
    set_stat("cstat")
    set_method('simplex')
    load_image(infile)
    load_psf("mypsf", psffile)
    set_psf(mypsf)
    set_model(gauss2d.q1+const2d.c0)
    guess(q1)
    fit()
    results = get_fit_results()
    save('core_source_fit.save', clobber=True)
    save_source('null_q1_c1.fits', clobber=True)
    covar()
    normgauss1d.g1
    g1.pos=q1.fwhm
    g1.fwhm = get_covar_results().parmaxes[0]
    set_prior(q1.fwhm,g1)
    set_sampler_opt('defaultprior', False)
    set_sampler_opt('priorshape', [True, False, False, False, False])
    set_sampler_opt('originalscale', [True, True, True, True, True])
    if mcmciter < num_sims*100:
        mcmciter = num_sims*100
    stats, accept, params = get_draws(1,niter=mcmciter)
    print('Simulating the null files')
    for i in range(num_sims):
        set_par(q1.fwhm,params[0][(i+1)*100])
        fake()
        save_image("sim_null_{}.fits".format(i), clobber=True)
    clean()


with open("lira_input_prep.yaml", 'r') as stream:
    try:
        params = create_the_inputs(
            prepare_the_params(yaml.safe_load(stream))
        )
            
            
    except yaml.YAMLError as exc:
        print(exc)
