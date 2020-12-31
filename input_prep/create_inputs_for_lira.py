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
    param_list = ['evt_file','binsize','core_reg','bkg_reg','n_psf_sims','n_null_sims','inp_size','psf_size','nH','add_gal','redshift','group','blur','center','no_core']
    param_types = ['str','float','str','str','int','int','int','int','float','int','float','int','float','list','int']
    param_req = [True,False,True,False,True,False,False,True,False,False,False,False,True,False,True]
    default_params = [None, 0.5,None,None,50,50,64,32,None,0.0,0.0,10,0.25,None,0]

    param_definition = {}
    for i in range(0,len(param_list)):
        param_definition[param_list[i]] = {'type':param_types[i],'required':param_req[i],'value':default_params[i]} 
        
    return param_definition
 
 
def create_the_inputs(params):
    
    #create the image
    image_file,core_cen,core_cen_adj = create_image(params['evt_file']['value'],
                 params['core_reg']['value'],
                 params['binsize']['value'],
                 params['inp_size']['value'],params['center']['value'])
    
    print('Image center: {0}'.format(core_cen))
    print('Image center-adjusted: {0}'.format(core_cen_adj))
    #raise()
    max_ra_dec = get_max_ra_dec(image_file,params['core_reg']['value'])
    #simulate the psf             
    simulated_psf_file = sim_core_psf(params['evt_file']['value'],
                 params['n_psf_sims']['value'],
                 params['core_reg']['value'],
                 params['bkg_reg']['value'],
                 params['binsize']['value'],
                 min(params['binsize']['value'],params['psf_size']['value']),
                max_ra_dec[0],max_ra_dec[1],
                 params['nH']['value'],params['add_gal']['value'],params['redshift']['value']
                 ,params['group']['value'],
                 params['blur']['value'])

    psf_file,psf_cen,psf_cen_adj = create_image(simulated_psf_file
        ,params['core_reg']['value']
        ,params['binsize']['value']
        ,params['psf_size']['value']
        ,None
        ,outfile='core_psf.fits')


    print('PSF centroid-adjusted: {0}'.format(psf_cen_adj))
    

    #create the baseline image and simulate images from it
    simulate_null_images(image_file,'core_psf.fits',params['n_null_sims']['value'],params['no_core']['value'])
        
def get_max_ra_dec(img_file,reg_file):
    dmstat.punlearn()
    dmstat("{0}[sky=region({1})]".format(img_file,reg_file))
    
    vals = [float(x) for x in dmstat.out_max_loc.split(',')]
    xval=vals[0]
    yval=vals[1]

    dmcoords.punlearn()
    dmcoords(img_file,option='sky',x=vals[0],y=vals[1],celfmt='deg')
    return (dmcoords.ra,dmcoords.dec)

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
 
def create_image(evt_file,reg_file,binsize,nsize,center,outfile=None):

    if center is None or len(center)==0:
        dmstat.punlearn()
        #get the centroid of the image
        dmstat("{0}[sky=region({1})][cols sky]".format(evt_file,reg_file),clip=True)
        vals = [float(x) for x in dmstat.out_mean.split(',')]
        xval=vals[0]
        yval=vals[1]
    else:
        xval,yval=center

    if(outfile is None):outfile = "img_{0}x{1}_{2}.fits".format(nsize,nsize,binsize)
    
    print("Creating the input image file")
    centroid = rebin_img(infile="{0}[energy=500:7000]".format(evt_file),outfile=outfile,
              binsize=binsize, nsize=nsize, xcen=xval,ycen=yval)
    
    return (outfile,center or vals,centroid)
    
#def create_psf_image(evt_file,reg_file,binsize,nsize,outfile=None):
#    dmstat.punlearn()
#    #get the centroid of the image
#    dmstat("{0}[sky=region({1})][cols sky]".format(evt_file,reg_file))
#    
#    vals = [float(x) for x in dmstat.out_mean.split(',')]
#    xval=vals[0]
#    yval=vals[1]
#
#    if(outfile is None):outfile = "img_{0}x{1}_{2}.fits".format(nsize,nsize,binsize)
#    
#    print("Creating the input image file")
#    centroid = rebin_psf(infile="{0}[energy=500:7000]".format(evt_file),outfile=outfile,
#              binsize=binsize, nsize=nsize, xcen=xval,ycen=yval)
#    
#    return (outfile,centroid)

def sim_core_psf(evt_file,npsf_sims,core_reg,bkg_reg,binsize,psf_size,xcen,ycen,nH=None,add_gal=0,redshift=0.0,group=10.0,blur=0.25):
    #extract the spectrum
    specextract.punlearn()
    print('Extracting the spectrum')
    specextract(infile="{0}[sky=region({1})]".format(evt_file,core_reg),outroot="core_spectrum"
    #,bkgfile="{0}[sky=region({1})]".format(evt_file,bkg_reg)
    ,clobber=True)  
    
    #fit the spectrum
    clean()
    zFlag=False
    print('Fitting the spectrum')
    load_pha("core_spectrum.pi")
    if (nH is not None) and (nH > 0.0):
        if(add_gal==1):
            set_source(xsphabs.gal*xszphabs.abs1*powlaw1d.srcp1)
            gal.nH=nH
            freeze(gal.nH)
            zFlag=True

        else:
            set_source(xsphabs.abs1*powlaw1d.srcp1)
            abs1.nH=nH
            freeze(abs1.nH)
    else:
        set_source(xszphabs.abs1*powlaw1d.srcp1)
        zFlag=True


    if zFlag is True:
        #print('REDSHIFT',redshift)
        abs1.redshift=redshift
        freeze(abs1.redshift)


    

    group_counts(1,group)
    ignore(":0.5")
    ignore("7:")
    show_filter()
    show_model()

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
    ,blur=blur  
    #,minsize=psf_size
    )
    save_all(outfile='lira_input_psfsim.log',clobber=True)
    return 'core_psf_sim_projrays.fits'

#slightly modified version of kmc's code from astrostat/LIRA            
def simulate_null_images(infile,psffile,num_sims,no_core,mcmciter=5000):
    print('Creating the null file')
    clean()
    set_stat("cstat")
    set_method('simplex')
    load_image(infile)
    load_psf("mypsf", psffile)
    set_psf(mypsf)

    if no_core>0:
        set_model(const2d.c0)
        set_par(c0.c0,min=0)
    else:
        set_model(gauss2d.q1+const2d.c0)
        set_par(c0.c0,min=0)
    #set_par(q1.fwhm,max=0.5)
        guess(q1)
    fit()
    results = get_fit_results()
    save('core_source_fit.save', clobber=True)
    save_source('null_q1_c1.fits', clobber=True)
    covar()

    if no_core:
        for i in range(num_sims):
            fake()
            save_image("sim_null_{}.fits".format(i), clobber=True)
        clean()
        return 0

    normgauss1d.g1
    g1.pos=q1.fwhm
    g1.fwhm = get_covar_results().parmaxes[0]

    #for i in range(num_sims):
    #        fake()
    #        save_image("sim_null_{}.fits".format(i), clobber=True)
    #clean()
    #return 0

    #check if there is a valid upper bound. 
    print(get_covar_results())
    if get_covar_results().parmaxes[0] is None or get_covar_results().parmins[1] is None or get_covar_results().parmins[0] is None:
        for i in range(num_sims):
            fake()
            save_image("sim_null_{}.fits".format(i), clobber=True)
        clean()
        return 0
    #if not go for the regular  
    set_prior(q1.fwhm,g1)
    set_sampler_opt('defaultprior', False)
    set_sampler_opt('priorshape', [True, False, False, False, False])
    set_sampler_opt('originalscale', [True, True, True, True, True])
    if mcmciter < num_sims*100:
        mcmciter = num_sims*100
    stats, accept, params = get_draws(1,niter=mcmciter)
    #print('Simulating the null files')
    for i in range(num_sims):
        set_par(q1.fwhm,params[0][(i+1)*100])
        fake()
        save_image("sim_null_{}.fits".format(i), clobber=True)
    save_all(outfile='lira_input_baseline_sim.log',clobber=True)
    clean()


with open("lira_input_prep.yaml", 'r') as stream:
    try:
        params = create_the_inputs(
            prepare_the_params(yaml.safe_load(stream))
        )
            
            
    except yaml.YAMLError as exc:
        print(exc)
