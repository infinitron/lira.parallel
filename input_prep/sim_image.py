#!/usr/bin/env python

from ciao_contrib.runtool import *
from sherpa.astro.ui import *
import numpy as np
import os
import sys
import yaml


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
    param_list = ['evt_file','binsize','core_reg','bkg_reg','n_psf_sims','n_null_sims','inp_size','psf_size','nH','add_gal','redshift','group','blur','center','no_core','sim_baselines']
    param_types = ['str','float','str','str','int','int','int','int','float','int','float','int','float','list','bool','bool']
    param_req = [True,False,True,False,True,False,False,True,False,False,False,False,True,False,True,False]
    default_params = [None, 0.5,None,None,50,50,64,32,None,0.0,0.0,10,0.25,None,False,True]

    param_definition = {}
    for i in range(0,len(param_list)):
        param_definition[param_list[i]] = {'type':param_types[i],'required':param_req[i],'value':default_params[i]} 
        
    return param_definition

def create_the_inputs(params):
    nsize=params['inp_size']['value']
    binsize=params['binsize']['value']
    img_file="img_{0}x{1}_{2}.fits".format(nsize,nsize,binsize)
    simulate_null_images(img_file,'core_psf.fits',params['n_null_sims']['value'],params['no_core']['value'])

def simulate_null_images(infile,psffile,num_sims,no_core,mcmciter=5000):
    print('Creating the null file')
    clean()
    set_stat("cstat")
    set_method('simplex')
    load_image(infile)
    load_psf("mypsf", psffile)
    set_psf(mypsf)

    if no_core:
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

    #the following code throws an error sometimes #bug 
    try:
        stats, accept, params = get_draws(1,niter=mcmciter)
    except:
        params = [np.repeat(q1.fwhm.val,mcmciter)]

    #print('Simulating the null files')
    for i in range(num_sims):
        set_par(q1.fwhm,params[0][(i+1)*100-1])
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



