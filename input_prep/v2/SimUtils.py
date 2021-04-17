from collections import namedtuple
from typing import Any, Dict, List, Match, NoReturn, Optional, Tuple, Union
import numpy as np
from scipy.special import factorial
from ciao_contrib.runtool import (
    simulate_psf
)
from sherpa.astro.ui import (
    set_stat,
    set_method,
    load_image,
    load_psf,
    set_psf,
    clean,
    set_model,
    set_par,
    guess,
    fit,
    save,
    save_source,
    covar,
    fake,
    save_image,
    normgauss1d,
    get_covar_results,
    set_prior,
    set_sampler_opt,
    get_draws,
    save_all,
    const2d,
    gauss2d,
    get_fit_results
    
)
from sherpa_contrib.chart import *
import os
from math import floor
from TypeDefs import FluxObj, PointCel,PointPhys, LiraInputConfig
from CoordUtils import CoordUtils


class SimUtils():
    def __init__(self) -> None:
        pass

    @staticmethod
    def simulate_null_images(img_file:str,psf_file:str,n_null_sims:int,no_core:bool=False,mcmciter:int=5000,**kwargs)->None:

        """
        Simulates a specified number of baseline images for a given input observation and a psf file

        :param img_file: Path to the input image file
        :param psf_file: Path to the psf image file
        :param n_null_sims: Number of baseline replicates to be simulated
        :param no_core: Setting this to True will only generate baseline replicates with a flat background while the default value includes a point source at the location of the core
        :param mcmciter: The number of MCMC samples to draw for simulating the baselines
        """
        print("Creating the null file")
        clean()
        set_stat("cstat")
        set_method("simplex")
        load_image(img_file)
        load_psf("mypsf", psf_file)
        set_psf(mypsf)

        if no_core:
            set_model(const2d.c0)
            set_par(c0.c0, min=0)
        else:
            set_model(gauss2d.q1 + const2d.c0)
            set_par(c0.c0, min=0)
            # set_par(q1.fwhm,max=0.5)
            guess(q1)
        fit()
        results = get_fit_results()
        save("core_source_fit.save", clobber=True)
        save_source("null_q1_c1.fits", clobber=True)
        covar()

        if no_core:
            for i in range(n_null_sims):
                fake()
                save_image("sim_null_{}.fits".format(i), clobber=True)
            clean()
            return

        normgauss1d.g1
        g1.pos = q1.fwhm
        g1.fwhm = get_covar_results().parmaxes[0]

        # check if there is a valid upper bound.
        print(get_covar_results())
        if (
            get_covar_results().parmaxes[0] is None
            or get_covar_results().parmins[1] is None
            or get_covar_results().parmins[0] is None
        ):
            for i in range(n_null_sims):
                fake()
                save_image("sim_null_{}.fits".format(i), clobber=True)
            clean()
            return
        # if not go for the regular
        set_prior(q1.fwhm, g1)
        set_sampler_opt("defaultprior", False)
        set_sampler_opt("priorshape", [True, False, False, False, False])
        set_sampler_opt("originalscale", [True, True, True, True, True])
        if mcmciter < n_null_sims * 100:
            mcmciter = n_null_sims * 100

        # the following code throws an error sometimes #bug
        try:
            stats, accept, params = get_draws(1, niter=mcmciter)
        except:
            params = [np.repeat(q1.fwhm.val, mcmciter)]

        # print('Simulating the null files')
        for i in range(n_null_sims):
            set_par(q1.fwhm, params[0][(i + 1) * 100 - 1])
            fake()
            save_image("sim_null_{}.fits".format(i), clobber=True)
        save_all(outfile="lira_input_baseline_sim.log", clobber=True)
        clean()


    @staticmethod
    def simulate_chandra_psf(evt_file:str,reg_file:str,n_psf_sims:int,blur:float,binsize:Union[int,float],flux_obj:Optional[FluxObj]=None,simulator:str='saotrace',**kwargs)->str:
        """
        Simulate a PSF of a point source for a given chandra observation using SAOTrace or MARX

        :param evt_file: Path for the Chandra events file
        :param regfile: Path to the region file enclosing the core
        :param n_psf_sims: Number of PSFs to be simulated for averaging
        :param photon_flux: The simulation will be performed using an integrated flux value. ``monoenergy`` must be specified for this to work
        :param monoenergy: The representative energy for the specified photon flux
        ...
        :return: Return the path to the simulated events file
        """

        coords = CoordUtils.get_centroid(evt_file,reg_file,'cel')

        if simulator not in ('saotrace','marx'):
            raise ValueError("Simulator must be either one of saotrace or marx")

        if simulator == 'saotrace':
            rayfiles = SAOTraceUtils.run_sao_raytrace_parallel(evt_file,coords,n_psf_sims,flux_obj)

            simulate_psf(
                infile=evt_file,
                outroot="core_psf_sim",
                simulator="file",
                rayfile=rayfiles,
                projector="marx",
                blur=blur,
                pixadj="none",
                binsize=binsize,
                ra=coords.ra,
                dec=coords.dec,
            )
        else:
            spectrum_file:Optional[str]='core_flux_chart.dat'
            if flux_obj is not None:
                spectrum_file=None
            else:
                flux_obj=FluxObj(None,None)

            simulate_psf(
                infile=evt_file,
                outroot="core_psf_sim",
                monoenergy=flux_obj.monoenergy,
                flux=flux_obj.photon_flux,
                numiter=n_psf_sims,
                ra=coords.ra,
                spectrum=spectrum_file,
                dec=coords.dec,
                binsize=binsize,
                blur=blur,
                pixadj="none"
            )
        return "core_psf_sim_projrays.fits"
        
        




