#!/usr/bin/env python


from prepare_spectra import prepare_spectra
import yaml
from rebin_img import rebin_img
from rebin_psf import rebin_psf
from ciao_contrib.runtool import *
from sherpa.astro.ui import *
from sherpa_contrib.chart import *

# from sim_image import sim_image
from saotrace_helpers import (
    run_sao_raytrace,
    run_sao_raytrace_parallel,
    run_sao_raytrace_parallel_flux,
)
import logging
import numpy as np
from datetime import datetime
from pathlib import Path
import matplotlib.pyplot as plt
import shutil

# initialize the logger
logger = logging.getLogger("sherpa")
logging.basicConfig(
    level=logging.INFO, filename="sherpa_session_for_lira.log", filemode="a"
)
logger.info("\n\n\n")
logger.info(datetime.now().strftime("%c"))
logger.info("\n\n\n")


def prepare_the_params(params):
    # parse the params and return the
    param_def = get_param_definition()

    # go through all the params and edit the if a user supplies it
    for key, defn in param_def.items():
        if (defn["required"]) and not key in params:
            raise Exception("{0} is required".format(key))
        if defn["required"] and defn["type"] != type(params[key]).__name__:
            raise Exception("Incorrect value for {0}".format(key))
        if key in params:
            param_def[key]["value"] = params[key]

    return param_def


def get_param_definition():
    param_list = [
        "evt_file",
        "binsize",
        "core_reg",
        "bkg_reg",
        "n_psf_sims",
        "n_null_sims",
        "inp_size",
        "psf_size",
        "nH",
        "add_gal",
        "redshift",
        "group",
        "blur",
        "center",
        "no_core",
        "sim_baselines",
        "extract_spectra",
        "fit_spectra",
    ]
    param_types = [
        "str",
        "float",
        "str",
        "str",
        "int",
        "int",
        "int",
        "int",
        "float",
        "int",
        "float",
        "int",
        "float",
        "list",
        "bool",
        "bool",
        "bool",
        "bool",
    ]
    param_req = [
        True,
        False,
        True,
        False,
        True,
        False,
        False,
        True,
        False,
        False,
        False,
        False,
        True,
        False,
        True,
        False,
        False,
        False,
    ]
    default_params = [
        None,
        0.5,
        None,
        None,
        50,
        50,
        64,
        32,
        None,
        0.0,
        0.0,
        10,
        0.25,
        None,
        False,
        True,
        True,
        True,
    ]

    param_definition = {}
    for i in range(0, len(param_list)):
        param_definition[param_list[i]] = {
            "type": param_types[i],
            "required": param_req[i],
            "value": default_params[i],
        }

    return param_definition


def create_the_inputs(params):

    # create the image
    image_file, core_cen, core_cen_adj = create_image(
        params["evt_file"]["value"],
        params["core_reg"]["value"],
        params["binsize"]["value"],
        params["inp_size"]["value"],
        params["center"]["value"],
    )

    print("Image center: {0}".format(core_cen))
    print("Image center-adjusted: {0}".format(core_cen_adj))
    print(f"Image_file_name {image_file}")
    # raise()
    max_ra_dec = get_max_ra_dec(image_file, params["core_reg"]["value"])
    # simulate the psf
    simulated_psf_file = sim_core_psf(
        params["evt_file"]["value"],
        params["n_psf_sims"]["value"],
        params["core_reg"]["value"],
        params["bkg_reg"]["value"],
        params["binsize"]["value"],
        min(params["binsize"]["value"], params["psf_size"]["value"]),
        max_ra_dec[0],
        max_ra_dec[1],
        image_file,
        params["nH"]["value"],
        params["add_gal"]["value"],
        params["redshift"]["value"],
        params["group"]["value"],
        params["blur"]["value"],
        params["extract_spectra"]["value"],
        params["fit_spectra"]["value"],
    )

    psf_file, psf_cen, psf_cen_adj = create_image(
        simulated_psf_file,
        params["core_reg"]["value"],
        params["binsize"]["value"],
        params["psf_size"]["value"],
        None,
        outfile="core_psf.fits",
        psf=True,
    )

    # print creating ecf profiles and plotting them
    # print('Creating ecf profiles')
    # compare_ecf_profiles(
    #    params['evt_file']['value'],
    #   simulated_psf_file,
    #   psf_cen_adj[0],psf_cen_adj[1]
    # )
    print("PSF centroid-adjusted: {0}".format(psf_cen_adj))
    print(f"psf_cntr_adj {psf_cen_adj[0]}_{psf_cen_adj[1]}")

    # create the baseline image and simulate images from it
    if params["sim_baselines"]["value"]:
        simulate_null_images(
            image_file,
            "core_psf.fits",
            params["n_null_sims"]["value"],
            params["no_core"]["value"],
        )


def compare_ecf_profiles(
    evt_file,
    psf_file,
    xcen,
    ycen,
    radius=20,
    binsize=0.2,
    ecf_fraction=[
        0.01,
        0.025,
        0.05,
        0.075,
        0.1,
        0.2,
        0.3,
        0.4,
        0.5,
        0.6,
        0.7,
        0.8,
        0.85,
        0.9,
        0.95,
        0.99,
    ],
):
    ecf_calc.punlearn()
    ecf_calc(
        infile=evt_file,
        outfile="core_ecf.fits",
        xpos=xcen,
        ypos=ycen,
        radius=radius,
        binsize=binsize,
        clobber=True,
        fraction=ecf_fraction,
    )
    ecf_calc(
        infile=psf_file,
        outfile="psf_ecf.fits",
        xpos=xcen,
        ypos=ycen,
        radius=radius,
        binsize=binsize,
        clobber=True,
        fraction=ecf_fraction,
    )

    # load and plot them
    load_data(1, "core_ecf.fits", colkeys=["r_mid", "frac_in"])
    load_data(2, "psf_ecf.fits", colkeys=["r_mid", "frac_in"])

    plot_data(1, yerrorbars=False, linestyle="solid", xlog=True, ylog=True)
    plot_data(2, yerrorbars=False, linestyle="solid", overplot=True)

    plt.savefig("ecf_profiles.eps")


def get_max_ra_dec(img_file, reg_file):
    dmstat.punlearn()
    dmstat("{0}[sky=region({1})]".format(img_file, reg_file))

    vals = [float(x) for x in dmstat.out_max_loc.split(",")]
    xval = vals[0]
    yval = vals[1]

    dmcoords.punlearn()
    dmcoords(img_file, option="sky", x=vals[0], y=vals[1], celfmt="deg")
    return (dmcoords.ra, dmcoords.dec)


def get_centroid_physical(evt_file, reg_file):
    dmstat.punlearn()
    dmstat("{0}[sky=region({1})][cols sky]".format(evt_file, reg_file))

    vals = [float(x) for x in dmstat.out_mean.split(",")]
    xval = vals[0]
    yval = vals[1]

    return (xval, yval)


def get_centroid_ra_dec(evt_file, reg_file):
    dmcoords.punlearn()
    cntrd_phys = get_centroid_physical(evt_file, reg_file)
    dmcoords(evt_file, option="sky", x=cntrd_phys[0], y=cntrd_phys[1], celfmt="deg")
    return (dmcoords.ra, dmcoords.dec)


def create_image(evt_file, reg_file, binsize, nsize, center, outfile=None, psf=False):

    if center is None or len(center) == 0:
        dmstat.punlearn()
        # get the centroid of the image
        dmstat("{0}[sky=region({1})][cols sky]".format(evt_file, reg_file), clip=True)
        vals = [float(x) for x in dmstat.out_mean.split(",")]
        xval = vals[0]
        yval = vals[1]
    else:
        xval, yval = center

    if outfile is None:
        outfile = "img_{0}x{1}_{2}.fits".format(nsize, nsize, binsize)

    print("Creating the input image file")
    if not psf:
        centroid = rebin_img(
            infile="{0}[energy=500:7000]".format(evt_file),
            outfile=outfile,
            binsize=binsize,
            nsize=nsize,
            xcen=xval,
            ycen=yval,
        )
    else:
        centroid = rebin_psf(
            infile="{0}[energy=500:7000]".format(evt_file),
            outfile=outfile,
            binsize=binsize,
            nsize=nsize,
            xcen=xval,
            ycen=yval,
        )

        # test if the maximum pixel is at the center. If not, rebin it again
        dmstat.punlearn()
        # get the centroid of the image
        dmstat(f"{outfile}[sky=region({reg_file})]", clip=True)
        vals = [float(x) for x in dmstat.out_max_loc.split(",")]
        if not np.all(centroid == vals):
            print(f"The maximum of the PSF is not at the center")
            print(f"Max location: {vals}")
            print(f"Center of the image: {centroid}")
            print("Re-binning...")
            centroid = rebin_psf(
                infile="{0}[energy=500:7000]".format(evt_file),
                outfile=outfile,
                binsize=binsize,
                nsize=nsize,
                xcen=vals[0],
                ycen=vals[1],
            )

    return (outfile, center or vals, centroid)


def sim_core_psf(
    evt_file,
    npsf_sims,
    core_reg,
    bkg_reg,
    binsize,
    psf_size,
    xcen,
    ycen,
    image_file,
    nH=None,
    add_gal=0,
    redshift=0.0,
    group=10.0,
    blur=0.25,
    extract_spectra=True,
    fit_spectra=True,
):

    # test for the number of counts. If the number of counts are less than 40,
    # don't use a spectrum file. Simply specify the flux and the monoenergy.
    dmstat.punlearn()
    dmstat(f"{image_file}[sky=region({core_reg})]", centroid=False)
    use_flux = False
    photon_flux = None
    monoenergy = None
    if float(dmstat.out_sum) < 40:
        print(
            f"The counts in the core {dmstat.out_sum} are too small for spectral fitting"
        )
        print(f"Flux at a monoenergy will be used for the PSF simulation")
        use_flux = True

    # extract the spectrum
    if extract_spectra:
        specextract.punlearn()
        print("Extracting the spectrum")

        specextract(
            infile="{0}[sky=region({1})]".format(evt_file, core_reg),
            outroot="core_spectrum",
            bkgfile="{0}[sky=region({1})]".format(evt_file, bkg_reg),
            clobber=True,
        )

    ##fit the spectrum
    if fit_spectra and not use_flux:
        prepare_spectra(group, nH, add_gal, redshift)

    # get the ra dec from sky coords
    ra_dec = get_centroid_ra_dec(evt_file, core_reg)

    if use_flux:
        # calculate the monoenergy
        # Using case 2 in https://cxc.cfa.harvard.edu/ciao/why/monochromatic_energy.html
        dmtcalc.punlearn()
        dmtcalc(
            "core_spectrum.arf",
            "arf_weights",
            "mid_energy=(energ_lo+energ_hi)/2.0;weights=(mid_energy*specresp)",
            clobber=True,
        )

        dmstat.punlearn()
        dmstat("arf_weights[mid_energy=2.0:7.0][cols weights,specresp]", verbose="0")
        out_sum = dmstat.out_sum
        out_sum = [float(val) for val in out_sum.split(",")]

        monoenergy = out_sum[0] / out_sum[1]

        srcflux.punlearn()
        srcflux(
            evt_file,
            ",".join(ra_dec),
            "flux",
            bands=f"0.5:7.0:{monoenergy}",
            psfmethod="quick",
            verbose="0",
            clobber=True,
        )

        dmkeypar.punlearn()
        dmkeypar("flux_0.5-7.0.flux", "net_photflux_aper")
        photon_flux = dmkeypar.rval

        dmkeypar.punlearn()
        dmkeypar("flux_0.5-7.0.flux", "net_rate")
        print(f"Net count rate: {dmkeypar.rval}")

        print(
            f"Monoenergy: {monoenergy} keV, photon flux: {photon_flux} photons/s/cm^2"
        )
        scale_fac = input("Scaling for photon flux (1): ")
        if scale_fac == "":
            scale_fac = 1
        photon_flux *= float(scale_fac)

    # simulate the psf
    print("Simulating the psf")
    saotrace = True
    simulate_psf.punlearn()
    if saotrace:
        print("Raytracing using SAOTrace")
        if use_flux:
            rayfiles = run_sao_raytrace_parallel_flux(
                evt_file, ra_dec[0], ra_dec[1], npsf_sims, photon_flux, monoenergy
            )
        else:
            rayfiles = run_sao_raytrace_parallel(
                evt_file, ra_dec[0], ra_dec[1], npsf_sims
            )
        print("Projecting on to the detector")
        simulate_psf(
            infile=evt_file,
            outroot="core_psf_sim",
            simulator="file",
            rayfile=rayfiles,
            projector="marx",
            blur=blur,
            pixadj="none",
            binsize=binsize,
            ra=ra_dec[0],
            dec=ra_dec[1],
        )
    else:
        if use_flux:
            simulate_psf(
                infile=evt_file,
                outroot="core_psf_sim",
                monoenergy=monoenergy,
                flux=photon_flux,
                numiter=npsf_sims,
                ra=ra_dec[0],
                dec=ra_dec[1],
                binsize=binsize,
                blur=blur,
                pixadj="none"
                # ,minsize=psf_size
            )
        else:
            simulate_psf(
                infile=evt_file,
                outroot="core_psf_sim",
                spectrum="core_flux_chart.dat",
                numiter=npsf_sims,
                ra=ra_dec[0],
                dec=ra_dec[1],
                binsize=binsize,
                blur=blur,
                pixadj="none"
                # ,minsize=psf_size
            )

    save_all(outfile="lira_input_psfsim.log", clobber=True)
    return "core_psf_sim_projrays.fits"


# slightly modified version of kmc's code from astrostat/LIRA
def simulate_null_images(infile, psffile, num_sims, no_core, mcmciter=5000):
    print("Creating the null file")
    clean()
    set_stat("cstat")
    set_method("simplex")
    load_image(infile)
    load_psf("mypsf", psffile)
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
        for i in range(num_sims):
            fake()
            save_image("sim_null_{}.fits".format(i), clobber=True)
        clean()
        return 0

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
        for i in range(num_sims):
            fake()
            save_image("sim_null_{}.fits".format(i), clobber=True)
        clean()
        return 0
    # if not go for the regular
    set_prior(q1.fwhm, g1)
    set_sampler_opt("defaultprior", False)
    set_sampler_opt("priorshape", [True, False, False, False, False])
    set_sampler_opt("originalscale", [True, True, True, True, True])
    if mcmciter < num_sims * 100:
        mcmciter = num_sims * 100

    # the following code throws an error sometimes #bug
    try:
        stats, accept, params = get_draws(1, niter=mcmciter)
    except:
        params = [np.repeat(q1.fwhm.val, mcmciter)]

    # print('Simulating the null files')
    for i in range(num_sims):
        set_par(q1.fwhm, params[0][(i + 1) * 100 - 1])
        fake()
        save_image("sim_null_{}.fits".format(i), clobber=True)
    save_all(outfile="lira_input_baseline_sim.log", clobber=True)
    clean()


with open("lira_input_prep.yaml", "r") as stream:
    try:
        params = create_the_inputs(prepare_the_params(yaml.safe_load(stream)))
        shutil.rmtree("raysdir")

    except yaml.YAMLError as exc:
        print(exc)
