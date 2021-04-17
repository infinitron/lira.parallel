#!/usr/bin/env python
import sys

if len(sys.argv) != 6:
    print(
        "Usage: blur_estimator.py <evt_file> <region_file> <blurs> <n_psf_sims> <radius>"
    )
    exit(-1)

from ciao_contrib.runtool import *
from sherpa.astro.ui import *
from sherpa_contrib.chart import *

import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from pathlib import Path
import similaritymeasures


def get_best_blur(
    evt_file,
    numiter,
    region_file,
    blur_vals,
    radius_100p,
    binsize=0.2,
    spectrum="core_flux_chart.dat",
    out_root_common="core_psf",
):

    xcen, ycen, ra, dec = get_centroids(evt_file, region_file)
    print(xcen, ycen)

    ecf_fractions = [
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
    ]
    n_fractions = len(ecf_fractions)
    ecf_profile_values = np.zeros((blur_vals.shape[0], n_fractions), dtype=float)
    ecf_radii_values = np.zeros((blur_vals.shape[0], n_fractions), dtype=float)

    # generate the profile of the core
    obs_ecf_file = "obs_ecf.fits"
    ecf_calc.punlearn()
    ecf_calc(
        infile=evt_file,
        outfile="obs_ecf.fits",
        xpos=xcen,
        ypos=ycen,
        radius=radius_100p,
        binsize=binsize,
        clobber=True,
        fraction=ecf_fractions,
    )

    for idx, blur in enumerate(blur_vals):
        print(f"Blur={blur}")
        outroot = out_root_common + "_blur_" + str(blur)
        if not Path(outroot + "_projrays.fits").is_file():
            print(f"Simulating PSF")

            simulate_psf(
                infile=evt_file,
                outroot=outroot,
                spectrum=spectrum,
                numiter=numiter,
                ra=ra,
                dec=dec,
                binsize=binsize,
                blur=blur,
            )

            print("Generating the ecf profile")
        out_ecf = generate_ecf_profile(
            outroot + "_projrays.fits",
            blur,
            xcen,
            ycen,
            binsize=binsize,
            ecf_fractions=ecf_fractions,
            radius=radius_100p,
        )
        load_table(idx, out_ecf, colkeys=["r_mid", "fraction"])
        ecf_profile = get_data(idx)
        ecf_profile_values[idx, :] = ecf_profile.y
        ecf_radii_values[idx, :] = ecf_profile.x
        plt.plot(ecf_profile.x, ecf_profile.y, label=f"blur={blur}")
        print("Done")

    load_table(idx + 1, obs_ecf_file, colkeys=["r_mid", "fraction"])
    obs_ecf = get_data(idx + 1)
    plt.plot(obs_ecf.x, obs_ecf.y, label=f"Observation", ls="--", lw=2, color="black")

    plt.legend()
    plt.xlabel("Radius (ACIS pixels)")
    plt.ylabel("ECF")
    plt.show()

    dtw_vals = np.zeros(idx + 1)
    pcm_vals = np.zeros(idx + 1)
    clm_vals = pcm_vals.copy()
    abtc_vals = pcm_vals.copy()
    for i in range(idx + 1):
        n1 = np.zeros((n_fractions, 2))
        n2 = n1.copy()
        n1[:, 0] = ecf_radii_values[i, :]
        n2[:, 0] = obs_ecf.x
        n1[:, 1] = ecf_profile_values[i, :]
        n2[:, 1] = obs_ecf.y
        dtw_vals[i], _ = similaritymeasures.dtw(n1, n2)
        pcm_vals[i] = similaritymeasures.pcm(n1, n2)
        clm_vals[i] = similaritymeasures.curve_length_measure(n1, n2)
        abtc_vals[i] = similaritymeasures.area_between_two_curves(n1, n2)

    print(f"Blur (Dynamic time warping) {blur_vals[np.argmin(dtw_vals)]}")
    print(f"Blur (Partial Curve Mapping): {blur_vals[np.argmin(pcm_vals)]}")
    print(f"Blur (Curve Length Measure): {blur_vals[np.argmin(clm_vals)]}")
    print(f"Blur (Area Curve Measure): {blur_vals[np.argmin(abtc_vals)]}")


def get_centroids(evt_file, reg_file):
    dmstat.punlearn()
    dmstat("{0}[sky=region({1})][cols sky]".format(evt_file, reg_file))

    vals = [float(x) for x in dmstat.out_mean.split(",")]
    xval = vals[0]
    yval = vals[1]

    dmcoords.punlearn()
    dmcoords(
        f"{evt_file}[sky=region({reg_file})][cols sky]",
        option="sky",
        x=xval,
        y=yval,
        celfmt="deg",
    )

    ra, dec = dmcoords.ra, dmcoords.dec

    return (xval, yval, ra, dec)


def generate_ecf_profile(
    evt_file,
    blur,
    xcen,
    ycen,
    radius=8.2,
    binsize=0.2,
    ecf_fractions=[
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
    out_file = "core_ecf_blur_" + str(blur) + ".fits"
    if Path(out_file).is_file():
        return out_file
    ecf_calc(
        infile=evt_file,
        outfile=out_file,
        xpos=xcen,
        ypos=ycen,
        radius=radius,
        binsize=binsize,
        clobber=True,
        fraction=ecf_fractions,
    )

    return out_file


evt_file = sys.argv[1]
region_file = sys.argv[2]
blurs = sys.argv[3]
n_psf_sims = sys.argv[4]
radius_100p = sys.argv[5]

blurs = np.asfarray(blurs.split(","), dtype=float)

get_best_blur(evt_file, n_psf_sims, region_file, blurs, radius_100p)
