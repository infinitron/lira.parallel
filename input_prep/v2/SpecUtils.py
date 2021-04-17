from collections import namedtuple
from typing import Any, Dict, List, NoReturn, Optional, Tuple, Union
from CoordUtils import CoordUtils
from ciao_contrib.runtool import specextract,dmtcalc,srcflux
import numpy as np
from scipy.special import factorial
from ciao_contrib.runtool import (
    dmstat,
    dmcoords,
    dmcopy,
    dmkeypar,
    dmlist
)
from SAOTraceUtils import SAOTraceUtils
from sherpa.astro.io import (
    read_pha,
)
from sherpa_contrib.chart import save_chart_spectrum
from sherpa.astro.ui import (
    xsphabs,
    powlaw1d,
    freeze,
    set_data,
    set_model,
    
    xszphabs
)
from scipy import special

from sherpa.optmethods import MonCar
from sherpa.stats import WStat, CStat
from sherpa.astro.instrument import RSPModelPHA
from sherpa.plot import ModelPlot
from sherpa.fit import Fit
from sherpa.astro.ui import *
from sherpa_contrib.chart import *
import os
from math import floor
from TypeDefs import FluxObj

class SpecUtils():
    def __init__(self) -> None:
        pass

    @staticmethod
    def estimate_gof_cstat(miu:Union[np.array],obs:Union[np.array],**kwargs)->Tuple[float,float,float]:
        """
        Estimate the goodness of fit for cstat statistic. Adapted from https://github.com/sherpa/sherpa/issues/594 
        
        :param miu: Model SED values
        :param obs: Observed SED values
        ...
        :return: Return the Co, Ce, Cv values that can be used to quantify the gof. See Kaastra 2017 (A&A 605, A51) for more details
        :rtype: Tuple[float,float,float]
        """
        ce = np.zeros(np.size(miu))
        cv = np.zeros(np.size(miu))

        lnnm = np.log(obs / miu)
        where_are_nan = np.isnan(lnnm)
        where_are_inf = np.isinf(lnnm)
        lnnm[where_are_nan] = 0
        lnnm[where_are_inf] = 0
        co = miu - obs + obs * lnnm
        Co = 2 * sum(co)

        for i in range(np.size(miu)):
            if miu[i] == 0:
                ce[i] = -0.25 * miu[i] ** 3 + 1.38 * miu[i] ** 2
            elif 0 < miu[i] and miu[i] <= 0.5:
                ce[i] = (
                    -0.25 * miu[i] ** 3 + 1.38 * miu[i] ** 2 - 2 * miu[i] * np.log(miu[i])
                )
            elif 0.5 < miu[i] and miu[i] <= 2:
                ce[i] = (
                    -0.00335 * miu[i] ** 5
                    + 0.04259 * miu[i] ** 4
                    - 0.27331 * miu[i] ** 3
                    + 1.381 * miu[i] ** 2
                    - 2 * miu[i] * np.log(miu[i])
                )
            elif 2 < miu[i] and miu[i] <= 5:
                ce[i] = 1.019275 + 0.1345 * miu[i] ** (0.461 - 0.9 * np.log(miu[i]))
            elif 5 < miu[i] and miu[i] <= 10:
                ce[i] = 1.00624 + 0.604 / miu[i] ** 1.68
            else:
                ce[i] = 1 + 0.1649 / miu[i] + 0.226 / (miu[i] ** 2)

        where_are_nan = np.isnan(ce)
        where_are_inf = np.isinf(ce)
        ce[where_are_nan] = 0
        ce[where_are_inf] = 0

        k = np.arange(1, 5)
        for i in range(np.size(miu)):
            if 0 <= miu[i] and miu[i] <= 0.1:
                cv[i] = (
                    4
                    * sum(
                        np.exp(-miu[i])
                        * miu[i] ** k
                        / factorial(k)
                        * (miu[i] - k + k * np.log(k / miu[i])) ** 2
                    )
                    + 4 * (np.exp(-miu[i]) * miu[i] ** 0 / factorial(0) * (miu[i]) ** 2)
                    - ce[i] ** 2
                )
            elif 0.1 < miu[i] and miu[i] <= 0.2:
                cv[i] = (
                    -262 * miu[i] ** 4
                    + 195 * miu[i] ** 3
                    - 51.24 * miu[i] ** 2
                    + 4.34 * miu[i]
                    + 0.77005
                )
            elif 0.2 < miu[i] and miu[i] <= 0.3:
                cv[i] = 4.23 * miu[i] ** 2 - 2.8254 * miu[i] + 1.12522
            elif 0.3 < miu[i] and miu[i] <= 0.5:
                cv[i] = -3.7 * miu[i] ** 3 + 7.328 * miu[i] ** 2 - 3.6926 * miu[i] + 1.20641
            elif 0.5 < miu[i] and miu[i] <= 1:
                cv[i] = (
                    1.28 * miu[i] ** 4
                    - 5.191 * miu[i] ** 3
                    + 7.666 * miu[i] ** 2
                    - 3.5446 * miu[i]
                    + 1.15431
                )
            elif 1 < miu[i] and miu[i] <= 2:
                cv[i] = (
                    0.1125 * miu[i] ** 4
                    - 0.641 * miu[i] ** 3
                    + 0.859 * miu[i] ** 2
                    + 1.0914 * miu[i]
                    - 0.05748
                )
            elif 2 < miu[i] and miu[i] <= 3:
                cv[i] = (
                    0.089 * miu[i] ** 3 - 0.872 * miu[i] ** 2 + 2.8422 * miu[i] - 0.67539
                )
            elif 3 < miu[i] and miu[i] <= 5:
                cv[i] = 2.12336 + 0.012202 * miu[i] ** (5.717 - 2.6 * np.log(miu[i]))
            elif 5 < miu[i] and miu[i] <= 10:
                cv[i] = 2.05159 + 0.331 * miu[i] ** (1.343 - np.log(miu[i]))
            else:
                cv[i] = 12 / (miu[i] ** 3) + 0.79 / (miu[i] ** 2) + 0.6747 / miu[i] + 2

        where_are_nan = np.isnan(cv)
        where_are_inf = np.isinf(cv)
        cv[where_are_nan] = 0
        cv[where_are_inf] = 0
        Ce = sum(ce)
        Cv = sum(cv)

        return (Co, Ce, Cv)

    @staticmethod
    def get_gof_cstat_sherpa(data:Any,model:Any,exposure:float,**kwargs)->Tuple[float,float,float]:
        """
        A wrapper for the ```estimate_gof_cstat`` function to use SHERPA data

        :param data: Observed SHERPA data object
        :param mdel: Model SHERPA data object
        :param exposure: The total exposure of the observation
        ...
        :return: Same as ``estimate_gof_cstat``
        :rtype: Tuple[float,float,float]
        """

        dat_y = data.y

        xlo = model.xlo
        xhi = model.xhi

        ymod = model.y

        miu = ymod * exposure * (xhi - xlo)
        obs = dat_y * exposure * (xhi - xlo)

        return SpecUtils.estimate_gof_cstat(miu, obs)

    @staticmethod
    def prepare_spectra(nH:float,group:int=1,add_gal:bool=False,redshift:Optional[float]=None,**kwargs)->float:
        """
        Fit the spectra using an absorbed powerlaw model using the Wstat statistic. The function also returns a p-value for the gof.
        :param nH: The galactic absorption column density in units of 10^22 /cm3
        :param group: The number of counts per energy bin
        :param add_gal: Setting this to True would add an intrinsic abrosption column density along side the galactic one
        :param redshift: The redshift to use in the fit. Only takes effect if add_gal is set to True
        ...
        :return: Returns the p-value of the gof. The null hypothesis states that the model and the observation differ while alternate says that the model explains the data
        """

        pha = read_pha("core_spectrum.pi")
        pha.set_analysis("energy")
        pha.notice(0.5, 7.0)
        tabs = ~pha.mask
        pha.group_counts(group, tabStops=tabs)
        x = pha.get_x()
        x = pha.apply_filter(x, pha._middle)
        y = pha.get_y(filter=True)
        pha.set_analysis("energy")

        model = xsphabs.abs1 * powlaw1d.srcp1
        print("Fitting the spectrum")

        zFlag = False
        if (nH is not None) and (nH > 0.0):
            if add_gal == 1:
                model = xsphabs.gal * xszphabs.abs1 * powlaw1d.srcp
                gal.nH = nH
                freeze(gal.nH)
                zFlag = True

            else:
                model = xsphabs.abs1 * powlaw1d.srcp1
                abs1.nH = nH
                freeze(abs1.nH)
        else:
            model = xszphabs.abs1 * powlaw1d.srcp1
            zFlag = True

        if zFlag is True and add_gal == 1:
            # print('REDSHIFT',redshift)
            abs1.redshift = redshift
            freeze(abs1.redshift)

        full_model = RSPModelPHA(pha.get_arf(), pha.get_rmf(), pha, pha.exposure * model)

        print(full_model)

        fit = Fit(pha, full_model, method=MonCar(), stat=WStat())
        res = fit.fit()

        print(res.format())
        print(fit.est_errors())

        # calculate the p-value for wstat
        mplot2 = ModelPlot()
        mplot2.prepare(pha, full_model)

        miu = mplot2.y * pha.exposure * 0.0146
        obs = y * pha.exposure * 0.0146

        c, ce, cv = SpecUtils.estimate_gof_cstat(miu, obs)

        #print(f"C0={c},C_e={ce},C_v={cv}")

        zval = (fit.calc_stat() - ce) / np.sqrt(cv)

        if zval > 0:
            pval = special.erfc(zval / np.sqrt(2))
        else:
            pval = special.erf(abs(zval) / np.sqrt(2))

        print(f"p-value for wstat = {pval}")

        set_data(pha)
        set_model(model)
        save_chart_spectrum("core_flux_chart.dat", elow=0.5, ehigh=7.0)
        # save_chart_spectrum("core_flux_chart.rdb",format='text/tsv', elow=0.5, ehigh=7.0)
        SAOTraceUtils.save_spectrum_rdb("core_flux_chart.dat")

        return pval

    @staticmethod
    def extract_and_fit_spectra(evt_file:str,img_file:str,reg_file:str,bkg_file:str,extract_spectra:bool=True,fit_spectra:bool=True,**kwargs)->Optional[FluxObj]:
        """
        Extract the spectrum from a specified region and fit it
        :param evt_file: Path to the events file
        :param img_file: Path to the image file
        :param reg_file: Path to the region file that will be used to extract the spectrum
        :param bkg_file: Path to the background region file
        :param extract_spectra: Set this to True to extract the spectrum
        :param fit_spectra: Set this to fit and save the model spectrum
        ...
        :return: If the number of counts inside is greater than 40, an absorbed powerlaw will be fit. If not, a representative flux and energy will be chosen. The flux can be interactively  scaled run time
        :rtype: Tuple[float,float] or None 
        """
        dmstat.punlearn()
        dmstat(f"{img_file}[sky=region({reg_file})]", centroid=False)
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

            specextract(
                infile="{0}[sky=region({1})]".format(evt_file, reg_file),
                outroot="core_spectrum",
                bkgfile="{0}[sky=region({1})]".format(evt_file, bkg_file),
                clobber=True,
            )

        if fit_spectra and not use_flux:
            SpecUtils.prepare_spectra(**kwargs)

            return None


        # get the ra dec from sky coords
        ra_dec = CoordUtils.get_centroid(evt_file, reg_file)

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
            scale_fac : Union[str,float] = input("Scaling for photon flux (1): ")
            if scale_fac == "":
                scale_fac = 1
            photon_flux *= float(scale_fac)
            return FluxObj(photon_flux,monoenergy)
        return None

