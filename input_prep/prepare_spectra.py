from sherpa.astro.io import read_pha
from sherpa.optmethods import MonCar
from sherpa.stats import WStat,CStat
from sherpa.astro.instrument import RSPModelPHA
from sherpa.plot import ModelPlot
from sherpa.fit import Fit
from cstat_gof import gof_cstat
from sherpa_contrib.chart import  save_chart_spectrum
from scipy import special
import numpy as np
from cstat_gof import gof_cstat
from sherpa.astro.ui import *


def prepare_spectra(group,nH,add_gal):
    pha = read_pha('core_spectrum.pi')
    pha.set_analysis('energy')
    pha.notice(0.5, 7.0)
    tabs = ~pha.mask
    pha.group_counts(group, tabStops=tabs)
    x = pha.get_x()
    x = pha.apply_filter(x, pha._middle)
    y = pha.get_y(filter=True)
    pha.set_analysis('energy')

    model = xsphabs.abs1*powlaw1d.srcp1
    print('Fitting the spectrum')

    zFlag=False
    if (nH is not None) and (nH > 0.0):
        if(add_gal==1):
            model = xsphabs.gal*xszphabs.abs1*powlaw1d.srcp
            gal.nH=nH
            freeze(gal.nH)
            zFlag=True

        else:
            model = xsphabs.abs1*powlaw1d.srcp1
            abs1.nH=nH
            freeze(abs1.nH)
    else:
        model = xszphabs.abs1*powlaw1d.srcp1
        zFlag=True


    if zFlag is True and add_gal==1:
        #print('REDSHIFT',redshift)
        abs1.redshift=redshift
        freeze(abs1.redshift)

    full_model = RSPModelPHA(pha.get_arf(), pha.get_rmf(), pha, pha.exposure * model)
    

    fit = Fit(pha, full_model,method=MonCar(),stat=WStat())
    res = fit.fit()

    print(res.format())
    print(fit.est_errors())


    #calculate the p-value for wstat
    mplot2 = ModelPlot()
    mplot2.prepare(pha, full_model)

    miu=mplot2.y * pha.exposure * 0.0146
    obs=y * pha.exposure * 0.0146

    c,ce,cv = gof_cstat(miu,obs)

    print(f'C0={c},C_e={ce},C_v={cv}')

    zval=(fit.calc_stat()-ce)/np.sqrt(cv)


    if zval>0:
        pval = special.erfc(zval/np.sqrt(2))
    else:
        pval = special.erf(abs(zval)/np.sqrt(2))
        
    print(f'p-value for wstat = {pval}')


    set_data(pha)
    set_model(model)
    save_chart_spectrum("core_flux_chart.dat", elow=0.5, ehigh=7.0)