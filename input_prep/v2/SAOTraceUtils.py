from collections import namedtuple

from typing import Any, Dict, List, NoReturn, Optional, Tuple, Union
import numpy as np
from scipy.special import factorial
from ciao_contrib.runtool import (
    dmstat,
    dmcoords,
    dmcopy,
    dmkeypar,
    dmlist
)
from sherpa.astro.ui import *
from sherpa_contrib.chart import *
import os
from math import floor
from TypeDefs import FluxObj, PointCel,PointPhys
import subprocess
from pathlib import Path
import multiprocessing as mp


class SAOTraceUtils():
    def __init__(self) -> None:
        pass

    @staticmethod
    def get_saotrace_keywords(evt_file:str,keys:Union[List,Tuple,np.array,None]=None,**kwargs)->Dict:
        """
        Get keywords from the events file that are required for simulating a PSF using SAOTrace

        :param evt_file: Path to the events file
        :param keys: The keys for which values need to be extracted. If None, ``"ROLL_PNT", "RA_PNT", "DEC_PNT", "ASOLFILE", "EXPOSURE"`` will be used
        ...
        :return: The dictionary containing key-value pairs
        :rtype: Dict
        """

        if keys is None:
            keys=["ROLL_PNT", "RA_PNT", "DEC_PNT", "ASOLFILE", "EXPOSURE"]

        out=dict()
        for key in keys:
            dmkeypar.punlearn()
            dmkeypar(evt_file, key=key)
            out[key] = dmkeypar.value

        dmlist.punlearn()
        data = dmlist(f"{evt_file}[GTI7][#row=1]", "data,raw")
        out["TSTART"] = data.splitlines()[1].split()[0]

        return out
    
    @staticmethod
    def create_saotrace_runconf_point_flux(keyvals:dict, coords:Union[PointCel,Tuple[float,float]], flux:float, monoenergy:float,**kwargs)->None:
        """
        Create a lua config script to simulate a point source with a specified flux at a specified energy

        :param keyvals: The keyvalue pairs that to configure the script, which can be obtained using ``get_saotrace_keywords`` function
        :param coords: ra,dec for the point
        :param flux: flux in ergs/cm2/s
        :param monoenergy: Energy for the given flux in keV
        """
        with open("saotrace_conf.lua", "w") as f:
            ra,dec=coords
            f.write(
                f"""
    ra_pnt={keyvals['RA_PNT']}
    dec_pnt={keyvals['DEC_PNT']}
    roll_pnt={keyvals['ROLL_PNT']}
    dither_asol{{
            file="{keyvals['ASOLFILE']}",
            ra={keyvals['RA_PNT']},
            dec={keyvals['DEC_PNT']},
            roll={keyvals['ROLL_PNT']}
        }}
    point{{
        position = {{ ra= tostring({ra}),
        dec = tostring({dec}),ra_aimpt=tostring(ra_pnt),dec_aimpt=tostring(dec_pnt) }},
        spectrum = {{ {{  {monoenergy},{flux}     }} }}
        }}
            """
            ) 

    @staticmethod
    def create_saotrace_runconf_spectrum(keyvals:dict, coords:Union[PointCel,Tuple[float,float]], spectrum_file,**kwargs)->None:
        """
        Similar to ``create_saotrace_runconf_point_flux`` but the point source will now have a specified spectrum

        :param keyvals: The keyvalue pairs that to configure the script, which can be obtained using ``get_saotrace_keywords`` function
        :param coords: ra,dec for the point
        :param spectrum_file: Path to the spectrum table. The columns are energy_low, energy_high, photon flux density
        """
        with open("saotrace_conf.lua", "w") as f:

            ra,dec=coords
            f.write(
                f"""
    ra_pnt={keyvals['RA_PNT']}
    dec_pnt={keyvals['DEC_PNT']}
    roll_pnt={keyvals['ROLL_PNT']}
    dither_asol{{
            file="{keyvals['ASOLFILE']}",
            ra={keyvals['RA_PNT']},
            dec={keyvals['DEC_PNT']},
            roll={keyvals['ROLL_PNT']}
        }}
    point{{
        position = {{ ra= tostring({ra}),
        dec = tostring({dec}),ra_aimpt=tostring(ra_pnt),dec_aimpt=tostring(dec_pnt) }},
        spectrum = {{ {{ type = 'file',
        format = 'rdb',
        file = "{spectrum_file}",
        units = 'photons/s/cm2',
        emin='elo',emax='ehi',flux='spectrum'
        }} }}
        }}
            """
            )

    @staticmethod
    def save_spectrum_rdb(file:str,**kwargs)->None:
        """
        Save the SHERPA generated spectum in rdb format, which is compatible with SAOTrace

        :param file: path to the input spectrum 
        """
        outfile = "core_flux_saotrace.rdb"

        with open(outfile, "w") as f:
            print("elo\tehi\tspectrum", file=f)
            print("N\tN\tN", file=f)
            with open(file, "r") as f2:
                rows = f2.readlines()[3:]
                for row in rows:
                    print("\t".join(row.strip().split(" ")), file=f)

    @staticmethod
    def generate_rayfile_saotrace(arg)->str:
        """
        Generate a ray file using SAOTrace with the specified parameters. Needs SAOTrace to be initialized before running this function

        :param arg: Should contain in order, the path to the directory where the rayfiles should be written to, a dictionary containing the start time and the exposure for the observation, t  he iteration number for the simulation
        ...
        :return: Path to the raytrace file
        """
        tempdir,saotrace_keyvals, iter_no = arg
        output = subprocess.check_output(
        [
            "trace-nest",
            f"tag={tempdir}/HRMS_rayfile_{iter}",
            "srcpars=saotrace_conf.lua",
            f'tstart={saotrace_keyvals["TSTART"]}',
            f'limit={saotrace_keyvals["EXPOSURE"]}',
            "limit_type=sec",
            f"seed1={int(np.random.random()*100000)}",
        ]
        )
        # logger.info(output.decode("utf-8") )
        return f"{tempdir}/HRMS_rayfile_{iter_no}.fits"


    @staticmethod
    def setup_saotrace_sim(evt_file:str, coords:PointCel, spectrum_file:str,flux_obj:Optional[FluxObj]=None)->Dict:

        keyvals = SAOTraceUtils.get_saotrace_keywords(evt_file)
        
        if spectrum_file is not None and flux_obj is None:
            SAOTraceUtils.create_saotrace_runconf_spectrum(keyvals, coords, spectrum_file)
        elif flux_obj is not None:
            SAOTraceUtils.create_saotrace_runconf_point_flux(keyvals,coords,flux_obj.photon_flux,flux_obj.monoenergy)
        elif spectrum_file is None and flux_obj is None:
            raise ValueError('Either a spectrum file or a flux must be specified')
        else:
            ...
        return keyvals 

    @staticmethod
    def run_sao_raytrace_parallel(evt_file:str,coords:PointCel,n_psf_sims:int,flux_obj:Optional[FluxObj]=None,tempdir:str='raysdir',spectrum_file:str='core_flux_saotrace.rdb')->List:
        """
        Simulate rayfiles for a point source using SAOTrace in parallel. A spectrum file or a flux can be specified
        :param evt_file: Path to the events file
        :param coords: WCS coordinates of the point source to be simulated
        :param n_psf_sims: Number of rayfiles to generate
        :param flux_obj: Set this to use a representative flux instead of a spectrum
        :param tempdir: Directory for the rayfiles
        :param spectrum_file: SED for the point source in rdb format
        ...
        :return: Returns a list of Paths to the simuated rayfiles
        """
        

        pool = mp.Pool(mp.cpu_count())

        keyvals = SAOTraceUtils.setup_saotrace_sim(evt_file,coords, spectrum_file,flux_obj=flux_obj)
        Path(tempdir).mkdir(exist_ok=True)

        results = pool.map(
            SAOTraceUtils.generate_rayfile_saotrace, [(tempdir, keyvals, iter) for iter in range(n_psf_sims)]
        )

        pool.close()

        return results
