import subprocess
import logging
logger = logging.getLogger("sherpa")
import numpy as np
from pathlib import Path 
from ciao_contrib.runtool import *
import multiprocessing as mp

def get_saotrace_keywords(evt_file):
    keys=['ROLL_PNT','RA_PNT','DEC_PNT','ASOLFILE','EXPOSURE']

    out = dict()

    for key in keys:
        dmkeypar.punlearn()
        dmkeypar(evt_file,key=key)
        out[key] = dmkeypar.value

    dmlist.punlearn()
    data=dmlist(f"{evt_file}[GTI7][#row=1]",'data,raw')
    out['TSTART'] = data.splitlines()[1].split()[0]

    return out        

def create_saotrace_runconf(keyvals,ra,dec,spectrum_file):

    #with open('saotrace.config','w') as f:
    #    print(f'ss_evt_file:"{evt_file}"',file=f)
    #    print('ss_src_list_file="srcs.rdb"',file=f)
    #    print(f'ss_n_simulations={npsf_sims}',file=f)
    #    print('ss_marx_output_dir="marx_outdir"',file=f)
    #    print('ss_rays_output_dir="saotrace_outdir"',file=f)
#
    with open('saotrace_conf.lua','w') as f:
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

def setup_saotrace_sim(evt_file,ra,dec,spectrum_file):
    keyvals = get_saotrace_keywords(evt_file)
    create_saotrace_runconf(keyvals,ra,dec,spectrum_file)
    return keyvals

def run_sao_raytrace_parallel(evt_file,ra,dec,npsf_sims,tempdir='raysdir',spectrum_file='core_flux_saotrace.rdb'):
    keyvals = setup_saotrace_sim(evt_file,ra,dec,spectrum_file)
    Path(tempdir).mkdir(exist_ok=True)
    print(f'Performing {npsf_sims} SAOTrace simulations')

    pool = mp.Pool(mp.cpu_count())

    results = pool.map(generate_rayfile,[(tempdir,keyvals,iter) for iter in range(npsf_sims)])
    #r.wait()
    pool.close()

    return results 


def generate_rayfile(arg):
    tempdir,keyvals,iter = arg
    output = subprocess.check_output(
            [
                'trace-nest',
                f'tag={tempdir}/HRMS_rayfile_{iter}',
                'srcpars=saotrace_conf.lua',
                f'tstart={keyvals["TSTART"]}',
                f'limit={keyvals["EXPOSURE"]}',
                'limit_type=sec',
                f'seed1={int(np.random.random()*100000)}'

            ]
        )
    logger.info(output.decode("utf-8") )
    return f'{tempdir}/HRMS_rayfile_{iter}.fits'


def run_sao_raytrace(evt_file,ra,dec,npsf_sims,tempdir='raysdir',spectrum_file='core_flux_saotrace.rdb'):
    keyvals = setup_saotrace_sim(evt_file,ra,dec,npsf_sims,spectrum_file)

    Path(tempdir).mkdir(exist_ok=True)
    rayfiles = []

    print(f'Performing {npsf_sims} SAOTrace simulations')
    for i in range(npsf_sims):
        output = subprocess.check_output(
            [
                'trace-nest',
                f'tag={tempdir}/HRMS_rayfile_{i}',
                'srcpars=saotrace_conf.lua',
                f'tstart={keyvals["TSTART"]}',
                f'limit={keyvals["EXPOSURE"]}',
                'limit_type=sec',
                f'seed1={int(np.random.random()*100000)}'

            ]
        )
        logger.info(output.decode("utf-8") )
        rayfiles = rayfiles + [f'{tempdir}/HRMS_rayfile_{i}.fits']

    return ','.join(rayfiles)

