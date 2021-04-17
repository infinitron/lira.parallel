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
from TypeDefs import PointCel,PointPhys

class CoordUtils():
    def __init__(self) -> None:
        pass

    @staticmethod
    def get_ra_dec_of_max(img_file: str,reg_file: str,**kwargs)->PointCel:
        """
        Return maxima's ra and dec in decimal degrees

        :param img_file: Path to the events/image file
        :param reg_file: Path to the region file with CIAO/PHYS format that covers the desired region
        ...
        :rtype: PointCel
        """
        dmstat.punlearn()
        dmstat(f"{img_file}[sky=region({reg_file})]")

        vals = [float(x) for x in dmstat.out_max_loc.split(",")]

        dmcoords.punlearn()
        dmcoords(img_file, option="sky", x=vals[0], y=vals[1], celfmt="deg")
        return PointCel(dmcoords.ra, dmcoords.dec)

    @staticmethod
    def get_centroid(evt_file:str,reg_file:str,fmt:str,**kwargs)->Union[PointPhys,PointCel]:
        """
        Return the centroid of a region

        :param img_file: Path to the events file
        :param reg_file: Path to the region file in CIAO/PHYS format
        :param fmt: Should be either of "phys" or "cel" to return the coordinates in physical detector or WCS, respectively
        ...
        :raises ValueError: Raised when an invalid format is specified
        ...
        :rtype: PointPhys or PointCel
        """

        if fmt not in ('phys','cel'):
            raise ValueError("fmt must be one of phys or cel")
        dmstat.punlearn()
        dmstat(f"{evt_file}[sky=region({reg_file})][cols sky]")

        vals = [float(x) for x in dmstat.out_mean.split(",")]

        if fmt == 'phys':
            return PointPhys(*vals)
        else:
            dmcoords.punlear()
            dmcoords(evt_file, option="sky", x=vals[0], y=vals[1], celfmt="deg")
            return PointCel(dmcoords.ra,dmcoords.dec)
