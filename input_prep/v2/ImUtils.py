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


class ImUtils():
    def __init__(self) -> None:
        pass

    @staticmethod
    def bin_image(evt_file:str,reg_file:Optional[str],binsize:Union[float,int],nsize:int,center:Union[PointPhys,List,None]=None,psf:Optional[bool]=False,out:Optional[str]=None,**kwargs)->Tuple[str,PointPhys,PointPhys]:
        """
        Generic wrapper for the rebin_img and rebin_psf functions. If a region file is specified, it will be used to compute the center of the output image or else the center will be used. If the psf flag is set to true, the maxima of the image will be centered.

        :param evt_file: Path to the events file
        :param reg_file: Path to the region file in CIAO/PHYS format
        :param binsize: The binning factor for the image
        :param nsize: Size of the output image. Currently only square images are supported
        :param center: If a region file is not specified, this will be used to center the image
        :param out: Name of the output file. If None, the outfile will be named ``img_<nsize>_<nsize>_<binsize>.fits``
        :param psf: If True, the maximum of the image will be at the center. nsize needs to be an odd number for this to take effect.
        ...
        :return: Returns a tuple with the output name,  center of the maxima if psf is set to True, and center of the image snapped to grid.
        """
        if center is None or len(center)==0:
            dmstat.punlearn()
            # get the centroid of the image
            dmstat(f"{evt_file}[sky=region({reg_file})][cols sky]", clip=True)
            vals = [float(x) for x in dmstat.out_mean.split(",")]
            xval = vals[0]
            yval = vals[1]
        else:
            xval,yval = center

        outfile = out or f'img_{nsize}x{nsize}_{binsize}.fits'

        if not psf:
            centroid = ImUtils.rebin_img(
                infile=f"{evt_file}[energy=500:7000]",
                outfile=outfile,
                binsize=binsize,
                nsize=nsize,
                xcen=xval,
                ycen=yval,
            )
        else:
            centroid = ImUtils.rebin_img(
                infile=f"{evt_file}[energy=500:7000]",
                outfile=outfile,
                binsize=binsize,
                nsize=nsize,
                xcen=xval,
                ycen=yval,
                psf=psf
            )

            # test if the maximum pixel is at the center. If not, rebin it again
            dmstat.punlearn()
            # get the centroid of the image
            dmstat(f"{outfile}[sky=region({reg_file})]", clip=True)
            vals = [float(x) for x in dmstat.out_max_loc.split(",")]
            if not np.all(centroid == vals):
                centroid = ImUtils.rebin_img(
                    infile=f"{evt_file}[energy=500:7000]",
                    outfile=outfile,
                    binsize=binsize,
                    nsize=nsize,
                    xcen=vals[0],
                    ycen=vals[1],
                    psf=psf
                )

        return (outfile,center or PointPhys(*vals),PointPhys(*centroid))

    @staticmethod
    def rebin_img(infile:str,outfile:str,binsize:Union[float,int]=0.5,nsize:int=64,xcen:Union[float,int]=4060.,ycen:Union[float,int]=4090.,psf:bool=False,**kwargs)->PointPhys:
        """
         Bins an event file with the specified size and binning factor The binning will be performed on a standard grid based on the specified factor.

         :param infile: Path to the events file
         :param outfile: Path to the output file
         :param binsize: Binning factor for for the image
         :param nsize: Size of the image
         :param xcen: X-coordinate of the center in PHYS coordinates
         :param ycen: Y-coordinate of the center in PHYS coordinates
         :param psf: Set this to True to output real values instead of integers in the output image
         ...
         :rtype: PointPhys
         """

        N = nsize
        B = binsize

        xcen3 = xcen
        ycen3 = ycen

        px = floor((xcen - floor(xcen)) / binsize) * binsize
        xcen3 = floor(xcen3) + px
        #
        py = floor((ycen - floor(ycen)) / binsize) * binsize
        ycen3 = floor(ycen3) + py


        if N % 2 == 0:
            xmin = round(xcen3 - (N * B) * 0.5, 2)
            xmax = round(xcen3 + (N * B) * 0.5, 2)
            ymin = round(ycen3 - (N * B) * 0.5, 2)
            ymax = round(ycen3 + (N * B) * 0.5, 2)
        else:
            xmin = round(xcen3 - ((N - 1) * B) * 0.5, 2)
            xmax = round(xcen3 + ((N + 1) * B) * 0.5, 2)
            ymin = round(ycen3 - ((N - 1) * B) * 0.5, 2)
            ymax = round(ycen3 + ((N + 1) * B) * 0.5, 2)
        
        outformat = 'r4' if psf else 'i4'
        newfile = f"{infile}[bin x={xmin}:{xmax}:#{N},y={ymin}:{ymax}:#{N}][opt type={outformat}]"

        dmcopy.punlearn()
        dmcopy(newfile, outfile, clobber=True)
        return PointPhys(xcen3, ycen3)
