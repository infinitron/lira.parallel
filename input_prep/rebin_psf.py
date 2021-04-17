#
# Required CIAO 4.5/contrib package
#
from ciao_contrib.runtool import *
import os
from math import floor, ceil, modf

# rebin_img function to make an image with the center pixel at the region centroid
# infile - input psf event file from SAOtrace
# outfile - output file, binned psf image
# binsize - binning of the image
# nsize - size of the image in image pixels
# xcenter, ycenter - initial centroid


def rebin_psf(
    infile="evt2.fits",
    outfile="test.fits",
    binsize=0.25,
    nsize=64,
    xcen=4060.0,
    ycen=4090.0,
):

    N = nsize
    B = binsize

    xcen3 = xcen
    ycen3 = ycen
    # nearest half-pixel vertex
    # xcen3 = floor(xcen-B/2.0) + round(modf(xcen- B/2.0)[0]*2.0/B)*B/2.0
    # ycen3 = floor(ycen-B/2.0) + round(modf(ycen- B/2.0)[0]*2.0/B)*B/2.0

    px = floor((xcen - floor(xcen)) / binsize) * binsize
    xcen3 = floor(xcen3) + px
    #
    py = floor((ycen - floor(ycen)) / binsize) * binsize
    ycen3 = floor(ycen3) + py

    ##get the nearest vertex
    # if abs(xcen3-xcen)>abs(xcen3-xcen+B/2):
    #    xcen3=xcen3+B/2
    #
    # if abs(ycen3-ycen)>abs(ycen3-ycen+B):
    #    ycen3=ycen3+B/2
    # print(xcen3,ycen3,xcen4,ycen4,xcen,ycen)

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
    newfile = "{}[bin x={}:{}:#{},y={}:{}:#{}][opt type=r4]".format(
        infile, xmin, xmax, N, ymin, ymax, N
    )
    # print('core_psf_bin_string ',"x={}:{}:#{},y={}:{}:#{}".format(xmin, xmax, N, ymin, ymax, N))

    dmcopy.punlearn()
    dmcopy(newfile, outfile, clobber=True)
    return (xcen3, ycen3)
