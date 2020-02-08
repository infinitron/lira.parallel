#
# Required CIAO 4.5/contrib package
#
from ciao_contrib.runtool import *
import os
from math import floor

# rebin_img function to make an image with the center pixel at the region centroid
# infile - input psf event file from SAOtrace
# outfile - output file, binned psf image 
# binsize - binning of the image
# nsize - size of the image in image pixels
# xcenter, ycenter - initial centroid

    
def rebin_img(infile="evt2.fits", outfile="test.fits",
              binsize=0.25, nsize=64, xcen=4060.0, ycen=4090.0):
    

    #for i in range(5):
    #    filename = "{}[sky=circle({},{},10)][cols sky]".format(infile,xcen, ycen)
    #    dmstat(filename, verbose=0)
    #    vals = [float(x) for x in dmstat.out_mean.split(',')]
    #    xcen=vals[0]
    #    ycen=vals[1]
    #    #print vals

    # calculate the size
    N = nsize
    B = binsize
    #xcen3 = vals[0]
    #ycen3 = vals[1]
    xcen3=xcen
    ycen3=ycen
    #xmin = xcen3 - (N*B+B/2)*0.5
    #xmax = xcen3 + (N*B+B/2)*0.5
    #ymin = ycen3 - (N*B+B/2)*0.5
    #ymax = ycen3 + (N*B+B/2)*0.5

    xmin = floor(xcen3)+0.5 - (N*B)*0.5
    xmax = floor(xcen3)+0.5 + (N*B)*0.5
    ymin = floor(ycen3)+0.5 - (N*B)*0.5
    ymax = floor(ycen3)+0.5 + (N*B)*0.5

    newfile = "{}[bin x={}:{}:{},y={}:{}:{}][opt type=i4]".format(infile, xmin, xmax, B, ymin, ymax, B)

    dmcopy(newfile, outfile, clobber=True)
    
    

