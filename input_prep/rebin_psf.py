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

    
def rebin_psf(infile="evt2.fits", outfile="test.fits",
              binsize=0.25, nsize=64, xcen=4060.0, ycen=4090.0):
    

    #for i in range(5):
    #    filename = "{}[sky=circle({},{},10)][cols sky]".format(infile,xcen, ycen)
    #    dmstat(filename, verbose=0)
    #    vals = [float(x) for x in dmstat.out_mean.split(',')]
    #    xcen=vals[0]
    #    ycen=vals[1]
    #    #print vals
    
    #get the location of the maxima of the psf
    #filename = "{}[sky=circle({},{},10)][cols sky]".format(infile,xcen, ycen)
    #dmstat(filename,v=0,centroid=True)
    
    #centroid = [float(x) for x in dmstat.out_cntrd_phys.split(',')]
    #maxloc = [float(x) for x in dmstat.out_max_loc.split(',')]
    

    # calculate the size
    N = nsize
    B = binsize
    
    xcen3=xcen
    ycen3=ycen
    
    if N%2==0:
        xmin = xcen3 - (N*B)*0.5
        xmax = xcen3 + (N*B)*0.5
        ymin = ycen3 - (N*B)*0.5
        ymax = ycen3 + (N*B)*0.5
    else:
        xmin = xcen3 - ((N-1)*B)*0.5
        xmax = xcen3 + ((N+1)*B)*0.5
        ymin = ycen3 - ((N-1)*B)*0.5
        ymax = ycen3 + ((N+1)*B)*0.5

    newfile = "{}[bin x={}:{}:{},y={}:{}:{}][opt type=i4]".format(infile, xmin, xmax, B, ymin, ymax, B)
    print(newfile)

    dmcopy(newfile, outfile, clobber=True)
    
    

