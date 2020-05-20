#
# Required CIAO 4.5/contrib package
#
from ciao_contrib.runtool import *
import os
from math import floor,ceil

# rebin_img function to make an image with the center pixel at the region centroid
# infile - input psf event file from SAOtrace
# outfile - output file, binned psf image 
# binsize - binning of the image
# nsize - size of the image in image pixels
# xcenter, ycenter - initial centroid

    
def rebin_img(infile="evt2.fits", outfile="test.fits",
              binsize=0.25, nsize=64, xcen=4060.0, ycen=4090.0):

    N = nsize
    B = binsize

    xcen3=xcen
    ycen3=ycen

    #get the nearest natural pixel
    px = floor((xcen-floor(xcen))/binsize)*binsize
    xcen3 = floor(xcen3) + px

    py = floor((ycen-floor(ycen))/binsize)*binsize
    ycen3 = floor(ycen3) + py
    
    
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
    return (xcen3,ycen3)

#make sure that the centroid lies on the edge of the physical pixel
# nearest_x = xcen3/binsize-0.5
# nearest_y = ycen3/binsize-0.5
# 
# dx = nearest_x-floor(nearest_x)
# dy = nearest_y-floor(nearest_y)
# 
# dx_bin = dx/binsize
# dy_bin = dy/binsize
# 
# d_dx_bin = dx_bin - floor(dx_bin)
# d_dy_bin = dy_bin - floor(dy_bin)   
# 
# if d_dx_bin>=0.5:
#     xcen3 =  (ceil(d_dx_bin) + floor(nearest_x)) * binsize
# else:
#     xcen3 =  (floor(d_dx_bin) + floor(nearest_x)) * binsize
# 
# if d_dy_bin>=0.5:
#     ycen3 = (ceil(d_dy_bin) + floor(nearest_y)) * binsize
# else:
#     ycen3 =  (floor(d_dy_bin) + floor(nearest_y)) * binsize 
      
