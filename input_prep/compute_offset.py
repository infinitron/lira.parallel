#!/usr/bin/env python

import sys

if len(sys.argv)!=5:
    print('Usage: compute_offset <out_file> <region_file> <r_xcen> <r_ycen>')
    sys.exit(-1)

import numpy as np
from numpy.lib.function_base import copy
#import region
import matplotlib.pyplot as plt
from pycrates import read_file,copy_piximgvals
from ciao_contrib.runtool import  dmcoords

import region

def centroid(mat,region_file):
    idx = np.nonzero(mat)

    imsize = mat.shape[0]

    idx_x = idx[1]+1
    idx_y = imsize-idx[0]

    cx = (mat[imsize-idx_y,idx_x-1] * idx_x).sum()/mat[imsize-idx_y,idx_x-1].sum()
    cy = (mat[imsize-idx_y,idx_x-1] * idx_y).sum()/mat[imsize-idx_y,idx_x-1].sum()

    #dmcoords(infile=region_file,option='logical',logicalx=cx,logicaly=cy)

    return cx,cy

def get_bin_size(reg_file):
    dmcoords.punlearn()
    dmcoords(infile=region_file,option='logical',logicalx=1,logicaly=1)
    x = dmcoords.x

    dmcoords.punlearn()
    dmcoords(infile=region_file,option='logical',logicalx=2,logicaly=1)
    x2 = dmcoords.x

    return np.abs(x-x2)

def get_centroids(out_file,region_file):

    region = np.flipud(copy_piximgvals(read_file(region_file)))
    output = np.loadtxt(out_file)

    imsize = output.shape[1]
    niter = int(output.shape[0]/imsize)

    avg_image = np.zeros((imsize,imsize))

    xvals = np.zeros(niter)
    yvals = np.zeros(niter)

    dmcoords.punlearn()
    for i in range(niter):
        avg_image += output[i*imsize:(i+1)*imsize,:]
        xvals[i], yvals[i] = centroid(np.flipud(output[i*imsize:(i+1)*imsize,:]) * region,region_file)
    
    return xvals, yvals

def calculate_offsets(x_xcen,x_ycen,r_xcen,r_ycen,region_file,binsize):

    #dmcoords(infile=region_file,option='sky',x=r_xcen,y=r_ycen)
    dmcoords(infile=region_file,option='cel',celfmt='hms',ra=r_xcen,dec=r_ycen)

    print(x_xcen.mean(),x_ycen.mean())
    
    print(dmcoords.logicalx,dmcoords.logicaly)

    x = x_xcen - dmcoords.logicalx
    y = x_ycen - dmcoords.logicaly

    offsets = np.sqrt(np.add(x**2,y**2))

    plt.hist(offsets*binsize*0.492);plt.show()

    return offsets.mean(), offsets.std()




outfile = sys.argv[1]
region_file = sys.argv[2]
r_xcen = sys.argv[3]
r_ycen = sys.argv[4]
bin_size = get_bin_size(region_file)

x_xcens, y_ycens = get_centroids(outfile,region_file)
offset_mean, offset_std = calculate_offsets(x_xcens, y_ycens, r_xcen, r_ycen,region_file,bin_size)

print(offset_mean,offset_std)

print(f'The offset is: {offset_mean * bin_size * 0.492} +/- {offset_std * bin_size * 0.492}')
