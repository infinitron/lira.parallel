#!/usr/bin/env python

import sys

if len(sys.argv) != 5:
    print("Usage: compute_offset <out_file> <region_file> <r_xcen> <r_ycen>")
    sys.exit(-1)

from pointpats import PointPattern
from pointpats.centrography import std_distance, ellipse
import numpy as np
from numpy.lib.function_base import copy

# import region
import matplotlib.pyplot as plt
from pycrates import read_file, copy_piximgvals
from ciao_contrib.runtool import dmcoords
import region
from scipy.stats import skew, mode, kurtosis


def centroid(mat, region_file):
    idx = np.nonzero(mat)

    imsize = mat.shape[0]

    idx_x = idx[1] + 1
    idx_y = imsize - idx[0]

    if np.isnan(mat.any()):
        print("Nans in the mat!")

    cx = (mat[imsize - idx_y, idx_x - 1] * idx_x).sum() / mat[
        imsize - idx_y, idx_x - 1
    ].sum()
    cy = (mat[imsize - idx_y, idx_x - 1] * idx_y).sum() / mat[
        imsize - idx_y, idx_x - 1
    ].sum()

    # dmcoords(infile=region_file,option='logical',logicalx=cx,logicaly=cy)

    return cx, cy


def get_bin_size(reg_file):
    dmcoords.punlearn()
    dmcoords(infile=region_file, option="logical", logicalx=1, logicaly=1)
    x = dmcoords.x

    dmcoords.punlearn()
    dmcoords(infile=region_file, option="logical", logicalx=2, logicaly=1)
    x2 = dmcoords.x

    return np.abs(x - x2)


def get_centroids(out_file, region_file):

    region = np.flipud(copy_piximgvals(read_file(region_file)))
    output = np.loadtxt(out_file)

    imsize = output.shape[1]
    niter = int(output.shape[0] / imsize)

    avg_image = np.zeros((imsize, imsize))

    xvals = np.zeros(niter)
    yvals = np.zeros(niter)

    count_sums = np.zeros(niter)

    dmcoords.punlearn()
    for i in range(niter):
        avg_image += output[i * imsize : (i + 1) * imsize, :]
        xvals[i], yvals[i] = centroid(
            np.flipud(output[i * imsize : (i + 1) * imsize, :]) * region, region_file
        )
        count_sums[i] = ((output[i * imsize : (i + 1) * imsize, :]) * region).sum()

    return xvals, yvals, count_sums


def calculate_offsets(x_xcen, x_ycen, r_xcen, r_ycen, region_file, binsize):

    # dmcoords(infile=region_file,option='sky',x=r_xcen,y=r_ycen)
    dmcoords(infile=region_file, option="cel", celfmt="hms", ra=r_xcen, dec=r_ycen)

    # print(x_xcen.mean(),x_ycen.mean())

    # print(dmcoords.logicalx,dmcoords.logicaly)

    x = x_xcen - dmcoords.logicalx
    y = x_ycen - dmcoords.logicaly

    nan_indices = np.logical_and(np.isnan(x_xcen), np.isnan(x_ycen))
    x = x[~nan_indices]
    y = y[~nan_indices]

    offsets = np.sqrt(np.add(x ** 2, y ** 2))

    plt.hist(offsets * binsize * 0.492)
    plt.show()

    return offsets.mean(), offsets.std(), offsets


outfile = sys.argv[1]
region_file = sys.argv[2]
r_xcen = sys.argv[3]
r_ycen = sys.argv[4]
bin_size = get_bin_size(region_file)

x_xcens, y_ycens, count_sums = get_centroids(outfile, region_file)
nan_indices = np.logical_and(np.isnan(x_xcens), np.isnan(y_ycens))
x_xcens = x_xcens[~nan_indices]
y_ycens = y_ycens[~nan_indices]
# calculate the standard distance deviation (spatial dispersion)
pp = PointPattern(list(zip(x_xcens, y_ycens)))

print(
    f"Standard distance deviation: {std_distance(pp.points) * bin_size * 0.492:.3f} arcsec"
)

# calculate the standard deviational ellipse
sx, sy, theta = ellipse(pp.points)

print(
    f"Standard deviational ellipse: X->{sx*bin_size*0.492:.3f}, Y->{sy*bin_size*0.492:.3f}, Theta->{np.degrees(theta):.2f}"
)

offset_mean, offset_std, offsets = calculate_offsets(
    x_xcens, y_ycens, r_xcen, r_ycen, region_file, bin_size
)

# print(offset_mean,offset_std)

print(
    f"The offset is: {offset_mean * bin_size * 0.492:.3f} +/- {offset_std * bin_size * 0.492:.3f} arcsec"
)

print(
    f"Skewness: {skew(offsets)}, kurtosis: {kurtosis(offsets)} mode: {mode(offsets).mode[0]* bin_size * 0.492:.3f}, median: {np.median(offsets)* bin_size * 0.492:.3f}"
)

np.savetxt("offsets.txt", offsets)
