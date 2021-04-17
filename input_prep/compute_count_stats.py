#!/usr/bin/env python

import sys

if len(sys.argv) != 3:
    print("Usage: compute_count_stats <out_file> <region_file>")
    sys.exit(-1)

# from pointpats import PointPattern
# from pointpats.centrography import std_distance,ellipse
import numpy as np
from numpy.lib.function_base import copy

# import region
import matplotlib.pyplot as plt
from pycrates import read_file, copy_piximgvals
from ciao_contrib.runtool import dmcoords

##import region
from scipy.stats import skew, mode, kurtosis


def get_count_stats(out_file, region_file):

    region = np.flipud(copy_piximgvals(read_file(region_file)))
    output = np.loadtxt(out_file)

    imsize = output.shape[1]
    niter = int(output.shape[0] / imsize)

    count_sums = np.zeros(niter)

    for i in range(niter):
        count_sums[i] = (
            np.flipud(output[i * imsize : (i + 1) * imsize, :]) * region
        ).sum()

    return count_sums


outfile = sys.argv[1]
region_file = sys.argv[2]

count_sums = get_count_stats(outfile, region_file)
plt.hist(count_sums)
plt.show()

print(
    f"Mean: {count_sums.mean():.3f}, Std: {count_sums.std():.3f}, Fraction: {count_sums.std()/count_sums.mean()}"
)
print(f"Quantiles (16,50,84): {np.quantile(count_sums,[ .16,.50,.84])}")
