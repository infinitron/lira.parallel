from astropy.io import fits
import numpy as np
import sys

posterior_name = sys.argv[1]
wcs_file = sys.argv[2]
if len(sys.argv)==4):
    out_name = sys.argv[3]
else:
    out_name = "post.fits"

post = np.loadtxt(posterior_name)
with fits.open(wcs_file) as hdul:
    hdul[0].data = post
    hdul.writeto(out_name)
