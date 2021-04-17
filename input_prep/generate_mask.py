#!/usr/bin/env python

from ciao_contrib.runtool import dmimgcalc, dmcopy
from os.path import exists, isfile
from pathlib import Path
from os import remove
import sys
import tempfile


def generate_pixel_mask_from_region_file(image_file, reg_file):
    """Generate a pixel mask from a region file that overlaps with the image file.
    The values of the pixels withing the region will be assigned one and zero to
    the ones outside"""

    if not (exists(image_file) and isfile(image_file)):
        return {
            "status_code": -1,
            "err_msg": "File {0} does not exist ".format(image_file),
        }

    if not (exists(reg_file) and isfile(reg_file)):
        return {
            "status_code": -1,
            "err_msg": "File {0} does not exist ".format(reg_file),
        }

    img_basename = Path(image_file).stem
    img_ext = Path(image_file).suffix
    reg_basename = Path(reg_file).stem

    # generate the temp and outfile names
    temp_file = next(tempfile._get_candidate_names()) + "." + img_ext
    mask_file = img_basename + "." + reg_basename + img_ext

    dmimgcalc.punlearn()
    # create an 'ones' image with the same size
    dmimgcalc(
        infile=image_file,
        infile2=None,
        outfile=temp_file,
        operation="imgout=(1+(img1-img1))",
        clobber=True,
    )

    dmcopy.punlearn()
    # apply the region and generate the pixel mask
    dmcopy(
        infile="{0}[sky=region({1})][opt full]".format(temp_file, reg_file),
        outfile=mask_file,
        clobber=True,
    )

    remove(temp_file)


if len(sys.argv) != 3:
    print("Usage: generate_mask.py <image_file> <region file>")
    sys.exit(-1)

image_file = sys.argv[1]
region_file = sys.argv[2]


generate_pixel_mask_from_region_file(image_file, region_file)
