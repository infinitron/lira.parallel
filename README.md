# LIRA parallel

This R module is a simple parallel processing wrapper for [LIRA](https://github.com/astrostat/LIRA), which runs LIRA on individual images on individual cores in parallel.

## What is LIRA?
The Low-counts Image Reconstruction and Analysis (LIRA) is a robust statistical tool that is designed to estimate the statistical significance of emission from low-count images (e.g., shallow *Chandra* observations) while accounting for emission from surrounding bright objects and Poisson background. See its github [page](https://github.com/astrostat/LIRA) for more details and references. 

## What does LIRA parallel do?
LIRA was written as an R module which processes one image at a time. LIRA processor can be used to run LIRA on multiple images at once and automatically estimate the significance of emission from chosen ROIs in the observation.

## Prerequisites
LIRA processor depends on the following:
* Latest version of CIAO (conda installation [instructions](https://cxc.harvard.edu/ciao/download/conda.html))
* R packages (will be automatically installed if using devtools)
    * devtools (optional but recommended)
    * lira
    * Fitsio
    * yaml
    * parallel
    * spatstat
    * tools
    * reticulate
    


## Installation 
After cloning/downloading the code, execute the following commands in the root directory to install the package.

### Using devtools (recommended)

```bash
$ cd lira.processor
$ R
R> library(devtools)
R> devtools::install() 
```

### Using default commands
```bash
$> R CMD build lira.processor
$> R
R> packages.install("lira.parallel_0.1.tar.gz",source=NULL)
```

The above commands install a packaged called ```lira.parallel``` to the location specified by the ```R_LIBS``` environment variable.

## Execution
The following input files are required before running LIRA:
* Input observation (size should be a power of two)
* PSF image (size should be $\leq$ the size of the image)
* Baseline image (usually point source+flat background; same size as the input image)
* Baseline replicates (*faked* images of the baseline model; same size as the input image)
* ROI/region files (CIAO format) to test for significance
* ```config.yaml```

A sample ```config.yaml``` file is located in the [input_prep](input_prep/) directory




