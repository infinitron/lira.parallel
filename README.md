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
    * knitr
    * sm
    


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
* PSF image (size should be <img src="https://render.githubusercontent.com/render/math?math=\leq"> the size of the image)
* Baseline image (usually point source+flat background; same size as the input image)
* Baseline replicates (*faked* images of the baseline model; same size as the input image)
* ROI/region files (CIAO format) to test for significance
* ```config.yaml```

Execute the following commands to run LIRA on the input files
```bash
$> R
R> library(lira.parallel)
R> detect_sources_LIRA() -> yes.please

Generating pixel masks...Done
Running LIRA on all the images...Done
Generating the distributions of Xi...Done


|region.names |p.values |p.values2 |
|:------------|:--------|:---------|
|core_lira    |0.005    |0.015     |
|knot_lira    |0.005    |0.005     |
|C            |1        |1         |
Done

Yeah, you are welcome!
R>
```
Two p-values (upperlimits) are reported for each region. For the first one, the <img src="https://render.githubusercontent.com/render/math?math=1-\gamma"> quantile is estimated directly on the distribution while a kernel density estimate is used for the second (not to be used for scientific purposes).

## How does it work?
```detect_sources_LIRA()``` executes the following steps:
* Creates mask images from the region files.
* Runs LIRA on both the observed image and the baseline replicates. 
* Post processes the output images from LIRA to extract the posterior distributions of <img src="https://render.githubusercontent.com/render/math?math=\xi">

## Input preparation
Python helper scripts to prepare inputs are located in the [input_prep](input_prep/) folder. These helper scripts are tailored for use with *Chandra* observations.

 After configuring ```lira_input_prep.yaml``` in the folder with inputs (event files, region files), exceute ```create_inputs_for_lira.py```. This will:
* Create a binned image from the events file
* Extract the spectrum of the core
* Fit an absorbed powerlaw model to the spectrum
* Create a PSF of the core with matched pixel boundaries
* Create a baseline image and the specified number of replicates.

These inputs can be directly used with the ```detect_sources_LIRA()``` function.










