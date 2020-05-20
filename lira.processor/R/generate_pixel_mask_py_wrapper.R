# Generate a pixel mask from a region file that overlaps with the image file.
# The values of the pixels withing the region will be assigned one and zero to
# the ones outside
generate_pixel_mask_from_region_file = function(image_file,reg_file){

    # do the necessary imports
    
    #NOTE: for some reason reticulate doesn't find the contrib scripts
    # It will be added manually for now until the root cause is known
    #syspy = reticulate::import("sys",convert=FALSE)
    ascds_contrib_path = Sys.getenv("ASCDS_CONTRIB")
    syspy$path$append(file.path(ascds_contrib_path,"lib/python3.5/site-packages")) #temporary workaround
    ascds_lib_path=Sys.getenv("ASCDS_LIB") #another temporary workaround
    syspy$path$append(file.path(ascds_lib_path,"python3.5/site-packages"))
    syspy$path$append(file.path(ascds_lib_path,"python3.5/site-packages/paramio"))
    #ciao.runtool_py = reticulate::import("ciao_contrib.runtool",convert=FALSE)
    
    os_py = reticulate::import("os")
    pathlib_py = reticulate::import("pathlib")
    tempfile_py = reticulate::import("tempfile")

    if  (!(os_py$path$exists(image_file) && os_py$path$isfile(image_file))){
        return(generate.status(-1,'File %s does not exist ' %--% c(image_file)))
    }
    
    if (!(os_py$path$exists(reg_file) && os_py$path$isfile(reg_file))){
        return(generate.status(-1,'File %s does not exist ' %--% c(reg_file)))
    }

    img_basename = pathlib_py$Path(image_file)$stem
    img_ext = pathlib_py$Path(image_file)$suffix
    reg_basename = pathlib_py$Path(reg_file)$stem

    #generate the temp and outfile names
    temp_file = '%s%s' %--% c(reticulate::iter_next(tempfile_py[["_get_candidate_names"]]()), img_ext)
    mask_file = '%s%s%s%s' %--% c(img_basename, '_', reg_basename, img_ext)

    ciao.runtool_py$dmimgcalc$punlearn()
    #create an 'ones' image with the same size
    ciao.runtool_py$dmimgcalc(infile=image_file
            ,infile2=NULL
            ,outfile=temp_file
            ,operation="imgout=(1+(img1-img1))"
            ,clobber=TRUE)
    
    ciao.runtool_py$dmcopy$punlearn()
    #apply the region and generate the pixel mask
    ciao.runtool_py$dmcopy(infile='%s[sky=region(%s)][opt full]' %--% c(temp_file,reg_file)
            ,outfile=mask_file,
            clobber=TRUE)

    unlink(temp_file)

    return(
        list(
            mask_files=mask_file,
            status=generate.status(0,'')
        )
    )

}
#generate_pixel_mask_from_region_file('img_64x64_0.5.fits','core.reg')