#source('utilities.R',chdir=TRUE)

#set the cwd
setwd(getwd())

#read the config.yaml file


#' Run LIRA to detect sources using the configuration from config.yaml
#' Following are the required parameters in config.yaml
#' obs_file: This is the X-ray observation file with square dimensions
#' psf_file: PSF file. It is assumed that this file represents the psf at all points in the image
#' null_file: The null model against which the observation will be compared against
#' replica_im_template: template name of the images replicated from the null/baseline model.
#'                      For example, if the replicated images are named sim_64x64_0.5_{0...n}.fits, 
#'                      the template will be sim_64x64_0.5_%s.fits
#' n_replicas: The number of replicated images. The count should start from 0
#' 
#' Optional parameters
#' max_iter: Total number of draws from the LIRA posterior. Defaults to 2000
#' burn_in: Total number of draws to be ignored from the beginning. Defaults to 1000
#' alpha.init: initial values of the smoothing parameters. Defaults to [0.3,0.4,0.5,0.6,0.7,0.8]
#' thin: Number of {thin}th draws to store. Defaults to 1
#' output_dir: Output directory to store all the outputs. Defaults to LIRA_outputs
#' n_cores: LIRA will be run simultaneously on n_cores each processing a different image.
#'          Defaults to all the available cores. 
#' @export
detect_sources_LIRA <- function(){
    if(!file.exists('config.yaml')){
        stop('config.yaml does not exist')
    }

    config.file <- 'config.yaml'
    config <- get_config(config.file)
    config <- process_config(config)

    #generate the pixel masks
    cat('Generating pixel masks...')
    masks <- generate.mask_files(config)
    if(masks$status$status_code !=0){
        return(list(output=NULL,status=masks$status$status))
    }
    cat('Done\n')
<<<<<<< Updated upstream
<<<<<<<< Updated upstream:lira.parallel/R/run_lira_on_images.R
    payloads <- generate.payloads(config)
========
    payloads = generate.payloads(config)
>>>>>>>> Stashed changes:lira.processor/R/run_lira_on_images.R
=======
    payloads <- generate.payloads(config)
>>>>>>> Stashed changes

    cat('Running LIRA on all the images...')
    cluster <- parallel::makeCluster(config$n_cores
                #,outfile=""
                )
    parallel::clusterEvalQ(cluster,{
        library(FITSio)
        library(lira)
    })
    
<<<<<<< Updated upstream
<<<<<<<< Updated upstream:lira.parallel/R/run_lira_on_images.R
    results <- parallel::parLapply(cluster,payloads,run_LIRA)
========
    results = parallel::parLapply(cluster,payloads,run_LIRA)
>>>>>>>> Stashed changes:lira.processor/R/run_lira_on_images.R
=======
    results <- parallel::parLapply(cluster,payloads,run_LIRA)
>>>>>>> Stashed changes

    parallel::stopCluster(cluster)
    cat('Done\n')
    #compute_and_save_xi_and_p.ul(results,config,masks)
    compute_and_save_xi_and_p_ul_t(results,config,masks)
    
    cat('Done')
    cat('\n\n\n\nYeah, you\'re welcome!\n')
    return("")

}

