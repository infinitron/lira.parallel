#library(yaml)
#library(parallel)
#source('./lira_processor.R',chdir=TRUE)


# An operator similar to .format in python
`%--%` <- function(x, y) {

  do.call(sprintf, c(list(x), y))

}

# Reads the YAML config file and return the config
get_config <- function(config_file){
    if(!file.exists(config_file)){
        stop("File '%s' doesnot exist: " %--% c(config_file))
    }

    return(yaml::read_yaml(config_file))
}

# Tarnish the configuration object. Perform existeny checks on the required params and provide
# defaults for the optional params if not specified
process_config <- function(config){
    #check for required objects
    required_params <- c('regions','obs_file','psf_file','null_file','replica_im_template','n_replicas')
    for (param in required_params){
        if(is.null(config[[param]])){
            stop("Parameter %s is required!" %--% c(param))
        }
    }

    #generate a vector of replicated images
    replicated_images <- sapply(0:(config$n_replicas-1)
            ,function(im){
                return(config$replica_im_template %--% c(im))
            })

    #all images for LIRA to be run on
    config$all_obs_files <- c(config$obs_file,replicated_images)

    #if the number of cores are not set, use all cores by default
    n_cores <- parallel::detectCores()


    #set the defaults if the optional params aren't provided
    optional_params <- list(max_iter=2000
        ,burn_in=1000
        ,alpha.init=c(0.3,0.4,0.5,0.6,0.7,0.8)
        ,thin=1
        ,output_dir='LIRA_outputs'
        ,n_cores=n_cores
        ,post_only=FALSE)

    for (param in names(optional_params)){
        if(is.null(config[[param]])){
            config[[param]] <- optional_params[[param]]
        }
    }

    if(config[['post_only']]==TRUE){
        cat('Skipping LIRA runs. Post-processing the data...\n')
    }

    if (!dir.exists(config$output_dir)){
        dir.create(config$output_dir)
    }

    return(config)
}

# Generate a status list object that functions can return to notify of success/failure
# Obj keys: status_code, err_msg
generate.status <- function(status_code,err_msg){
    return(
            list(status_code=status_code,err_msg=err_msg)
        )
}

# The processed config object will contain all the images that are to be compared 
# against the baseline model. An array of payloads each which can be individually  
# processed by LIRA are generated
generate.payloads <- function(processed_config){
    all_obs_files <- processed_config[['all_obs_files']]

    payloads <- vector("list",length(all_obs_files))

    i <-1

    for (obs_file in all_obs_files){
        payload <- processed_config
        payload$obs_file <- obs_file
        payloads[[i]] <- payload
        i <- i+1
    }

    return(payloads)
}

#Generate a pixel mask (image) for each mask file
generate.mask_files = function(config){
    obs_file <- config$obs_file
    regions <- config$regions
    masks <- c()
    region_names <- c()

    for(region in regions){
        region_name <- tools::file_path_sans_ext(region)
        region_names <- c(region_names,region_name)
        output <- generate_pixel_mask_from_region_file(obs_file,region)
        if(output$status$status_code !=0){
            return(output$status)
        }
        masks <- c(masks,output$mask_file)
    }

    return(
        list(
            mask_files = masks,
            status = generate.status(0,''),
            region_names = region_names
        )
    )
}

# Check if all the files required by LIRA exist on the disk
payload.consistency_check <- function(payload){
    #check if all the files in the payload exist
    required_files <- c('obs_file','psf_file','null_file')
    for (file in required_files){
        
        if(!file.exists(payload[[file]])){
            #print(payload[[file]])
            #return(generate.status(-1,"File %s doesn't exist" %--% c(payload[[file]])))
            stop("File %s doesn't exist" %--% c(payload[[file]]))
        }
    }

    for (region in payload$regions){
        if(!file.exists(region)){
            return(generate.status(-1,"File %s doesn't exist" %--% c(region)))
        }
    }

    return(generate.status(0,''))
}

# This function is a wrapper around a slightly modified version of KMc's liraOutput function
# Accepts a payload and runs LIRA on it. 
run_LIRA = function(payload){
    #consistency checks of the files
    consistency = payload.consistency_check(payload)
    if(consistency[['status_code']]!=0){
        return(list(output=list(),status=consistency))
    }

    #Launch!
    lira_output <- launchLIRA(payload$obs_file
            ,startFile=payload$null_file
            ,bkgFile=payload$null_file
            ,psfFile=payload$psf_file
            ,outDir=payload$output_dir
            ,maxIter=payload$max_iter
            ,alpha.init=sample(payload$alpha_init)
            ,thin=payload$thin
            ,burn=payload$burn_in
            ,postOnly=payload$post_only)
    
    return(
        list(
            output=lira_output,
            status=generate.status(0,'')
        )
    )
}

get_output_images <- function(results){
    output_images <- c()
    for(result in results){
        if(result$status$status_code != 0){
            stop(result$status$err_msg)
        }
        output_images <- c(output_images,result$output$out_images_file)
    }

    return(output_images)
}

