#library(yaml)
#library(parallel)
#source('./lira_processor.R',chdir=TRUE)

`%--%` <- function(x, y) {

  do.call(sprintf, c(list(x), y))

}

get_config = function(config_file){
    #read the YAML config file and return the config
    #sprintf()
    if(!file.exists(config_file)){
        stop("File '%s' doesnot exist: " %--% c(config_file))
    }

    return(yaml::read_yaml(config_file))
}

process_config = function(config){
    #Use a configuration objet to setup LIRA to run on multiple cores

    #check for required objects
    required_params = c('regions','obs_file','psf_file','null_file','replica_im_template','n_replicas')
    for (param in required_params){
        if(is.null(config[[param]])){
            stop("Parameter %s is required!" %--% c(param))
        }
    }

    #generate a vector of replicated images
    replicated_images = sapply(0:(config$n_replicas-1)
            ,function(im){
                return(config$replica_im_template %--% c(im))
            })

    #all images for LIRA to be run on
    config$all_obs_files = c(config$obs_file,replicated_images)

    #if the number of cores are not set, use all cores by default
    n_cores = parallel::detectCores()


    #set the defaults if the optional params aren't provided
    optional_params = list(max_iter=2000
        ,burn_in=1000
        ,alpha.init=c(0.3,0.4,0.5,0.6,0.7,0.8)
        ,thin=1
        ,output_dir='LIRA_outputs'
        ,n_cores=n_cores)

    for (param in names(optional_params)){
        if(is.null(config[[param]])){
            config[[param]] = optional_params[[param]]
        }
    }

    return(config)
}

generate.status = function(status_code,err_msg){
    return(
            list(status_code=status_code,err_msg=err_msg)
        )
}

generate.payloads = function(processed_config){
    all_obs_files = processed_config[['all_obs_files']]

    payloads=vector("list",length(all_obs_files))

    i=1

    for (obs_file in all_obs_files){
        payload = processed_config
        payload$obs_file = obs_file
        payloads[[i]] = payload
        i = i+1
    }

    return(payloads)
}

payload.consistency_check = function(payload){
    #check if all the files in the payload exist
    required_files = c('obs_file','psf_file','null_file')
    for (file in required_files){
        
        if(!file.exists(payload[[file]])){
        print(payload[[file]])
            return(generate.status(-1,"File %s doesn't exist" %--% c(payload[[file]])))
        }
    }

    for (region in payload$regions){
        if(!file.exists(region)){
            return(generate.status(-1,"File %s doesn't exist" %--% c(region)))
        }
    }

    return(generate.status(0,''))
}

run_LIRA = function(payload){
    #consistency checks of the files
    consistency = payload.consistency_check(payload)
    if(consistency['status_code']!=0){
        return(list(output=list(),status=consistency))
    }

    #Launch!
    lira_output = launchLIRA(payload$obs_file
            ,startFile=payload$null_file
            ,bkgFile=payload$null_file
            ,psfFile=payload$psf_file
            ,outDir=payload$output_dir
            ,maxIter=payload$max_iter
            ,alpha.init=payload$alpha_init
            ,thin=payload$thin
            ,burn=payload$burn_in)
    
    return(
        list(
            output=lira_output,
            status=generate.status(0,'')
        )
    )
}