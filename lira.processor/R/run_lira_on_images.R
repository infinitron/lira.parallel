#source('utilities.R',chdir=TRUE)

#set the cwd
setwd(getwd())

#read the config.yaml file


#' Run LIRA to detect sources using the configuration from config.yaml
#' @export
detect_sources_LIRA = function(){
    if(!file.exists('config.yaml')){
        stop('config.yaml does not exist')
    }

    config.file = 'config.yaml'
    config = get_config(config.file)
    payloads = generate.payloads(
    process_config(
            config
        )
    )

    cluster = parallel::makeCluster(config$n_cores)
    parallel::clusterEvalQ(cluster,{
    #source('utilities.R')
    #    source('lira_processor.R')
    library(FITSio)
    library(lira)
    })
    
    #run_LIRA(payloads[[1]])
    #print('ok')
    results = parallel::parSapply(cluster,payloads,run_LIRA)
    parallel::stopCluster(cluster)
    return(results)

}

