# Generate distributions for Xi using the output images and masks
get_matrix_from_fits = function(file_name){
    if(!file.exists(file_name)) return(list(data_mat=NULL,nrows=NULL,ncols=NULL))

    #readfits prints the filename. It will be suppressed
    invisible(capture.output(data <- FITSio::readFITS(file_name)))


    data_mat = matrix(data=data$imDat
        ,nrow=data$axDat$len[1]
	    ,ncol=data$axDat$len[2])

    return(list(data_mat=data_mat,nrows=data$axDat$len[1],ncols=data$axDat$len[2]))
}

get_matrix_from_file = function(file_name){
    if(!file.exists(file_name)) return(list(data_mat=NULL,nrows=NULL,ncols=NULL))
}


#' Compute the xi for all the iterations in each image
#' @export
generate.distribution_xi = function(images,mask_file,null_file,n_iter,thin=1){
    


    #print(images)
    #print('in')
    #read each mask in to a matrix
    region_mat_obj = get_matrix_from_fits(mask_file)
    #print('out')
    #print(names(region_mat_obj))
    if(is.null(region_mat_obj$data_mat)){
        print('File %s does not exist' %--% mask_file )
        return(list(output=NULL,generate.status(-1,'File %s does not exist' %--% mask_file)))
    }

    np = reticulate::import("numpy")
    null_mat_obj = get_matrix_from_fits(null_file)

    
    #compute the null counts
    nrows = null_mat_obj$nrows
    ncols = null_mat_obj$ncols
    mask_counts.frac = sum(region_mat_obj$data_mat)/(nrows*ncols)
    total_null_counts = sum(null_mat_obj$data_mat) * mask_counts.frac
    #print(total_null_counts)
    null_counts = sum(null_mat_obj$data_mat * region_mat_obj$data_mat)
    null_nrows = null_mat_obj$nrows

    # For each image, compute the distribution of Xi
    xi.distribution.region=list()
    n_images = length(images)
    xi.all_iter=matrix(0,nrow=n_images,ncol=as.integer(n_iter/thin))

    for(im_number in 1:n_images){

        #read the image
        #This image contains all the draws from LIRA
        image_file = images[[im_number]]

        if(!file.exists(image_file)){
            print('File %s does not exist' %--% image_file)
            next
        }

        im_mat = np$loadtxt(image_file)

        im_mat_obj = list(data_mat=im_mat,nrows=dim(im_mat)[[1]],ncols=dim(im_mat)[[2]])


        xi.iteration_wise = c()

        for(i in seq(1,im_mat_obj$nrows,null_nrows)){

            image_iter = im_mat_obj$data_mat[i:(i+null_nrows-1),]
            #transpose the matrix read from numpy since the LIRA output images are joined side ways and not top bottom
            im_counts = sum(t(image_iter) * region_mat_obj$data_mat)
            if(im_counts==0) {
                log_xi = -8.5 #some far away value so as to not to interfere with the rest of the distribution
            }else{
                log_xi = log(im_counts/(im_counts+null_counts),10)
            }
            

            xi.iteration_wise = c(xi.iteration_wise,log_xi)
        }

        xi.all_iter[im_number,] = xi.iteration_wise

    }

    return(list(output=xi.all_iter,status=generate.status(0,''),total_null_counts=null_counts))
}

generate.distribution_xi.wrapper = function(payload){
    return(
        generate.distribution_xi(payload$images
                                ,payload$mask_file
                                ,payload$null_file
                                ,payload$n_iter
                                ,payload$thin)
    )
}

#' Take the distibutions of xi, save them as pdf/mat files, and compute the upper bound on p value
#' @export
save_distributions_xi = function(distributions,region_name,out_dir){

    #print('in')
    #first row contains the obs vs baseline
    #the rest would contain replicas vs baseline

    dims = dim(distributions)
    nrows = dims[[1]]
    ncols = dims[[2]]

    # transform the vector to have all the distributions on a single row
    all_dist_vec = matrix(t(distributions),nrow=1)

    #compute the upper bound
    gamma = 0.005
    c = quantile(10^(all_dist_vec[1,(ncols+1):(nrows*ncols)]),1-gamma)
    c2=spatstat::quantile.density(density(10^(all_dist_vec[1,(ncols+1):(nrows*ncols)])),1-gamma)
    
    t_c.yobs = 1/(ncols)*sum(10^(all_dist_vec[,1:ncols])>=c)
    t_c2.yobs = 1/(ncols)*sum(10^(all_dist_vec[,1:ncols])>=c2)

    #cat("c: ",c,", c2: ",c2,", tc: ",t_c.yobs,", tc2: ",t_c2.yobs,"\n")
    if(t_c.yobs==0) t_c.yobs=gamma
    if(t_c2.yobs==0) t_c2.yobs=gamma
    p.value.upper_lim = gamma/t_c.yobs
    p.value.upper_lim2 = gamma/t_c2.yobs

    if(p.value.upper_lim>1)p.value.upper_lim=1
    if(p.value.upper_lim2>1)p.value.upper_lim2=1


    groups = rep(1:nrows,rep(ncols,nrows))
    groups = c(groups,rep(nrows+1,(nrows-1)*ncols))
    colors = rep("gray",nrows)
    colors[[1]] = "cyan"
    colors=c(colors,"black") #for the mean curve

    all_dist_vec = c(all_dist_vec[1,],all_dist_vec[1,(ncols+1):(nrows*ncols)])
    line_types = rep(2,(nrows+1))
    line_types[[1]] = 1
    line_types[(nrows+1)]=1


    line_widths = rep(1,(nrows+1))
    line_widths[(nrows+1)] = 2

    #graphing options
    #xlim = c(-10,0)
    #ylim=c(0,1)
    ngrid=200


    
    #save the distribution as a pdf
    pdf(file.path(out_dir,paste(region_name,".pdf",sep="")), width=12,height=4.25)
        sm::sm.density.compare(all_dist_vec,groups,col=colors,lty=line_types,lwd=line_widths
                #,xlim=xlim
                #,ylim=ylim
                ,ngrid=ngrid)
    dev.off()

    #also dump the data to a file
    cat(all_dist_vec,file=file.path(out_dir,paste(region_name,"_all_dist.out",sep="")))
    cat(groups,file=file.path(out_dir,paste(region_name,"_groups.out",sep="")))

    return(list(
        p_value=p.value.upper_lim
        ,p_value2 = p.value.upper_lim2
        ,status=generate.status(0,'')    
    ))
}

empty_datastructure_xi_distribution <- function(images,n_iter){
    ds <- list()
    for(image in images){
        ds[[image]] <- rep(0.0,n_iter)
    }
    return(ds)
}

#This function takes list of output images, mask files, a null file, total number of iterations and the thinning factor
#For each output image and for each mask file, it computes a distribution of Xi
#Each output image is read once and computations are made on it
generate_distribution_xi_t = function(images,mask_files,null_file,n_iter){
    #read all the masks into an array
    masks <- list()
    #hosts the distribution of Xi for each mask
    xi.distribution.region <- list()

    n_images <- length(images)

    #read the null file.
    null.model <- get_matrix_from_fits(null_file)

    complimentary.mask = matrix(0,
        nrow=null.model$nrows,
        ncol=null.model$ncols
    )

    #compute the null counts for each region
    null.counts <- list()

    for(mask in mask_files){
        
        masks[[mask]] <- get_matrix_from_fits(mask)
        complimentary.mask <- complimentary.mask + masks[[mask]]$data_mat
        xi.distribution.region[[mask]] <- empty_datastructure_xi_distribution(images,n_iter)
    }

    #add a mask that excludes all the regions
    complimentary.mask <- (1- - complimentary.mask/complimentary.mask)
    complimentary.mask[is.nan(complimentary.mask)] <- 1

    masks[["complimentary"]] <- list(
        data_mat=complimentary.mask,
        nrows=null.model$nrows,
        ncols=null.model$ncols
    )

    mask_files <- c(mask_files,"complimentary")

    xi.distribution.region[['complimentary']] <- empty_datastructure_xi_distribution(images,n_iter)


    for(mask in mask_files){
        null.counts[[mask]] <- sum(null.model$data_mat * masks[[mask]]$data_mat)
    }

    
    #dimension of the image 
    imsize <- null.model$nrows
    log.xi = 0

    #iterate through each image and gather the distribution of Xi for all the masks
    for(out.imagename in images){

        #read the image
        if(!file.exists(out.imagename)){
            print('File %s does not exist' %--% image_file)
            next
        }
        #read the output image; TODO: figure out an efficient way to get output from LIRA
        image.df <- data.matrix(read.table(out.imagename))

        #for each LIRA draw in this image, compute Xi for each mask

        for(i in 1:n_iter){

            #loop over the masks
            for(mask.name in mask_files){

                #the images are joined sideways in the LIRA draws,so transpose it
                im.counts <- sum(t(image.df[((i-1)*imsize+1):(i*imsize),]
                ) * masks[[mask.name]]$data_mat)
                
                if(im.counts==0) {
                    log.xi <- -8.5
                }
                else {
                    log.xi <- log10(im.counts/(im.counts + null.counts[[mask]]))
                }
                xi.distribution.region[[mask.name]][[out.imagename]][[i]] <- log.xi
            }
        }

    }
    return(list(output=xi.distribution.region,status=generate.status(0,'')))

}

save_distribution_xi_t <- function(xi.distribution,region.name,out.dir){
    #take the xi distributions of each image and compute an upper limit on the p value
    #for the sake of consistency the qunatile will be estimated using the quantile function and also the quantile.density function

    #The first output image the observed one and the rest are simulated images

    total.images <- length(xi.distribution)
    total.lira.draws <- length(xi.distribution[[1]])

    #print(total.images)
    #print(total.lira.draws)
    observed.xi.distribution <- xi.distribution[[1]]
    #this would be the mean distribution of Xi for the simulated images
    simulated.xi.distribution.all <- unlist(xi.distribution[2:total.images], use.names=F)

    gamma=0.005

    #compute c using quantile
    c <- quantile(10^simulated.xi.distribution.all, 1-gamma)
    
    #compute c using quantile.density
    c2 <- spatstat::quantile.density(density(10^(simulated.xi.distribution.all)),1-gamma)

    t_c.yobs <- 1/total.lira.draws * sum(10^observed.xi.distribution>=c)

    t_c2.yobs <- 1/total.lira.draws * sum(10^observed.xi.distribution>=c2)

    if(t_c.yobs<=gamma) t_c.yobs=gamma
    if(t_c2.yobs<=gamma) t_c2.yobs=gamma

    p.upper_lim <- gamma/t_c.yobs
    p.upper_lim2 <- gamma/t_c2.yobs

    #generate a quantile density plot comparing the observed and simulated image xi distributions.

    xi.distribution.groups <- rep(1:total.images,rep(total.lira.draws,total.images))
    xi.distribution.groups <- c(xi.distribution.groups,rep((total.images+1),(total.images-1)*total.lira.draws))

    colors <- rep("gray",total.images)
    colors[[1]] <- "cyan" #for the observed image
    colors <- c(colors,"black") #for the mean distribution

    line.types <- rep(2,total.images+1)
    line.types[[1]] <- 1 # for the observed image
    line.types[[total.images+1]] <- 1 #for the mean distribution

    line.widths <- rep(1,total.images+1)
    line.widths[[total.images+1]] <- 2 #for the mean distribution

    ngrid=200

    #TODO: Let the user add the customization
    pdf(file.path(out.dir,paste(region.name,".pdf",sep="")),width=12,height=4.25)
        sm::sm.density.compare(
            unlist(c(observed.xi.distribution,simulated.xi.distribution.all,simulated.xi.distribution.all)),
            xi.distribution.groups,
            col=colors,
            lty=line.types,
            lwd=line.widths,
            ngrid=ngrid,
            xlab="Log10(Xi)"
        )
    dev.off()


    #also dump the data to a file
    cat(c(observed.xi.distribution,simulated.xi.distribution.all,simulated.xi.distribution.all),file=file.path(out.dir,paste(region.name,"_all_dist.out",sep="")))

    cat(xi.distribution.groups,file=file.path(out.dir,paste(region.name,"_groups.out",sep="")))

    return(list(
        p_value=p.upper_lim
        ,p_value2 = p.upper_lim2
        ,status=generate.status(0,'')    
    ))


}

compute_and_save_xi_and_p_ul_t = function(results,config,masks){

    region.names <- masks$region_names
    masks <-masks$mask_files

    #generate the distributions of Xi for each region
    xi.distributions <- generate_distribution_xi_t(
        get_output_images(results),
        masks,
        config$null_file,
        config$max_iter
    )

    if(xi.distributions$status$status_code != 0){
        return(xi.distributions$status)
    }

    region.names = append(region.names,"C")
    masks <- append(masks,"complimentary")
    i <- 0

    p.values = c()
    p.values2 = c()

    for(mask in masks){
        i <- i+1
        xi.processed <- save_distribution_xi_t(
            xi.distributions$output[[mask]],
            region.names[i],
            config$output_dir
        )

        p.values <- round(c(p.values,xi.processed$p_value),digits=3)
        p.values2 <- round(c(p.values2,xi.processed$p_value2),digits=3)
    }
    cat('Done\n')
    print(knitr::kable(cbind(
        region.names,
        p.values,
        p.values2
    )))
}

compute_and_save_xi_and_p.ul= function(results,config,masks){
    output_images = get_output_images(results)
    null_file = config$null_file
    mask_files = masks$mask_files
    region_names = masks$region_names
    n_iter = config$max_iter
    thin = config$thin

    payloads = list()

    #weird, the type shows up as char and length(mask_files) gives a vector!!
    n_masks = 0
    for(mask in mask_files) n_masks = n_masks + 1

    for(i in 1:n_masks){
        payloads[[i]]=list(images=output_images
            ,mask_file=mask_files[[i]]
            ,null_file=null_file
            ,n_iter=n_iter
            ,thin=thin)
    }

    cluster = parallel::makeCluster(1
        #config$n_cores
                    #,outfile=""
                    )
    parallel::clusterEvalQ(cluster,{
        library(FITSio)
        library(spatstat)
        library(sm)
        library(reticulate)
    })
    
    
    results = parallel::parLapply(cluster,payloads,generate.distribution_xi.wrapper)


    cat('Computing the p-values\n')
    cat('---------------------------------','\n')
    cat("Region","\t","p1","\t","p2","\n")
    for(i in 1:n_masks){
        xi.processed = save_distributions_xi(results[[i]]$output,region_names[[i]],config$output_dir)

        cat(region_names[[i]]
        ,"\t",round(xi.processed$p_value,3)
        ,"\t",round(xi.processed$p_value2,3)
        ,"\n")
    }
    cat('---------------------------------','\n')

    parallel::stopCluster(cluster)
}

