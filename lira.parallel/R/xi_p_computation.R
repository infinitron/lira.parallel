# Generate distributions for Xi using the output images and masks
get_matrix_from_fits <- function(file_name,fill.na=0){
    if(!file.exists(file_name)) return(list(data_mat=NULL,nrows=NULL,ncols=NULL))

    #readfits prints the filename. It will be suppressed
    invisible(utils::capture.output(data <- FITSio::readFITS(file_name)))

    #set all non-finite values to zero
    data$imDat[!is.finite(data$imDat)] <- fill.na

    data_mat = matrix(data=data$imDat
        ,nrow=data$axDat$len[1]
	    ,ncol=data$axDat$len[2])

    return(list(data_mat=data_mat,nrows=data$axDat$len[1],ncols=data$axDat$len[2]))
}

get_matrix_from_file <- function(file_name){
    if(!file.exists(file_name)) return(list(data_mat=NULL,nrows=NULL,ncols=NULL))
}

empty_datastructure_xi_distribution <- function(images,n_iter){
    ds <- list()
    for(image in images){
        ds[[image]] <- rep(0.0,n_iter)
    }
    return(ds)
}


# Compute the xi for all the iterations in each image
# This function takes list of output images, mask files, a null file, total number of iterations and the thinning factor
#  For each output image and for each mask file, it computes a distribution of Xi
#  Each output image is read once and computations are made on it

generate_distribution_xi_t <- function(images,mask_files,null_file,n_iter){
    #read all the masks into an array
    masks <- list()
    #hosts the distribution of Xi for each mask
    xi.distribution.region <- list()

    n_images <- length(images)

    #read the null file.
    null.model <- get_matrix_from_fits(null_file)

    complimentary.mask <- matrix(0,
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
    complimentary.mask <- (1- complimentary.mask/complimentary.mask)
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
    log.xi <- 0
    obs.file.flag <- 0
    avg.obs.LIRA.image <- matrix(0,
        nrow=null.model$nrows,
        ncol=null.model$ncols)

    #iterate through each image and gather the distribution of Xi for all the masks
    for(out.imagename in images){

        obs.file.flag <- obs.file.flag + 1

        #read the image
        if(!file.exists(out.imagename)){
            print('File %s does not exist' %--% out.imagename)
            next
        }
        #read the output image; TODO: figure out an efficient way to get output from LIRA
        image.df <- data.matrix(utils::read.table(out.imagename))

        #treat the NANs
        image.df[!is.finite(image.df)] <- 0

        #for each LIRA draw in this image, compute Xi for each mask

        for(i in 1:n_iter){

            image.this.iteration <- image.df[((i-1)*imsize+1):(i*imsize),]
            if(obs.file.flag==1){
                avg.obs.LIRA.image = avg.obs.LIRA.image + image.this.iteration
            }

            #loop over the masks
            for(mask.name in mask_files){

                #the images are joined sideways in the LIRA draws,so transpose it

                im.counts <- sum(t(image.this.iteration) * masks[[mask.name]]$data_mat)
                
                if(!is.finite(im.counts)){
                    im.counts <- 0
                }

                if(im.counts==0) {
                    log.xi <- -8.5
                }
                else {
                    log.xi <- log10(im.counts/(im.counts + null.counts[[mask]]))
                }
                xi.distribution.region[[mask.name]][[out.imagename]][[i]] <- log.xi
            }
        }

        if(obs.file.flag==1){
            avg.obs.LIRA.image <- avg.obs.LIRA.image/n_iter
            write_FITS_image_with_wcs(avg.obs.LIRA.image,null_file)
        }

    }
    return(list(output=xi.distribution.region,status=generate.status(0,'')))

}

#Write a matrix to a FITS file using wcs information from another file
write_FITS_image_with_wcs <- function(data,wcsfile.name){
    hdul <- astropy$io$fits$open(wcsfile.name)
    reticulate::py_set_attr(hdul[0]$header,'BITPIX',-32L)
    reticulate::py_set_attr(hdul[0],'data',data)
    hdul$writeto("avg_LIRA.fits",overwrite=T)
    hdul$close()
}

# Take the distibutions of xi, save them as pdf/mat files, and compute the upper bound on p value

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
    c <- stats::quantile(10^simulated.xi.distribution.all, 1-gamma)
    
    #compute c using quantile.density
    c2 <- spatstat::quantile.density(stats::density(10^(simulated.xi.distribution.all)),1-gamma)

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

    ngrid <- 200


    #also dump the data to a file
    cat(c(observed.xi.distribution,simulated.xi.distribution.all,simulated.xi.distribution.all),file=file.path(out.dir,paste(region.name,"_all_dist.out",sep="")))

    cat(xi.distribution.groups,file=file.path(out.dir,paste(region.name,"_groups.out",sep="")))

    #TODO: Let the user add the customization
    tryCatch({
        grDevices::pdf(file.path(out.dir,paste(region.name,".pdf",sep="")),width=12,height=4.25)
            sm::sm.density.compare(
                unlist(c(observed.xi.distribution,simulated.xi.distribution.all,simulated.xi.distribution.all)),
                xi.distribution.groups,
                col=colors,
                lty=line.types,
                lwd=line.widths,
                ngrid=ngrid,
                xlab="Log10(Xi)"
            )
        grDevices::dev.off()
    }, warning=function(warn){
        #do nothing
    },error=function(err){
        
        cat('\nXi distribution plots will not be generated\n')
        
    },finally={
        
    })



    return(list(
        p_value=p.upper_lim
        ,p_value2 = p.upper_lim2
        ,status=generate.status(0,'')    
    ))


}

compute_and_save_xi_and_p_ul_t <- function(results,config,masks){

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

    region.names <- append(region.names,"C")
    masks <- append(masks,"complimentary")
    i <- 0

    p.values = c()
    p.values2 = c()

    cat('Generating the distributions of Xi...')
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
    cat(knitr::kable(cbind(
        region.names,
        p.values,
        p.values2
    )),file=file.path(config$output_dir,'p_ul.values.txt'))
}


