# Generate distributions for Xi using the output images and masks
get_matrix_from_fits = function(file_name){
    if(!file.exists(file_name)) return(list(data_mat=NULL,nrows=NULL,ncols=NULL))

    data = FITSio::readFITS(file_name)
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
    print('in')
    #read each mask in to a matrix
    region_mat_obj = get_matrix_from_fits(mask_file)
    print('out')
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
        #print(dim(im_mat))
        im_mat_obj = list(data_mat=im_mat,nrows=dim(im_mat)[[1]],ncols=dim(im_mat)[[2]])


        xi.iteration_wise = c()
        #null_counts=1
        #print(im_mat_obj$nrows)
        #print(null_nrows)
        for(i in seq(1,im_mat_obj$nrows,null_nrows)){

            image_iter = im_mat_obj$data_mat[i:(i+null_nrows-1),]
            im_counts = sum(image_iter * region_mat_obj$data_mat)
            if(im_counts==0) {
                log_xi = -8.5 #some far away value so as to not to interfere with the rest of the distribution
            }else{
                log_xi = log(im_counts/(im_counts+null_counts),10)
            }
            

            xi.iteration_wise = c(xi.iteration_wise,log_xi)
        }
        #print(length(xi.iteration_wise))
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
    #the rest would contain replicas vs baselines

    dims = dim(distributions)
    nrows = dims[[1]]
    ncols = dims[[2]]

    # log transform the vector
    all_dist_vec = matrix(t(distributions),nrow=1)

    #compute the upper bound
    gamma = 0.005
    c = quantile(10^(all_dist_vec[1,(ncols+1):(nrows*ncols)]),1-gamma)

    #print('c1')
    #print(c)

    #print('c2')
    c2=spatstat::quantile.density(density(10^(all_dist_vec[1,(ncols+1):(nrows*ncols)])),1-gamma)
    #print(c2)
    

    t_c.yobs = 1/(ncols)*sum(10^(all_dist_vec[,1:ncols])>=c)
    t_c2.yobs = 1/(ncols)*sum(10^(all_dist_vec[,1:ncols])>=c2)

    #print(t_c.yobs)
    if(t_c.yobs==0) t_c.yobs=gamma
    if(t_c2.yobs==0) t_c2.yobs=gamma
    p.value.upper_lim = gamma/t_c.yobs
    p.value.upper_lim2 = gamma/t_c2.yobs


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
    #print(all_dist_vec)
    #print(length(groups))
    #print(length(all_dist_vec))


    
    #save the distribution as a pdf
    #print('saving dist')
    pdf(file.path(out_dir,paste(region_name,".pdf",sep="")), width=12,height=4.25)
        sm::sm.density.compare(all_dist_vec,groups,col=colors,lty=line_types,lwd=line_widths
                #,xlim=xlim
                #,ylim=ylim
                ,ngrid=ngrid)
    dev.off()

    #also dump the data to a file
    write(all_dist_vec,file.path(out_dir,paste(region_name,"_all_dist.out",sep="")))
    write(groups,file.path(out_dir,paste(region_name,"_groups.out",sep="")))

    return(list(
        p_value=p.value.upper_lim
        ,p_value2 = p.value.upper_lim2
        ,status=generate.status(0,'')    
    ))
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

    cluster = parallel::makeCluster(config$n_cores
                    #,outfile=""
                    )
    parallel::clusterEvalQ(cluster,{
        library(FITSio)
        library(spatstat)
        library(sm)
        library(reticulate)
    })
    
    
    results = parallel::parLapply(cluster,payloads,generate.distribution_xi.wrapper)


    #iterate over the results and get the upper limits
    cat('Computing the p-values\n')
    #print('Displaying upper bounds on p values')
    cat('---------------------------------','\n')
    cat("Region","\t","p1","\t","p2","\n")
    for(i in 1:n_masks){
        #print('total null counts')
        #print(results[[i]]$total_null_counts)
        #print(dim(results[[i]]$output))
        xi.processed = save_distributions_xi(results[[i]]$output,region_names[[i]],config$output_dir)

        cat(region_names[[i]]
        ,"\t",round(xi.processed$p_value,3)
        ,"\t",round(xi.processed$p_value2,3)
        ,"\n")
    }
    cat('---------------------------------','\n')

    parallel::stopCluster(cluster)
}

