#require(lira)
#require(FITSio)

launchLIRA<-function(obsFile,startFile=FALSE,mapFile=F, bkgFile=FALSE,psfFile=FALSE,fit.bkg.scale=T,outDir='/data/reu/kmckeough/KM_lira/outputs/',maxIter,alpha.init,thin=1,burn=0,mcmc=TRUE){

#INPUTS:
	#obsFile   - string;a 2^n x 2^n matrix which you would like to
		   #analyze; put this file in the current directory
	#startFile --initial values for residuals
	#mapFile   -string; exposure map
	#bkgFile   -null model
	#psfFile   -psf
	#outDir	   - output directory specific to machine
	#maxIter   - lira() input; the maximum number of iterations
	#alpha.init- set smoothing parameters (higher=more smoothing)
	#WARNING: Other inputs specific to lira must be coded into script below

	#NOTE: set thin=1 to avoid floating point exception and subsequent crash on Mac OS X

#OUTPUTS:
	#.out
	#.param
	#.pdf

#Extract data and format to array/matrix

#Read in FITS files
	obs <- FITSio::readFITS(obsFile)
	obsmat <- matrix(data=obs$imDat,nrow=obs$axDat$len[1],
	ncol=obs$axDat$len[2])
		
	if(startFile!= F){
		   strt <- FITSio::readFITS(startFile)
		   strtmat <- matrix(data=strt$imDat,
		   nrow=strt$axDat$len[1],
		   ncol=strt$axDat$len[2])
		   }else{
		   strtmat<-matrix(1, nrow(obsmat),
		   ncol(obsmat))}
	if(mapFile!= F){
		   map <- FITSio::readFITS(mapFile)
		   mapmat<- matrix(data=map$imDat,
		   nrow=map$axDat$len[1],
		   ncol=strt$axDat$len[2])
		   }else{
		   mapmat<-matrix(1, nrow(obsmat),
		   ncol(obsmat))}
	if(bkgFile != F){
		   bkgd <- FITSio::readFITS(bkgFile)
		   bkgdmat <- matrix(data=bkgd$imDat,
		   nrow=bkgd$axDat$len[1],
		   ncol=bkgd$axDat$len[2])
		   }else{
		   bkgdmat<-matrix(0, nrow(obsmat),
		   ncol(obsmat))}
	if(psfFile != F){
		   psf <- FITSio::readFITS(psfFile)
		   psfmat <- matrix(data=psf$imDat,
		   nrow=psf$axDat$len[1],
		   ncol=psf$axDat$len[2])
		   #psfmat <- matrix(data=psf$imDat,
		   #nrow=psf$axDat$len[1])
		   }else{
		   psfmat<-matrix(1, 1, 1)}
#Create output file names		   		 
	basename<-tools::file_path_sans_ext(obsFile)
	outsave<-file.path(outDir,paste(basename,'.out',sep=''))
	paramsave<-file.path(outDir,paste(basename,'.param',sep=''))
	pdfsave<-file.path(outDir,paste(basename,'.pdf',sep=''))
    posteriorsave<-file.path(outDir,paste(basename,'.posterior',sep=''))

#Run lira (see lira documentation for help)
	img<-lira::lira(obs.matrix=obsmat, start.matrix=strtmat,map.matrix=mapmat,
	  bkg.matrix=bkgdmat, psf.matrix=psfmat, out.file=outsave,
	  fit.bkg.scale=fit.bkg.scale,thin=thin,burn=burn,
	  param.file=paramsave,max.iter=maxIter, alpha.init=alpha.init,mcmc=mcmc)

#Write PDF ofimages
	pdf(pdfsave, width=12,height=4.25)
	  par(mfrow=c(1,3))
	  image(obsmat, xaxt="n", yaxt="n", main="Observed Data")
	  image(psfmat, xaxt="n", yaxt="n", main="Null/Best-Fit/Background Model")
	  image(img$final, xaxt="n", yaxt="n", main="Mean MultiScale of Data/Model MisMatch")
	dev.off()
	
	write(img$final,file=posteriorsave,ncolumns=obs$axDat$len[2])
	#return(img)

    return(list(out_images_file=outsave,params_file=paramsave,pdf_file=pdfsave,posterior_file=posteriorsave))

}
