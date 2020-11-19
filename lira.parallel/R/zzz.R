syspy <- NULL
astropy <- NULL
ciao.runtool_py <- NULL

.onLoad <- function(libname, pkgname) {
    #env.ascds_contrib_path <- Sys.getenv("ASCDS_CONTRIB")
    #env.path <- Sys.getenv("PATH")
    #env.ascds_lib_path <- Sys.getenv("ASCDS_LIB") 
    #Sys.setenv(
    #  PATH=paste(
    #    env.path,
    #    file.path(env.ascds_contrib_path,"lib/python3.5/site-packages"),
    #    file.path(env.ascds_lib_path,"python3.5/site-packages"),
    #    file.path(env.ascds_lib_path,"python3.5/site-packages/paramio"),
    #    sep=":"
    #    )
    #  )
    #print(Sys.getenv("PATH"))
  #reticulate::use_python("/home/batman/software/ciao/ciao-4.11/bin/python")
  syspy <<- reticulate::import("sys", delay_load = TRUE,convert=FALSE)
  astropy <<- reticulate::import("astropy",delay_load=TRUE,convert=FALSE)
  ciao.runtool_py <<- reticulate::import("ciao_contrib.runtool",delay_load=TRUE,convert=FALSE)
}