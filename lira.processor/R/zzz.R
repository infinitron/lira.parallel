numpy <- NULL
syspy <- NULL
astropy <- NULL
ciao.runtool_py <- NULL
ospy <- NULL
.onLoad <- function(libname, pkgname) {
  numpy <<- reticulate::import("numpy", delay_load = TRUE,convert=FALSE)
  syspy <<- reticulate::import("sys", delay_load = TRUE,convert=FALSE)
  astropy <<- reticulate::import("astropy",delay_load=TRUE,convert=FALSE)
  ciao.runtool_py <<- reticulate::import("ciao_contrib.runtool",delay_load=TRUE,convert=FALSE)
  os_py <<- reticulate::import("os",delay_load=TRUE,convert=FALSE)
  pathlib_py <<- reticulate::import("pathlib",delay_load=TRUE,convert=FALSE)
  tempfile_py <<- reticulate::import("tempfile",delay_load=TRUE,convert=FALSE)

}