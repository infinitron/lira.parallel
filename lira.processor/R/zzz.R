numpy <- NULL
syspy <- NULL
.onLoad <- function(libname, pkgname) {
  numpy <<- reticulate::import("numpy", delay_load = TRUE,convert=FALSE)
  syspy <<- reticulate::import("sys", delay_load = TRUE,convert=FALSE)
}