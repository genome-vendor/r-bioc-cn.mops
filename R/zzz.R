# TODO: Add comment
# 
# Author: klambaue
###############################################################################


.onLoad <- function(libname, pkgname) {
	library.dynam("cn.mops", pkgname, libname)
}

.onUnload <- function(libpath)
{
	library.dynam.unload("cn.mops", libpath)
}

