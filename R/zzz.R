#--------------------------------------------------------------------
#   zzz.R (npsp package)
#--------------------------------------------------------------------
#   .onAttach(libname, pkgname)
#
#   (c) R. Fernandez-Casal
#--------------------------------------------------------------------


#--------------------------------------------------------------------
.onAttach <- function(libname, pkgname){
#--------------------------------------------------------------------
#   pkg.info <- utils::packageDescription(pkgname, libname, fields = c("Title", "Version", "Date"))
    pkg.info <- drop( read.dcf( file = system.file("DESCRIPTION", package = "npsp"),
                      fields = c("Title", "Version", "Date") ))
    packageStartupMessage( 
      paste0(" Package npsp: ", pkg.info["Title"], ",\n"),
      paste0(" version ", pkg.info["Version"], " (built on ", pkg.info["Date"], ").\n"),
      " Copyright (C) R. Fernandez-Casal 2012-2021.\n",
      " Type `help(npsp)` for an overview of the package or\n",
      ' visit https://rubenfcasal.github.io/npsp.\n')
}
