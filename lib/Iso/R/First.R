.onAttach <- function(lib, pkg) {
	ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
        packageStartupMessage(paste(pkg, ver))
        msg <- paste("\n     An \"infelicity\" in the function ufit() (whereby",
                     "\n     it was all too easy to conflate the location of",
                     "\n     the mode with its index in the entries of the",
                     "\n     \"x\" argument) has been corrected.  To this end,",
                     "\n     ufit() now has arguments \"lmode\" (the location",
                     "\n     of the mode), and \"imode\" (its index).  At most",
                     "\n     one of these arguments should be specified.  See",
                     "\n     the help for ufit().")
       packageStartupMessage(msg)
}
