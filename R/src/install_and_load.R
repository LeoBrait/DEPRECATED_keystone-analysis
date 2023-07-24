#' @title Install and load libraries
#' @description This function verifies the installed libs them installs
#' the absent ones, load them and warns you if the installed version is
#' different from the specified.
#' @param libs is a named chr of libraries to install and load. The names
#' are the library names and the values are the ideal version.
#' @param loc is the location where the libraries will be installed.
#' set to NULL to install in the default R library directory.
#' @warning The function DOESN'T install the libraries in the specified version
#' it only verifies if the installed version is the same as the specified and
#' warns you.
#' @author Brait D. Mage
install_and_load <- function(libs, loc = NULL) {

  installed_packages <- installed.packages(lib.loc = loc)

  #install
  for (lib in seq_along(libs)) {
    if (!names(libs)[lib] %in% installed.packages(lib.loc = loc)) {
      install.packages(names(libs)[lib], dependencies = TRUE, lib = loc)
    }
  }

  #load
   for (lib in seq_along(libs)) {
    require(names(libs)[lib], character.only = TRUE, lib.loc = loc)
    }

# Version Warnings
installed_packages <- installed.packages(lib.loc = loc)
unmatched_libs <- libs[!libs %in% installed_packages[, "Version"]]

  if (length(unmatched_libs) > 0) {
  print("The following packages don't match the ideal version")

  for(i in seq_along(unmatched_libs)){

    print("current version:")
    print(
      paste(
        names(unmatched_libs)[i],
        packageVersion(names(unmatched_libs)[i], lib.loc = loc),
        sep = " "))

    print("ideal version:")
    print(
      paste(
        names(unmatched_libs)[i],
        unmatched_libs[i],
        sep = " "))
    }
  print("This is just a warning. Don't mean that the script will not work.")
  }
}