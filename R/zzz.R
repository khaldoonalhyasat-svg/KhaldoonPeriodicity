# Package Startup and Load Messages

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("KhaldoonPeriodicity v1.0.0 loaded successfully.")
  packageStartupMessage("Use khaldoon_test() for periodicity detection in time series.")
}

.onLoad <- function(libname, pkgname) {
  # Package initialization code if needed
  invisible()
}
