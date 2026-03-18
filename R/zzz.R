# # INLA install flag.
# .onLoad <- function(libname, pkgname) {
#     if (!requireNamespace("INLA", quietly = TRUE)) {
#         stop(
#             "Package '", pkgname, "' requires INLA.\n\n",
#             "Install INLA with:\n",
#             "install.packages('INLA', repos = c(getOption('repos'), INLA = 'https://inla.r-inla-download.org/R/stable'), dependencies = TRUE)\n",
#             call. = FALSE
#         )
#     }
# }
