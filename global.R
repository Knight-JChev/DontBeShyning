#' *------------------------------------------------------*
#' Author: Julien Chevreau
#' Mail: julien.chevreau(at)univ-rouen.fr
#' Date: 08/10/25
#' Affiliation: University of Rouen Normandie
#' Code: This script lists and install all dependencies needed 
#' to run the associated shiny application
#' NB: From a fresh install, you might need to update and upgrade and
#' install curl from libcurl4-openssl-dev and openSSL from libssl-dev
#' *------------------------------------------------------*
list_of_dependencies = c("shiny", "bslib", "bs4Dash", "fresh", 
                         "plotly", "DT")

for (package in list_of_dependencies){
  if (!require(package, character.only = T)){
    install.packages(package)
    library(package, character.only = T)
  }
}
