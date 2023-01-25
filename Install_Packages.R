#Install packages if not installed
installReqs <- function(package_name, bioc) {
  if (requireNamespace(package_name, quietly = TRUE) == FALSE) {
    if (bioc == FALSE)
      install.packages(package_name)
    else if (bioc == TRUE) #install using Bioconductor package manager
      BiocManager::install(package_name)
  }
}

#Install pacakges:
installReqs("BiocManager", bioc = FALSE)
installReqs("ggplot2", bioc = TRUE)
installReqs("Seurat", bioc = TRUE)
installReqs("stringr", bioc = TRUE)
installReqs("dplyr", bioc = TRUE)
installReqs("ggridges", bioc = TRUE)
installReqs("ggrepel", bioc = TRUE)
installReqs("ComplexHeatmap", bioc = TRUE)

