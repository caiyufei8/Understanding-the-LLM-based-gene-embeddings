rm(list = ls())
set.seed(7)

library(dplyr)
# install.packages("this.path",lib="~/myRlibs",repos='https://cran.r-project.org')
library(this.path, lib = '~/myRlibs')
script_dir <- this.dir()
setwd(script_dir)
# install.packages("purrr",lib="~/myRlibs",repos='https://cran.r-project.org')
library(purrr, lib = '~/myRlibs')
# install.packages("BiocManager",lib="~/myRlibs",repos='https://cran.r-project.org')
# library(BiocManager,lib='~/myRlibs')
# BiocManager::install("msigdb",lib="~/myRlibs")
# BiocManager::install("GSEABase",lib="~/myRlibs")
# BiocManager::install("ExperimentHub",lib="~/myRlibs")
# BiocManager::install("qusage",lib="~/myRlibs")
library(qusage, lib = '~/myRlibs')
library(msigdb, lib = '~/myRlibs')
library(GSEABase, lib = '~/myRlibs')
library(ExperimentHub, lib = '~/myRlibs')
count_times <- function(dictionary, vector_sea) {
  sum(dictionary %in% vector_sea)
}
gene_sets <- c('h', 'c2')
file_names <- c('h.all.v2024.1.Hs.symbols.gmt',
                'c2.cp.v2024.1.Hs.symbols.gmt')
subset_file_names <- c(
  'c2.cp.biocarta.v2024.1.Hs.symbols.gmt',
  'c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt',
  'c2.cp.pid.v2024.1.Hs.symbols.gmt',
  'c2.cp.reactome.v2024.1.Hs.symbols.gmt',
  'c2.cp.wikipathways.v2024.1.Hs.symbols.gmt',
  'c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt'
)
original_gene_symbols <- read.csv('../original_text-embedding-3-small_filtered.csv')$X
for (i in 1:length(gene_sets))
{
  if (gene_sets[i] == 'c2') {
    gene_set <- list()
    for (subset_file in subset_file_names) {
      subset_gene_set <- read.gmt(paste('c2_subsets/', subset_file, sep = ''))
      gene_set <- c(gene_set, subset_gene_set)
    }
  }
  else {
    gene_set <- read.gmt(file_names[i])
  }
  new_gene_set <- list()
  removed <- 0
  for (key in names(gene_set)) {
    value <- gene_set[[key]]
    in_count <- count_times(value, original_gene_symbols)
    # if (in_count / length(value) >= 0.96)
    if (in_count >= 10)
    {
      new_gene_set[[key]] <- value
    } else {
      # cat("Removing", key, "with count", in_count, "\n")
      removed <- removed + 1
    }
  }
  print(paste("Removed", removed, "pathways from", gene_sets[i], "gene set."))
  save_path <- paste0("./data/pathways/", gene_sets[i], "_test.rds")
  saveRDS(new_gene_set, save_path)
}
