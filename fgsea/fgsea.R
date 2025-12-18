rm(list = ls())
set.seed(7)
# install.packages("dplyr",lib="~/myRlibs",repos='https://cran.r-project.org')
library(dplyr, lib = '~/myRlibs')
# install.packages("this.path",lib="~/myRlibs",repos='https://cran.r-project.org')
library(this.path, lib = '~/myRlibs')
script_dir <- this.dir()
setwd(script_dir)
# install.packages("devtools",lib="~/myRlibs",repos='https://cran.r-project.org')
# library(devtools,lib='~/myRlibs')
# install.packages("farver",lib="~/myRlibs",repos='https://cran.r-project.org')
library(farver, lib = '~/myRlibs')
# install.packages("RColorBrewer",lib="~/myRlibs",repos='https://cran.r-project.org')
library(RColorBrewer, lib = '~/myRlibs')
# install.packages("desc",lib="~/myRlibs",repos='https://cran.r-project.org')
# library(desc,lib='~/myRlibs')
# install.packages("ps",lib="~/myRlibs",repos='https://cran.r-project.org')
# library(ps,lib='~/myRlibs')
# install.packages("callr",lib="~/myRlibs",repos='https://cran.r-project.org')
# library(callr,lib='~/myRlibs')
# install_github("ctlab/fgsea",lib="~/myRlibs")
library(fgsea, lib = '~/myRlibs')
# install.packages("data.table",lib="~/myRlibs",repos='https://cran.r-project.org')
library(data.table, lib = "~/myRlibs")
args <- commandArgs(trailingOnly = T)
U_protein_coding <- read.csv(args[2], row.names = 1)
names(U_protein_coding) <- colnames(U_protein_coding)
pwy_collections <- c(paste(args[1], "_test", sep = '')) # "h", "c2", "c6", "c7", "c8")
top_all_embed_pwy <- c()
# tail_all_embed_pwy <- c()
set.seed(1027)
for (j in 1:length(pwy_collections)) {
  pwy <- pwy_collections[j]
  signature <- readRDS(paste0("./data/pathways/", pwy, ".rds"))
  for (i in 1:ncol(U_protein_coding)) {
    embed_idx <- i
    chosen_embed <- U_protein_coding[, embed_idx]
    names(chosen_embed) <- rownames(U_protein_coding)
    top <- sort(chosen_embed, decreasing = TRUE)#[1:num_genes]
    # tail <- sort(chosen_embed, decreasing = FALSE)#[1:num_genes]
    # tail <- -tail
    # Run fGSEA----
    top_res <- fgsea(signature, top, minSize = 10, maxSize = 500)
    top_res$embed_idx <- i
    top_all_embed_pwy <- rbind(top_all_embed_pwy, top_res)
    # tail_res <- fgsea(signature, tail, minSize = 10, maxSize = 500)
    # tail_res$embed_idx <- i
    # tail_all_embed_pwy <- rbind(tail_all_embed_pwy, tail_res)
  }
}
top_all_embed_pwy_sig <- top_all_embed_pwy#[top_all_embed_pwy$padj < 0.05, ]
# top_all_embed_pwy_sig$type <- "top"
# tail_all_embed_pwy_sig <- tail_all_embed_pwy[tail_all_embed_pwy$padj < 0.05, ]
# tail_all_embed_pwy_sig$type <- "tail"
all_embed_pwy_sig <- top_all_embed_pwy_sig %>% as.data.frame()
all_embed_pwy_sig$leadingEdge <- as.character(all_embed_pwy_sig$leadingEdge)
write.csv(all_embed_pwy_sig, paste("./", args[3], ".csv", sep = ''))
