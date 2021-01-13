# Analyse partitioned PRS
# Author: Michael Geaghan
# Date: 2021-01-13

# imports
library(reshape2)
library(ggplot2)
library(factoextra)

# function definitions
# loading data
name_integer_cols <- function(df, column, integers, labels) {
  # replace integer values with labels in a specified column of a data.frame
  if(length(integers) != length(labels)) {
    warning("'integers' argument of different length to 'labels' argument")
    return(NA)
  }
  if(length(integers) != length(unique(integers))) {
    warning("'integers' argument contains non-unique values")
    return(NA)
  }
  if(!is.numeric(df[[column]])) {
    warning(paste("column '", column, "' of data frame is not numeric", sep = ""))
    return(NA)
  }
  if(!is.numeric(integers)) {
    warning("'integers' argument must be numeric")
    return(NA)
  }
  new_df <- df
  for (i in 1:length(integers)) {
    new_df[column][new_df[column] == integers[i]] <- labels[i]
  }
  return(new_df)
}

relabel_pheno_df <- function(pheno, character_sex = TRUE, character_phenotype = TRUE, phenotype_labels = c('control', 'case')) {
  new_pheno <- pheno
  if(character_sex) {
    new_pheno <- name_integer_cols(new_pheno, 'sex', c(1, 2), c('male', 'female'))
  }
  if(character_phenotype) {
    new_pheno <- name_integer_cols(new_pheno, 'phenotype', c(1, 2), phenotype_labels)
  }
  return(new_pheno)
}

load_phenotypes <- function(fam_file, sep = " ", relabel_pheno = TRUE, ...) {
  # requires a PLINK fam file
  # produces a data frame with individual IDs as the row names, and two columns: sex and phenotype
  # sex: 1 = male; 2 = female
  # phenotype: 1 = control; 2 = case
  # ... = arguments to be passed to relabel_pheno_df()
  pheno <- read.delim(fam_file, sep = sep, row.names = 2, header = FALSE, stringsAsFactors = FALSE)[c(4, 5)]
  colnames(pheno) <- c("sex", "phenotype")
  pheno$sex[pheno$sex == 0] <- NA
  pheno$phenotype[pheno$phenotype %in% c(0, -9)] <- NA
  if (relabel_pheno) {
    pheno <- relabel_pheno_df(pheno, ...)
  }
  return(pheno)
}

load_prs <- function(prs_file, pheno, sep = "\t", row.names = 1) {
  # expects a tab-delimited file containing individual PRS scores, one individual per line
  # each column should correspond to a different partitioned PRS score for that individual
  # the first column is expected to contain the individual IDs
  # the pheno argument expects a phenotype data.frame, as loaded by load_phenotypes()
  # individual IDs should correspond to those in pheno
  prs <- read.delim(prs_file, sep = sep, row.names = row.names, stringsAsFactors = FALSE)
  prs_long <- prs
  prs_long$ID <- rownames(prs_long)
  prs_long$Pheno <- pheno[prs_long$ID, "phenotype"]
  prs_long$Sex <- pheno[prs_long$ID, "sex"]
  new_df <- melt(prs_long, id.vars = c("ID", "Pheno", "Sex"))
  colnames(new_df) <- c("ID", "Phenotype", "Sex", "PRS.Partition", "PRS")
  return(new_df)
}

# plotting PRS histograms
prs_hist_cc <- function(df, title = "PRS Histogram - Case vs Control") {
  return(ggplot(df, aes(x = PRS, fill = Phenotype)) +
           geom_histogram(color = "black", position = "identity", alpha = 0.5) +
           theme(legend.position = "bottom") +
           ggtitle(title) +
           theme(plot.title = element_text(hjust = 0.5)))
}

prs_hist_cc_partition <- function(df, partition, title = "PRS Histogram - Case vs Control") {
  return(ggplot(df[df$PRS.Partition == partition,], aes(x = PRS, fill = Phenotype)) +
           geom_histogram(color = "black", position = "identity", alpha = 0.5) +
           theme(legend.position = "bottom") +
           ggtitle(title) +
           theme(plot.title = element_text(hjust = 0.5)))
}

prs_hist_partition_merge <- function(df, phenotype = NA, title = "PRS Histogram - Case vs Control") {
  if(!is.na(phenotype)) {
    new_df <- df[df$Phenotype == phenotype,]
  } else {
    new_df <- df
  }
  return(ggplot(new_df, aes(x = PRS, fill = PRS.Partition)) +
           geom_histogram(color = "black", position = "identity", alpha = 0.5) +
           theme(legend.position = "bottom") +
           ggtitle(title) +
           theme(plot.title = element_text(hjust = 0.5)))
}

prs_violin_cc_partition <- function(df, title = "PRS Violin Plots - Case vs Control") {
  return(ggplot(df) +
           geom_violin(aes(x = Phenotype, y = PRS, fill = Phenotype)) +
           facet_wrap(~ PRS.Partition, ncol = 4, scales = "free_y"))
}
