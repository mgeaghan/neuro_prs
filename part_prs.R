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
prs_hist <- function(g, title = "PRS Histogram", color = "black", position = "identity", alpha = 0.5, legend.position = "bottom", hjust = 0.5) {
  return(g + geom_histogram(color = color, position = position, alpha = alpha) +
           theme(legend.position = legend.position) +
           ggtitle(title) +
           theme(plot.title = element_text(hjust = hjust)))
}

prs_hist_cc <- function(df, title = "PRS Histogram - Case vs Control") {
  return(prs_hist(ggplot(df, aes(x = PRS, fill = Phenotype)), title = title))
}

prs_hist_cc_partition <- function(df, partition, title = "PRS Histogram - Case vs Control") {
  return(prs_hist(ggplot(df[df$PRS.Partition == partition,], aes(x = PRS, fill = Phenotype)), title = title))
}

prs_hist_partition_merge <- function(df, phenotype = NA, title = "PRS Histogram - Case vs Control") {
  if(!is.na(phenotype)) {
    new_df <- df[df$Phenotype == phenotype,]
  } else {
    new_df <- df
  }
  return(prs_hist(ggplot(new_df, aes(x = PRS, fill = PRS.Partition)), title = title))
}

prs_violin_cc_partition <- function(df, title = "PRS Violin Plots - Case vs Control") {
  return(ggplot(df) +
           geom_violin(aes(x = Phenotype, y = PRS, fill = Phenotype)) +
           facet_wrap(~ PRS.Partition, ncol = 4, scales = "free_y"))
}

# label and count extreme PRS values
prs_extremes <- function(df, frac, subsets = NA) {
  # given a fraction, identify and label each PRS as "TOP", "BOTTOM", or "MIDDLE" depending on whether it is in the top, bottom, or middle fraction of PRS
  # PRS quantiles are calculated relative to the groups in the given column
  # additionally, add three count columns, giving the total number of PRS within the top, bottom, and middle fractions of scores
  new_df <- df
  new_df$PRS.Partition <- as.character(new_df$PRS.Partition)
  if(!is.na(subsets)) {
    new_subsets <- subsets[!(subsets == "PRS.Partition")]
    if (length(new_subsets) == 0) {
      new_subsets <- NA
      quantile_col <- "PRS.Quantile.Partition"
      quantile_count_col <- "PRS.Quantile.Count.Paritition"
    } else {
      for(s in new_subsets) {
        new_df[[s]] <- as.character(new_df[[s]])
      }
      quantile_col <- paste("PRS.Quantile.Partition", paste(new_subsets, collapse = "."), sep = ".")
      quantile_count_col <- paste("PRS.Quantile.Count.Paritition", paste(new_subsets, collapse = "."), sep = ".")
    }
  } else {
    new_subsets <- NA
    quantile_col <- "PRS.Quantile.Partition"
    quantile_count_col <- "PRS.Quantile.Count.Paritition"
  }
  cols <- list()
  for(i in 1:length(colnames(new_df))) {
    cols[[colnames(new_df)[i]]] <- i
  }
  quantiles.counts <- do.call(rbind, apply(new_df, 1, function(x) {
    part <- as.character(x[cols$PRS.Partition])
    prs <- as.numeric(x[cols$PRS])
    tmp_df <- new_df[new_df$PRS.Partition == part,]
    if(!is.na(new_subsets)) {
      for(i in 1:length(new_subsets)) {
        s <- new_subsets[i]
        group <- as.character(x[cols[[s]]])
        tmp_df <- tmp_df[tmp_df[[s]] == group,]
      }
    }
    pctls <- quantile(tmp_df$PRS, c(frac, 1 - frac))
    pctl_bottom <- pctls[1]
    pctl_top <- pctls[2]
    group.count <- data.frame(group = character(1), count = numeric(1))
    if (prs <= pctl_bottom) {
      group.count$group <- "BOTTOM"
      group.count$count <- sum(tmp_df$PRS <= pctl_bottom)
    } else if (prs > pctl_top) {
      group.count$group <- "TOP"
      group.count$count <- sum(tmp_df$PRS > pctl_top)
    } else {
      group.count$group <- "MIDDLE"
      group.count$count <- sum(tmp_df$PRS > pctl_bottom & tmp_df$PRS <= pctl_top)
    }
    return(group.count)
  }))
  colnames(quantiles.counts) <- c(quantile_col, quantile_count_col)
  new_df <- cbind(new_df, quantiles.counts)
  return(new_df)
}
