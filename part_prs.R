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

# label extreme PRS values
prs_extremes <- function(df, frac, subsets = NA) {
  # given a fraction, identify and label each PRS as "TOP", "BOTTOM", or "MIDDLE" depending on whether it is in the top, bottom, or middle fraction of PRS
  # PRS quantiles are calculated relative to the PRS partition group and any additional groups in the given subset columns
  new_df <- df
  new_df$PRS.Partition <- as.character(new_df$PRS.Partition)
  if(!is.na(subsets)) {
    new_subsets <- subsets[!(subsets == "PRS.Partition")]
    if (length(new_subsets) == 0) {
      new_subsets <- NA
      quantile_col <- "PRS.Quantile.Partition"
    } else {
      for(s in new_subsets) {
        new_df[[s]] <- as.character(new_df[[s]])
      }
      quantile_col <- paste("PRS.Quantile.Partition", paste(new_subsets, collapse = "."), sep = ".")
    }
  } else {
    new_subsets <- NA
    quantile_col <- "PRS.Quantile.Partition"
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
    group.count <- data.frame(group = character(1))
    if (prs <= pctl_bottom) {
      group.count$group <- "BOTTOM"
    } else if (prs > pctl_top) {
      group.count$group <- "TOP"
    } else {
      group.count$group <- "MIDDLE"
    }
    return(group.count)
  }))
  colnames(quantiles.counts) <- c(quantile_col)
  new_df <- cbind(new_df, quantiles.counts)
  return(structure(new_df, subsets = new_subsets, quantile_col = quantile_col))
}

prs_count_extremes <- function(df, subsets, long = TRUE) {
  # takes a data.frame produced by prs_extremes()
  # calculates the proportion of PRS in the top, bottom, and middle ranges for each group in the given subset columns
  if(is.na(subsets)) {
    warning("Must supply a vector of subsets")
    return(NA)
  }
  if(!is.na(attributes(df)$subsets)) {
    new_subsets <- subsets[!(subsets %in% attributes(df)$subsets)]
    if(is.na(new_subsets) || length(new_subsets) == 0) {
      warning("Must supply a valid vector of subsets")
      return(NA)
    }
    extreme_cols <- unique(c("PRS.Partition", attributes(df)$subsets))
  } else {
    new_subsets <- subsets
    extreme_cols <- "PRS.Partition"
  }
  subset_groups <- unique(df[new_subsets])
  extreme_groups <- unique(df[extreme_cols])
  subset_groups_ncol <- dim(subset_groups)[2]
  extreme_groups_ncol <- dim(extreme_groups)[2]
  ret <- (do.call(rbind, apply(extreme_groups, 1, function(x) {
    e_groups <- list()
    for(i in 1:extreme_groups_ncol) {
      e_groups[[i]] <- as.character(x[i])
    }
    names(e_groups) <- colnames(extreme_groups)
    new_df_e <- df[df[[names(e_groups)[1]]] == e_groups[[1]],]
    if(extreme_groups_ncol > 1) {
      for(idx in 2:extreme_groups_ncol) {
        new_df_e <- new_df_e[new_df_e[[names(e_groups)[idx]]] == e_groups[[idx]],]
      }
    }
    return(do.call(rbind, apply(subset_groups, 1, function(y) {
      s_groups <- list()
      for(j in 1:subset_groups_ncol) {
        s_groups[[j]] <- as.character(y[j])
      }
      names(s_groups) <- colnames(subset_groups)
      new_df_s <- new_df_e[new_df_e[[names(s_groups)[1]]] == s_groups[[1]],]
      if(subset_groups_ncol > 1) {
        for(jdx in 2:subset_groups_ncol) {
          new_df_s <- new_df_s[new_df_s[[names(s_groups)[jdx]]] == s_groups[[jdx]],]
        }
      }
      # count top, bottom and middle in new_df_s
      count_top <- sum(new_df_s[[attributes(df)$quantile_col]] == "TOP")
      count_bottom <- sum(new_df_s[[attributes(df)$quantile_col]] == "BOTTOM")
      count_middle <- sum(new_df_s[[attributes(df)$quantile_col]] == "MIDDLE")
      # count total in new_df_e
      count_total <- dim(new_df_e)[1]
      prop_top <- count_top/count_total
      prop_bottom <- count_bottom/count_total
      prop_middle <- count_middle/count_total
      return(data.frame(e_groups, s_groups,
                        count_top = count_top,
                        count_bottom = count_bottom,
                        count_middle = count_middle,
                        count_total = count_total,
                        prop_top = prop_top,
                        prop_bottom = prop_bottom,
                        prop_middle = prop_middle))
    })))
  })))
  if(long) {
    df_count <- ret[!grepl("prop_", colnames(ret))]
    df_prop <- ret[!grepl("count_", colnames(ret))]
    df_count_id_vars <- colnames(df_count)[(!grepl("count_", colnames(df_count))) | grepl("count_total", colnames(df_count))]
    df_prop_id_vars <- colnames(df_prop)[!grepl("prop_", colnames(df_prop))]
    df_count_long <- melt(df_count, id.vars = df_count_id_vars)
    df_prop_long <- melt(df_prop, id.vars = df_prop_id_vars)
    df_count_long$variable <- as.character(df_count_long$variable)
    df_prop_long$variable <- as.character(df_prop_long$variable)
    df_count_long$variable[df_count_long$variable == "count_top"] <- "TOP"
    df_count_long$variable[df_count_long$variable == "count_bottom"] <- "BOTTOM"
    df_count_long$variable[df_count_long$variable == "count_middle"] <- "MIDDLE"
    df_prop_long$variable[df_prop_long$variable == "prop_top"] <- "TOP"
    df_prop_long$variable[df_prop_long$variable == "prop_bottom"] <- "BOTTOM"
    df_prop_long$variable[df_prop_long$variable == "prop_middle"] <- "MIDDLE"
    colnames(df_count_long)[colnames(df_count_long) %in% c("variable", "value")] <- c("PRS.Count.Range", "Count")
    colnames(df_prop_long)[colnames(df_prop_long) %in% c("variable", "value")] <- c("PRS.Proportion.Range", "Proportion")
    common_cols <- colnames(df_prop_long) %in% colnames(df_count_long)
    ret <- cbind(df_count_long, df_prop_long[!common_cols])
    ret <- ret[colnames(ret) != "PRS.Proportion.Range"]
    colnames(ret)[colnames(ret) == "PRS.Count.Range"] <- attributes(df)$quantile_col
  }
  return(structure(ret, subsets = new_subsets, extreme_subsets = extreme_cols))
}
