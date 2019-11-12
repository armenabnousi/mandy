run_limma <- function(df, groups) {
  f <- factor(groups)
  design <- model.matrix(~f)
  
  qnorm_columns <- grep("_qnorm", colnames(df), value = TRUE)
  df_select <- select(df, qnorm_columns)

  fit <- lmFit(df_select, design)
  
  result = tryCatch({
    message("running ebayes")
    fit <- eBayes(fit, trend = FALSE)
    result <- topTable(fit, number = nrow(df), coef=ncol(design), sort.by = "none")
  }, warning = function(e){
    message(e)
    message("Above warning received while running eBayes")
    fit <- eBayes(fit, trend = FALSE)
    result <- topTable(fit, number = nrow(df), coef=ncol(design), sort.by = "none")
  }, error = function(e) {
    message(e)
    message("Above error received while running eBayes")
    result <- data.frame(B = rep(NA, nrow(df)))
  })
  df <- cbind(df, result)
  df$fdr <- p.adjust(df$P.Value, method = 'fdr')
  return(df)
}

compute_qnorm <- function (df, rep_names) {
  lognorm_columns <- paste0(rep_names, "_log2norm")
  qnorms <- as.data.frame(normalize.quantiles(as.matrix(df[, lognorm_columns])))
  colnames(qnorms) <- paste0(rep_names, "_qnorm")
  df <- cbind(df, qnorms)
  return(df)
}

process_replicate <- function(bedpe_dir, chrom, rep_name, testables) {
  filenames <- system(paste0("ls ", bedpe_dir, "/reg_raw.chr", chrom, ".*| egrep and$\\|xor$"), intern = TRUE)
  if (length(filenames) != 2) {
    stop("for each chromosome there should be two files with .and and .or ending")
  }
  rep_reads <- rbind(cbind(fread(filenames[1], select = c("bin1_mid", "bin2_mid", "logl", "loggc", "logm", "logShortCount", "count")), 
                           data.frame(type = substr(filenames[1], nchar(filenames[1]) - 2, nchar(filenames[1])))), 
                     cbind(fread(filenames[2], select = c("bin1_mid", "bin2_mid", "logl", "loggc", "logm", "logShortCount", "count")),
                           data.frame(type = substr(filenames[2], nchar(filenames[2]) - 2, nchar(filenames[1])))))
  
  rep_reads$chr <- paste0("chr", chrom)
  rep_reads <- merge(rep_reads, chrom_testables, by.x = c("chr", "bin1_mid", "bin2_mid"), by.y = c("chr1", "start1", "start2"))
  rep_reads$expected <- 0
  rep_reads$pval <- 1
  ##for "and" and "xor" sets compute the regression and add the expected values to the dataframe
  for (type in unique(rep_reads$type)) {
    rep_reads_type <- rep_reads[rep_reads$type == type]
    fit <- vglm(count ~ logl + loggc + logm + logShortCount, family = pospoisson(), data = rep_reads_type)
    expected <- fitted(fit)
    pvals <- ppois(rep_reads_type$count, expected, lower.tail = FALSE, log.p = FALSE) / ppois(0, expected, lower.tail = FALSE, log.p = FALSE)
    rep_reads[rep_reads$type == type, "expected"] <- expected
    rep_reads[rep_reads$type == type, "pval"]$pval <- pvals
  }
  rep_reads$normalized <- rep_reads$count / rep_reads$expected
  rep_reads$log2norm <- log2(rep_reads$normalized + 1)
  
  #rename the columns by appending replicate name
  cols_to_rename <- c('type', 'count', 'expected', 'pval', 'normalized', 'log2norm')
  rep_reads <- select(rep_reads, c("chr", "bin1_mid", "bin2_mid", cols_to_rename))
  for (col in cols_to_rename) {
    colnames(rep_reads)[colnames(rep_reads) == col] <- paste0(rep_name, "_", col)
  }
  
  return(rep_reads)
}

find_testables <- function(all_interactions, chip_filenames) {
  if (length(unique(all_interactions$end1 - all_interactions$start1)) > 1 ) {
    stop("the input files have different binsize!")
  } else {
    binsize <- as.integer(all_interactions[1, "end1"] - all_interactions[1, "start1"])
  }
  
  for (chip_filename in chip_filenames) {
    message(paste("matching", nrow(all_interactions), "interactions with ChIP peaks in", chip_filename))
    chip_peaks <- bin_bed(read.csv(chip_filename, sep = "\t", header = FALSE), binsize)
    chip_peaks$V3 <- NULL
    chip_peaks$bin1_chip <- TRUE
    all_interactions <- merge(all_interactions, chip_peaks, all.x = TRUE, by.x = c("chr1", "start1"), by.y = c("V1", "V2"))
    colnames(chip_peaks)[3] <- "bin2_chip"
    all_interactions <- merge(all_interactions, chip_peaks, all.x = TRUE, by.x = c("chr1", "start2"), by.y = c("V1", "V2"))
    all_interactions[is.na(all_interactions)] <- FALSE
    all_interactions <- all_interactions[all_interactions$bin1_chip | all_interactions$bin2_chip,]
    all_interactions <- all_interactions[,c("chr1", "start1", "start2")]
    message(paste("remaining interactions:", nrow(all_interactions)))
  }
  return(all_interactions)
}


##this function binned bed file, if a line exceeds one bin, multiple lines replace it
bin_bed <- function(d, binsize) {
  d$V2 <- as.integer(as.character(d$V2)) %/% binsize * binsize
  d$V3 <- as.integer(as.character(d$V3)) %/% binsize * binsize
  mask <- which(d$V3 == d$V2)
  d_binsize <- d[mask,]
  d_binsize_exceeding <- d[-mask,]
  d_binsize_exceeding <- rbindlist(apply(d_binsize_exceeding, 1, function(entry, binsize) {
    start_points <- seq(as.integer(entry[2]), as.integer(entry[3]) + binsize, binsize)
    expanded_d <- data.frame(V1 = rep(entry[1], length(start_points)), V2 = start_points)
    return(expanded_d)
  }, binsize))
  
  d <- unique(rbind(select(d_binsize, c("V1", "V2")), select(d_binsize_exceeding, c("V1", "V2"))))
  d$V3 <- d$V2 + binsize
  return(d)
}
