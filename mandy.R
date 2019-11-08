suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(VGAM))) #for pospoisson regression
suppressWarnings(suppressMessages(library(preprocessCore))) #for quantile normalization
suppressWarnings(suppressMessages(library(limma)))
suppressWarnings(suppressMessages(library(timereg))) #for quantile-cutting by distance
suppressWarnings(suppressMessages(library(yaml)))
suppressWarnings(suppressMessages(library(argparse)))

parser <- ArgumentParser()
parser$add_argument("-c", "--configfile", type="character", 
                    help="path to the YAML file for all parameters")
args <- parser$parse_args()
#inputs <- read_yaml("~/mandy.yaml")

inputs <- read_yaml(args$configfile)

#parse the yaml arguments
output_filename <- inputs[['output_filename']]
chip_filenames <- trimws(unlist(strsplit(inputs[['chip_filenames']], ",")))
bedpe_dirs <- trimws(unlist(strsplit(inputs[['bedpe_dirs']], ",")))
maps_peaks_filenames <- trimws(unlist(strsplit(inputs[['maps_peaks_filenames']], ",")))
genome <- inputs[['genome']]
groups <- trimws(unlist(strsplit(inputs[['groups']], ",")))
utils_path <- inputs[['mandy_utils_path']]

source(utils_path)

##generating list of chromosomes
if (genome == "hg") chromosomes <- 1:22 else chromosomes <- 1:19

##combine two MAPS significant interaction calls
all_interactions <- unique(do.call(rbind, lapply(maps_peaks_filenames, fread, select = c("chr1", "start1", "end1", "start2"))))

##select testable and significant binpairs from all significant interactions
testables <- find_testables(all_interactions, chip_filenames)
rm(all_interactions)


##generating replicate names to later use for columns.
##for each entry in "i" in "groups" it creates a name like "seti_repj"
groups_df <- data.table(group = groups)[, rep := rowid(group)]
rep_names <- paste0("set", groups_df$group, "_rep", groups_df$rep)

##for each chromosome, read bedpe files for all replicates, run regression, 
##get expected and normalized values
##at the end keep only testable binpairs 
testable_counts <- data.frame()
for (chrom in chromosomes) {
  message(paste("processing chromosome", chrom))
  chrom_testables <- testables[testables$chr1 == paste0("chr", chrom),] #this is constant for each chrom to extract reads
  chrom_testables_counts <- chrom_testables #this gets two new columns for each replicate: logShortCount and count
  
  for (replicate_number in seq_along(bedpe_dirs)) {
    rep_name <- rep_names[replicate_number]
    bedpe_dir <- bedpe_dirs[replicate_number]
    
    #read counts and run regression
    rep_counts <- process_replicate(bedpe_dir, chrom, rep_name)
    
    #subset to only testables
    rep_counts <- merge(rep_counts, chrom_testables, by.x = c("chr", "bin1_mid", "bin2_mid"), by.y = c("chr1", "start1", "start2"))
    
    ##append the replica data to the data from previous replicates
    chrom_testables_counts <- merge(chrom_testables_counts, rep_counts, all.x = TRUE, by.x = c("chr1", "start1", "start2"), by.y = c("chr", "bin1_mid", "bin2_mid"))
  }
  message(paste("number of processed binpairs:", nrow(chrom_testables_counts)))
  
  ##add the chrom data to the genome-wide dataframe
  testable_counts <- rbind(testable_counts, chrom_testables_counts)
  message(paste("total number of binpairs processed:", nrow(testable_counts)))
}

#check whether data from all testable binpairs have been found, if not set them to zero and throw a warning
if (nrow(testable_counts) != nrow(testables)) {
  message(paste("number of testable binpairs:", nrow(testables), "; number of binpairs with counts:", nrow(testable_counts)))
  warning("some testable binpairs were not found in any of the replicates. Setting all values for them to zero.")
  testable_counts[is.na(testable_counts)] <- 0
}

distances <- testable_counts$start2 - testable_counts$start1
testable_counts$distance_stratum <- qcut(distances, cuts = 10)

##quantile normalization
merged_qnormed <- testable_counts %>% group_by(distance_stratum) %>% do(qnorms = compute_qnorm(., rep_names))
merged_qnormed <- rbindlist(merged_qnormed$qnorms, fill = TRUE)

##run limma
merged_tested <- merged_qnormed %>% group_by(distance_stratum) %>% do(test_results = run_limma(., groups))
merged_tested <- rbindlist(merged_tested$test_results, fill = TRUE)

#write output to file
write.table(merged_tested, output_filename, sep = "\t", quote = F, row.names = F)
