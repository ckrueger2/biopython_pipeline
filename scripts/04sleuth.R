library(sleuth)
library(data.table)
library(dplyr)
library(argparse)

#set up directory argument parsers
parser <- ArgumentParser()
parser$add_argument("--dir", help = "directory")
args <- parser$parse_args()

#read in table of samples, conditions, and kallisto file paths
input_file <- paste0(args$dir, "/sleuth_input.txt")

#assign input file
input <- read.table(input_file, header=TRUE)

#initialize slueth object
sleuth_obj <- sleuth_prep(input)

#fit a model comparing the two conditions
sleuth_obj <- sleuth_fit(sleuth_obj, ~condition, "full")

#fit the reduced model to compare in the likelihood ratio test
sleuth_obj <- sleuth_fit(sleuth_obj, ~1, "reduced")

#perform the likelihood ratio test for differential expression between conditions # nolint
sleuth_obj <- sleuth_lrt(sleuth_obj, "reduced", "full")

#extract the test results from the sleuth object
sleuth_table <- sleuth_results(sleuth_obj, "reduced:full", "lrt", show_all=FALSE)

#filter most significant results (FDR/qval < 0.05) and sort by pvalue
sleuth_significant <- sleuth_table %>% dplyr::filter(qval <= 0.05) %>% dplyr::arrange(pval)

#filter only desired columns
filtered_signif <- sleuth_significant %>% dplyr::select(target_id, test_stat, pval, qval)

#assign output directory for log file
log_file <- paste0(args$dir, "/PipelineProject.log")

#add spacer before adding to log file
write("", file = log_file, append = TRUE) 

#write FDR < 0.05 transcripts to file
write.table(filtered_signif, log_file, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t", append=TRUE)