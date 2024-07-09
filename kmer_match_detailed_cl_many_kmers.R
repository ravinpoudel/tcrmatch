#Command line for `kmer_match_detailed.`
#
#Rscript kmer_match_detailed_cl.R Set1_file Set2_file output_file mismatch_max
#
#It takes as input `Set1_file` and `Set2_file`, which contains a list of kmers for set1 and set2, respectively.
#It output matches to `output_file` and summary to `output_file + .summary`. 
#The number of mismatches to be considered is given by mismatch_max.
#
#
#Apr/05/2022
#
#library("rstudioapi") 
#library(here)
#cur_dir = dirname(getSourceEditorContext()$path )

#if(!require("this.path",quietly = T)){
# install.packages("this.path")
#}

#cur_dir = dirname( this.path::this.path())
cur_dir=getwd()
#cur_dir = dir((sys.frame(1)$ofile))
PROJECT_PATH = sprintf("%s/../",cur_dir)

source("code/kmer_match_detailed.R")
`%notin%` <- Negate(`%in%`)

args = commandArgs(trailingOnly = T);

set1_file = args[[1]]
set2_file = args[[2]]
output_file_for_crags = args[[3]]
output_file_for_summary = args[[4]]
mismatch_max = as.numeric(args[[5]])
output_file_for_crags_only = args[[6]]

# set1_file <-"/Users/ravinpoudel/Documents/FY2022_Q2_BIOP_pangenome_netmhcpan/input_data/coh_A_neoantigens/P29_neoans_unique.txt"
# set2_file <- "/Users/ravinpoudel/Documents/FY2022_Q2_BIOP_pangenome_netmhcpan/output/consolidated_mde_by_patients_and_bugs_at_genus/A_P29_Bacteroides_unique_mde.txt"
# output_file_for_crags <- "output/test.txt"
# output_file_for_summary <- "output/output_file_for_summary.txt"
# mismatch_max = 1
# output_file_for_crags_only ="output/only_crags.txt"

######################


message("Processing set1_files:", set1_file)
message("Processing set2_files:", set2_file)
message("Output file:", output_file_for_crags)
message("Output file:", output_file_for_summary)
message("output_file_for_crags_only: ", output_file_for_crags_only)
# set1_file = "~/temp1p1"
# set2_file = "~/temp2p1"

set1 = fread(set1_file,header = F)
set2 = fread(set2_file,header = F)



set1 = unique(set1[[1]])
set2 = unique(set2[[1]])

set1 <- set1[!set1 %in% c("Peptide")]
set2 <- set2[!set2 %in% c("Peptide")]



kmers_sizes_set1 <- sort(unique(nchar(set1)))
message("Kmer sizes in set1_file: ", paste(as.character(kmers_sizes_set1), sep="' '", collapse=","))
kmers_sizes_set2 <- sort(unique(nchar(set2)))
message("Kmer sizes in set1_file: ", paste(as.character(kmers_sizes_set2), sep="' '", collapse=","))


temp_dir <- tempdir(tempfile())

sprintf("Temp dir for storage: %s", temp_dir)


## write into temp dir
# set1
fwrite(data.table(set1[nchar(set1)==8]), paste0(temp_dir,"/","f1_8mers.csv"), col.names=FALSE)
fwrite(data.table(set1[nchar(set1)==9]), paste0(temp_dir,"/","f1_9mers.csv"), col.names=FALSE)
fwrite(data.table(set1[nchar(set1)==10]), paste0(temp_dir,"/","f1_10mers.csv"), col.names=FALSE)
fwrite(data.table(set1[nchar(set1)==11]), paste0(temp_dir,"/","f1_11mers.csv"), col.names=FALSE)

# set2
fwrite(data.table(set2[nchar(set2)==8]), paste0(temp_dir,"/","f2_8mers.csv"), col.names=FALSE)
fwrite(data.table(set2[nchar(set2)==9]), paste0(temp_dir,"/","f2_9mers.csv"), col.names=FALSE)
fwrite(data.table(set2[nchar(set2)==10]), paste0(temp_dir,"/","f2_10mers.csv"), col.names=FALSE)
fwrite(data.table(set2[nchar(set2)==11]), paste0(temp_dir,"/","f2_11mers.csv"), col.names=FALSE)




##########################
#print("test")

kmer_match_detailed(set1_file = paste0(temp_dir,"/","f1_11mers.csv"),
                    set2_file = paste0(temp_dir,"/","f2_11mers.csv"),
                    output_file =  paste0(temp_dir,"/","f1_f2_11mers_kmer_match.csv"),
                    mismatch_max = mismatch_max)

kmer_match_detailed(set1_file = paste0(temp_dir,"/","f1_8mers.csv"),
                        set2_file = paste0(temp_dir,"/","f2_8mers.csv"),
                        output_file =  paste0(temp_dir,"/","f1_f2_8mers_kmer_match.csv"),
                        mismatch_max = mismatch_max)

kmer_match_detailed(set1_file = paste0(temp_dir,"/","f1_9mers.csv"),
                         set2_file = paste0(temp_dir,"/","f2_9mers.csv"),
                         output_file =  paste0(temp_dir,"/","f1_f2_9mers_kmer_match.csv"),
                         mismatch_max = mismatch_max)
kmer_match_detailed(set1_file = paste0(temp_dir,"/","f1_10mers.csv"),
                         set2_file = paste0(temp_dir,"/","f2_10mers.csv"),
                         output_file =  paste0(temp_dir,"/","f1_f2_10mers_kmer_match.csv"),
                         mismatch_max = mismatch_max)



## read in the results and consolidate
all_kmer_match_files = unlist(list(paste0(temp_dir,"/","f1_f2_8mers_kmer_match.csv"), paste0(temp_dir,"/","f1_f2_9mers_kmer_match.csv"), 
                            paste0(temp_dir,"/","f1_f2_10mers_kmer_match.csv"), paste0(temp_dir,"/","f1_f2_11mers_kmer_match.csv")))

all_files_list <- list.files(paste0(temp_dir,"/"), pattern = "*_kmer_match.csv", full.names = TRUE)

all_kmer_match_files <- all_files_list[!grepl("N_summary", all_files_list)]

all_kmer_match <- lapply(all_kmer_match_files, function(x){ fread(x, header = TRUE)})
all_kmer_match_df <- do.call(rbind, all_kmer_match)
fwrite(all_kmer_match_df,output_file_for_crags)
# 
# all_kmer_match_summary_files <- unlist(list(paste0(temp_dir,"/","f1_f2_8mers_kmer_match.csv.N_summary"), paste0(temp_dir,"/","f1_f2_9mers_kmer_match.csv.N_summary"), 
#                                             paste0(temp_dir,"/","f1_f2_10mers_kmer_match.csv.N_summary"), paste0(temp_dir,"/","f1_f2_11mers_kmer_match.csv.N_summary")))

all_kmer_match_summary_files <- all_files_list[grepl("N_summary", all_files_list)]
all_kmer_match_summary <- lapply(all_kmer_match_summary_files, function(x){ fread(x, header = TRUE)})
all_kmer_match_summary_df <- do.call(rbind, all_kmer_match_summary)
fwrite(all_kmer_match_summary_df,output_file_for_summary)

library(tidyverse)
crags_only_df <- all_kmer_match_df %>% filter(mismatch_size==0 | mismatch_size==1)


# this will take both kmer.x -- microbome sides crags
crags_only <- unique(crags_only_df$kmer.x)

fwrite(data.table(crags_only),output_file_for_crags_only , col.names =FALSE)

# delete tempdiroutput_file_for_crags_only
unlink(temp_dir, recursive = TRUE)



