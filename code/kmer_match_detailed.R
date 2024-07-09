#Kmer match detailed algoritm
#Mar/31/2022
#
#This script is a detailed version of `kmer_match.R`. 
#It not only tells whether there is a match, but it also returns the pairwise match.
#
#
#Usage:
#Rscript kmer_match_cl.R set1_file set2_file output_file
#set1_file
#   - Single column file with list of kmers from set1
#set2_file
#   - Single column file with list of kmers from set2
#output_file
# - Output file
#   - The output file will have 5 columns:
#      - match_type (no_match; perfect_match; imperfect_match)
#      - k_minus1_mer (the `k-minus1_mer` used to match)
#      - kmer.x, kmer.y (kmer from set1 and set2)
#      - letter_skipped (the letter skipped in the `k-minus1_mer` match).
#   - It also outputs a N_summary at output_file + `.N_summary`.
#   - It also accepts `superimperfect` matches.
#
#Algorithm
#For fast comparison on imperfect match, I compare all "(k-1)"mer in each set by removing each letter position from kmer.
#If match exists for any `(k-1)`-mer, I consider that kmer to contain an imperfect match 
#
#It allows to find`superimperfect` (2 mismatches) matches.
#
#

if(!require("data.table",quietly = T)){
  install.packages("data.table")
}
library(data.table);

kmer_match_detailed <- function(set1_file, set2_file, output_file=NULL, input_as_vector=F, mismatch_max=1,warning_loop_step=T){
  
  if(((file.size(set1_file)>0 & file.size(set2_file)>0))==TRUE){
    if(input_as_vector==F){
      set1 = fread(set1_file,header = F)
      set2 = fread(set2_file,header = F)
      set1 = unique(set1[[1]])
      set2 = unique(set2[[1]])
      
    }else{
      set1 = unique(set1_file);
      set2 = unique(set2_file);
    }
    
    kmer_size = unique( c(nchar(set1),nchar(set2)))
    if(length(kmer_size) != 1){
      stop("All kmer must have the same size! Try function `kmer_match_detailed_many_mers()` instead")
    }
    
    
    #Identify perfect match set:
    perfect_match_set = set1[set1 %in% set2]
    1
    
    dt_set1_and_2_match = data.table(kmer.x=c(perfect_match_set,NA),
                                     kmer.y=c(perfect_match_set,NA),
                                     k_m_less_mer=NA,
                                     letter_skipped=NA,
                                     mismatch_size=0);
    
    N_summary = data.table(perfect_match=length(perfect_match_set))
    
    #Select kmers that are not perfect match
    for(m in 1:mismatch_max){
      set1_pending = setdiff(set1,dt_set1_and_2_match$kmer.x);
      set2_pending = setdiff(set2,dt_set1_and_2_match$kmer.y);
      
      if(length(set1_pending)==0 | length(set2_pending)==0){
        N_summary[[sprintf("N_mismatch_%d_set1", m)]] = 0
        N_summary[[sprintf("N_mismatch_%d_set2", m)]] = 0
        next;
      }
      d_set_imperfect_cur = kmer_imperfect_match_general(set1 = set1_pending,
                                                         set2 = set2_pending,
                                                         m_mismatch = m,
                                                         warning_loop_step=warning_loop_step)
      
      dt_set1_and_2_match = rbind(dt_set1_and_2_match,
                                  d_set_imperfect_cur)
      
      N_summary[[sprintf("mismatch_%d_set1", m)]] = length(unique(d_set_imperfect_cur$kmer.x))
      N_summary[[sprintf("mismatch_%d_set2", m)]] = length(unique(d_set_imperfect_cur$kmer.y))
      
    }
    set1_pending = setdiff(set1,dt_set1_and_2_match$kmer.x);
    set2_pending = setdiff(set2,dt_set1_and_2_match$kmer.y);
    
    N_summary[["no_match_set1"]] = length(unique(set1_pending))
    N_summary[["no_match_set2"]] = length(unique(set2_pending))
    
    dt_set1_and_2_no_match = rbind(data.table(kmer.x=c(set1_pending,NA),
                                              kmer.y=NA,
                                              k_m_less_mer=NA,
                                              letter_skipped=NA,
                                              mismatch_size=NA),
                                   data.table(kmer.x=NA,
                                              kmer.y=c(set2_pending,NA),
                                              k_m_less_mer=NA,
                                              letter_skipped=NA,
                                              mismatch_size=NA))
    
    dt_set1_and_2_match = rbind(dt_set1_and_2_match,
                                dt_set1_and_2_no_match)
    
    dt_set1_and_2_match = dt_set1_and_2_match[!(is.na(kmer.x) & is.na(kmer.y))]
    
    if(!is.null(output_file)){
      write.table(dt_set1_and_2_match, output_file,quote = F,sep = "\t",row.names = F);
      write.table(N_summary,sprintf("%s.N_summary", output_file),quote = F,sep = "\t",row.names = F);
    }
    
    
    return(list(d_set_match=dt_set1_and_2_match,
                N_summary = N_summary,
                N_summary_relative = N_summary/sum(N_summary)));
    
    
  }}
  
  kmer_imperfect_match_2 <- function(set1,set2){
    #It run a second order imperfect match.
    #Only runs it on 
    
    #  warning("prototyping script. Input must containt non-matching sets. I.e. sets without any perfect or 1 mismatch imperfect match");
    
    set1 = unique(set1);
    set2 = unique(set2);
    
    kmer_size = unique( c(nchar(set1),nchar(set2)))
    if(length(kmer_size) != 1){
      stop("All kmer must have the same size! Try function `kmer_match_detailed_many_mers()` instead")
    }
    
    # set1_as_vector_l = strsplit(set1,"")
    # set2_as_vector_l = strsplit(set2,"")
    
    for(skip_ind_i in 1:(kmer_size-1)){
      
      for(skip_ind_j in (skip_ind_i+1):(kmer_size)){
        
        write(sprintf(" query super-imperfect match: (%d,%d)/%d",skip_ind_i, skip_ind_j, kmer_size),stderr());
        
        # dt_set1_pending_cur = data.table(kmer = set1, k_minus2_mer = unlist(lapply(set1_as_vector_l,function(x) {return (paste(x[-c(skip_ind_i, skip_ind_j)],collapse = ""))})))
        # dt_set2_pending_cur = data.table(kmer = set2, k_minus2_mer = unlist(lapply(set2_as_vector_l,function(x) {return (paste(x[-c(skip_ind_i, skip_ind_j)],collapse = ""))})))
        # 
        
        set1_pending_one_less_mer =  paste(substr(set1,1,skip_ind_i-1),
                                           substr(set1,skip_ind_i+1,skip_ind_j-1),
                                           substr(set1,skip_ind_j+1,kmer_size),
                                           sep="")
        
        set2_pending_one_less_mer =  paste(substr(set2,1,skip_ind_i-1),
                                           substr(set2,skip_ind_i+1,skip_ind_j-1),
                                           substr(set2,skip_ind_j+1,kmer_size),
                                           sep="")
        
        dt_set1_pending_cur = data.table(kmer = set1, k_minus2_mer = set1_pending_one_less_mer)
        dt_set2_pending_cur = data.table(kmer = set2, k_minus2_mer = set2_pending_one_less_mer)
        
        dt_set1_and_2_superimperfect_match_cur = merge(dt_set1_pending_cur,dt_set2_pending_cur,
                                                       by="k_minus2_mer");
        
        
        dt_set1_and_2_superimperfect_match_cur$letter_skipped_i = skip_ind_i;
        dt_set1_and_2_superimperfect_match_cur$letter_skipped_j = skip_ind_j;
        
        if(skip_ind_i==1 && skip_ind_j==2){
          dt_set1_and_2_superimperfect_match = dt_set1_and_2_superimperfect_match_cur;
        }else{
          dt_set1_and_2_superimperfect_match = rbind(dt_set1_and_2_superimperfect_match,
                                                     dt_set1_and_2_superimperfect_match_cur);
        }
        
      }
    }
    
    return(dt_set1_and_2_superimperfect_match);
    
  }
  
  
kmer_imperfect_match_general <- function(set1,set2, m_mismatch=1,warning_loop_step=T){
    #It run all combinations of `m` mismatches.
    #It uses function combn
    
    #  warning("prototyping script. Input must containt non-matching sets. I.e. sets without any perfect or 1 mismatch imperfect match");
    
    set1 = unique(set1);
    set2 = unique(set2);
    
    kmer_size = unique( c(nchar(set1),nchar(set2)))
    if(length(kmer_size) != 1){
      stop("All kmer must have the same size! Try function `kmer_match_detailed_many_mers()` instead.")
    }
    
    # set1_as_vector_l = strsplit(set1,"")
    # set2_as_vector_l = strsplit(set2,"")
    comb_set = combn(kmer_size, m_mismatch);
    
    for(i in 1:ncol(comb_set)){
      
      skipping_points = comb_set[,i];
      
      if(warning_loop_step){
        write(sprintf("  query imperfect match: %d/%d; index_skipped=%s",
                      i,ncol(comb_set),
                      paste(skipping_points,collapse = ",")),
              stderr());
      }
      
      j=1
      set1_pending_m_less_mer = substr(set1,1,skipping_points[j]-1)
      set2_pending_m_less_mer = substr(set2,1,skipping_points[j]-1);
      if(m_mismatch>1){
        for(j in 2:m_mismatch){
          set1_pending_m_less_mer = paste(set1_pending_m_less_mer,
                                          substr(set1,skipping_points[j-1]+1,skipping_points[j]-1),
                                          sep="");
          set2_pending_m_less_mer = paste(set2_pending_m_less_mer,
                                          substr(set2,skipping_points[j-1]+1,skipping_points[j]-1),
                                          sep="");
        }
      }
      set1_pending_m_less_mer = paste(set1_pending_m_less_mer,
                                      substr(set1,skipping_points[j]+1,kmer_size),
                                      sep="");
      set2_pending_m_less_mer = paste(set2_pending_m_less_mer,
                                      substr(set2,skipping_points[j]+1,kmer_size),
                                      sep="");
      
      dt_set1_pending_cur = data.table(kmer = set1, k_m_less_mer = set1_pending_m_less_mer)
      dt_set2_pending_cur = data.table(kmer = set2, k_m_less_mer = set2_pending_m_less_mer)
      
      dt_set1_and_2_superimperfect_match_cur = merge(dt_set1_pending_cur,dt_set2_pending_cur,
                                                     by="k_m_less_mer");
      
      dt_set1_and_2_superimperfect_match_cur$letter_skipped = paste(skipping_points,collapse = ";");
      dt_set1_and_2_superimperfect_match_cur$mismatch_size=m_mismatch
      
      if(i==1){
        dt_set1_and_2_superimperfect_match = dt_set1_and_2_superimperfect_match_cur;
      }else{
        dt_set1_and_2_superimperfect_match = rbind(dt_set1_and_2_superimperfect_match,
                                                   dt_set1_and_2_superimperfect_match_cur);
      }
      
    }
    
    
    return(dt_set1_and_2_superimperfect_match);
    
  }
  
  
  
kmer_Evalue_binomial_estimation <- function(N1, N2, alphabet_size=20, kmer_size=8, max_mismatch=2){
    #This is the rough* calculation of expected hits for a set of k-mer with `m` mismatches.
    #Let N1 and N2 be the number of k-mers in set1 and set2. The expected number of hits with `m` mismatches is computed as:
    #  Expected_hits_null_set1 = N1*N2*choose(k,m)/alphabet_size^(k-m)
    #
    #Expected hits set1:
    #N1*(1 - (1-N2/alphabet_size^(k-m))^choose(k,m))
    #Expected hits set2:
    #N2*(1 - (1-N1/alphabet_size^(k-m))^choose(k,m))
    #
    #*This is not fully accurate because it assumes all `k-m`-mer from a set is unique, which is not true. (edited) 
    
    #kmer_size = unique( c(nchar(set1),nchar(set2)) );
    
    if(length(kmer_size) != 1){
      stop("All kmer must have the same size! Try function `kmer_match_detailed_many_mers()` instead")
    }
    
    # N1 = length(set1);
    # N2 = length(set2);
    
    for(m in 0:max_mismatch){
      N_m_mismatch_total = alphabet_size^(kmer_size-m)
      # N_m_mismatch_set1 = N1*choose(kmer_size,m)
      # N_m_mismatch_set2 = N2*choose(kmer_size,m)
      
      # p_hits_null = N1*N2*choose(kmer_size,m)/(N_m_mismatch_total^2);
      # E_hit_null = N1*N2*p_hits_null;
      
      #p_hits_null = N1*N2*choose(kmer_size,m)/(N_m_mismatch_total^2);
      #E_hit_null = as.double(N1)*N2*choose(kmer_size,m)/(N_m_mismatch_total)
      
      #Set1
      E_hit_null_N1 = N1*(1 - (1-(N2/(N_m_mismatch_total)))^choose(kmer_size,m))
      E_hit_null_N2 = N2*(1 - (1-(N1/(N_m_mismatch_total)))^choose(kmer_size,m))
      
      dt_E_hits_set1_cur = data.table(N_mismatch = m,
                                      Set1_size = N1,
                                      Set2_size = N2,
                                      E_match=E_hit_null_N1,
                                      E_match_se=sqrt(E_hit_null_N1),
                                      E_match_ci_low = max(E_hit_null_N1 - 1.96*sqrt(E_hit_null_N1),0),
                                      E_match_ci_high = E_hit_null_N1 + 1.96*sqrt(E_hit_null_N1),
                                      set="set1",
                                      model="binomial")
      
      dt_E_hits_set2_cur = data.table(N_mismatch = m,
                                      Set1_size = N1,
                                      Set2_size = N2,
                                      E_match=E_hit_null_N2,
                                      E_match_se=sqrt(E_hit_null_N2),
                                      E_match_ci_low = max(E_hit_null_N2 - 1.96*sqrt(E_hit_null_N2),0),
                                      E_match_ci_high = E_hit_null_N2 + 1.96*sqrt(E_hit_null_N2),
                                      set="set2",
                                      model="binomial")
      
      dt_E_hits_cur = rbind(dt_E_hits_set1_cur,
                            dt_E_hits_set2_cur)#p_null = p_hits_null);
      
      # colnames(dt_E_hits_cur) = c(sprintf("E_null_%d_mismatched",m),
      #                             sprintf("E_null_E_null_ci_low_%d_mismatched",m))
      
      if(m==0){
        dt_E_hits = dt_E_hits_cur
      }else{
        dt_E_hits = rbind(dt_E_hits,
                          dt_E_hits_cur)
      }
    }
    
    return(dt_E_hits);
  }
  
  
kmer_Evalue_binomial_and_correction_estimation <- function(N1, N2, alphabet_size=20, kmer_size=8, max_mismatch=2){
    #This is an attempt to correct for the bias in E-value estimation by binomial prediction.
    #
    #The binomial estimation of E_value is slightlty higher than what is actually being observed.
    #I postulate the smaller E_value in simulation occurs because chances of match were removed from at prior matching stage.
    #
    #Thus, I propose `alpha` the following correction of E_value
    #
    #alpha1(N1,N2,kmer_size,m ) = alpha1(N1,N2,kmer_size,m -1)*(E_value_binomial(k,m-1))/(N2/(N_m_mismatch_total)*choose(kmer_size,m))
    #
    #E_value_corrected ~ E_value_binomial(k,m)*alpha
    #
    #and alpha is correct by the average expected match the relative number of kmer divided by 
    #where alpha1 = (E_value_binomial(k,m))/(N2/(N_m_mismatch_total)*choose(kmer_size,m))
    #and alpha2 = (E_value_binomial(k,m))/(N1/(N_m_mismatch_total)*choose(kmer_size,m))
    #
    
    stop("This was not finished yet! I need to spend more time thinking about how to make this correction.")
    
    for(m in 0:max_mismatch){
      N_m_mismatch_total = alphabet_size^(kmer_size-m)
      
      #Set1
      E_hit_null_N1 = N1*(1 - (1-(N2/(N_m_mismatch_total)))^choose(kmer_size,m))
      E_hit_null_N2 = N2*(1 - (1-(N1/(N_m_mismatch_total)))^choose(kmer_size,m))
      
      alpha1 = (1 - (1-(N2/(N_m_mismatch_total)))^choose(kmer_size,m)) / (choose(kmer_size,m)*N2/(N_m_mismatch_total))
      alpha2 = (1 - (1-(N1/(N_m_mismatch_total)))^choose(kmer_size,m)) / (choose(kmer_size,m)*N1/(N_m_mismatch_total))
      
      dt_E_hits_set1_cur = data.table(N_mismatch = m,
                                      Set1_size = N1,
                                      Set2_size = N2,
                                      E_match=E_hit_null_N1*alpha1,
                                      E_match_se=sqrt(E_hit_null_N1*alpha1),
                                      E_match_ci_low = max(E_hit_null_N1*alpha1 - 1.96*sqrt(E_hit_null_N1*alpha1),0),
                                      E_match_ci_high = E_hit_null_N1*alpha1 + 1.96*sqrt(E_hit_null_N1*alpha1),
                                      alpha=alpha1,
                                      set="set1",
                                      model="binomial_and_correction")
      
      dt_E_hits_set2_cur = data.table(N_mismatch = m,
                                      Set1_size = N1,
                                      Set2_size = N2,
                                      E_match=E_hit_null_N2*alpha2,
                                      E_match_se=sqrt(E_hit_null_N2*alpha2),
                                      E_match_ci_low = max(E_hit_null_N2*alpha2 - 1.96*sqrt(E_hit_null_N2*alpha2),0),
                                      E_match_ci_high = E_hit_null_N2*alpha2 + 1.96*sqrt(E_hit_null_N2*alpha2),
                                      alpha=alpha2,
                                      set="set2",
                                      model="binomial_and_correction")
      
      dt_E_hits_cur = rbind(dt_E_hits_set1_cur,
                            dt_E_hits_set2_cur)#p_null = p_hits_null);
      
      # colnames(dt_E_hits_cur) = c(sprintf("E_null_%d_mismatched",m),
      #                             sprintf("E_null_E_null_ci_low_%d_mismatched",m))
      
      if(m==0){
        dt_E_hits = dt_E_hits_cur
      }else{
        dt_E_hits = rbind(dt_E_hits,
                          dt_E_hits_cur)
      }
    }
    
    return(dt_E_hits);
  }
  
kmer_Evalue_simulation_estimation <- function(Set1_size, Set2_size, alphabet_size=20, kmer_size=8, max_mismatch=2, simulation_repeats=5){
    #This simulate Evalue of kmers accoring to a simulation
    
    alphabet = LETTERS[1:alphabet_size];
    
    for(j in 1:simulation_repeats){
      
      M_set1 = matrix(sample(alphabet, 8*Set1_size,replace = T), ncol = 8)
      S1 = unique(apply(M_set1,1,function(x) {paste(x,collapse = "")}))
      
      M_set2 = matrix(sample(alphabet, 8*Set2_size,replace = T), ncol = 8)
      S2 = unique(apply(M_set2,1,function(x) {paste(x,collapse = "")}))
      
      print(sprintf("simulation: %d/%d; Set1_size:%d, Set2_size:%d",j,simulation_repeats,Set1_size,Set2_size))
      
      o_null=kmer_match_detailed(S1, S2,
                                 #output_file = "~/temp_out",
                                 mismatch_max = max_mismatch,
                                 input_as_vector = T);
      
      dt_Evalue_sim_cur = cbind(data.table(Set1_size=Set1_size,Set2_size=Set2_size, simulation_index=j),
                                o_null$N_summary);
      
      if(j==1){
        dt_Evalue_sim = dt_Evalue_sim_cur
      }else{
        dt_Evalue_sim = rbind(dt_Evalue_sim,
                              dt_Evalue_sim_cur);
      }
    }
    
    dt_Evalue_sim_long = melt.data.table(dt_Evalue_sim,id.vars = c("Set1_size","Set2_size","simulation_index"),
                                         value.name = "N_kmers",variable.name = "match_type")
    
    dt_Evalue_sim_long_mean_and_se = dt_Evalue_sim_long[,.(E_match=mean(N_kmers),E_match_se=sqrt(var(N_kmers)/.N)),
                                                        .(Set1_size,Set2_size,match_type)]
    
    dt_Evalue_sim_long_mean_and_se[,c("E_match_ci_low","E_match_ci_high") := .(E_match-1.96*E_match_se,E_match+1.96*E_match_se)]
    dt_Evalue_sim_long_mean_and_se[E_match_ci_low<0,E_match_ci_low:=0]
    
    dt_Evalue_sim_long_mean_and_se$N_mismatch = gsub("_set.","",gsub("mismatch_","",dt_Evalue_sim_long_mean_and_se$match_type))
    dt_Evalue_sim_long_mean_and_se[N_mismatch=="perfect_match",N_mismatch:=0]
    dt_Evalue_sim_long_mean_and_se[N_mismatch=="no_match",N_mismatch:=.(sprintf("%d+",max_mismatch+1))]
    dt_Evalue_sim_long_mean_and_se[, set :=.(gsub("mismatch_._","",match_type))]
    dt_Evalue_sim_long_mean_and_se[grep("no_match",set), set :=.(gsub("no_match_","",match_type))]
    dt_Evalue_sim_long_mean_and_se[grep("perfect_match",set), set :="set1"]
    
    dt_Evalue_sim_long_mean_and_se$model="simulation";
    dt_Evalue_sim_long_mean_and_se$simulation_repeats = simulation_repeats
    return(dt_Evalue_sim_long_mean_and_se);
  }

  
 