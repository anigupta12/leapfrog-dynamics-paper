options(stringsAsFactors = F)
library(stringr)
library(reshape2)
library(pracma)
library(glmnet)
library(XML)
library(ggplot2)
library(randomForest)
library(caret)
library(Metrics)
#library(xgboost)
library(pbapply)

source("app_toolbox_6.R")
source("additional_app_toolbox_6.R")

set.seed(101)

################ bacteria and virus mutation profile

# updated code and files to create bacteria phylogeny is in bac_tree

# bac_table_doc = htmlParse("bac_anno_50.html", encoding = "UTF-8")
# bac_table = data.frame(readHTMLTable(bac_table_doc, stringsAsFactors = FALSE, check.names = F))
# colnames(bac_table) = gsub("-", "_", bac_table[1, ])
# bac_table = bac_table[which(bac_table[, 1] != "position"), ]
# rownames(bac_table) = 1:nrow(bac_table)
# bac_table = data.frame(bac_table[, 1:2], an_b = rep("", nrow(bac_table)), bac_table[, 3:ncol(bac_table)])
# 
# # deal with 8th row only
# #bac_table_8 = bac_table
# #bac_table_8[9, which(bac_table[8, ] == "100%")] = "100%"
# #bac_table_8 = bac_table_8[-8, ]
# #bac_table = bac_table_8
# # deal with the 1st and 8th row
# #bac_dealed = bac_table
# #bac_dealed[9, which(bac_dealed[8, ] == "100%")] = "100%"
# #bac_dealed = bac_dealed[c(-1, -2, -8), ]
# #bac_dealed = bac_dealed[c(-1), ]
# bac_table = bac_table[-1, ]
# 
# bac_feature_raw = bac_table[, 3:53]
# bac_feature_binary = ifelse(bac_feature_raw == "" | bac_feature_raw == "?", 0, 1)
# #bac_feature_binary = ifelse(bac_feature_raw == "", 0, 1)
# bac_table = data.frame(bac_table, init_appear_day = get_initial_appear(bac_feature_binary), nonsyn_mutation = check_nonsyn(bac_table[, 54]))
# bac_feature_binary_add_an = bac_feature_binary    #only add genome but not mutation/ancestor indicator
# 
# no_mutation_bac_strains = which(colSums(bac_feature_binary_add_an) == 0)
# 
# bac_duplicated_genomes = which(duplicated(t(bac_feature_binary_add_an)))
# bac_non_duplicated_genomes = which(!duplicated(t(bac_feature_binary_add_an)))
# bac_feature_binary_add_an_uniq = bac_feature_binary_add_an[, bac_non_duplicated_genomes]
# bac_duplicated_genomes_mapping = get_dup_genome_mapping(bac_non_duplicated_genomes, bac_feature_binary_add_an)
# bac_feature_binary_add_an_uniq_nonsyn_only = bac_feature_binary_add_an_uniq[which(bac_table$nonsyn_mutation == 1), ]
# 
# bac_uniq_genomes_name = colnames(bac_feature_binary)[bac_non_duplicated_genomes]
# #bac_uniq_genome_table = bac_table_8[, c(1, 2, bac_non_duplicated_genomes+2, 54:56)]
# 
# #write.table(bac_uniq_genome_table, file = "bac_uniq_genome_table.tsv", sep = "\t", quote = F, row.names = F)

# phage mutation profile to create phylogeny
phage_table_doc = htmlParse("phage_all_annotated_cut50.html", encoding = "UTF-8")
phage_table = data.frame(readHTMLTable(phage_table_doc, stringsAsFactors = FALSE, check_names = F))
colnames(phage_table) = gsub("-", "_", phage_table[1, ])
phage_table = phage_table[which(phage_table[, 1] != "position"), ]
rownames(phage_table) = 1:nrow(phage_table)
phage_table = data.frame(phage_table[, 1:2], an_p = rep("", nrow(phage_table)), phage_table[, 3:ncol(phage_table)])

phage_feature_raw = phage_table[, 3:47]
phage_feature_binary = ifelse(phage_feature_raw == "" | phage_feature_raw == "?", 0, 1)
#phage_feature_binary = ifelse(phage_feature_raw == "", 0, 1)
phage_table = data.frame(phage_table, init_appear_day = get_initial_appear(phage_feature_binary), nonsyn_mutation = check_nonsyn(phage_table[, 48]))
phage_feature_binary_add_an = phage_feature_binary

phage_duplicated_genomes = which(duplicated(t(phage_feature_binary_add_an)))
phage_non_duplicated_genomes = which(!duplicated(t(phage_feature_binary_add_an)))
phage_feature_binary_add_an_uniq = phage_feature_binary_add_an[, phage_non_duplicated_genomes]
phage_duplicated_genomes_mapping = get_dup_genome_mapping(phage_non_duplicated_genomes, phage_feature_binary_add_an)

phage_uniq_genomes_name = colnames(phage_feature_binary)[phage_non_duplicated_genomes]
phage_uniq_genome_table = phage_table[, c(1, 2, phage_non_duplicated_genomes + 2, 48:52)]

#phage nonsynonymous mutations only
phage_feature_binary_add_an_uniq_nonsyn_only = phage_feature_binary_add_an_uniq[which(phage_table$nonsyn_mutation == 1), ]
phage_duplicated_genomes_after_nonsyn =  which(duplicated(t(phage_feature_binary_add_an_uniq_nonsyn_only)))
phage_non_duplicated_genomes_after_nonsyn = which(!duplicated(t(phage_feature_binary_add_an_uniq_nonsyn_only)))
phage_feature_binary_add_an_uniq_nonsyn_only_uniq_again = phage_feature_binary_add_an_uniq_nonsyn_only[, phage_non_duplicated_genomes_after_nonsyn]
phage_duplicated_genomes_mapping_after_nonsyn = get_dup_genome_mapping(phage_non_duplicated_genomes_after_nonsyn, phage_feature_binary_add_an_uniq_nonsyn_only)


################ phenotypes load

raw_list = read.csv("raw_infection_list_mod8-1.csv", header=T, check.names = F)
row_idx = which((raw_list$Host != "A") & !str_detect(raw_list$Host, "^606") & (raw_list$Phage != "A"))
raw_list_noan = raw_list[row_idx, ]
raw_list_noan[, 1] = paste0("P_D_", raw_list_noan[, 1])
raw_list_noan[, 1] = gsub("-", "_", raw_list_noan[, 1])
raw_list_noan[, 2] = paste0("B_D_", raw_list_noan[, 2])
raw_list_noan[, 2] = gsub("-", "_", raw_list_noan[, 2])

#using eop as phenotype
pheno_matrix = dcast(data = raw_list_noan, Host ~ Phage, value.var="EOP")

rownames(pheno_matrix) = pheno_matrix[,1]
pheno_matrix = pheno_matrix[,-1]
pheno_matrix_add_an = data.frame(an_p = rep(0, nrow(pheno_matrix)), pheno_matrix)
pheno_matrix_add_an = rbind(an_b = rep(1, ncol(pheno_matrix_add_an)), pheno_matrix_add_an)

pheno_matrix_add_an_ordered = pheno_matrix_add_an
pheno_matrix_add_an_ordered = pheno_matrix_add_an_ordered[colnames(bac_feature_binary), ]
pheno_matrix_add_an_ordered = pheno_matrix_add_an_ordered[, colnames(phage_feature_binary)]

#using abolute as phenotype
pheno_matrix_abs = dcast(data = raw_list_noan, Host ~ Phage, value.var="Plaques per ml")

rownames(pheno_matrix_abs) = pheno_matrix_abs[,1]
pheno_matrix_abs = pheno_matrix_abs[,-1]
# add ancestor phage vs all except for ancestor bacteria
pheno_matrix_abs_add_an = data.frame(an_p = rep(0, nrow(pheno_matrix_abs)), pheno_matrix_abs)
# add ancestor bacteria vs all phage MANUAL!
##### rows in orginal list (need to -1 to remove header): (2252, 2297, 2342), 2309, 2318, 2319, 2310:2317, 2275, 2284, 2285, 2276:2283, 2286, 2295, 2296, 2287:2294, 
##### an_bac_phe = c((41625000 + 67640625 + 286171875)/3, 249750, 208125, 1332, 3996, 49950, 6660, 1998, 13320, 49950, 6660, 416250, 41625, 49950, 16650, 9990,19980, 33300, 8325, 13320, 9990, 249750, 291375, 83250, 58275, 2289375, 49950, 49950, 333000, 66600, 1040625, 208125, 249750, 1040625)

an_bac_phe = c((raw_list[2251, 5] + raw_list[2296, 5] + raw_list[2341, 5])/3, (raw_list[c(2252:2295), 5] + raw_list[c(2297:2340), 5] + raw_list[c(2342:2385), 5])/3)

pheno_matrix_abs_add_an = rbind(an_b = an_bac_phe, pheno_matrix_abs_add_an)

pheno_matrix_abs_add_an_ordered = pheno_matrix_abs_add_an
pheno_matrix_abs_add_an_ordered = pheno_matrix_abs_add_an_ordered[colnames(bac_feature_binary), ]
pheno_matrix_abs_add_an_ordered = pheno_matrix_abs_add_an_ordered[, colnames(phage_feature_binary)]

#avg same genome phenotyp to get unique matrix

pheno_matrix_add_an_ordered_uniq_avg = get_unique_phenom_by_avg(pheno_matrix_add_an_ordered, bac_duplicated_genomes_mapping, phage_duplicated_genomes_mapping)

pheno_matrix_abs_add_an_ordered_uniq_avg = get_unique_phenom_by_avg(pheno_matrix_abs_add_an_ordered, bac_duplicated_genomes_mapping, phage_duplicated_genomes_mapping)


################ check same genome have more similar phenotypes

bac_all_names = rownames(pheno_matrix_add_an_ordered)
bac_unique_names = colnames(bac_feature_binary_add_an_uniq_nonsyn_only)
phage_all_names = colnames(pheno_matrix_add_an_ordered)
phage_unique_names = colnames(phage_feature_binary_add_an_uniq_nonsyn_only)

#using eop
bac_cor_res = get_phylo_cor_res(t(pheno_matrix_add_an_ordered), bac_duplicated_genomes_mapping)
phage_cor_res = get_phylo_cor_res(pheno_matrix_add_an_ordered, phage_duplicated_genomes_mapping)

# pdf("same_geno_diff_pheno2.pdf", height = 12, width = 12, useDingbats = F)
# 
# par(mfrow = c(2,2))
# 
# hist(bac_cor_res[["nondup_cor"]], breaks = 50, col = rgb(1, 0, 0, 0.5), main = "bacteria relative infectivity correlation (EOP)", xlab = "Pearson's correlation")
# hist(bac_cor_res[["dup_cor"]], breaks = 20, col = rgb(0, 0, 1, 0.5), add = T)
# legend("topright", c("Not same genome pairs", "Same genome pairs"), col=c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), lwd=5)
# 
# 
# hist(phage_cor_res[["nondup_cor"]], breaks = 50, col = rgb(1, 0, 0, 0.5), main = "phage relative infectivity correlation (EOP)", xlab = "Pearson's correlation")
# hist(phage_cor_res[["dup_cor"]], breaks = 20, col = rgb(0, 0, 1, 0.5), add = T)
# 
# #using abs
# bac_cor_res_abs = get_phylo_cor_res(t(pheno_matrix_abs_add_an_ordered), bac_duplicated_genomes_mapping)
# phage_cor_res_abs = get_phylo_cor_res(pheno_matrix_abs_add_an_ordered, phage_duplicated_genomes_mapping)
# 
# hist(bac_cor_res_abs[["nondup_cor"]], breaks = 50, col = rgb(1, 0, 0, 0.5), main = "bacteria absolute infectivity correlation\n(Plaque per ml)", xlab = "Pearson's correlation")
# hist(bac_cor_res_abs[["dup_cor"]], breaks = 20, col = rgb(0, 0, 1, 0.5), add = T)
# 
# hist(phage_cor_res_abs[["nondup_cor"]], breaks = 50, col = rgb(1, 0, 0, 0.5), main = "phage absolute infectivity correlation\n(Plaque per ml)", xlab = "Pearson's correlation")
# hist(phage_cor_res_abs[["dup_cor"]], breaks = 20, col = rgb(0, 0, 1, 0.5), add = T)
# 
# dev.off()

bac_dist_vec = dist(t(bac_feature_binary_add_an), method = 'manhattan')
bac_pheno_vec = get_phylo_cor_res_nosep(t(pheno_matrix_abs_add_an_ordered))

phage_dist_vec = dist(t(phage_feature_binary_add_an), method = 'manhattan')
phage_pheno_vec = get_phylo_cor_res_nosep(pheno_matrix_abs_add_an_ordered)

#plot(bac_pheno_vec, bac_dist_vec, xlab = "Phenotype Correlation", ylab = "Genome Hamming Distance", main = "Bacteria")
#plot(phage_pheno_vec, phage_dist_vec, xlab = "Phenotype Correlation", ylab = "Genome Hamming Distance", main = "Virus")

################ calc phage dn/ds value per clone

# phage_mut_types <- phage_table[, ncol(phage_table)]
# phage_table_nonsyn <- phage_table[phage_mut_types == 1, ]
# 
# phage_dnds_per_clone_uniq <- apply(phage_feature_binary_add_an_uniq, 2, function(x){
#   mut_rows <- which(x == 1)
#   mut_types <- phage_mut_types[mut_rows]
#   if(length(mut_types) - sum(mut_types) == 0){
#     dnds <- 999
#   }else{
#     dnds <- sum(mut_types)/(length(mut_types) - sum(mut_types))
#   }
#   return(dnds)
# })
# 
# phage_dnds_ratio_per_clone_uniq <- apply(phage_feature_binary_add_an_uniq, 2, function(x){
#   mut_rows <- which(x == 1)
#   mut_types <- phage_mut_types[mut_rows]
#   dnds <- paste0(sum(mut_types), " vs ", length(mut_types) - sum(mut_types))
#   return(dnds)
# })
# 
# phage_dnds_cnt_per_clone_uniq <- t(apply(phage_feature_binary_add_an_uniq, 2, function(x){
#   mut_rows <- which(x == 1)
#   mut_types <- phage_mut_types[mut_rows]
#   #dnds <- paste0(sum(mut_types), " vs ", length(mut_types) - sum(mut_types))
#   return(c("S"= length(mut_types) - sum(mut_types), "N" = sum(mut_types)))
# }))
# 
# phage_dnds_table = cbind(ratio = round(phage_dnds_per_clone_uniq, 2), value = phage_dnds_ratio_per_clone_uniq)
# phage_dnds_table[phage_dnds_table == 999] = "Inf"

#write.csv(phage_dnds_table, file = "phage_dnds_table.csv", quote = F)

#write.csv(phage_dnds_cnt_per_clone_uniq, file = "phage_strain_level_dnds.csv", quote = F)


################ modelling prepare

#bac nonsyn mutation and all genome 16 by 51
bac_feature_binary_add_an_m_allmut_allgen <- rbind(b_an_ind = c(1, rep(0, ncol(bac_feature_binary_add_an) - 1)), bac_feature_binary_add_an)
bac_feature_binary_add_an_m_allmut_allgen[1, no_mutation_bac_strains] <- 1    # set ancestor indicator to 1 for no mutation strains

#bac nonsyn mutation and uniq genome 16 by 16
bac_feature_binary_add_an_m <- rbind(b_an_ind = c(1, rep(0, ncol(bac_feature_binary_add_an_uniq_nonsyn_only) - 1)), bac_feature_binary_add_an_uniq_nonsyn_only)
bac_feature_binary_add_an_m[1, no_mutation_bac_strains] <- 1        # set ancestor indicator to 1 for no mutation strains


#phage all mutation and all genome 177 by 45
phage_feature_binary_add_an_m_allmut_allgen <- rbind(p_an_ind = c(1, rep(0, ncol(phage_feature_binary_add_an) - 1)), phage_feature_binary_add_an)

#phage nonsyn mutation and all genome 61 by 45
phage_feature_binary_add_an_nonsyn_only = phage_feature_binary_add_an[which(phage_table$nonsyn_mutation == 1), ]
phage_feature_binary_add_an_m_nonsyn_allgen <- rbind(p_an_ind = c(1, rep(0, ncol(phage_feature_binary_add_an_nonsyn_only) - 1)), phage_feature_binary_add_an_nonsyn_only)

#phage all mutation and uniq genome 147 by 35
phage_feature_binary_add_an_m_allmut_uniqgen <- rbind(p_an_ind = c(1, rep(0, ncol(phage_feature_binary_add_an_uniq) - 1)), phage_feature_binary_add_an_uniq)

#phage nonsyn mutation and uniq genome 61 by 35
phage_feature_binary_add_an_m <- rbind(p_an_ind = c(1, rep(0, ncol(phage_feature_binary_add_an_uniq_nonsyn_only) - 1)), phage_feature_binary_add_an_uniq_nonsyn_only)
#phage_feature_binary_add_an_m_rref <- rref(phage_feature_binary_add_an_m)
#phage_feature_binary_add_an_m_reduce_idx <- which(rowSums(phage_feature_binary_add_an_m_rref)!=0)
phage_feature_binary_add_an_m_reduce <- phage_feature_binary_add_an_m

pheno_matrix_add_an_ordered_uniq_avg_m <- as.matrix(pheno_matrix_add_an_ordered_uniq_avg)
#pheno_matrix_add_an_ordered_uniq_avg_m <- as.matrix(pheno_matrix_abs_add_an_ordered_uniq_avg)

# determine cut or not
all_labels_raw = as.vector(as.matrix(pheno_matrix_add_an_ordered))
#all_labels_raw = as.vector(pheno_matrix_add_an_ordered_uniq_avg_m)

pheno_cap = quantile(all_labels_raw, 0.99)

#exclude_pairs = which(all_labels_raw >= pheno_cap)
exclude_pairs = which(all_labels_raw >= pheno_cap | all_labels_raw == 0)

all_labels = all_labels_raw
all_labels = all_labels[-exclude_pairs]

pairwise_concatenation = as.matrix(merge(x = t(bac_feature_binary_add_an_m_allmut_allgen), y = t(phage_feature_binary_add_an_nonsyn_only), by = NULL))
#pairwise_concatenation = as.matrix(merge(x = t(bac_feature_binary_add_an_m), y = t(phage_feature_binary_add_an_m_reduce), by = NULL))
pairwise_concatenation = pairwise_concatenation[-exclude_pairs, ]

pairwise_convolution = kronecker(t(phage_feature_binary_add_an_nonsyn_only), t(bac_feature_binary_add_an_m_allmut_allgen))
#pairwise_convolution = kronecker(t(phage_feature_binary_add_an_m_reduce), t(bac_feature_binary_add_an_m))
pairwise_convolution = pairwise_convolution[-exclude_pairs, ]

pairwise_conc_p_conv = cbind(pairwise_concatenation, pairwise_convolution)
colnames(pairwise_conc_p_conv) = NULL

knowledge_gene_col_idx = 7
knowledge_super_resistant_pairs = which(pairwise_conc_p_conv[, knowledge_gene_col_idx] == 1)
knowledge_super_other = which(pairwise_conc_p_conv[, knowledge_gene_col_idx] == 0)


