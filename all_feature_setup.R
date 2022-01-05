phage_nonsyn_table = phage_table[which(phage_table$nonsyn_mutation == 1), ]
rownames(phage_nonsyn_table) = 1:nrow(phage_nonsyn_table)


#### 4 types of phage table
phage_only_eqtl_m_allmut_allgen = phage_feature_binary_add_an_m_allmut_allgen
rownames(phage_only_eqtl_m_allmut_allgen) = c("phage_an", paste0("phage_mut_", 1:(nrow(phage_only_eqtl_m_allmut_allgen)-1)))

phage_only_eqtl_m_nonsyn_allgen = phage_feature_binary_add_an_m_nonsyn_allgen
rownames(phage_only_eqtl_m_nonsyn_allgen) = c("phage_an", paste0("phage_mut_nonsyn_", 1:(nrow(phage_only_eqtl_m_nonsyn_allgen)-1)))

phage_only_eqtl_m_allmut_uniqgen = phage_feature_binary_add_an_m_allmut_uniqgen
rownames(phage_only_eqtl_m_allmut_uniqgen) = c("phage_an", paste0("phage_mut_nonsyn_", 1:(nrow(phage_only_eqtl_m_allmut_uniqgen)-1)))

phage_only_eqtl_m_nonsyn_uniqgen = phage_feature_binary_add_an_m
rownames(phage_only_eqtl_m_nonsyn_uniqgen) = c("phage_an", paste0("phage_mut_nonsyn_", 1:(nrow(phage_only_eqtl_m_nonsyn_uniqgen)-1)))

#### 2 types of bac table
bac_only_eqtl_m_allmut_allgen = bac_feature_binary_add_an_m_allmut_allgen
rownames(bac_only_eqtl_m_allmut_allgen) = c("bac_an", paste0("bac_mut_", 1:(nrow(bac_only_eqtl_m_allmut_allgen)-1)))

bac_only_eqtl_m_allmut_uniqgen = bac_feature_binary_add_an_m
rownames(bac_only_eqtl_m_allmut_uniqgen) = c("bac_an", paste0("bac_mut_", 1:(nrow(bac_only_eqtl_m_allmut_uniqgen)-1)))


#### geno-table construction 1st order 2nd order

# 1. phage all mut all gen bac all mut all gen
p_aa_b_aa_1st <- construct_eqtl_feature_table_1st_inworking(phage_only_eqtl_m_allmut_allgen, bac_only_eqtl_m_allmut_allgen)
p_aa_b_aa_2nd <- construct_eqtl_feature_table_2nd(phage_only_eqtl_m_allmut_allgen, bac_only_eqtl_m_allmut_allgen)

# 3. phage nonsyn mut all gen bac all mut all gen
p_na_b_aa_1st <- construct_eqtl_feature_table_1st_inworking(phage_only_eqtl_m_nonsyn_allgen, bac_only_eqtl_m_allmut_allgen)
p_na_b_aa_2nd <- construct_eqtl_feature_table_2nd(phage_only_eqtl_m_nonsyn_allgen, bac_only_eqtl_m_allmut_allgen)

# 6. phage all mut uniq gen bac all mut uniq gen
p_au_b_au_1st <- construct_eqtl_feature_table_1st_inworking(phage_only_eqtl_m_allmut_uniqgen, bac_only_eqtl_m_allmut_uniqgen)
p_au_b_au_2nd <- construct_eqtl_feature_table_2nd(phage_only_eqtl_m_allmut_uniqgen, bac_only_eqtl_m_allmut_uniqgen)

# 8. phage nonsyn mut uniq gen bac all mut uniq gen
p_nu_b_au_1st <- construct_eqtl_feature_table_1st_inworking(phage_only_eqtl_m_nonsyn_uniqgen, bac_only_eqtl_m_allmut_uniqgen)
p_nu_b_au_2nd <- construct_eqtl_feature_table_2nd(phage_only_eqtl_m_nonsyn_uniqgen, bac_only_eqtl_m_allmut_uniqgen)



#### pheno-table and label construction 

p_allgen_name_order = colnames(phage_only_eqtl_m_allmut_allgen)
p_uniqgen_name_order = colnames(phage_only_eqtl_m_allmut_uniqgen)
b_allgen_name_order = colnames(bac_only_eqtl_m_allmut_allgen)
b_uniqgen_name_order = colnames(bac_only_eqtl_m_allmut_uniqgen)

p_allgen_b_allgen_pheno_m = pheno_matrix_add_an[, p_allgen_name_order]
p_allgen_b_allgen_pheno_m = p_allgen_b_allgen_pheno_m[b_allgen_name_order, ]
p_uniqgen_b_uniqgen_pheno_m = pheno_matrix_add_an_ordered_uniq_avg[, p_uniqgen_name_order]
p_uniqgen_b_uniqgen_pheno_m = p_uniqgen_b_uniqgen_pheno_m[b_uniqgen_name_order, ]

p_allgen_b_allgen_pheno_name = as.vector(as.matrix(sapply(colnames(p_allgen_b_allgen_pheno_m), function(x) paste0(x, ":", rownames(p_allgen_b_allgen_pheno_m)))))
p_uniqgen_b_uniqgen_pheno_name = as.vector(as.matrix(sapply(colnames(p_uniqgen_b_uniqgen_pheno_m), function(x) paste0(x, ":", rownames(p_uniqgen_b_uniqgen_pheno_m)))))

p_allgen_b_allgen_pheno_l = as.vector(as.matrix(p_allgen_b_allgen_pheno_m))
names(p_allgen_b_allgen_pheno_l) = p_allgen_b_allgen_pheno_name
p_uniqgen_b_uniqgen_pheno_l = as.vector(as.matrix(p_uniqgen_b_uniqgen_pheno_m))
names(p_uniqgen_b_uniqgen_pheno_l) = p_uniqgen_b_uniqgen_pheno_name

label_allgen = p_allgen_b_allgen_pheno_l
label_uniqgen = p_uniqgen_b_uniqgen_pheno_l


p_phen_allgen = as.matrix(p_allgen_b_allgen_pheno_m)
b_phen_allgen = t(as.matrix(p_allgen_b_allgen_pheno_m))

p_phen_uniqgen = as.matrix(p_uniqgen_b_uniqgen_pheno_m)
b_phen_uniqgen = t(as.matrix(p_uniqgen_b_uniqgen_pheno_m))


# consider rand reduce
# bac_feature_binary_add_an_m = as.matrix(sapply(bac_feature_binary_add_an, as.numeric))
# phage_feature_binary_add_an_m = as.matrix(sapply(phage_feature_binary_add_an, as.numeric))
# phage_feature_binary_add_an_m_rref = rref(phage_feature_binary_add_an_m)
# phage_feature_binary_add_an_m_reduce_idx = which(rowSums(phage_feature_binary_add_an_m_rref)!=0)
# phage_feature_binary_add_an_m_reduce = phage_feature_binary_add_an_m_rref[phage_feature_binary_add_an_m_reduce_idx, ]
# pheno_matrix_ordered_add_an_m = as.matrix(pheno_matrix_ordered_add_an)

#write file for neural network

dnn = data.frame(y=label_allgen, p_na_b_aa_1st)
dnn2 = data.frame(y=label_allgen, p_na_b_aa_2nd)
dnn3 = data.frame(y=label_allgen, p_na_b_aa_1st, p_na_b_aa_2nd)
# write.csv(label_allgen, file = "./nn/full_label.csv", quote = F, row.names = F)
# write.csv(dnn, file = "./nn/first_order_input.csv", quote = F, row.names = F)
# write.csv(dnn2, file = "./nn/second_order_input.csv", quote = F, row.names = F)
# write.csv(dnn3, file = "./nn/comb_fs_input.csv", quote = F, row.names = F)

dnn_pos = dnn[which(dnn$y > 0), ]
dnn2_pos = dnn2[which(dnn2$y > 0), ]
dnn3_pos = dnn3[which(dnn3$y > 0), ]
# write.csv(label_allgen[label_allgen > 0], file = "./nn/positive_only_label.csv", quote = F, row.names = F)
# write.csv(dnn_pos, file = "./nn/first_order_input_positive_only.csv", quote = F, row.names = F)
# write.csv(dnn2_pos, file = "./nn/second_order_input_positive_only.csv", quote = F, row.names = F)
# write.csv(dnn3_pos, file = "./nn/comb_fs_input_positive_only.csv", quote = F, row.names = F)
