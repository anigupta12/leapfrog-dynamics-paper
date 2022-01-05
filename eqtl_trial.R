library(MatrixEQTL)

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

# # 2. phage all mut uniq gen bac all mut all gen
# p_au_b_aa_1st <- p_na_b_aa_1st(phage_only_eqtl_m_allmut_uniqgen, bac_only_eqtl_m_allmut_allgen)
# p_au_b_aa_2nd <- construct_eqtl_feature_table_2nd(phage_only_eqtl_m_allmut_uniqgen, bac_only_eqtl_m_allmut_allgen)

# 3. phage nonsyn mut all gen bac all mut all gen
p_na_b_aa_1st <- construct_eqtl_feature_table_1st_inworking(phage_only_eqtl_m_nonsyn_allgen, bac_only_eqtl_m_allmut_allgen)
p_na_b_aa_2nd <- construct_eqtl_feature_table_2nd(phage_only_eqtl_m_nonsyn_allgen, bac_only_eqtl_m_allmut_allgen)
  
# # 4. phage nonsyn mut uniq gen bac all mut all gen
# p_nu_b_aa_1st <- p_na_b_aa_1st(phage_only_eqtl_m_nonsyn_uniqgen, bac_only_eqtl_m_allmut_allgen)
# p_nu_b_aa_2nd <- construct_eqtl_feature_table_2nd(phage_only_eqtl_m_nonsyn_uniqgen, bac_only_eqtl_m_allmut_allgen)

# # 5. phage all mut all gen bac all mut uniq gen
# p_aa_b_au_1st <- p_na_b_aa_1st(phage_only_eqtl_m_allmut_allgen, bac_only_eqtl_m_allmut_uniqgen)
# p_aa_b_au_2nd <- construct_eqtl_feature_table_2nd(phage_only_eqtl_m_allmut_allgen, bac_only_eqtl_m_allmut_uniqgen)

# 6. phage all mut uniq gen bac all mut uniq gen
p_au_b_au_1st <- construct_eqtl_feature_table_1st_inworking(phage_only_eqtl_m_allmut_uniqgen, bac_only_eqtl_m_allmut_uniqgen)
p_au_b_au_2nd <- construct_eqtl_feature_table_2nd(phage_only_eqtl_m_allmut_uniqgen, bac_only_eqtl_m_allmut_uniqgen)

# # 7. phage nonsyn mut all gen bac all mut uniq gen
# p_na_b_au_1st <- p_na_b_aa_1st(phage_only_eqtl_m_nonsyn_allgen, bac_only_eqtl_m_allmut_uniqgen)
# p_na_b_au_2nd <- construct_eqtl_feature_table_2nd(phage_only_eqtl_m_nonsyn_allgen, bac_only_eqtl_m_allmut_uniqgen)

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
# 
# #write file for neural network
# 
# dnn = data.frame(y=label_allgen, p_na_b_aa_1st)
# dnn2 = data.frame(y=label_allgen, p_na_b_aa_2nd)
# dnn3 = data.frame(y=label_allgen, p_na_b_aa_1st, p_na_b_aa_2nd)
# write.csv(dnn, file = "./nn/first_order_input.csv", quote = F, row.names = F)
# write.csv(dnn2, file = "./nn/second_order_input.csv", quote = F, row.names = F)
# 
# dnn_pos = dnn[which(dnn$y > 0), ]
# dnn2_pos = dnn2[which(dnn2$y > 0), ]
# dnn3_pos = dnn3[which(dnn3$y > 0), ]
# write.csv(dnn, file = "./nn/first_order_input_positive_only.csv", quote = F, row.names = F)
# write.csv(dnn2, file = "./nn/second_order_input_positive_only.csv", quote = F, row.names = F)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # model 1 phage all genome all mut bac all genome all mut
# # label_allgen
# # p_phen_allgen
# # b_phen_allgen
# # phage_only_eqtl_m_allmut_allgen
# # bac_only_eqtl_m_allmut_allgen
# # p_aa_b_aa_2nd
# 
# p_me_model1 = prepRunMatrixEqtl(gen_data = phage_only_eqtl_m_allmut_allgen, phen_data = p_phen_allgen, outprefix = "model1_phage")
# b_me_model1 = prepRunMatrixEqtl(gen_data = bac_only_eqtl_m_allmut_allgen, phen_data = b_phen_allgen, outprefix = "model1_bac")
# p_b_2nd_me_model1 = prepRunMatrixEqtl(gen_data = t(p_aa_b_aa_2nd), phen_data = label_allgen, outprefix = "model1_secOrder")
# 
# p_me_res_model1 = summaryMatrixEqtl(p_me_model1, "phage")
# b_me_res_model1 = summaryMatrixEqtl(b_me_model1, "E. coli")
# b_p_2nd_me_res_model1 = addcol_2ndOrder_eqtlSummary_res(summaryMatrixEqtl(p_b_2nd_me_model1, "2nd order"))
# 
# 
# # model 2 phage all genome nonsyn mut bac all genome all mut
# # label_allgen
# # p_phen_allgen
# # b_phen_allgen
# # phage_only_eqtl_m_nonsyn_allgen
# # bac_only_eqtl_m_allmut_allgen
# # p_na_b_aa_2nd
# 
# 
# p_me_model2 = prepRunMatrixEqtl(gen_data = phage_only_eqtl_m_nonsyn_allgen, phen_data = p_phen_allgen, outprefix = "model2_phage")
# b_me_model2 = prepRunMatrixEqtl(gen_data = bac_only_eqtl_m_allmut_allgen, phen_data = b_phen_allgen, outprefix = "model2_bac")
# p_b_2nd_me_model2 = prepRunMatrixEqtl(gen_data = t(p_na_b_aa_2nd), phen_data = label_allgen, outprefix = "model2_secOrder")
# 
# p_me_res_model2 = summaryMatrixEqtl(p_me_model2, "phage")
# b_me_res_model2 = summaryMatrixEqtl(b_me_model2, "E. coli")
# b_p_2nd_me_res_model2 = addcol_2ndOrder_eqtlSummary_res(summaryMatrixEqtl(p_b_2nd_me_model2, "2nd order"))
# 
# # model 3 phage uniq genome all mut bac uniq genome all mut
# # label_uniqgen
# # p_au_b_au_1st
# # p_au_b_au_2nd
# 
# # model 4 phage uniq genome nonsyn mut bac uniq genome all mut
# # label_uniqgen
# # p_phen_uniqgen
# # b_phen_uniqgen
# # phage_only_eqtl_m_nonsyn_uniqgen
# # bac_only_eqtl_m_allmut_uniqgen
# # p_nu_b_au_2nd
# 
# p_me_model4 = prepRunMatrixEqtl(gen_data = phage_only_eqtl_m_nonsyn_uniqgen, phen_data = p_phen_uniqgen, outprefix = "model4_phage")
# b_me_model4 = prepRunMatrixEqtl(gen_data = bac_only_eqtl_m_allmut_uniqgen, phen_data = b_phen_uniqgen, outprefix = "model4_bac")
# p_b_2nd_me_model4 = prepRunMatrixEqtl(gen_data = t(p_nu_b_au_2nd), phen_data = label_uniqgen, outprefix = "model4_secOrder")
# 
# p_me_res_model4 = summaryMatrixEqtl(p_me_model4, "phage")
# b_me_res_model4 = summaryMatrixEqtl(b_me_model4, "E. coli")
# b_p_2nd_me_res_model4 = addcol_2ndOrder_eqtlSummary_res(summaryMatrixEqtl(p_b_2nd_me_model4, "2nd order"))