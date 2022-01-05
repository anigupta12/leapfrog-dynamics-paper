pheno_m_an = pheno_matrix_add_an_ordered
pheno_m_an = ifelse(pheno_m_an > 0, 1, 0)
pheno_m_nophagean = pheno_m_an[, -1]

write.table(pheno_m_nophagean, file = "bimat_input_matrix_no_phage_an.csv", sep = ",", quote = F, row.names = F, col.names = F)


pheno_m_noan_noempty = pheno_m_an[-1, -1]
pheno_m_noan_noempty = pheno_m_noan_noempty[-which(rowSums(pheno_m_noan_noempty) == 0), -which(colSums(pheno_m_noan_noempty) == 0)]

write.table(pheno_m_noan_noempty, file = "bimat_input_matrix_no_empty_rowcol.csv", sep = ",", quote = F, row.names = F, col.names = F)