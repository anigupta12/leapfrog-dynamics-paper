options(stringsAsFactors = F)

prefix = "ecoli"
dir.create(prefix)
input_name = "bac_mutation_profile_for_circos_unique.csv"
label_file_name = paste0(prefix, ".labels.txt")
genome_end_pos = 4629812
gap_size = 4000

ecoli = read.table(input_name, header = T, sep = ",")
print(unique(ecoli$gene))
ecoli_df = ecoli[, c(1,2,4:(length(ecoli)-3))]
ecoli_df$position = as.numeric(gsub(",", "", ecoli_df$position))
ecoli_df$mutation = gsub("→", "->", ecoli_df$mutation)
ecoli_df$mutation = gsub("Δ", "delta_", ecoli_df$mutation)
ecoli_df$mutation = gsub(".bp", "_bp", ecoli_df$mutation, fixed=F)

labels_df = data.frame(chr = "chr1", start = ecoli_df$position, end = ecoli_df$position + 1, label = ecoli_df$mutation)
write.table(labels_df, paste0(prefix, "/", label_file_name), row.names = F, col.names = F, quote = F)

# sort unique genomes
ecoli_df = ecoli_df[, -2]
ecoli_m = ecoli_df[, -1] == "100%"
ecoli_color = c("vvdgreen", "vdgreen", "dgreen", "green", "lgreen", "vlgreen", "dblue", "blue", "lblue", "vlblue", "240,0,240", "vvdpurple", "vdpurple", "dpurple", "purple", "lpurple", "vlpurple", "vvlpurple", "vdorange", "dorange", "orange", "lorange")
names(ecoli_color) = gsub("_[0-9]+$", "", colnames(ecoli_m), fixed = F)
##### CHECK that no. of colors below in ecoli_color should match the no. of unique genomes in ecoli_df and with correct color coding shades #############
num_unique_genome = ncol(ecoli_df) - 1
length(ecoli_color) == num_unique_genome

lapply(1:num_unique_genome, function(i){
  out_file_name = paste0(prefix, ".highlight.", i, ".txt")
  ends = ecoli_df$position[ecoli_m[, i]]
  diff_ends = diff(ends) # Difference between consecutive mutations found on the genome
  sel_ends = c(TRUE, diff_ends > gap_size)  
  ends = ends[sel_ends] # If the gap size > diff_ends then it will combine those mutations
  ends = c(ends, genome_end_pos)
  starts = c(1, (ends + gap_size)[-length(ends)])
  df = data.frame(chr = "chr1", starts = starts, ends = ends, col = paste0("fill_color=", ecoli_color[i]))
  write.table(df, paste0(prefix, "/", out_file_name), row.names = F, col.names = F, quote = F)
})

############### highlight file to mark genes' start and end positions in circos

input_name_gene = "bac_genes_for_circos_unique.csv"
gene_file_name = paste0(prefix, ".highlight.genes.txt")
#gap_size = 4000

ecoli_genes = read.table(input_name_gene, header = T, sep = ",")
ecoli_genes_df = ecoli_genes[,]
ecoli_genes_df$start_position = as.numeric(gsub(",", "", ecoli_genes_df$start_position))
ecoli_genes_df$end_position = as.numeric(gsub(",", "", ecoli_genes_df$end_position))

num_genes = nrow(ecoli_genes_df)

file.create(paste0(prefix, "/", gene_file_name))
lapply(1:num_genes, function(i){
  starts = ecoli_genes_df$start_position[i]
  ends = ecoli_genes_df$end_position[i]
  df = data.frame(chr = "chr1", starts = starts, ends = ends, col = paste0("fill_color=", "lgrey"))
  write.table(df, paste0(prefix, "/", gene_file_name), row.names = F, col.names = F, quote = F, append=TRUE)
})

################ highlight file to mark whole population mutations in circos

input_name_pop = "bac_whole_pop_for_circos.csv"
pop_file_name = paste0(prefix, ".highlight.ideogram.txt")
gap_size = 3000

pop_mut = read.table(input_name_pop, header = T, sep = ",")
pop_mut$mutation_position = as.numeric(gsub(",", "", pop_mut$mutation_position))
#num_pop_mut = nrow(pop_mut)

file.create(paste0(prefix, "/", pop_file_name))
starts = pop_mut$mutation_position
df = data.frame(chr = "chr1", starts = starts, ends = starts + gap_size, col = paste0("fill_color=gneg"))
# Marking the manZ and malT whole pop mutations that are present in the dominant lineage of Day 37
df$col[match(1882915,df$starts)]="fill_color=dorange"
df$col[match(3482802,df$starts)]="fill_color=dorange"
write.table(df, paste0(prefix, "/", pop_file_name), row.names = F, col.names = F, quote = F)

########################################################### phage



prefix = "phage"
input_name = "phage_uniq_genome_table.tsv"
label_file_name = paste0(prefix, ".labels.mutations.txt")
genome_end_pos = 42507
gap_size = 50

ecoli = read.table(input_name, header = T, sep = "\t")
print(unique(ecoli$gene))
ecoli_df = ecoli[, c(1,2,4:37)]
ecoli_df$position = as.numeric(gsub(",", "", ecoli_df$position))
ecoli_df$mutation = gsub("→", "->", ecoli_df$mutation)
ecoli_df$mutation = gsub("Δ", "delta_", ecoli_df$mutation)
ecoli_df$mutation = gsub(".bp", "_bp", ecoli_df$mutation, fixed=F)

ecoli_nonsyn_only_idx = ecoli$nonsyn_mutation == 1

labels_df = data.frame(chr = "chr1", start = ecoli_df$position, end = ecoli_df$position + gap_size, label = ecoli_df$mutation)
labels_df_sub = labels_df[ecoli_nonsyn_only_idx, ]
write.table(labels_df_sub, paste0(prefix, "/", label_file_name), row.names = F, col.names = F, quote = F)

# sort unique genomes
ecoli_df = ecoli_df[, -2]
ecoli_m = ecoli_df[, -1] == "100%"
ecoli_color = c("greens-8-seq-8","greens-8-seq-7","greens-8-seq-6","greens-8-seq-5","greens-8-seq-4","greens-8-seq-3","greens-8-seq-2","greens-8-seq-1","blues-8-seq-8","blues-8-seq-7","blues-8-seq-6","blues-8-seq-5","blues-8-seq-4","blues-8-seq-3","blues-8-seq-2","blues-8-seq-1","vvdpurple","vdpurple","dpurple","purple","lpurple","vlpurple","vvlpurple","oranges-9-seq-9","oranges-8-seq-8","oranges-9-seq-8","oranges-8-seq-7","oranges-9-seq-7","oranges-8-seq-6","oranges-9-seq-6","oranges-9-seq-5","oranges-9-seq-4","oranges-9-seq-3","oranges-9-seq-2")
names(ecoli_color) = gsub("_[0-9]+$", "", colnames(ecoli_m), fixed = F)


num_unique_genome = ncol(ecoli_df) - 1

lapply(1:num_unique_genome, function(i){
  out_file_name = paste0(prefix, ".highlight.", i, ".txt")
  ends = ecoli_df$position[ecoli_m[, i]]
  diff_ends = diff(ends)
  sel_ends = c(TRUE, diff_ends > gap_size)
  ends = ends[sel_ends]
  ends = c(ends, genome_end_pos)
  starts = c(1, (ends + gap_size)[-length(ends)])
  df = data.frame(chr = "chr1", starts = starts, ends = ends, col = paste0("fill_color=", ecoli_color[i]))
  write.table(df, paste0(prefix, "/", out_file_name), row.names = F, col.names = F, quote = F)
})

############### highlight file to mark genes' start and end positions in circos

input_name_gene = "phage_genes_for_circos_unique.csv"
gene_file_name = paste0(prefix, ".highlight.genes.txt")
#gap_size = 4000

ecoli_genes = read.table(input_name_gene, header = T, sep = ",")
ecoli_genes_df = ecoli_genes[,]
ecoli_genes_df$start_position = as.numeric(gsub(",", "", ecoli_genes_df$start_position))
ecoli_genes_df$end_position = as.numeric(gsub(",", "", ecoli_genes_df$end_position))

num_genes = nrow(ecoli_genes_df)

file.create(paste0(prefix, "/", gene_file_name))
lapply(1:num_genes, function(i){
  starts = ecoli_genes_df$start_position[i]
  ends = ecoli_genes_df$end_position[i]
  df = data.frame(chr = "chr1", starts = starts, ends = ends, col = paste0("fill_color=", "lgrey"))
  write.table(df, paste0(prefix, "/", gene_file_name), row.names = F, col.names = F, quote = F, append=TRUE)
})

################ highlight file to mark whole population mutations in circos

input_name_pop = "phage_whole_pop_for_circos.csv"
pop_file_name = paste0(prefix, ".highlight.ideogram.txt")
gap_size = 20

pop_mut = read.table(input_name_pop, header = T, sep = ",")
pop_mut$mutation_position = as.numeric(gsub(",", "", pop_mut$mutation_position))
#num_pop_mut = nrow(pop_mut)

file.create(paste0(prefix, "/", pop_file_name))
starts = pop_mut$mutation_position
df = data.frame(chr = "chr1", starts = starts, ends = starts + gap_size, col = paste0("fill_color=gneg"))
# Marking the H mutation in whole population sequencing that is present in the dominant lineage of Day 28
df$col[match(11451,df$starts)]="fill_color=red"
df$ends[match(11451,df$starts)]= df$ends[match(11451,df$starts)] + 10 # adding extra width to the H mutation because we want to highlight in the figure 
write.table(df, paste0(prefix, "/", pop_file_name), row.names = F, col.names = F, quote = F)
