options(stringsAsFactors = F)

library(ggplot2)
library(gridExtra)

pheno_raw_m = pheno_matrix
pheno_raw_m = data.frame(bac_name = rownames(pheno_raw_m), pheno_raw_m)
pheno_raw_l = melt(pheno_raw_m)
colnames(pheno_raw_l) = c("bac_name", "phage_name", "eop")
pheno_raw_l = data.frame(pheno_raw_l, bac_day = as.numeric(gsub("_", "", substr(pheno_raw_l$bac_name, 5, 6))), phage_day = as.numeric(gsub("_", "", substr(pheno_raw_l$phage_name, 5, 6))))

### bac
bac_days = sort(unique(pheno_raw_l$bac_day))
bac_plot_names = paste0("bp_", bac_days)

for(idx in 1:length(bac_days)){
  i = bac_days[idx]
  sub_bac_i = pheno_raw_l[pheno_raw_l$bac_day == i, ]
  sub_bac_i_mean_each = aggregate(sub_bac_i$eop, by = list(sub_bac_i$bac_name, sub_bac_i$phage_day), mean)
  colnames(sub_bac_i_mean_each) = c("bac_name", "phage_day", "meop")
  
  # i_bac_anova = aov(eop ~ phage_day, data = sub_bac_i)
  # i_bac_anova_res = summary(i_bac_anova)
  # i_bac_pv = format(i_bac_anova_res[[1]]$`Pr(>F)`[1], digits = 4)
  # i_bac_f = format(i_bac_anova_res[[1]]$`F value`[1], digits = 4)
  
  bac_phage_days_i = sort(unique(sub_bac_i$phage_day))
  sub_bac_i_list = split(sub_bac_i, f=sub_bac_i$phage_day)
  tt_p_list = sapply(1:(length(bac_phage_days_i)-1), function(j){
    tt = t.test(sub_bac_i_list[[j]]$eop, sub_bac_i_list[[length(sub_bac_i_list)]]$eop, alternative = "less", paired = TRUE)$p.value
  })
  tt_p_max = format(max(tt_p_list), digits = 4)
  
  sub_bac_i_mean_all = aggregate(sub_bac_i$eop, by = list(sub_bac_i$phage_day), mean)
  sub_bac_i_sd_all = aggregate(sub_bac_i$eop, by = list(sub_bac_i$phage_day), sd)
  sub_bac_i_mean_all = data.frame(sub_bac_i_mean_all, sub_bac_i_sd_all[, 2])
  colnames(sub_bac_i_mean_all) = c("phage_day", "meop", "sd")
  
  ploti = ggplot() + 
    geom_line(data = sub_bac_i_mean_each, aes(x = phage_day, y = meop, group = bac_name), linetype = 3, color = "gray50") + 
    geom_point(data = sub_bac_i_mean_each, aes(x = phage_day, y = meop, group = bac_name), color = "gray50") + 
    geom_vline(xintercept = i, linetype = 2) +
    geom_line(data = sub_bac_i_mean_all, aes(x = phage_day, y = meop), size = 1.1) + 
    geom_point(data = sub_bac_i_mean_all, aes(x = phage_day, y = meop), size = 2) +
    geom_errorbar(data = sub_bac_i_mean_all, aes(x = phage_day, ymin = meop-sd, ymax = meop+sd), width=0.5, position=position_dodge(0.05)) + 
    coord_cartesian(xlim = c(7, 29), ylim = c(0, 3.6)) + 
    scale_x_continuous(breaks = c(8, 15, 22, 28)) + 
    scale_y_continuous(breaks = c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5)) +
    labs(x = "Phage sample day", y = "Mean EOP value", title = paste0("Bacteria Day ", i)) +
    annotate("text", x = 9, y = 3.5, label = paste0("P-value<", tt_p_max), fontface = 1, hjust = 0, size = 2) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  assign(bac_plot_names[idx], ploti)
  
  if(i==22){
    bac_tmp = ggplot() + 
                geom_vline(xintercept = 18, linetype = 2) + 
                geom_line(data = sub_bac_i_mean_all, aes(x = phage_day, y = meop), size = 1.1) +
                coord_cartesian(xlim = c(7, 29), ylim = c(0, 3.6)) + 
                scale_x_continuous(breaks = c(8, 15, 22, 28)) + 
                labs(x = "Phage sample day", y = "Mean EOP value", title = "Schematic of time-shift analysis ") + 
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
  }
}

# pdf("fig3_bac_timeshift.pdf", height = 9, width = 6, useDingbats = F)
# grid.arrange(bac_tmp, bp_8, bp_15, bp_22, bp_28, bp_37, nrow = 3)
# dev.off()



### phage
phage_days = sort(unique(pheno_raw_l$phage_day))
phage_plot_names = paste0("pp_", phage_days)

for(idx in 1:length(phage_days)){
  i = phage_days[idx]
  sub_phage_i = pheno_raw_l[pheno_raw_l$phage_day == i, ]
  sub_phage_i_mean_each = aggregate(sub_phage_i$eop, by = list(sub_phage_i$phage_name, sub_phage_i$bac_day), mean)
  colnames(sub_phage_i_mean_each) = c("phage_name", "bac_day", "meop")
  
  #Anova test, but for day 8 phages, there is no variance so need to be treated separately
  # if(idx != 1){
  #   i_phage_anova = aov(eop ~ bac_day, data = sub_phage_i)
  #   i_phage_anova_res = summary(i_phage_anova)
  #   i_phage_pv = format(i_phage_anova_res[[1]]$`Pr(>F)`[1], digits = 4)
  #   i_phage_f = format(i_phage_anova_res[[1]]$`F value`[1], digits = 4)
  # }else{
  #   i_phage_pv = "NA"
  #   i_phage_f = "NA"
  # }
  
  if(idx != 1){
    phage_bac_days_i = sort(unique(sub_phage_i$bac_day))
    sub_phage_i_list = split(sub_phage_i, f=sub_phage_i$bac_day)
    tt_p_list = sapply(2:length(phage_bac_days_i), function(j){
      tt = t.test(sub_phage_i_list[[1]]$eop, sub_phage_i_list[[j]]$eop, alternative = "greater", paired = TRUE)$p.value
    })
    tt_p_max = format(max(tt_p_list), digits = 4)
  }else{
    tt_p_max = "NA"
  }
  
  sub_phage_i_mean_all = aggregate(sub_phage_i$eop, by = list(sub_phage_i$bac_day), mean)
  sub_phage_i_sd_all = aggregate(sub_phage_i$eop, by = list(sub_phage_i$bac_day), sd)
  sub_phage_i_mean_all = data.frame(sub_phage_i_mean_all, sub_phage_i_sd_all[, 2])
  colnames(sub_phage_i_mean_all) = c("bac_day", "meop", "sd")
  
  plotphage = ggplot() + 
    geom_line(data = sub_phage_i_mean_each, aes(x = bac_day, y = meop, group = phage_name), linetype = 3, color = "gray50") + 
    geom_point(data = sub_phage_i_mean_each, aes(x = bac_day, y = meop, group = phage_name), color = "gray50") + 
    geom_vline(xintercept = i, linetype = 2) +
    geom_line(data = sub_phage_i_mean_all, aes(x = bac_day, y = meop), size = 1.1) + 
    geom_point(data = sub_phage_i_mean_all, aes(x = bac_day, y = meop), size = 2) +
    geom_errorbar(data = sub_phage_i_mean_all, aes(x = bac_day, ymin = meop-sd, ymax = meop+sd), width=0.5, position=position_dodge(0.05)) + 
    coord_cartesian(xlim = c(7, 38), ylim = c(0, 2.5)) + 
    scale_x_continuous(breaks = c(8, 15, 22, 28, 37)) + 
    annotate("text", x = 29, y = 2.5, label = paste0("P-value<", tt_p_max), fontface = 1, hjust = 0, size = 2) + 
    labs(x = "Bacteria sample day", y = "Mean EOP value", title = paste0("Phage Day ", i)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  assign(phage_plot_names[idx], plotphage)
}

# pdf("fig3_phage_timeshift.pdf", height = 6, width = 6, useDingbats = F)
# grid.arrange(pp_8, pp_15, pp_22, pp_28, nrow = 2)
# dev.off()

pdf("fig3_bp_comb_timeshift.pdf", height = 6, width = 15, useDingbats = F)
grid.arrange(bac_tmp, pp_8, pp_15, pp_22, pp_28, bp_8, bp_15, bp_22, bp_28, bp_37, nrow = 2)
dev.off()









# pheno_raw_l_diff_added = data.frame(pheno_raw_l, bac_diff = pheno_raw_l$bac_day - pheno_raw_l$phage_day, phage_diff = pheno_raw_l$phage_day - pheno_raw_l$bac_day)
# 
# bac_view_ts = aggregate(pheno_raw_l_diff_added$eop, by = list(pheno_raw_l_diff_added$bac_name, pheno_raw_l_diff_added$bac_diff), mean)
# colnames(bac_view_ts) = c("bac_name", "bac_diff", "meop")
# bac_mean_all = aggregate(pheno_raw_l_diff_added$eop, by = list(pheno_raw_l_diff_added$bac_diff), mean)
# bac_sd_all = aggregate(pheno_raw_l_diff_added$eop, by = list(pheno_raw_l_diff_added$bac_diff), sd)
# bac_mean_all = data.frame(bac_mean_all, bac_sd_all[, 2])
# colnames(bac_mean_all) = c("bac_diff", "meop", "sd")
# 
# bac_meta_plot = ggplot() + 
#   geom_line(data = bac_view_ts, aes(x = bac_diff, y = meop, group = bac_name), linetype = 3, color = "gray50") + 
#   geom_point(data = bac_view_ts, aes(x = bac_diff, y = meop, group = bac_name), color = "gray50") + 
#   geom_vline(xintercept = 0, linetype = 2) +
#   geom_line(data = bac_mean_all, aes(x = bac_diff, y = meop), size = 1.1) + 
#   geom_point(data = bac_mean_all, aes(x = bac_diff, y = meop), size = 2) +
#   geom_errorbar(data = bac_mean_all, aes(x = bac_diff, ymin = meop-sd, ymax = meop+sd), width=0.5, position=position_dodge(0.05)) + 
#   coord_cartesian(xlim = c(-21, 30)) + 
#   scale_x_continuous(breaks = c(-20, -14, -7, 0, 7, 14, 22, 29)) + 
#   labs(x = "Phage sampling day difference", y = "Mean EOP value") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
# 
# 
# phage_view_ts = aggregate(pheno_raw_l_diff_added$eop, by = list(pheno_raw_l_diff_added$phage_name, pheno_raw_l_diff_added$phage_diff), mean)
# colnames(phage_view_ts) = c("phage_name", "phage_diff", "meop")
# phage_mean_all = aggregate(pheno_raw_l_diff_added$eop, by = list(pheno_raw_l_diff_added$phage_diff), mean)
# phage_sd_all = aggregate(pheno_raw_l_diff_added$eop, by = list(pheno_raw_l_diff_added$phage_diff), sd)
# phage_mean_all = data.frame(phage_mean_all, phage_sd_all[, 2])
# colnames(phage_mean_all) = c("phage_diff", "meop", "sd")
# 
# phage_meta_plot = ggplot() + 
#   geom_line(data = phage_view_ts, aes(x = phage_diff, y = meop, group = phage_name), linetype = 3, color = "gray50") + 
#   geom_point(data = phage_view_ts, aes(x = phage_diff, y = meop, group = phage_name), color = "gray50") + 
#   geom_vline(xintercept = 0, linetype = 2) +
#   geom_line(data = phage_mean_all, aes(x = phage_diff, y = meop), size = 1.1) + 
#   geom_point(data = phage_mean_all, aes(x = phage_diff, y = meop), size = 2) +
#   geom_errorbar(data = phage_mean_all, aes(x = phage_diff, ymin = meop-sd, ymax = meop+sd), width=0.5, position=position_dodge(0.05)) + 
#   coord_cartesian(xlim = c(-30, 21)) + 
#   scale_x_continuous(breaks = c(-29, -22, -14, -7, 0, 7, 14, 20)) + 
#   labs(x = "Host sampling day difference", y = "Mean EOP value") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
# 
# pdf("fig3_supp.pdf", height = 3, width = 6, useDingbats = F)
# grid.arrange(bac_meta_plot, phage_meta_plot, nrow = 1)
# dev.off()

