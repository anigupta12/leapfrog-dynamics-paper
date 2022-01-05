#run approximation6_uniq_nonsyn and eqtl_trial first
library(plyr)

raw_list_noan_nodash = raw_list_noan
raw_list_noan_nodash[, 1] = gsub("-", "_", raw_list_noan_nodash[, 1])
raw_list_noan_nodash[, 2] = gsub("-", "_", raw_list_noan_nodash[, 2])

ptb = phage_only_eqtl_m_allmut_allgen[-1, -1]
btb = bac_only_eqtl_m_allmut_allgen[-1, -1]

phage_mut = data.frame(as.numeric(gsub("_", "", substr(colnames(ptb), 5, 6))), colSums(ptb))
phage_mut[, ncol(phage_mut) + 1] = rownames(phage_mut)
colnames(phage_mut) = c("time", "num_mut", "name")

bac_mut = data.frame(as.numeric(gsub("_", "", substr(colnames(btb), 5, 6))), colSums(btb))
bac_mut[, ncol(bac_mut) + 1] = rownames(bac_mut)
colnames(bac_mut) = c("time", "num_mut", "name")


phage_aggregate <- ddply(raw_list_noan_nodash, ~Phage, summarise, num_eq_0 = sum(EOP == 0), num_gt_0 = sum(EOP > 0))
colnames(phage_aggregate) = c("name", "eq0", "gt0")

bac_aggregate <- ddply(raw_list_noan_nodash, ~Host, summarise, num_eq_0 = sum(EOP == 0), num_gt_0 = sum(EOP > 0))
colnames(bac_aggregate) = c("name", "num_eq_0", "num_gt_0")


phage_comb <- merge(phage_mut, phage_aggregate, by = "name")

bac_comb <- merge(bac_mut, bac_aggregate, by = "name")


phage_host_range <- ddply(phage_comb, ~ time, summarise, host_range = sum(gt0)/(11*50))
phage_host_range <- rbind(c(0, 0), phage_host_range)
bac_resistance <- ddply(bac_comb, ~ time, summarise, resistance = sum(num_eq_0)/(10*44))
bac_resistance <- rbind(c(0,0), bac_resistance)

# plot(phage_host_range[, 1], phage_host_range[, 2], xlab = "Time (days)", ylab = "Host Range", type = "l", lty = 2, xlim = c(0, 40), ylim = c(0, 400), main = "phage host range pairs by date", axes = F) + axis(1, at = c(0, 8, 15, 22, 28, 37)) + axis(2, at = c(0, 100, 200, 300, 400))
# 
# plot(bac_host_range[, 1], bac_host_range[, 2], xlab = "Time (days)", ylab = "Resistance", type = "l", xlim = c(0, 40), ylim = c(0, 400), main = "bacteria resistant pairs by date", axes = F) + axis(1, at = c(0, 8, 15, 22, 28, 37)) + axis(2, at = c(0, 100, 200, 300, 400))

pdf("fig1b1.pdf", height = 6.1, width = 5, useDingbats = F)
par(mar = c(5.1, 2.1, 4.1, 5.1))
plot(phage_host_range[, 1], phage_host_range[, 2], xlab = "Time (days)", ylab = "Host Range Percentage", type = "l", lty = 2, xlim = c(0, 40), ylim = c(0, 1), main = "phage host range pairs by date", axes = F) + axis(1, at = c(0, 8, 15, 22, 28, 37)) + axis(4, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), las = 1) + text(par("usr")[2]*1.2, mean(par("usr")[3:4]*1.3), "Host Range Percentage", srt = -90, xpd = TRUE, pos = 4)
# + mtext("Resistance Percentage", 4, line = 3)
dev.off()

pdf("fig1b2.pdf", height = 6.1, width = 5, useDingbats = F)
par(mar = c(5.1, 4.1, 4.1, 2.1))
plot(bac_resistance[, 1], bac_resistance[, 2], xlab = "Time (days)", ylab = "Resistance Percentage", type = "l", xlim = c(0, 40), ylim = c(0, 1), main = "bacteria resistant pairs by date", axes = F) + axis(1, at = c(0, 8, 15, 22, 28, 37)) + axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), las = 1)
dev.off()



