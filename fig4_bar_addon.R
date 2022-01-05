options(stringsAsFactors = FALSE)

library(ggplot2)
library(grDevices)

#bac
bactree_time_comp = read.delim("bac_tree_compo.tsv", header = T)
bactree_time_comp$branch = factor(bactree_time_comp$branch, levels = c("black", "red", "orange"))
bactree_time_comp$count = bactree_time_comp$count/10

pdf("fig4_bactree_comp.pdf", height = 4, width = 4, useDingbats = F)
ggplot() + 
  geom_bar(data = bactree_time_comp, aes(y = count, x = day, fill = branch), stat = "identity", show.legend = F) +
  scale_fill_manual(values = c("black", "red", "orange")) + 
  xlim(0, 40) + 
  labs(x = "", y = "", title = "") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()


#phage
phagetree_time_comp = read.delim("phage_tree_compo.tsv", header = T)
phagetree_time_comp$branch = factor(phagetree_time_comp$branch, levels = c("black", "green", "blue"))
phagetree_time_comp$count = phagetree_time_comp$count/11

pdf("fig4_phagetree_comp.pdf", height = 4, width = 4, useDingbats = F)
ggplot() + 
  geom_bar(data = phagetree_time_comp, aes(y = count, x = day, fill = branch), stat = "identity", show.legend = F) +
  scale_fill_manual(values = c("black", "forestgreen", "deepskyblue3")) + 
  xlim(0, 40) + 
  labs(x = "", y = "", title = "") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()
