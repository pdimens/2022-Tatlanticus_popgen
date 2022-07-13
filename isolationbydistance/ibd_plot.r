library(ggplot2)
library(data.table)

ibd <- read.csv("../Desktop/ibd.csv")
distnames <- as.character(ibd[2,])
ibd <- ibd[-c(1:2),]
cn <- ibd[,1]
ibd <- ibd[,-1]
ibd[] <- lapply(ibd, function(x) as.numeric(x))
ibd <- transpose(ibd)
colnames(ibd) <- cn
dn <- distnames[2:length(distnames)]
ibd$distance <- factor(dn, ordered = T, levels = dn)
colors <- c("estimate" = "dodgerblue", "95% null CI" = "#aa6699", "95% bootstrap CI" = "grey60")

ggplot(data = ibd) +
  geom_hline(yintercept = 0, color = "grey80") +
  geom_point(aes(x = 2, y = r, color = "estimate"), size = 2) +
  geom_segment(aes(x = 1, xend=1, y = U , yend = L, color = "95% null CI"), size = 1.5) +
  geom_segment(aes(x = 3, xend=3, y = bootlow , yend = boothi, color = "95% bootstrap CI"), size = 1.5) +
  facet_grid(~distance) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()
    ) +
  labs(
    x = "",
    y = "r",
    color = "Legend") +
  coord_cartesian(c(-2,6),c(-0.0015, 0.0035), expand = F) +
  labs(title = "Combined r values at different distance (km) class sizes") +
  scale_color_manual(values = colors)
  
ggsave(filename = "ibdplot.png", height = 3, width = 12, units = "in")
