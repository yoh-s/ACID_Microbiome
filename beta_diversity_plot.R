library(ggplot2)
library(ape)
library(tidyverse)
library(multcomp)
library(extrafont)  ##intalled 09/13/2024
library(patchwork)
library(vegan)
library(ggtern)
library(reshape2)
library(ggsci)
library(scales)

##subset Baseline level
base_study8_rare <- subset_samples(ps_study8_rare, Time == "Baseline") ##321

base_study8_rare <- prune_taxa(taxa_sums(base_study8_rare) > 0, base_study8_rare)  ##296

#Bray-Curtis at the ASV level
base_genus_bray

pcoa <- pcoa(base_genus_bray,correction = "none", rn = NULL) 

groups <- base_ps_genus_rare ##sample data
PC1 <- pcoa$vectors[,1]
PC2 <- pcoa$vectors[,2]

groups<-groups[match(rownames(pcoa$vectors),rownames(groups)),]
pcoadata <- data.frame(rownames(pcoa$vectors),
                       PC1,PC2,groups$Age_quartile)
colnames(pcoadata) <- c("sample","PC1","PC2","group")

pcoadata$group <- factor(pcoadata$group,levels = c("60-64_years","65-68_years","68-73_years","73-80_years"))

#Adonis test
set.seed(1988)
otu.adonis <- adonis2(base_rare_bray~group,data = pcoadata,permutations = 9999)
write.table(otu.adonis,'beta_bray-curtis_adonis.tsv',sep = '\t',quote = F)

#Boxplot
df <- pcoadata
df1 <- df %>% group_by(group) %>% summarise(Max = max(PC1))
df2 <- df %>% group_by(group) %>% summarise(Max = max(PC2))
df1$Max <- df1$Max + max(df1$Max)*0.1
df2$Max <- df2$Max + max(df2$Max)*0.1

tuk1 <- aov(PC1~group,data = pcoadata) %>% 
  glht(alternative = 'two.sided',linfct = mcp(group = "Tukey")) %>% cld(alpah = 0.05)
#mod<-aov(PC1~group,data = pcoadata)
#glt<-glht(mod, alternative = 'two.sided', linfct = mcp(group = 'Tukey'))
#summary(glt)
tuk2 <- aov(PC2~group,data = pcoadata) %>% 
  glht(alternative = 'two.sided',linfct = mcp(group = "Tukey")) %>% cld(alpah = 0.05)

res <- data.frame(PC1 = tuk1$mcletters$Letters,PC2 = tuk2$mcletters$Letters,
                  df1 = df1$Max,df2 = df2$Max,group = df1$group)

p1 <- ggplot(pcoadata, aes(PC1, PC2,colour = group)) +
  geom_point(aes(colour = group),size = 2)+
  labs(x = (floor(pcoa$values$Relative_eig[1]*100)) %>% 
         paste0("PCoA1 ( ", ., "%", " )"),
       y = (floor(pcoa$values$Relative_eig[2]*100)) %>% 
         paste0("PCoA2 ( ", ., "%", " )")) +
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD")) +
  theme(text = element_text(size = 15))+
  geom_vline(aes(xintercept = 0),linetype = "dotted")+
  geom_hline(aes(yintercept = 0),linetype = "dotted")+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.title.x = element_text(colour = 'black', size = 10),
        axis.title.y = element_text(colour = 'black', size = 10),
        axis.text = element_text(colour = 'black',size = 10),
        legend.title = element_blank(),
        legend.key.height = unit(0.4,"cm"),
        legend.position = 'none',
        panel.grid = element_blank())

# Rename factor levels for better readability
levels(pcoadata$group) <- gsub("_", " ", levels(pcoadata$group))

p2 <- ggplot(pcoadata,aes(group,PC1)) +
  geom_boxplot(aes(fill = group),outlier.colour = NA)+
  scale_fill_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD")) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_blank(),
        axis.line = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(colour = 'black',size = 10,face = "plain"),
        axis.text.x = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())+coord_flip()
p3 <- ggplot(pcoadata,aes(group,PC2)) +
  geom_boxplot(aes(fill = group),outlier.colour = NA) +
  scale_fill_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD")) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_blank(),
        axis.line = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = 'black',size = 12,angle = 90,vjust = 0.5,hjust = 0.5,face = "plain"),
        axis.text.y = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())
p4 <- ggplot(pcoadata,aes(PC1, PC2))+
  geom_text(aes(x = -0.5,y = 0.6,
                label = paste("PERMANOVA:\ndf = ",braygenus_agequartile$Df[1],"\nR2 = ",
                              round(braygenus_agequartile$R2[1],4),"\np-value = ",
                              braygenus_agequartile$`Pr(>F)`[1],
                              sep = "")),size = 2.8,family = "sans",fontface = 1)+
  theme_bw() + xlab(NULL) + ylab(NULL) +
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

p <- p2 + plot_spacer() + p1 + p3 + 
  plot_layout(heights = c(1, 4), widths = c(4, 1), ncol = 2, nrow = 2)
ggsave(plot = p,"area_plot/finaltop10/bray_curtis_agequartile.tiff",dpi = 600,width = 6,height = 5)

base_rare_wunifrac

pcoa <- pcoa(base_rare_wunifrac,correction = "none", rn = NULL) 

groups <- base_ps_genus_rare ##sample data
PC1 <- pcoa$vectors[,1]
PC2 <- pcoa$vectors[,2]

groups<-groups[match(rownames(pcoa$vectors),rownames(groups)),]
pcoadata <- data.frame(rownames(pcoa$vectors),
                       PC1,PC2,groups$Age_quartile)
colnames(pcoadata) <- c("sample","PC1","PC2","group")

pcoadata$group <- factor(pcoadata$group,levels = c("60-64_years","65-68_years","68-73_years","73-80_years"))

#Adonis test
set.seed(1988)
otu.adonis <- adonis2(base_rare_wunifrac~group,data = pcoadata,permutations = 999)
write.table(otu.adonis,'beta_bray-curtis_adonis.tsv',sep = '\t',quote = F)

#Boxplot
df <- pcoadata
df1 <- df %>% group_by(group) %>% summarise(Max = max(PC1))
df2 <- df %>% group_by(group) %>% summarise(Max = max(PC2))
df1$Max <- df1$Max + max(df1$Max)*0.1
df2$Max <- df2$Max + max(df2$Max)*0.1

tuk1 <- aov(PC1~group,data = pcoadata) %>% 
  glht(alternative = 'two.sided',linfct = mcp(group = "Tukey")) %>% cld(alpah = 0.05)
#mod<-aov(PC1~group,data = pcoadata)
#glt<-glht(mod, alternative = 'two.sided', linfct = mcp(group = 'Tukey'))
#summary(glt)
tuk2 <- aov(PC2~group,data = pcoadata) %>% 
  glht(alternative = 'two.sided',linfct = mcp(group = "Tukey")) %>% cld(alpah = 0.05)

res <- data.frame(PC1 = tuk1$mcletters$Letters,PC2 = tuk2$mcletters$Letters,
                  df1 = df1$Max,df2 = df2$Max,group = df1$group)

p1 <- ggplot(pcoadata, aes(PC1, PC2,colour = group)) +
  geom_point(aes(colour = group),size = 2)+
  labs(x = (floor(pcoa$values$Relative_eig[1]*100)) %>% 
         paste0("PCoA1 ( ", ., "%", " )"),
       y = (floor(pcoa$values$Relative_eig[2]*100)) %>% 
         paste0("PCoA2 ( ", ., "%", " )")) +
  scale_colour_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD")) +
  theme(text = element_text(size = 15))+
  geom_vline(aes(xintercept = 0),linetype = "dotted")+
  geom_hline(aes(yintercept = 0),linetype = "dotted")+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.title.x = element_text(colour = 'black', size = 10),
        axis.title.y = element_text(colour = 'black', size = 10),
        axis.text = element_text(colour = 'black',size = 10),
        legend.title = element_blank(),
        legend.key.height = unit(0.4,"cm"),
        legend.position = 'none',
        panel.grid = element_blank())

# Rename factor levels for better readability
levels(pcoadata$group) <- gsub("_", " ", levels(pcoadata$group))

p2 <- ggplot(pcoadata,aes(group,PC1)) +
  geom_boxplot(aes(fill = group),outlier.colour = NA)+
  scale_fill_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD")) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_blank(),
        axis.line = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(colour = 'black',size = 10,face = "plain"),
        axis.text.x = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())+coord_flip()
p3 <- ggplot(pcoadata,aes(group,PC2)) +
  geom_boxplot(aes(fill = group),outlier.colour = NA) +
  scale_fill_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD")) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_blank(),
        axis.line = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = 'black',size = 12,angle = 90,vjust = 0.5,hjust = 0.5,face = "plain"),
        axis.text.y = element_blank(),
        legend.position = "none",
        panel.grid = element_blank())
p4 <- ggplot(pcoadata,aes(PC1, PC2))+
  geom_text(aes(x = -0.5,y = 0.6,
                label = paste("PERMANOVA:\ndf = ",braygenus_agequartile$Df[1],"\nR2 = ",
                              round(braygenus_agequartile$R2[1],4),"\np-value = ",
                              braygenus_agequartile$`Pr(>F)`[1],
                              sep = "")),size = 2.8,family = "sans",fontface = 1)+
  theme_bw() + xlab(NULL) + ylab(NULL) +
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

p <- p2 + plot_spacer() + p1 + p3 + 
  plot_layout(heights = c(1, 4), widths = c(4, 1), ncol = 2, nrow = 2)
ggsave(plot = p,"area_plot/finaltop10/weightedunifrac_agequartile.tiff",dpi = 600,width = 6,height = 5)

