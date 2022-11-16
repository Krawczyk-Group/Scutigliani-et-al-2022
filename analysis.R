#### INSTALL PACKAGES ####

# install.packages

install.packages("ggplot2")
install.packages("tidyverse")
install.packages("msigdbr")
install.packages("VennDiagram")
install.packages("cowplot")
install.packages("DOSE")
install.packages("ggupset")
install.packages("pheatmap")
install.packages("RColorBrewer")
install.packages("matrixTests")
install.packages("ggfortify")
install.packages("ggrepel")
install.packages("UpSetR")

# bioconductor packages

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("genefilter")

# activate libraries

library(ggplot2)
library(tidyverse)
library(clusterProfiler)
library(msigdbr)
library(VennDiagram)
library(cowplot)
library(pheatmap)
library(RColorBrewer)
library(org.Hs.eg.db)
library(genefilter)
library(ggfortify)
library(ggrepel)

rm(list=ls())

# set working directory ...

#### IMPORT EXPRESSION DATA ####

df_ratio_long_filter <- read.csv("expression_all.csv", sep = ",", row.names = 1) %>%
  dplyr::rename("Publication" = publication, "Gene" = symbol, "Ratio" = FC)

# annotate long data with experimental info

df_ratio_long_filter_Tabuchi2008 <- df_ratio_long_filter %>%
  dplyr::filter(Publication=="Tabuchi2008") %>%
  dplyr::mutate(temperature_C = "41",
                duration_m = "30",
                timepoint_h = "3",
                CEM43 = "1.875") %>%
  distinct(Gene, .keep_all = T)

df_ratio_long_filter_Tabuchi2011a <- df_ratio_long_filter %>%
  dplyr::filter(Publication=="Tabuchi2011a") %>%
  dplyr::mutate(temperature_C = "42",
                duration_m = "15",
                timepoint_h = "0",
                CEM43 = "3.75") %>%
  distinct(Gene, .keep_all = T)

df_ratio_long_filter_Tabuchi2011b <- df_ratio_long_filter %>%
  dplyr::filter(Publication=="Tabuchi2011b") %>%
  dplyr::mutate(temperature_C = "42",
                duration_m = "15",
                timepoint_h = "1",
                CEM43 = "3.75") %>%
  distinct(Gene, .keep_all = T)

df_ratio_long_filter_Tabuchi2011c <- df_ratio_long_filter %>%
  dplyr::filter(Publication=="Tabuchi2011c") %>%
  dplyr::mutate(temperature_C = "42",
                duration_m = "15",
                timepoint_h = "3",
                CEM43 = "3.75") %>%
  distinct(Gene, .keep_all = T)

df_ratio_long_filter_Tabuchi2011d <- df_ratio_long_filter %>%
  dplyr::filter(Publication=="Tabuchi2011d") %>%
  dplyr::mutate(temperature_C = "42",
                duration_m = "15",
                timepoint_h = "6",
                CEM43 = "3.75") %>%
  distinct(Gene, .keep_all = T)

df_ratio_long_filter_Tabuchi2011e <- df_ratio_long_filter %>%
  dplyr::filter(Publication=="Tabuchi2011e") %>%
  dplyr::mutate(temperature_C = "44",
                duration_m = "15",
                timepoint_h = "0",
                CEM43 = "30") %>%
  distinct(Gene, .keep_all = T)

df_ratio_long_filter_Tabuchi2011f <- df_ratio_long_filter %>%
  dplyr::filter(Publication=="Tabuchi2011f") %>%
  dplyr::mutate(temperature_C = "44",
                duration_m = "15",
                timepoint_h = "1",
                CEM43 = "30") %>%
  distinct(Gene, .keep_all = T)

df_ratio_long_filter_Tabuchi2011g <- df_ratio_long_filter %>%
  dplyr::filter(Publication=="Tabuchi2011g") %>%
  dplyr::mutate(temperature_C = "44",
                duration_m = "15",
                timepoint_h = "3",
                CEM43 = "30") %>%
  distinct(Gene, .keep_all = T)

df_ratio_long_filter_Tabuchi2011h <- df_ratio_long_filter %>%
  dplyr::filter(Publication=="Tabuchi2011h") %>%
  dplyr::mutate(temperature_C = "44",
                duration_m = "15",
                timepoint_h = "6",
                CEM43 = "30") %>%
  distinct(Gene, .keep_all = T)

df_ratio_long_filter_Amaya2014a <- df_ratio_long_filter %>%
  dplyr::filter(Publication=="Amaya2014a") %>%
  dplyr::mutate(temperature_C = "45",
                duration_m = "30",
                timepoint_h = "4",
                CEM43 = "120") %>%
  distinct(Gene, .keep_all = T)

df_ratio_long_filter_Amaya2014b <- df_ratio_long_filter %>%
  dplyr::filter(Publication=="Amaya2014b") %>%
  dplyr::mutate(temperature_C = "45",
                duration_m = "30",
                timepoint_h = "4",
                CEM43 = "120") %>%
  distinct(Gene, .keep_all = T)

df_ratio_long_filter_Amaya2014c <- df_ratio_long_filter %>%
  dplyr::filter(Publication=="Amaya2014c") %>%
  dplyr::mutate(temperature_C = "45",
                duration_m = "30",
                timepoint_h = "4",
                CEM43 = "120") %>%
  distinct(Gene, .keep_all = T)

df_ratio_long_filter_Court2017 <- df_ratio_long_filter %>%
  dplyr::filter(Publication=="Court2017") %>%
  dplyr::mutate(temperature_C = "43",
                duration_m = "30",
                timepoint_h = "0",
                CEM43 = "30") %>%
  distinct(Gene, .keep_all = T)

df_ratio_long_filter_Andocs2015a <- df_ratio_long_filter %>%
  dplyr::filter(Publication=="Andocs2015a") %>%
  dplyr::mutate(temperature_C = "42",
                duration_m = "30",
                timepoint_h = "0",
                CEM43 = "7.5") %>%
  distinct(Gene, .keep_all = T)

df_ratio_long_filter_Andocs2015b <- df_ratio_long_filter %>%
  dplyr::filter(Publication=="Andocs2015b") %>%
  dplyr::mutate(temperature_C = "42",
                duration_m = "30",
                timepoint_h = "0",
                CEM43 = "7.5") %>%
  distinct(Gene, .keep_all = T)

df_ratio_long_filter_Yunoki2016 <- df_ratio_long_filter %>%
  dplyr::filter(Publication=="Yunoki2016") %>%
  dplyr::mutate(temperature_C = "44",
                duration_m = "90",
                timepoint_h = "24",
                CEM43 = "180") %>%
  distinct(Gene, .keep_all = T)

df_ratio_long_filter_Scutigliani2022a <- df_ratio_long_filter %>%
  dplyr::filter(Publication=="Scutigliani2022a") %>%
  dplyr::mutate(temperature_C = "42",
                duration_m = "60",
                timepoint_h = "6",
                CEM43 = "15") %>%
  distinct(Gene, .keep_all = T)

df_ratio_long_filter_Scutigliani2022b <- df_ratio_long_filter %>%
  dplyr::filter(Publication=="Scutigliani2022b") %>%
  dplyr::mutate(temperature_C = "42",
                duration_m = "60",
                timepoint_h = "24",
                CEM43 = "15") %>%
  distinct(Gene, .keep_all = T)

df_all <- rbind(df_ratio_long_filter_Tabuchi2008,
                df_ratio_long_filter_Tabuchi2011a,
                df_ratio_long_filter_Tabuchi2011b,
                df_ratio_long_filter_Tabuchi2011c,
                df_ratio_long_filter_Tabuchi2011d,
                df_ratio_long_filter_Tabuchi2011e,
                df_ratio_long_filter_Tabuchi2011f,
                df_ratio_long_filter_Tabuchi2011g,
                df_ratio_long_filter_Tabuchi2011h,
                df_ratio_long_filter_Amaya2014a,
                df_ratio_long_filter_Amaya2014b,
                df_ratio_long_filter_Amaya2014c,
                df_ratio_long_filter_Court2017,
                df_ratio_long_filter_Andocs2015a,
                df_ratio_long_filter_Andocs2015b,
                df_ratio_long_filter_Yunoki2016,
                df_ratio_long_filter_Scutigliani2022a,
                df_ratio_long_filter_Scutigliani2022b)

#### GLOBAL OVERVIEW - ALL DATASETS ####

# fetch DEGs per publication

df_all %>%
  group_by(Publication) %>%
  ggplot(., aes(x=Publication, y=log10(Ratio))) +
  geom_boxplot(fill="#E5E5E5", outlier.size = 0.5) +
  labs(title="All ratios",subtitle = "",xlab = "", ylab = "") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(NA),
        panel.grid.minor = element_line(NA),
        title = element_text(colour = "black", size=15),
        axis.title.x = element_text(colour = "black", size=15),
        axis.title.y = element_text(colour = "black", size=15),
        axis.line.x = element_line(colour = 'black', size = .5),
        axis.line.y = element_line(colour = 'black', size = .5),
        axis.ticks.x = element_line(colour = "black", size = 0),
        axis.ticks.y = element_line(colour = "black", size = 0),
        axis.text.x = element_text(colour = "black", size=9, angle = 90, hjust=1, vjust=0.5),
        axis.text.y = element_text(colour = "black", size=9),
        aspect.ratio = .5) +
  geom_hline(yintercept = 0, linetype = "dashed")

ggsave("Ratio_all_genes_all_studies.pdf", dpi=1000)


df_DEG <- df_all %>%
  group_by(Publication) %>%
  filter(Ratio >= 1.5 | Ratio <= 1/1.5)

for (i in 1:nrow(df_DEG))
  if (df_DEG$Ratio[i] > 1.5) df_DEG$Change[i] = "up" else df_DEG$Change[i] = "down"

ggplot(df_DEG, aes(x=Publication, fill=Change)) +
  geom_bar(position = "dodge") +
  labs(title="Differential genes", subtitle = "(1/1.5 > FC > 1.5)", xlab(""), ylab = "Count") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(NA),
        panel.grid.minor = element_line(NA),
        title = element_text(colour = "black", size=15),
        axis.title.x = element_text(colour = "black", size=15),
        axis.title.y = element_text(colour = "black", size=15),
        axis.line.x = element_line(colour = 'black', size = .5),
        axis.line.y = element_line(colour = 'black', size = .5),
        axis.ticks.x = element_line(colour = "black", size = 0),
        axis.ticks.y = element_line(colour = "black", size = 0),
        axis.text.x = element_text(colour = "black", size=9, angle = 90, hjust=1, vjust=0.5),
        axis.text.y = element_text(colour = "black", size=9),
        aspect.ratio = .5) +
  scale_fill_manual(values = c("#4169E1" ,"#FFA500"))

test <- df_DEG %>%
  ungroup() %>%
  group_by(Publication, Change) %>%
  summarize(n = n())

write.csv(test, "DEG_UP_DN_ALL.csv")

ggsave("all_DEG_barplot.pdf", dpi=1000)

dev.off()

#### TOP OVER/UNDEREXPRESSED GENES - ALL DATASETS ####

df_top <- df_all %>% ungroup() %>% group_by(Publication) %>% top_n(20, Ratio) %>% dplyr::arrange(Gene, Ratio)
df_bot <- df_all %>% ungroup() %>% group_by(Publication) %>% top_n(20, -Ratio) %>% dplyr::arrange(Gene, Ratio)

df_topbot <- rbind(df_top, df_bot)

textsize <- 8
theme <-   theme(plot.title = element_text(size = textsize, face = "bold", hjust = 0.5),
                 plot.subtitle = element_text(size = textsize, hjust = 0.5),
                 panel.background = element_rect(fill="white"),
                 panel.grid.major = element_line(NA),
                 panel.grid.minor = element_line(NA),
                 axis.title.x = element_text(colour = "black", size=textsize),
                 axis.title.y = element_text(colour = "black", size=textsize),
                 axis.line.x = element_line(colour = 'black', size = 0.5),
                 axis.line.y = element_line(colour = 'black', size = 0.5),
                 axis.ticks.x = element_line(colour = "black", size = 0.5),
                 axis.ticks.y = element_line(colour = "black", size = 0.5),
                 axis.text.x = element_text(colour = "black", size=textsize),
                 axis.text.y = element_text(colour = "black", size=textsize),
                 legend.justification = c("right", "top"),
                 legend.background = element_blank(),
                 legend.title = element_text(face = "bold", size = textsize),
                 legend.text = element_text(size = textsize, colour = "black"),
                 legend.key = element_blank(),
                 aspect.ratio = 1)

ggplot(df_topbot, aes(x = log10(Ratio), y = fct_reorder(Gene, -Ratio))) +
  geom_bar(stat = "identity") +
  facet_wrap(~Publication, scales = "free", ncol = 5) +
  theme

ggsave("top20_over_under_all.pdf", dpi = 500)

#### TOP OVER/UNDEREXPRESSED GENES - OVERLAP PER CLUSTER ####

df_topbot_antitabuchi <- df_topbot %>% dplyr::filter(Publication %in% c("Amaya2014",
                                                                        "Amaya2016",
                                                                        "Court2017",
                                                                        "Andocs2015a",
                                                                        "Andocs2015b",
                                                                        "Yokoni2016"))

df_topbot_tabuchi <- df_topbot %>% dplyr::filter(!Publication %in% c("Amaya2014",
                                                                     "Amaya2016",
                                                                     "Court2017",
                                                                     "Andocs2015a",
                                                                     "Andocs2015b",
                                                                     "Yokoni2016"))


df_topbot_antitabuchi_overlap <- df_topbot_antitabuchi %>%
  dplyr::select(Gene, Publication) %>%
  dplyr::mutate(Presence = 1) %>%
  tidyr::pivot_wider(Gene, names_from = Publication, values_from = Presence) %>%
  tibble::column_to_rownames("Gene") %>%
  replace(is.na(.), 0) %>%
  dplyr::mutate(sum = rowSums(.)) %>%
  rownames_to_column(var = "Gene")

df_topbot_tabuchi_overlap <- df_topbot_tabuchi %>%
  dplyr::select(Gene, Publication) %>%
  dplyr::mutate(Presence = 1) %>%
  tidyr::pivot_wider(Gene, names_from = Publication, values_from = Presence) %>%
  tibble::column_to_rownames("Gene") %>%
  replace(is.na(.), 0) %>%
  dplyr::mutate(sum = rowSums(.)) %>%
  rownames_to_column(var = "Gene")

# VENN of both clusters at least mentioned once

df_venn <- list(df_topbot_antitabuchi_overlap = df_topbot_antitabuchi_overlap$Gene,
                df_topbot_tabuchi_overlap = df_topbot_tabuchi_overlap$Gene)
partitions <- get.venn.partitions(df_venn)

temp <- venn.diagram(x = df_venn,
                     category.names = c("antitabuchi" ,
                                        "tabuchi"),
                     filename = NULL,
                     height = 5000,
                     width = 5000,
                     resolution =
                       2000,
                     units = "px",
                     compression = "lzw",
                     na = "stop",
                     sub = NULL,
                     col = c("orange", 'lightblue'),
                     fill = c(alpha("orange",0.3), alpha('lightblue',0.3)),
                     main = "",
                     main.pos = c(0.5, 1.05),
                     main.fontface = "plain",
                     main.fontfamily = "serif",
                     main.col = "black",
                     main.cex = 3,
                     main.just = c(0.5, 1),
                     cex = 3,
                     cat.cex = 3,
                     sub.fontface = "plain",
                     sub.fontfamily = "serif",
                     sub.col = "black",
                     sub.cex = 3,
                     sub.just = c(0.5, 1),
                     force.unique = TRUE,
                     print.mode = "raw",
                     sigdigs = 3,
                     direct.area = FALSE,
                     area.vector = 0,
                     hyper.test = FALSE,
                     total.population = NULL,
                     lower.tail = TRUE)
pdf(file="overlapping_top_genes_antitabuchi_vs_tabuchi_atleastONCE_VENN.pdf", width = 10, height = 10)
grid.draw(temp)
dev.off()

# VENN of both clusters at least mentioned wive

df_topbot_antitabuchi_overlap_filter <- df_topbot_antitabuchi_overlap %>% filter(sum >1)
df_topbot_tabuchi_overlap_filter <- df_topbot_tabuchi_overlap %>% filter(sum >2)

df_venn <- list(df_topbot_antitabuchi_overlap = df_topbot_antitabuchi_overlap_filter$Gene,
                df_topbot_tabuchi_overlap = df_topbot_tabuchi_overlap_filter$Gene)
partitions <- get.venn.partitions(df_venn)

temp <- venn.diagram(x = df_venn,
                     category.names = c("antitabuchi" ,
                                        "tabuchi"),
                     filename = NULL,
                     height = 5000,
                     width = 5000,
                     resolution =
                       2000,
                     units = "px",
                     compression = "lzw",
                     na = "stop",
                     sub = NULL,
                     col = c("orange", 'lightblue'),
                     fill = c(alpha("orange",0.3), alpha('lightblue',0.3)),
                     main = "",
                     main.pos = c(0.5, 1.05),
                     main.fontface = "plain",
                     main.fontfamily = "serif",
                     main.col = "black",
                     main.cex = 3,
                     main.just = c(0.5, 1),
                     cex = 3,
                     cat.cex = 3,
                     sub.fontface = "plain",
                     sub.fontfamily = "serif",
                     sub.col = "black",
                     sub.cex = 3,
                     sub.just = c(0.5, 1),
                     force.unique = TRUE,
                     print.mode = "raw",
                     sigdigs = 3,
                     direct.area = FALSE,
                     area.vector = 0,
                     hyper.test = FALSE,
                     total.population = NULL,
                     lower.tail = TRUE)
pdf(file="overlapping_top_genes_antitabuchi_vs_tabuchi_filter_VENN.pdf", width = 10, height = 10)
grid.draw(temp)
dev.off()

#### PCA - OPTIMAL DATASET ####

df_all_wide <- df_all %>%
  dplyr::select(-temperature_C, -duration_m, -timepoint_h, -CEM43) %>%
  pivot_wider(names_from = Publication, values_from = Ratio)

df_PCA <- df_all_wide %>%
  dplyr::select(-Gene) %>%
  na.omit()

df_PCA_input <- t(df_PCA)

samples.pca <- prcomp(df_PCA_input, center = TRUE, scale. = TRUE)

# plot eigenvalues

summary_PCA <- summary(samples.pca)
summary_PCA <- summary_PCA[["importance"]]
summary_PCA <- as.data.frame(summary_PCA)
summary_PCA <- summary_PCA %>% rownames_to_column(var = "factor")

summary_PCA_long <- summary_PCA %>%
  pivot_longer(!factor, names_to = "component", values_to = "value") %>%
  filter(factor %in% c("Cumulative Proportion")) %>%
  arrange(value)

summary_PCA_long$component <- factor(summary_PCA_long$component, levels = c("PC1", "PC2", "PC3", "PC4",
                                                                            "PC5", "PC6", "PC7", "PC8",
                                                                            "PC9", "PC10", "PC11", "PC12", "PC13",
                                                                            "PC14", "PC15", "PC16"))
summary_PCA_long_filter <- summary_PCA_long %>% filter(component %in% c("PC1", "PC2", "PC3", "PC4",
                                                                        "PC5", "PC6", "PC7", "PC8",
                                                                        "PC9", "PC10"))

write.csv(summary_PCA_long_filter, "eigenvalues.csv")

textsize <- 12
theme <-   theme(plot.title = element_text(size = textsize, face = "bold", hjust = 0.5),
                 plot.subtitle = element_text(size = textsize, hjust = 0.5),
                 panel.background = element_rect(fill="white"),
                 panel.grid.major = element_line(NA),
                 panel.grid.minor = element_line(NA),
                 axis.title.x = element_text(colour = "black", size=textsize),
                 axis.title.y = element_text(colour = "black", size=textsize),
                 axis.line.x = element_line(colour = 'black', size = 0.5),
                 axis.line.y = element_line(colour = 'black', size = 0.5),
                 axis.ticks.x = element_line(colour = "black", size = 0.5),
                 axis.ticks.y = element_line(colour = "black", size = 0.5),
                 axis.text.x = element_text(colour = "black", size=textsize),
                 axis.text.y = element_text(colour = "black", size=textsize),
                 legend.justification = c("right", "top"),
                 legend.background = element_blank(),
                 legend.title = element_text(face = "bold", size = textsize),
                 legend.text = element_text(size = textsize, colour = "black"),
                 legend.key = element_blank(),
                 aspect.ratio = 1)

ggplot(summary_PCA_long_filter, aes(x = component, y = value)) +
  geom_bar(stat = "identity") + theme +
  labs(title = "Eigenvalues", y = "Cumulative proportion of variance")

ggsave("eigenvalues.pdf", dpi = 1000)

# generate PCA plots

df_all_metrics <- df_all %>%
  dplyr::select(Publication, temperature_C, duration_m, timepoint_h, CEM43) %>%
  distinct()

publications <- df_all_metrics %>%
  dplyr::select(Publication)
publications <- publications$Publication

temp <- df_all_metrics %>%
  dplyr::select(temperature_C)
temp <- temp$temperature_C

time <- df_all_metrics %>%
  dplyr::select(timepoint_h)
time <- time$timepoint_h

duration <- df_all_metrics %>%
  dplyr::select(duration_m)
duration <- duration$duration_m

CEM43 <- df_all_metrics %>%
  dplyr::select(CEM43)
CEM43 <- CEM43$CEM43

time <- factor(time, levels = c("0", "1", "3", "4", "6", "24"))
CEM43 <- as.numeric(CEM43)
CEM43 <- factor(CEM43, levels = c("1.875", "3.75", "7.5", "15", "30", "120", "180", "480"))

textsize <- 12
theme <-   theme(plot.title = element_text(size = textsize, face = "bold", hjust = 0.5),
                 plot.subtitle = element_text(size = textsize, hjust = 0.5),
                 panel.background = element_rect(fill="white"),
                 panel.grid.major = element_line(NA),
                 panel.grid.minor = element_line(NA),
                 axis.title.x = element_text(colour = "black", size=textsize),
                 axis.title.y = element_text(colour = "black", size=textsize),
                 axis.line.x = element_line(colour = 'black', size = 0.5),
                 axis.line.y = element_line(colour = 'black', size = 0.5),
                 axis.ticks.x = element_line(colour = "black", size = 0.5),
                 axis.ticks.y = element_line(colour = "black", size = 0.5),
                 axis.text.x = element_text(colour = "black", size=textsize),
                 axis.text.y = element_text(colour = "black", size=textsize),
                 legend.justification = c("right", "top"),
                 legend.background = element_blank(),
                 legend.title = element_text(face = "bold", size = textsize),
                 legend.text = element_text(size = textsize, colour = "black"),
                 legend.key = element_blank(),
                 aspect.ratio = 1)
scale <- scale_color_brewer(palette = "Set3")

p1 <- ggplot(samples.pca, aes(x=PC1, y=PC2)) +
  geom_point() + geom_text_repel(label=publications,
                                                  color = "black",
                                                  max.overlaps = 27,
                                                  size = 4) + theme + scale
p2 <- ggplot(samples.pca, aes(x=PC2, y=PC3)) +
  geom_point() + geom_text_repel(label=publications,
                                                  color = "black",
                                                  max.overlaps = 26,
                                                  size = 4) + theme + scale
cowplot::plot_grid(p1,p2,ncol=2, align = "vh")

ggsave("PCA_labels.pdf", dpi=1000)

p3 <- ggplot(samples.pca, aes(x=PC1, y=PC2, color = time)) +
  geom_point(aes(size = CEM43)) + theme + scale

p4 <- ggplot(samples.pca, aes(x=PC2, y=PC3, color = time)) +
  geom_point(aes(size = CEM43)) + theme + scale

cowplot::plot_grid(p3,p4,ncol=2, align = "vh")

ggsave("PCA_param.pdf", dpi=1000)

testp1 <- layer_data(p1) %>%
  mutate(publication = publications, CEM = CEM43, time = time, duration = duration)

write.csv(testp1, "PCA1vs2.csv")

testp2 <- layer_data(p2) %>%
  mutate(publication = publications, CEM = CEM43, time = time, duration = duration)

write.csv(testp2, "PCA2vs3.csv")

#### HEATMAP - OVER/UNDEREPXRESSED GENES - ALL DATASETS ####

df_all_publication_filter <- df_all

df_genes <- df_all_publication_filter %>%
  dplyr::filter(Ratio >= 1.5 | Ratio <= 1/1.5) %>%
  dplyr::filter(!Ratio == 0)
df_genes <- unique(df_genes$Gene)

df_all_wide_test <- df_all_publication_filter %>%
  dplyr::filter(Gene %in% df_genes) %>%
  dplyr::select(-temperature_C, -duration_m, -timepoint_h, -CEM43) %>%
  tidyr::pivot_wider(names_from = Publication, values_from = Ratio) %>%
  na.omit()

df_all_wide_heatmap <- column_to_rownames(df_all_wide_test, var = "Gene") %>% na.omit()

df_all_metrics <- df_all %>%
  dplyr::select(Publication, temperature_C, duration_m, timepoint_h, CEM43) %>%
  distinct()

publications <- df_all_metrics %>%
  dplyr::select(Publication)
publications <- publications$Publication

temp <- df_all_metrics %>%
  dplyr::select(temperature_C)
temp <- temp$temperature_C

time <- df_all_metrics %>%
  dplyr::select(timepoint_h)
time <- time$timepoint_h

duration <- df_all_metrics %>%
  dplyr::select(duration_m)
duration <- duration$duration_m

CEM43 <- df_all_metrics %>%
  dplyr::select(CEM43)
CEM43 <- CEM43$CEM43

cat_df = data.frame("time" = time,
                    "CEM43" = CEM43)

rownames(cat_df) = colnames(df_all_wide_heatmap)

ann_colors = list(
  temp = c("41" = "#EEDBBD", "42" = "#E8B878", "43" = "#CF8552", "44" = "#C16840", "45" = "#9F3632"),
  time = c("0" = "#FEFFD9", "1" = "#EDF8BC", "3" = "#BFE8B2", "4" = "#6BC6BE", "6" = "#63C3BF", "24" = "#41B7C4"),
  CEM43 = c("1.875" = "#FFC685", "3.75" = "#F8B15B", "7.5" = "#F5A040","15" = "#b56f60","30" = "#ED7222", "120" = "#E26320", "180" = "#BA4C23", "480" = "#9E3D22"))

heatmap <- pheatmap(df_all_wide_heatmap,
         scale = "row",
         annotation_col = cat_df,
         annotation_color = ann_colors,
         cutree_cols = 3,
         cellwidth = 20,
         show_rownames = F,
         annotation_names_col = F,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100))
dev.off()

# export for figlinq

rownames <- rownames(df_all_wide_heatmap[heatmap$tree_row[["order"]],])

df_all_wide_heatmap.sorted <- df_all_wide_heatmap[match(rownames, rownames(df_all_wide_heatmap)),]

scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}

df_heatmap_num.sorted.scaled <- scale_rows(df_all_wide_heatmap.sorted)

write.csv(df_heatmap_num.sorted.scaled, "heatmap_DEG_all_studies_figlinq.csv")





pheatmap(df_all_wide_heatmap,
         scale = "row",
         annotation_col = cat_df,
         annotation_color = ann_colors,
         cutree_cols = 3,
         cutree_rows = 2,
         cellwidth = 20,
         cellheight = 0.5,
         show_rownames = F,
         annotation_names_col = F,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
         filename = "heatmap_overexp_underexp.pdf")
dev.off()



#### UPSET PLOT - OVER/UNDEREXPRESSED GENES - ALL ####

df_diff <- df_all %>%
  dplyr::filter(Ratio >= 1.5| Ratio <= 1/1.5) %>%
  dplyr::select(Gene, Publication, Ratio) %>%
  tidyr::pivot_wider(names_from = "Publication", values_from = "Ratio") %>%
  replace(is.na(.), 0)

df_test <- df_diff[2:ncol(df_diff)]
df_test[!df_test == 0] <- 1

df_diff_binary <- cbind(df_diff[1], df_test)

df_binary_all <- df_diff_binary %>%
  tibble::column_to_rownames(., var = "Gene") %>%
  dplyr::mutate(sum = rowSums(.)) %>%
  rownames_to_column(var = "Gene")

write.csv(df_binary_all, "gene_count_over_under_1point5_all.csv")

p <- ggplot(df_binary_all, aes(x=sum)) +
  geom_bar(stat="count")

layer <- layer_data(p, 1)

write.csv(layer, "gene_count_over_under_1point5_all_plot.csv")



#### UPSET PLOT - OVER/UNDEREXPRESSED GENES - ANTITABUCHI ####

df_diff <- df_all %>%
  dplyr::filter(Ratio >= 1.5| Ratio <= 1/1.5) %>%
  dplyr::select(Gene, Publication, Ratio) %>%
  tidyr::pivot_wider(names_from = "Publication", values_from = "Ratio") %>%
  replace(is.na(.), 0)

colnames(df_diff)

df_test <- df_diff[2:ncol(df_diff)]
df_test[!df_test == 0] <- 1

df_diff_binary <- cbind(df_diff[1], df_test)

df_binary_antitabuchi <- df_diff_binary %>%
  tibble::column_to_rownames(., var = "Gene") %>%
  dplyr::select(Amaya2014a,
                Amaya2014b,
                Amaya2014c,
                Court2017,
                Andocs2015a,
                Yunoki2016,
                Scutigliani2022a,
                Scutigliani2022b,
                Tabuchi2008) %>%
  dplyr::mutate(sum = rowSums(.)) %>%
  rownames_to_column(var = "Gene")

write.csv(df_binary_antitabuchi, "gene_count_over_under_1point5_antitabuchi.csv")

p <- ggplot(df_binary_antitabuchi, aes(x=sum)) +
  geom_bar(stat="count")

layer <- layer_data(p, 1)

write.csv(layer, "gene_count_over_under_1point5_cluster1_plot.csv")

#### UPSET PLOT - OVER/UNDEREXPRESSED GENES - TABUCHI ####

df_diff <- df_all %>%
  dplyr::filter(Ratio >= 1.5| Ratio <= 1/1.5) %>%
  dplyr::select(Gene, Publication, Ratio) %>%
  tidyr::pivot_wider(names_from = "Publication", values_from = "Ratio") %>%
  replace(is.na(.), 0)

df_test <- df_diff[2:ncol(df_diff)]
df_test[!df_test == 0] <- 1

df_diff_binary <- cbind(df_diff[1], df_test)

df_binary_tabuchi <- df_diff_binary %>%
  tibble::column_to_rownames(., var = "Gene") %>%
  dplyr::select(-Amaya2014a,
                -Amaya2014b,
                -Amaya2014c,
                -Court2017,
                -Andocs2015a,
                -Yunoki2016,
                -Scutigliani2022a,
                -Scutigliani2022b,
                -Tabuchi2008) %>%
  dplyr::mutate(sum = rowSums(.)) %>%
  rownames_to_column(var = "Gene")

write.csv(df_binary_tabuchi, "gene_count_over_under_1point5_tabuchi.csv")

p <- ggplot(df_binary_tabuchi, aes(x=sum)) +
  geom_bar(stat="count")

layer <- layer_data(p, 1)

write.csv(layer, "gene_count_over_under_1point5_cluster2_plot.csv")

#### UPSET PLOT - OVER/UNDEREXPRESSED GENES - OVERLAP BETWEEN CLUSTERS ####

df_intersections_antitabuchi <- read.csv("gene_count_over_under_1point5_antitabuchi.csv")
df_intersections_tabuchi <- read.csv("gene_count_over_under_1point5_tabuchi.csv")

df_intersect_filter_all <- df_intersections %>%
  dplyr::filter(sum > 1) # intersecting genes between datasets - ALL

df_intersect_filter_antitabuchi <- df_intersections_antitabuchi %>%
  dplyr::filter(sum > 3) # intersecting genes between datasets - antitabuchi

df_intersect_filter_tabuchi <- df_intersections_tabuchi %>%
  dplyr::filter(sum > 7) # intersecting genes between datasets - tabuchi

df_venn <- list(df_intersect_filter_antitabuchi = df_intersect_filter_antitabuchi$Gene,
                df_intersect_filter_tabuchi = df_intersect_filter_tabuchi$Gene)
partitions <- get.venn.partitions(df_venn)



#### UPSET PLOTS - ORA OF CLUSTER OVERLAP ####

## ORA ##

# read upset plot dataframes and perform inclusion

df_intersections_antitabuchi <- read.csv("gene_count_over_under_1point5_antitabuchi.csv") %>%
  dplyr::filter(sum > 3)

df_intersections_tabuchi <- read.csv("gene_count_over_under_1point5_tabuchi.csv") %>%
  dplyr::filter(sum > 7)

df_overlap <- list(df_intersections_antitabuchi = df_intersections_antitabuchi$Gene,
                   df_intersections_tabuchi = df_intersections_tabuchi$Gene)
partitions <- get.venn.partitions(df_overlap)
df_overlap <- unlist(partitions$..values..$`1`)

## ORA - H ##

# select background

m_t2g <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>% #selecting a species and category of terms
  dplyr::select(gs_name, gene_symbol) #filtering gene name per process as input for GSEA

# ORA

df_enrich_antitabuchi <- clusterProfiler::enricher(gene = df_intersections_antitabuchi$Gene, TERM2GENE = m_t2g)
df_enrich_tabuchi <- clusterProfiler::enricher(gene = df_intersections_tabuchi$Gene, TERM2GENE = m_t2g)
df_enrich_overlap <- clusterProfiler::enricher(gene = df_overlap, TERM2GENE = m_t2g)

df_enrich_antitabuchi <- df_enrich_antitabuchi@result
df_enrich_tabuchi <- df_enrich_tabuchi@result
df_enrich_overlap <- df_enrich_overlap@result

# export

write.csv(df_enrich_antitabuchi, "ORA_SumOver3_cluster1.csv")
write.csv(df_enrich_tabuchi, "ORA_SumOver7_cluster2.csv")
write.csv(df_enrich_overlap, "ORA_SumOver2And7_overlap.csv")

## ORA - GOBP ##

# select background

m_t2g <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>% #selecting a species and category of terms
  dplyr::select(gs_name, gene_symbol) #filtering gene name per process as input for GSEA

# ORA

df_enrich_antitabuchi <- clusterProfiler::enricher(gene = df_intersections_antitabuchi$Gene, TERM2GENE = m_t2g)
df_enrich_tabuchi <- clusterProfiler::enricher(gene = df_intersections_tabuchi$Gene, TERM2GENE = m_t2g)
df_enrich_overlap <- clusterProfiler::enricher(gene = df_overlap, TERM2GENE = m_t2g)

df_enrich_antitabuchi <- df_enrich_antitabuchi@result
df_enrich_tabuchi <- df_enrich_tabuchi@result
df_enrich_overlap <- df_enrich_overlap@result

# export

write.csv(df_enrich_antitabuchi, "ORA_SumOver3_cluster1.csv")
write.csv(df_enrich_tabuchi, "ORA_SumOver7_cluster2.csv")
write.csv(df_enrich_overlap, "ORA_GOBP_SumOver2And7_overlap.csv")

## heatmap visualization ##

# select genes of interest 

genes_of_interest <- df_enrich_overlap %>% 
  dplyr::select(geneID) %>% 
  tidyr::separate(geneID, sep = "/", paste("target", 1:10, sep = "_")) %>%
  tidyr::pivot_longer(paste("target", 1:10, sep = "_"), names_to = "Target_n", values_to = "Target") %>%
  dplyr::filter(!is.na(Target)) %>% 
  dplyr::select(Target) %>%
  unique()

df_exp <- df_all %>% 
  dplyr::filter(Gene %in% genes_of_interest$Target) %>% 
  dplyr::mutate(log10FC = log10(Ratio))

df_heatmap <- df_exp %>% 
  dplyr::select(Gene, log10FC, Publication) %>% 
  tidyr::pivot_wider(names_from = Gene, values_from = log10FC) 
df_heatmap[is.na(df_heatmap)] <- 0

df_heatmap_num <- df_heatmap[2:ncol(df_heatmap)]
df_heatmap_num <- as.data.frame(df_heatmap_num)
rownames(df_heatmap_num) = df_heatmap$Publication

df_heatmap_input <- t(as.matrix(df_heatmap_num))

df_all_metrics <- df_all %>%
  dplyr::select(Publication, temperature_C, duration_m, timepoint_h, CEM43) %>%
  distinct()

publications <- df_all_metrics %>%
  dplyr::filter(Publication %in% colnames(df_heatmap_input)) %>%
  dplyr::select(Publication)
publications <- publications$Publication

temp <- df_all_metrics %>%
  dplyr::filter(Publication %in% colnames(df_heatmap_input)) %>%
  dplyr::select(temperature_C)
temp <- temp$temperature_C

time <- df_all_metrics %>%
  dplyr::filter(Publication %in% colnames(df_heatmap_input)) %>%
  dplyr::select(timepoint_h)
time <- time$timepoint_h

CEM43 <- df_all_metrics %>%
  dplyr::filter(Publication %in% colnames(df_heatmap_input)) %>%
  dplyr::select(CEM43)
CEM43 <- CEM43$CEM43

cat_df = data.frame("time" = time,
                    "CEM43" = CEM43)

rownames(cat_df) = colnames(df_heatmap_input)

ann_colors = list(
  time = c("0" = "#FEFFD9", "1" = "#EDF8BC", "3" = "#BFE8B2", "4" = "#6BC6BE", "6" = "#63C3BF", "24" = "#41B7C4"),
  CEM43 = c("1.875" = "#FFC685", "3.75" = "#F8B15B", "7.5" = "#F5A040","15" = "#b56f60","30" = "#ED7222", "120" = "#E26320", "180" = "#BA4C23", "480" = "#9E3D22"))

heatmap <- pheatmap(df_heatmap_input,
                    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                    annotation_col = cat_df,
                    annotation_color = ann_colors,
                    annotation_names_col = T,
                    annotation_names_row = T,
                    legend = T,
                    scale = "none")
dev.off()

# export for figlinq

rownames <- rownames(df_heatmap_input[heatmap$tree_row[["order"]],])
colnames <- colnames(df_heatmap_input[,heatmap$tree_col[["order"]]])

df_heatmap_input.sorted <- df_heatmap_input[match(rownames, rownames(df_heatmap_input)),]
df_heatmap_input.sorted <- df_heatmap_input.sorted[,match(colnames, colnames(df_heatmap_input.sorted))]

write.csv(df_heatmap_input.sorted, "heatmap_gene_overlap_after_ORA_cluster1a2_figlinq.csv")

#### GSEA - DATA PREPARATION - ALL DATASETS ####

# start df

df_all_GSEA <- df_all %>%
  dplyr::select(-temperature_C, -duration_m, -timepoint_h) %>%
  tidyr::pivot_wider(names_from = Publication, values_from = Ratio)

# create input gene list

Tabuchi2008 <- df_all_GSEA %>% dplyr::select(Tabuchi2008, Gene) %>% na.omit()

Tabuchi2008_ranked <- Tabuchi2008 %>%
  mutate(rank = rank(-Tabuchi2008,  ties.method = "random")) %>%
  arrange(rank)

geneList.Tabuchi2008 <- Tabuchi2008_ranked$rank
names(geneList.Tabuchi2008) <- Tabuchi2008_ranked$Gene #result is a named vector
geneList.Tabuchi2008 <- sort(geneList.Tabuchi2008, decreasing = TRUE) #for ranking the vector, obligatory for input



Tabuchi2011a <- df_all_GSEA %>% dplyr::select(Tabuchi2011a, Gene) %>% na.omit()

Tabuchi2011a_ranked <- Tabuchi2011a %>%
  mutate(rank = rank(-Tabuchi2011a,  ties.method = "random")) %>%
  arrange(rank)

geneList.Tabuchi2011a <- Tabuchi2011a_ranked$rank
names(geneList.Tabuchi2011a) <- Tabuchi2011a_ranked$Gene #result is a named vector
geneList.Tabuchi2011a <- sort(geneList.Tabuchi2011a, decreasing = TRUE) #for ranking the vector, obligatory for input



Tabuchi2011b <- df_all_GSEA %>% dplyr::select(Tabuchi2011b, Gene) %>% na.omit()

Tabuchi2011b_ranked <- Tabuchi2011b %>%
  mutate(rank = rank(-Tabuchi2011b,  ties.method = "random")) %>%
  arrange(rank)

geneList.Tabuchi2011b <- Tabuchi2011b_ranked$rank
names(geneList.Tabuchi2011b) <- Tabuchi2011b_ranked$Gene #result is a named vector
geneList.Tabuchi2011b <- sort(geneList.Tabuchi2011b, decreasing = TRUE) #for ranking the vector, obligatory for input



Tabuchi2011c <- df_all_GSEA %>% dplyr::select(Tabuchi2011c, Gene) %>% na.omit()

Tabuchi2011c_ranked <- Tabuchi2011c %>%
  mutate(rank = rank(-Tabuchi2011c,  ties.method = "random")) %>%
  arrange(rank)

geneList.Tabuchi2011c <- Tabuchi2011c_ranked$rank
names(geneList.Tabuchi2011c) <- Tabuchi2011c_ranked$Gene #result is a named vector
geneList.Tabuchi2011c <- sort(geneList.Tabuchi2011c, decreasing = TRUE) #for ranking the vector, obligatory for input


Tabuchi2011d <- df_all_GSEA %>% dplyr::select(Tabuchi2011d, Gene) %>% na.omit()

Tabuchi2011d_ranked <- Tabuchi2011d %>%
  mutate(rank = rank(-Tabuchi2011d,  ties.method = "random")) %>%
  arrange(rank)

geneList.Tabuchi2011d <- Tabuchi2011d_ranked$rank
names(geneList.Tabuchi2011d) <- Tabuchi2011d_ranked$Gene #result is a named vector
geneList.Tabuchi2011d <- sort(geneList.Tabuchi2011d, decreasing = TRUE) #for ranking the vector, obligatory for input



Tabuchi2011e <- df_all_GSEA %>% dplyr::select(Tabuchi2011e, Gene) %>% na.omit()

Tabuchi2011e_ranked <- Tabuchi2011e %>%
  mutate(rank = rank(-Tabuchi2011e,  ties.method = "random")) %>%
  arrange(rank)

geneList.Tabuchi2011e <- Tabuchi2011e_ranked$rank
names(geneList.Tabuchi2011e) <- Tabuchi2011e_ranked$Gene #result is a named vector
geneList.Tabuchi2011e <- sort(geneList.Tabuchi2011e, decreasing = TRUE) #for ranking the vector, obligatory for input



Tabuchi2011f <- df_all_GSEA %>% dplyr::select(Tabuchi2011f, Gene) %>% na.omit()

Tabuchi2011f_ranked <- Tabuchi2011f %>%
  mutate(rank = rank(-Tabuchi2011f,  ties.method = "random")) %>%
  arrange(rank)

geneList.Tabuchi2011f <- Tabuchi2011f_ranked$rank
names(geneList.Tabuchi2011f) <- Tabuchi2011f_ranked$Gene #result is a named vector
geneList.Tabuchi2011f <- sort(geneList.Tabuchi2011f, decreasing = TRUE) #for ranking the vector, obligatory for input



Tabuchi2011g <- df_all_GSEA %>% dplyr::select(Tabuchi2011g, Gene) %>% na.omit()

Tabuchi2011g_ranked <- Tabuchi2011g %>%
  mutate(rank = rank(-Tabuchi2011g,  ties.method = "random")) %>%
  arrange(rank)

geneList.Tabuchi2011g <- Tabuchi2011g_ranked$rank
names(geneList.Tabuchi2011g) <- Tabuchi2011g_ranked$Gene #result is a named vector
geneList.Tabuchi2011g <- sort(geneList.Tabuchi2011g, decreasing = TRUE) #for ranking the vector, obligatory for input




Tabuchi2011h <- df_all_GSEA %>% dplyr::select(Tabuchi2011h, Gene) %>% na.omit()

Tabuchi2011h_ranked <- Tabuchi2011h %>%
  mutate(rank = rank(-Tabuchi2011h,  ties.method = "random")) %>%
  arrange(rank)

geneList.Tabuchi2011h <- Tabuchi2011h_ranked$rank
names(geneList.Tabuchi2011h) <- Tabuchi2011h_ranked$Gene #result is a named vector
geneList.Tabuchi2011h <- sort(geneList.Tabuchi2011h, decreasing = TRUE) #for ranking the vector, obligatory for input


Amaya2014a <- df_all_GSEA %>% dplyr::select(Amaya2014a, Gene) %>% na.omit()

Amaya2014a_ranked <- Amaya2014a %>%
  mutate(rank = rank(-Amaya2014a,  ties.method = "random")) %>%
  arrange(rank)

geneList.Amaya2014a <- Amaya2014a_ranked$rank
names(geneList.Amaya2014a) <- Amaya2014a_ranked$Gene #result is a named vector
geneList.Amaya2014a <- sort(geneList.Amaya2014a, decreasing = TRUE) #for ranking the vector, obligatory for input


Amaya2014b <- df_all_GSEA %>% dplyr::select(Amaya2014b, Gene) %>% na.omit()

Amaya2014b_ranked <- Amaya2014b %>%
  mutate(rank = rank(-Amaya2014b,  ties.method = "random")) %>%
  arrange(rank)

geneList.Amaya2014b <- Amaya2014b_ranked$rank
names(geneList.Amaya2014b) <- Amaya2014b_ranked$Gene #result is a named vector
geneList.Amaya2014b <- sort(geneList.Amaya2014b, decreasing = TRUE) #for ranking the vector, obligatory for input



Amaya2014c <- df_all_GSEA %>% dplyr::select(Amaya2014c, Gene) %>% na.omit()

Amaya2014c_ranked <- Amaya2014c %>%
  mutate(rank = rank(-Amaya2014c,  ties.method = "random")) %>%
  arrange(rank)

geneList.Amaya2014c <- Amaya2014c_ranked$rank
names(geneList.Amaya2014c) <- Amaya2014c_ranked$Gene #result is a named vector
geneList.Amaya2014c <- sort(geneList.Amaya2014c, decreasing = TRUE) #for ranking the vector, obligatory for input



Court2017 <- df_all_GSEA %>% dplyr::select(Court2017, Gene) %>% na.omit()

Court2017_ranked <- Court2017 %>%
  mutate(rank = rank(-Court2017,  ties.method = "random")) %>%
  arrange(rank)

geneList.Court2017 <- Court2017_ranked$rank
names(geneList.Court2017) <- Court2017_ranked$Gene #result is a named vector
geneList.Court2017 <- sort(geneList.Court2017, decreasing = TRUE) #for ranking the vector, obligatory for input



Andocs2015a <- df_all_GSEA %>% dplyr::select(Andocs2015a, Gene) %>% na.omit()

Andocs2015a_ranked <- Andocs2015a %>%
  mutate(rank = rank(-Andocs2015a,  ties.method = "random")) %>%
  arrange(rank)

geneList.Andocs2015a <- Andocs2015a_ranked$rank
names(geneList.Andocs2015a) <- Andocs2015a_ranked$Gene #result is a named vector
geneList.Andocs2015a <- sort(geneList.Andocs2015a, decreasing = TRUE) #for ranking the vector, obligatory for input



Andocs2015b <- df_all_GSEA %>% dplyr::select(Andocs2015b, Gene) %>% na.omit()

Andocs2015b_ranked <- Andocs2015b %>%
  mutate(rank = rank(-Andocs2015b,  ties.method = "random")) %>%
  arrange(rank)

geneList.Andocs2015b <- Andocs2015b_ranked$rank
names(geneList.Andocs2015b) <- Andocs2015b_ranked$Gene #result is a named vector
geneList.Andocs2015b <- sort(geneList.Andocs2015b, decreasing = TRUE) #for ranking the vector, obligatory for input


Yunoki2016 <- df_all_GSEA %>% dplyr::select(Yunoki2016, Gene) %>% na.omit()

Yunoki2016_ranked <- Yunoki2016 %>%
  mutate(rank = rank(-Yunoki2016,  ties.method = "random")) %>%
  arrange(rank)

geneList.Yunoki2016 <- Yunoki2016_ranked$rank
names(geneList.Yunoki2016) <- Yunoki2016_ranked$Gene #result is a named vector
geneList.Yunoki2016 <- sort(geneList.Yunoki2016, decreasing = TRUE) #for ranking the vector, obligatory for input

Scutigliani2022a <- df_all_GSEA %>% dplyr::select(Scutigliani2022a, Gene) %>% na.omit()

Scutigliani2022a_ranked <- Scutigliani2022a %>%
  mutate(rank = rank(-Scutigliani2022a,  ties.method = "random")) %>%
  arrange(rank)

geneList.Scutigliani2022a <- Scutigliani2022a_ranked$rank
names(geneList.Scutigliani2022a) <- Scutigliani2022a_ranked$Gene #result is a named vector
geneList.Scutigliani2022a <- sort(geneList.Scutigliani2022a, decreasing = TRUE) #for ranking the vector, obligatory for input

Scutigliani2022b <- df_all_GSEA %>% dplyr::select(Scutigliani2022b, Gene) %>% na.omit()

Scutigliani2022b_ranked <- Scutigliani2022b %>%
  mutate(rank = rank(-Scutigliani2022b,  ties.method = "random")) %>%
  arrange(rank)

geneList.Scutigliani2022b <- Scutigliani2022b_ranked$rank
names(geneList.Scutigliani2022b) <- Scutigliani2022b_ranked$Gene #result is a named vector
geneList.Scutigliani2022b <- sort(geneList.Scutigliani2022b, decreasing = TRUE) #for ranking the vector, obligatory for input



#### GSEA - HALLMARKS - ALL DATASETS ####

# load molecular signatures of interest using the msigdbr package

library("msigdbr")

m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% #selecting a species and category of terms
  dplyr::select(gs_name, gene_symbol) #filtering gene name per process as input for GSEA

# perform GSEA using the ClusterProfiler package

library("clusterProfiler")


gseaRes.Tabuchi2008 <- GSEA(
  geneList.Tabuchi2008,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = m_t2g,
  TERM2NAME = NA,
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea",
)

df_gseaRes.Tabuchi2008 <- gseaRes.Tabuchi2008@result

gseaRes.Tabuchi2011a <- GSEA(
  geneList.Tabuchi2011a,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = m_t2g,
  TERM2NAME = NA,
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea",
)

df_gseaRes.Tabuchi2011a <- gseaRes.Tabuchi2011a@result

gseaRes.Tabuchi2011b <- GSEA(
  geneList.Tabuchi2011b,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = m_t2g,
  TERM2NAME = NA,
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea",
)

df_gseaRes.Tabuchi2011b <- gseaRes.Tabuchi2011b@result

gseaRes.Tabuchi2011c <- GSEA(
  geneList.Tabuchi2011c,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = m_t2g,
  TERM2NAME = NA,
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea",
)

df_gseaRes.Tabuchi2011c <- gseaRes.Tabuchi2011c@result

gseaRes.Tabuchi2011d <- GSEA(
  geneList.Tabuchi2011d,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = m_t2g,
  TERM2NAME = NA,
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea",
)

df_gseaRes.Tabuchi2011d <- gseaRes.Tabuchi2011d@result

gseaRes.Tabuchi2011e <- GSEA(
  geneList.Tabuchi2011e,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = m_t2g,
  TERM2NAME = NA,
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea",
)

df_gseaRes.Tabuchi2011e <- gseaRes.Tabuchi2011e@result

gseaRes.Tabuchi2011f <- GSEA(
  geneList.Tabuchi2011f,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = m_t2g,
  TERM2NAME = NA,
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea",
)

df_gseaRes.Tabuchi2011f <- gseaRes.Tabuchi2011f@result

gseaRes.Tabuchi2011g <- GSEA(
  geneList.Tabuchi2011g,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = m_t2g,
  TERM2NAME = NA,
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea",
)

df_gseaRes.Tabuchi2011g <- gseaRes.Tabuchi2011g@result

gseaRes.Tabuchi2011h <- GSEA(
  geneList.Tabuchi2011h,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = m_t2g,
  TERM2NAME = NA,
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea",
)

df_gseaRes.Tabuchi2011h <- gseaRes.Tabuchi2011h@result

gseaRes.Amaya2014a <- GSEA(
  geneList.Amaya2014a,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = m_t2g,
  TERM2NAME = NA,
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea",
)

df_gseaRes.Amaya2014a <- gseaRes.Amaya2014a@result

gseaRes.Amaya2014b <- GSEA(
  geneList.Amaya2014b,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = m_t2g,
  TERM2NAME = NA,
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea",
)

df_gseaRes.Amaya2014b <- gseaRes.Amaya2014b@result

gseaRes.Amaya2014c <- GSEA(
  geneList.Amaya2014c,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = m_t2g,
  TERM2NAME = NA,
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea",
)

df_gseaRes.Amaya2014c <- gseaRes.Amaya2014c@result

gseaRes.Court2017 <- GSEA(
  geneList.Court2017,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = m_t2g,
  TERM2NAME = NA,
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea",
)

df_gseaRes.Court2017 <- gseaRes.Court2017@result

gseaRes.Andocs2015a <- GSEA(
  geneList.Andocs2015a,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = m_t2g,
  TERM2NAME = NA,
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea",
)

df_gseaRes.Andocs2015a <- gseaRes.Andocs2015a@result

gseaRes.Andocs2015b <- GSEA(
  geneList.Andocs2015b,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = m_t2g,
  TERM2NAME = NA,
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea",
)

df_gseaRes.Andocs2015b <- gseaRes.Andocs2015b@result

gseaRes.Yunoki2016 <- GSEA(
  geneList.Yunoki2016,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = m_t2g,
  TERM2NAME = NA,
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea",
)

df_gseaRes.Yunoki2016 <- gseaRes.Yunoki2016@result

gseaRes.Scutigliani2022a <- GSEA(
  geneList.Scutigliani2022a,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = m_t2g,
  TERM2NAME = NA,
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea",
)

df_gseaRes.Scutigliani2022a <- gseaRes.Scutigliani2022a@result

gseaRes.Scutigliani2022b <- GSEA(
  geneList.Scutigliani2022b,
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 1e-10,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = m_t2g,
  TERM2NAME = NA,
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea",
)

df_gseaRes.Scutigliani2022b <- gseaRes.Scutigliani2022b@result

# visualization in heatmap

# filter signifcantly enriched terms

df_ENRICH_Tabuchi2008 <- dplyr::filter(df_gseaRes.Tabuchi2008, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2008")
df_ENRICH_Tabuchi2011a <- dplyr::filter(df_gseaRes.Tabuchi2011a, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011a")
df_ENRICH_Tabuchi2011b <- dplyr::filter(df_gseaRes.Tabuchi2011b, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011b")
df_ENRICH_Tabuchi2011c <- dplyr::filter(df_gseaRes.Tabuchi2011c, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011c")
df_ENRICH_Tabuchi2011d <- dplyr::filter(df_gseaRes.Tabuchi2011d, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011d")
df_ENRICH_Tabuchi2011e <- dplyr::filter(df_gseaRes.Tabuchi2011e, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011e")
df_ENRICH_Tabuchi2011f <- dplyr::filter(df_gseaRes.Tabuchi2011f, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011f")
df_ENRICH_Tabuchi2011g <- dplyr::filter(df_gseaRes.Tabuchi2011g, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011g")
df_ENRICH_Tabuchi2011h <- dplyr::filter(df_gseaRes.Tabuchi2011h, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust)%>%
  dplyr::mutate(condition = "Tabuchi2011h")
df_ENRICH_Amaya2014a <- dplyr::filter(df_gseaRes.Amaya2014a, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Amaya2014a")
df_ENRICH_Amaya2014b <- dplyr::filter(df_gseaRes.Amaya2014b, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Amaya2014b")
df_ENRICH_Amaya2014c <- dplyr::filter(df_gseaRes.Amaya2014c, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Amaya2014c")
df_ENRICH_Court2017<- dplyr::filter(df_gseaRes.Court2017, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Court2017")
df_ENRICH_Andocs2015a <- dplyr::filter(df_gseaRes.Andocs2015a, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Andocs2015a")
df_ENRICH_Andocs2015b <- dplyr::filter(df_gseaRes.Andocs2015b, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Andocs2015b")
df_ENRICH_Yunoki2016 <- dplyr::filter(df_gseaRes.Yunoki2016, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Yunoki2016")

df_ENRICH_Scutigliani2022a <- dplyr::filter(df_gseaRes.Scutigliani2022a, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Scutigliani2022a")

df_ENRICH_Scutigliani2022b <- dplyr::filter(df_gseaRes.Scutigliani2022b, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Scutigliani2022b")

df_ENRICH_all <- rbind(df_ENRICH_Tabuchi2008,
                       df_ENRICH_Tabuchi2011a,
                       df_ENRICH_Tabuchi2011b,
                       df_ENRICH_Tabuchi2011c,
                       df_ENRICH_Tabuchi2011d,
                       df_ENRICH_Tabuchi2011e,
                       df_ENRICH_Tabuchi2011f,
                       df_ENRICH_Tabuchi2011g,
                       df_ENRICH_Tabuchi2011h,
                       df_ENRICH_Amaya2014a,
                       df_ENRICH_Amaya2014b,
                       df_ENRICH_Amaya2014c,
                       df_ENRICH_Court2017,
                       df_ENRICH_Andocs2015a,
                       df_ENRICH_Andocs2015b,
                       df_ENRICH_Yunoki2016,
                       df_ENRICH_Scutigliani2022a,
                       df_ENRICH_Scutigliani2022b)
list_all_terms <- unique(df_ENRICH_all$Description)

df_ENRICH_Tabuchi2008_sign <- dplyr::filter(df_gseaRes.Tabuchi2008, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2008")
df_ENRICH_Tabuchi2011a_sign <- dplyr::filter(df_gseaRes.Tabuchi2011a, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011a")
df_ENRICH_Tabuchi2011b_sign <- dplyr::filter(df_gseaRes.Tabuchi2011b, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011b")
df_ENRICH_Tabuchi2011c_sign <- dplyr::filter(df_gseaRes.Tabuchi2011c, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011c")
df_ENRICH_Tabuchi2011d_sign <- dplyr::filter(df_gseaRes.Tabuchi2011d, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011d")
df_ENRICH_Tabuchi2011e_sign <- dplyr::filter(df_gseaRes.Tabuchi2011e, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011e")
df_ENRICH_Tabuchi2011f_sign <- dplyr::filter(df_gseaRes.Tabuchi2011f, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011f")
df_ENRICH_Tabuchi2011g_sign <- dplyr::filter(df_gseaRes.Tabuchi2011g, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011g")
df_ENRICH_Tabuchi2011h_sign <- dplyr::filter(df_gseaRes.Tabuchi2011h, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust)%>%
  dplyr::mutate(condition = "Tabuchi2011h")
df_ENRICH_Amaya2014a_sign <- dplyr::filter(df_gseaRes.Amaya2014a, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Amaya2014a")
df_ENRICH_Amaya2014b_sign <- dplyr::filter(df_gseaRes.Amaya2014b, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Amaya2014b")
df_ENRICH_Amaya2014c_sign <- dplyr::filter(df_gseaRes.Amaya2014c, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Amaya2014c")
df_ENRICH_Court2017_sign <- dplyr::filter(df_gseaRes.Court2017, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Court2017")
df_ENRICH_Andocs2015a_sign <- dplyr::filter(df_gseaRes.Andocs2015a, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Andocs2015a")
df_ENRICH_Andocs2015b_sign <- dplyr::filter(df_gseaRes.Andocs2015b, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Andocs2015b")
df_ENRICH_Yunoki2016_sign <- dplyr::filter(df_gseaRes.Yunoki2016, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Yunoki2016")
df_ENRICH_Scutigliani2022a_sign <- dplyr::filter(df_gseaRes.Scutigliani2022a, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Scutigliani2022a")
df_ENRICH_Scutigliani2022b_sign <- dplyr::filter(df_gseaRes.Scutigliani2022b, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Scutigliani2022b")

df_ENRICH_sign <- rbind(df_ENRICH_Tabuchi2008_sign,
                        df_ENRICH_Tabuchi2011a_sign,
                        df_ENRICH_Tabuchi2011b_sign,
                        df_ENRICH_Tabuchi2011c_sign,
                        df_ENRICH_Tabuchi2011d_sign,
                        df_ENRICH_Tabuchi2011e_sign,
                        df_ENRICH_Tabuchi2011f_sign,
                        df_ENRICH_Tabuchi2011g_sign,
                        df_ENRICH_Tabuchi2011h_sign,
                        df_ENRICH_Amaya2014a_sign,
                        df_ENRICH_Amaya2014b_sign,
                        df_ENRICH_Amaya2014c_sign,
                        df_ENRICH_Court2017_sign,
                        df_ENRICH_Andocs2015a_sign,
                        df_ENRICH_Andocs2015b_sign,
                        df_ENRICH_Yunoki2016_sign,
                        df_ENRICH_Scutigliani2022a_sign,
                        df_ENRICH_Scutigliani2022b_sign)

for (i in 1:nrow(df_ENRICH_sign))
  if (df_ENRICH_sign$p.adjust[i] < 0.05) df_ENRICH_sign$sign[i] = 1 else df_ENRICH_sign$sign[i] = 0

df_ENRICH_binary <- dplyr::select(df_ENRICH_sign, -p.adjust)
df_ENRICH_shaped <- pivot_wider(df_ENRICH_binary, names_from = condition, values_from = sign) %>%
  replace(is.na(.), 0)

df_heatmap_num <- df_ENRICH_shaped[2:ncol(df_ENRICH_shaped)]
rownames(df_heatmap_num) = df_ENRICH_shaped$Description
df_heatmap_input <- as.matrix(df_heatmap_num)

df_all_metrics <- df_all %>%
  dplyr::select(Publication, temperature_C, duration_m, timepoint_h, CEM43) %>%
  distinct()

publications <- df_all_metrics %>%
  dplyr::filter(Publication %in% colnames(df_heatmap_input)) %>%
  dplyr::select(Publication)
publications <- publications$Publication

temp <- df_all_metrics %>%
  dplyr::filter(Publication %in% colnames(df_heatmap_input)) %>%
  dplyr::select(temperature_C)
temp <- temp$temperature_C

time <- df_all_metrics %>%
  dplyr::filter(Publication %in% colnames(df_heatmap_input)) %>%
  dplyr::select(timepoint_h)
time <- time$timepoint_h

CEM43 <- df_all_metrics %>%
  dplyr::filter(Publication %in% colnames(df_heatmap_input)) %>%
  dplyr::select(CEM43)
CEM43 <- CEM43$CEM43

cat_df = data.frame("time" = time,
                    "CEM43" = CEM43)

rownames(cat_df) = colnames(df_heatmap_input)

ann_colors = list(
  time = c("0" = "#FEFFD9", "1" = "#EDF8BC", "3" = "#BFE8B2", "4" = "#6BC6BE", "6" = "#63C3BF", "24" = "#41B7C4"),
  CEM43 = c("1.875" = "#FFC685", "3.75" = "#F8B15B", "7.5" = "#F5A040","15" = "#b56f60","30" = "#ED7222", "120" = "#E26320", "180" = "#BA4C23", "480" = "#9E3D22"))

heatmap <- pheatmap(df_heatmap_input,
         color = c("white", "lightblue"),
         annotation_col = cat_df,
         annotation_color = ann_colors,
         cutree_cols = 2,
         cutree_rows = 2,
         cellheight = 10,
         cellwidth = 20,
         annotation_names_col = F,
         legend = F,
         file = "heatmap_hallmarks_GSEA.pdf")
dev.off()

rownames <- rownames(df_heatmap_input[heatmap$tree_row[["order"]],])
colnames <- colnames(df_heatmap_input[,heatmap$tree_col[["order"]]])

df_heatmap_num.sorted <- df_heatmap_input[match(rownames, rownames(df_heatmap_input)),]
df_heatmap_num.sorted <- df_heatmap_num.sorted[,match(colnames, colnames(df_heatmap_num.sorted))]

write.csv(df_heatmap_num.sorted, "heatmap_GSEA_H.csv")

#### GSEA - HALLMARKS - barplot NES/P ####

# filter signifcantly enriched terms

df_ENRICH_Tabuchi2008 <- dplyr::filter(df_gseaRes.Tabuchi2008, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust, NES) %>%
  dplyr::mutate(condition = "Tabuchi2008")
df_ENRICH_Tabuchi2011a <- dplyr::filter(df_gseaRes.Tabuchi2011a, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust, NES) %>%
  dplyr::mutate(condition = "Tabuchi2011a")
df_ENRICH_Tabuchi2011b <- dplyr::filter(df_gseaRes.Tabuchi2011b, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust, NES) %>%
  dplyr::mutate(condition = "Tabuchi2011b")
df_ENRICH_Tabuchi2011c <- dplyr::filter(df_gseaRes.Tabuchi2011c, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust, NES) %>%
  dplyr::mutate(condition = "Tabuchi2011c")
df_ENRICH_Tabuchi2011d <- dplyr::filter(df_gseaRes.Tabuchi2011d, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust, NES) %>%
  dplyr::mutate(condition = "Tabuchi2011d")
df_ENRICH_Tabuchi2011e <- dplyr::filter(df_gseaRes.Tabuchi2011e, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust, NES) %>%
  dplyr::mutate(condition = "Tabuchi2011e")
df_ENRICH_Tabuchi2011f <- dplyr::filter(df_gseaRes.Tabuchi2011f, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust, NES) %>%
  dplyr::mutate(condition = "Tabuchi2011f")
df_ENRICH_Tabuchi2011g <- dplyr::filter(df_gseaRes.Tabuchi2011g, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust, NES) %>%
  dplyr::mutate(condition = "Tabuchi2011g")
df_ENRICH_Tabuchi2011h <- dplyr::filter(df_gseaRes.Tabuchi2011h, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust, NES)%>%
  dplyr::mutate(condition = "Tabuchi2011h")
df_ENRICH_Amaya2014a <- dplyr::filter(df_gseaRes.Amaya2014a, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust, NES) %>%
  dplyr::mutate(condition = "Amaya2014a")
df_ENRICH_Amaya2014b <- dplyr::filter(df_gseaRes.Amaya2014b, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust, NES) %>%
  dplyr::mutate(condition = "Amaya2014b")
df_ENRICH_Amaya2014c <- dplyr::filter(df_gseaRes.Amaya2014c, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust, NES) %>%
  dplyr::mutate(condition = "Amaya2014c")
df_ENRICH_Court2017 <- dplyr::filter(df_gseaRes.Court2017, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust, NES) %>%
  dplyr::mutate(condition = "Court2017")
df_ENRICH_Andocs2015a <- dplyr::filter(df_gseaRes.Andocs2015a, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust, NES) %>%
  dplyr::mutate(condition = "Andocs2015a")
df_ENRICH_Andocs2015b <- dplyr::filter(df_gseaRes.Andocs2015b, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust, NES) %>%
  dplyr::mutate(condition = "Andocs2015b")
df_ENRICH_Yunoki2016 <- dplyr::filter(df_gseaRes.Yunoki2016, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust, NES) %>%
  dplyr::mutate(condition = "Yunoki2016")
df_ENRICH_Scutigliani2022a <- dplyr::filter(df_gseaRes.Scutigliani2022a, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust, NES) %>%
  dplyr::mutate(condition = "Scutigliani2022a")
df_ENRICH_Scutigliani2022b <- dplyr::filter(df_gseaRes.Scutigliani2022b, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust, NES) %>%
  dplyr::mutate(condition = "Scutigliani2022b")

df_ENRICH_all <- rbind(df_ENRICH_Tabuchi2008,
                       df_ENRICH_Tabuchi2011a,
                       df_ENRICH_Tabuchi2011b,
                       df_ENRICH_Tabuchi2011c,
                       df_ENRICH_Tabuchi2011d,
                       df_ENRICH_Tabuchi2011e,
                       df_ENRICH_Tabuchi2011f,
                       df_ENRICH_Tabuchi2011g,
                       df_ENRICH_Tabuchi2011h,
                       df_ENRICH_Amaya2014a,
                       df_ENRICH_Amaya2014b,
                       df_ENRICH_Amaya2014c,
                       df_ENRICH_Court2017,
                       df_ENRICH_Andocs2015a,
                       df_ENRICH_Andocs2015b,
                       df_ENRICH_Yunoki2016,
                       df_ENRICH_Scutigliani2022a,
                       df_ENRICH_Scutigliani2022b)

df_ENRICH_all <- df_ENRICH_all %>% group_by(Description)
df_ENRICH_all$Description <- as.factor(df_ENRICH_all$Description)

#write.csv(df_ENRICH_all, "GSEA_H_barplot_all.csv")



#### LEADING EDGE ANALYSIS - HALLMARKS - ALL DATASETS ####

df_ENRICH_Tabuchi2008 <- dplyr::filter(df_gseaRes.Tabuchi2008, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2008")
df_ENRICH_Tabuchi2011a <- dplyr::filter(df_gseaRes.Tabuchi2011a, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011a")
df_ENRICH_Tabuchi2011b <- dplyr::filter(df_gseaRes.Tabuchi2011b, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011b")
df_ENRICH_Tabuchi2011c <- dplyr::filter(df_gseaRes.Tabuchi2011c, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011c")
df_ENRICH_Tabuchi2011d <- dplyr::filter(df_gseaRes.Tabuchi2011d, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011d")
df_ENRICH_Tabuchi2011e <- dplyr::filter(df_gseaRes.Tabuchi2011e, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011e")
df_ENRICH_Tabuchi2011f <- dplyr::filter(df_gseaRes.Tabuchi2011f, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011f")
df_ENRICH_Tabuchi2011g <- dplyr::filter(df_gseaRes.Tabuchi2011g, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011g")
df_ENRICH_Tabuchi2011h <- dplyr::filter(df_gseaRes.Tabuchi2011h, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust)%>%
  dplyr::mutate(condition = "Tabuchi2011h")
df_ENRICH_Amaya2014a <- dplyr::filter(df_gseaRes.Amaya2014a, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Amaya2014a")
df_ENRICH_Amaya2014b <- dplyr::filter(df_gseaRes.Amaya2014b, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Amaya2014b")
df_ENRICH_Amaya2014c <- dplyr::filter(df_gseaRes.Amaya2014c, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Amaya2014c")
df_ENRICH_Court2017 <- dplyr::filter(df_gseaRes.Court2017, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Court2017")
df_ENRICH_Andocs2015a <- dplyr::filter(df_gseaRes.Andocs2015a, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Andocs2015a")
df_ENRICH_Andocs2015b <- dplyr::filter(df_gseaRes.Andocs2015b, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Andocs2015b")
df_ENRICH_Yunoki2016 <- dplyr::filter(df_gseaRes.Yunoki2016, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Yunoki2016")
df_ENRICH_Scutigliani2022a <- dplyr::filter(df_gseaRes.Scutigliani2022a, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Scutigliani2022a")
df_ENRICH_Scutigliani2022b <- dplyr::filter(df_gseaRes.Scutigliani2022b, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Scutigliani2022b")

df_ENRICH_all <- rbind(df_ENRICH_Tabuchi2008,
                       df_ENRICH_Tabuchi2011a,
                       df_ENRICH_Tabuchi2011b,
                       df_ENRICH_Tabuchi2011c,
                       df_ENRICH_Tabuchi2011d,
                       df_ENRICH_Tabuchi2011e,
                       df_ENRICH_Tabuchi2011f,
                       df_ENRICH_Tabuchi2011g,
                       df_ENRICH_Tabuchi2011h,
                       df_ENRICH_Amaya2014a,
                       df_ENRICH_Amaya2014b,
                       df_ENRICH_Amaya2014c,
                       df_ENRICH_Court2017,
                       df_ENRICH_Andocs2015a,
                       df_ENRICH_Andocs2015b,
                       df_ENRICH_Yunoki2016,
                       df_ENRICH_Scutigliani2022a,
                       df_ENRICH_Scutigliani2022b)

list_all_terms <- unique(df_ENRICH_all$Description)

df_ENRICH_Tabuchi2008_sign <- dplyr::filter(df_gseaRes.Tabuchi2008, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, core_enrichment) %>%
  dplyr::mutate(condition = "Tabuchi2008")
df_ENRICH_Tabuchi2011a_sign <- dplyr::filter(df_gseaRes.Tabuchi2011a, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, core_enrichment) %>%
  dplyr::mutate(condition = "Tabuchi2011a")
df_ENRICH_Tabuchi2011b_sign <- dplyr::filter(df_gseaRes.Tabuchi2011b, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, core_enrichment) %>%
  dplyr::mutate(condition = "Tabuchi2011b")
df_ENRICH_Tabuchi2011c_sign <- dplyr::filter(df_gseaRes.Tabuchi2011c, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, core_enrichment) %>%
  dplyr::mutate(condition = "Tabuchi2011c")
df_ENRICH_Tabuchi2011d_sign <- dplyr::filter(df_gseaRes.Tabuchi2011d, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, core_enrichment) %>%
  dplyr::mutate(condition = "Tabuchi2011d")
df_ENRICH_Tabuchi2011e_sign <- dplyr::filter(df_gseaRes.Tabuchi2011e, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, core_enrichment) %>%
  dplyr::mutate(condition = "Tabuchi2011e")
df_ENRICH_Tabuchi2011f_sign <- dplyr::filter(df_gseaRes.Tabuchi2011f, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, core_enrichment) %>%
  dplyr::mutate(condition = "Tabuchi2011f")
df_ENRICH_Tabuchi2011g_sign <- dplyr::filter(df_gseaRes.Tabuchi2011g, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, core_enrichment) %>%
  dplyr::mutate(condition = "Tabuchi2011g")
df_ENRICH_Tabuchi2011h_sign <- dplyr::filter(df_gseaRes.Tabuchi2011h, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, core_enrichment)%>%
  dplyr::mutate(condition = "Tabuchi2011h")
df_ENRICH_Amaya2014a_sign <- dplyr::filter(df_gseaRes.Amaya2014a, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, core_enrichment) %>%
  dplyr::mutate(condition = "Amaya2014a")
df_ENRICH_Amaya2014b_sign <- dplyr::filter(df_gseaRes.Amaya2014b, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, core_enrichment) %>%
  dplyr::mutate(condition = "Amaya2014b")
df_ENRICH_Amaya2014c_sign <- dplyr::filter(df_gseaRes.Amaya2014c, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, core_enrichment) %>%
  dplyr::mutate(condition = "Amaya2014c")
df_ENRICH_Court2017_sign <- dplyr::filter(df_gseaRes.Court2017, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, core_enrichment) %>%
  dplyr::mutate(condition = "Court2017")
df_ENRICH_Andocs2015a_sign <- dplyr::filter(df_gseaRes.Andocs2015a, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, core_enrichment) %>%
  dplyr::mutate(condition = "Andocs2015a")
df_ENRICH_Andocs2015b_sign <- dplyr::filter(df_gseaRes.Andocs2015b, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, core_enrichment) %>%
  dplyr::mutate(condition = "Andocs2015b")
df_ENRICH_Yunoki2016_sign <- dplyr::filter(df_gseaRes.Yunoki2016, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, core_enrichment) %>%
  dplyr::mutate(condition = "Yunoki2016")
df_ENRICH_Scutigliani2022a_sign <- dplyr::filter(df_gseaRes.Scutigliani2022a, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, core_enrichment) %>%
  dplyr::mutate(condition = "Scutigliani2022a")
df_ENRICH_Scutigliani2022b_sign <- dplyr::filter(df_gseaRes.Scutigliani2022b, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, core_enrichment) %>%
  dplyr::mutate(condition = "Scutigliani2022b")

df_ENRICH_sign <- rbind(df_ENRICH_Tabuchi2008_sign,
                        df_ENRICH_Tabuchi2011a_sign,
                        df_ENRICH_Tabuchi2011b_sign,
                        df_ENRICH_Tabuchi2011c_sign,
                        df_ENRICH_Tabuchi2011d_sign,
                        df_ENRICH_Tabuchi2011e_sign,
                        df_ENRICH_Tabuchi2011f_sign,
                        df_ENRICH_Tabuchi2011g_sign,
                        df_ENRICH_Tabuchi2011h_sign,
                        df_ENRICH_Amaya2014a_sign,
                        df_ENRICH_Amaya2014b_sign,
                        df_ENRICH_Amaya2014c_sign,
                        df_ENRICH_Court2017_sign,
                        df_ENRICH_Andocs2015a_sign,
                        df_ENRICH_Andocs2015b_sign,
                        df_ENRICH_Yunoki2016_sign,
                        df_ENRICH_Scutigliani2022a_sign,
                        df_ENRICH_Scutigliani2022b_sign)

df_ENRICH_leading_edge <- df_ENRICH_sign %>%
  dplyr::select(Description, core_enrichment) %>%
  tibble::remove_rownames() %>%
  tidyr::separate(core_enrichment, sep = "/", paste("target", 1:200, sep = "_")) %>%
  tidyr::pivot_longer(paste("target", 1:200, sep = "_"), names_to = "Target_n", values_to = "Target") %>%
  dplyr::filter(!is.na(Target)) %>%
  dplyr::select(Description, Target) %>%
  dplyr::mutate(Presence = "1") %>%
  unique() %>%
  tidyr::pivot_wider(Target, names_from = Description, values_from = Presence) %>%
  tibble::column_to_rownames("Target") %>%
  replace(is.na(.), 0)

# write.csv(df_ENRICH_leading_edge, "leading_edge_hallmarks_all_genes.csv")

df_ENRICH_leading_edge[1:ncol(df_ENRICH_leading_edge)] <- as.numeric(unlist(df_ENRICH_leading_edge[1:ncol(df_ENRICH_leading_edge)]))

df_ENRICH_leading_edge <- df_ENRICH_leading_edge %>% dplyr::filter(rowSums(df_ENRICH_leading_edge[1:ncol(df_ENRICH_leading_edge)]) > 5)

pheatmap(df_ENRICH_leading_edge,
         color = c("white", "lightblue"),
         cellwidth = 10,
         cellheight = 10,
         annotation_names_col = F,
         show_rownames = T,
         legend = F,
         filename = "leading_edge_hallmarks.pdf")
dev.off()

df_ENRICH_leading_edge_gene_count <- df_ENRICH_sign %>%
  dplyr::select(Description, core_enrichment) %>%
  tibble::remove_rownames() %>%
  tidyr::separate(core_enrichment, sep = "/", paste("target", 1:127, sep = "_")) %>%
  tidyr::pivot_longer(paste("target", 1:127, sep = "_"), names_to = "Target_n", values_to = "Target") %>%
  dplyr::filter(!is.na(Target)) %>%
  dplyr::select(Description, Target) %>%
  dplyr::mutate(Presence = "1") %>%
  unique() %>%
  ungroup() %>%
  group_by(Target) %>%
  summarize(count = n()) %>%
  arrange(-count)

# write.csv(df_ENRICH_leading_edge_gene_count, "leading_edge_hallmarks_all_genes_count_all_datasets.csv")

ggplot(df_ENRICH_leading_edge_gene_count, aes(x = fct_reorder(Target, count), y = count)) +
  geom_col() +
  labs(title="Leading edge analysis") + xlab("") + ylab("count") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(NA),
        panel.grid.minor = element_line(NA),
        title = element_text(colour = "black", size=19),
        axis.title.x = element_text(colour = "black", size=22),
        axis.title.y = element_text(colour = "black", size=22),
        axis.line.x = element_line(colour = 'black', size = .5),
        axis.line.y = element_line(colour = 'black', size = .5),
        axis.ticks.x = element_line(colour = "black", size = 0),
        axis.ticks.y = element_line(colour = "black", size = 0),
        axis.text.x = element_text(colour = "black", size=9),
        axis.text.y = element_text(colour = "black", size=9),
        aspect.ratio = 2) +
  coord_flip()

ggsave("leading_edge_genes_top20.pdf", dpi=1000)



#### LEADING EDGE ANALYSIS - HALLMARKS - ANTITABUCHI CLUSTER ####

df_ENRICH_sign <- rbind(df_ENRICH_Amaya2014a_sign,
                        df_ENRICH_Amaya2014c_sign,
                        df_ENRICH_Court2017_sign,
                        df_ENRICH_Andocs2015a_sign,
                        df_ENRICH_Andocs2015b_sign,
                        df_ENRICH_Yunoki2016_sign,
                        df_ENRICH_Scutigliani2022a_sign,
                        df_ENRICH_Scutigliani2022b_sign)

df_ENRICH_leading_edge <- df_ENRICH_sign %>%
  dplyr::select(Description, core_enrichment) %>%
  tibble::remove_rownames() %>%
  tidyr::separate(core_enrichment, sep = "/", paste("target", 1:200, sep = "_")) %>%
  tidyr::pivot_longer(paste("target", 1:200, sep = "_"), names_to = "Target_n", values_to = "Target") %>%
  dplyr::filter(!is.na(Target)) %>%
  dplyr::select(Description, Target) %>%
  dplyr::mutate(Presence = "1") %>%
  unique() %>%
  tidyr::pivot_wider(Target, names_from = Description, values_from = Presence) %>%
  tibble::column_to_rownames("Target") %>%
  replace(is.na(.), 0)

# write.csv(df_ENRICH_leading_edge, "leading_edge_hallmarks_all_genes_antitabuchi.csv")

df_ENRICH_leading_edge[1:ncol(df_ENRICH_leading_edge)] <- as.numeric(unlist(df_ENRICH_leading_edge[1:ncol(df_ENRICH_leading_edge)]))

df_ENRICH_leading_edge <- df_ENRICH_leading_edge %>%
  dplyr::filter(rowSums(.) > 4) %>%
  dplyr::select_if(colSums(.) > 0)

pheatmap(df_ENRICH_leading_edge,
         color = c("white", "lightblue"),
         cellwidth = 10,
         cellheight = 10,
         annotation_names_col = F,
         show_rownames = T,
         legend = F,
         filename = "leading_edge_hallmarks_antitabuchi_cluster.pdf")
dev.off()

df_ENRICH_leading_edge_gene_count <- df_ENRICH_sign %>%
  dplyr::select(Description, core_enrichment) %>%
  tibble::remove_rownames() %>%
  tidyr::separate(core_enrichment, sep = "/", paste("target", 1:200, sep = "_")) %>%
  tidyr::pivot_longer(paste("target", 1:200, sep = "_"), names_to = "Target_n", values_to = "Target") %>%
  dplyr::filter(!is.na(Target)) %>%
  dplyr::select(Description, Target) %>%
  dplyr::mutate(Presence = "1") %>%
  unique() %>%
  ungroup() %>%
  group_by(Target) %>%
  summarize(count = n()) %>%
  arrange(-count)

# write.csv(df_ENRICH_leading_edge_gene_count, "leading_edge_hallmarks_all_genes_count_antitabuchi.csv")

ggplot(df_ENRICH_leading_edge_gene_count, aes(x = fct_reorder(Target, count), y = count)) +
  geom_col() +
  labs(title="Leading edge analysis") + xlab("") + ylab("count") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(NA),
        panel.grid.minor = element_line(NA),
        title = element_text(colour = "black", size=19),
        axis.title.x = element_text(colour = "black", size=22),
        axis.title.y = element_text(colour = "black", size=22),
        axis.line.x = element_line(colour = 'black', size = .5),
        axis.line.y = element_line(colour = 'black', size = .5),
        axis.ticks.x = element_line(colour = "black", size = 0),
        axis.ticks.y = element_line(colour = "black", size = 0),
        axis.text.x = element_text(colour = "black", size=9),
        axis.text.y = element_text(colour = "black", size=9),
        aspect.ratio = 2) +
  coord_flip()

ggsave("leading_edge_genes_top20_antitabuchi_cluster.pdf", dpi=1000)



#### LEADING EDGE ANALYSIS - HALLMARKS - TABUCHI CLUSTER ####

df_ENRICH_sign <- rbind(df_ENRICH_Tabuchi2008_sign,
                        df_ENRICH_Tabuchi2011a_sign,
                        df_ENRICH_Tabuchi2011b_sign,
                        df_ENRICH_Tabuchi2011c_sign,
                        df_ENRICH_Tabuchi2011d_sign,
                        df_ENRICH_Tabuchi2011e_sign,
                        df_ENRICH_Tabuchi2011f_sign,
                        df_ENRICH_Tabuchi2011g_sign,
                        df_ENRICH_Tabuchi2011h_sign,
                        df_ENRICH_Amaya2015_sign)

df_ENRICH_leading_edge <- df_ENRICH_sign %>%
  dplyr::select(Description, core_enrichment) %>%
  tibble::remove_rownames() %>%
  tidyr::separate(core_enrichment, sep = "/", paste("target", 1:127, sep = "_")) %>%
  tidyr::pivot_longer(paste("target", 1:127, sep = "_"), names_to = "Target_n", values_to = "Target") %>%
  dplyr::filter(!is.na(Target)) %>%
  dplyr::select(Description, Target) %>%
  dplyr::mutate(Presence = "1") %>%
  unique() %>%
  tidyr::pivot_wider(Target, names_from = Description, values_from = Presence) %>%
  tibble::column_to_rownames("Target") %>%
  replace(is.na(.), 0)

# write.csv(df_ENRICH_leading_edge, "leading_edge_hallmarks_all_genes_tabuchi.csv")

df_ENRICH_leading_edge[1:ncol(df_ENRICH_leading_edge)] <- as.numeric(unlist(df_ENRICH_leading_edge[1:ncol(df_ENRICH_leading_edge)]))

df_ENRICH_leading_edge <- df_ENRICH_leading_edge %>%
  dplyr::filter(rowSums(.) > 3) %>%
  dplyr::select_if(colSums(.) > 0)

pheatmap(df_ENRICH_leading_edge,
         color = c("white", "lightblue"),
         cellwidth = 10,
         cellheight = 10,
         annotation_names_col = F,
         show_rownames = T,
         legend = F,
         filename = "leading_edge_hallmarks_tabuchi_cluster.pdf")
dev.off()

df_ENRICH_leading_edge_gene_count <- df_ENRICH_sign %>%
  dplyr::select(Description, core_enrichment) %>%
  tibble::remove_rownames() %>%
  tidyr::separate(core_enrichment, sep = "/", paste("target", 1:127, sep = "_")) %>%
  tidyr::pivot_longer(paste("target", 1:127, sep = "_"), names_to = "Target_n", values_to = "Target") %>%
  dplyr::filter(!is.na(Target)) %>%
  dplyr::select(Description, Target) %>%
  dplyr::mutate(Presence = "1") %>%
  unique() %>%
  ungroup() %>%
  group_by(Target) %>%
  summarize(count = n()) %>%
  arrange(-count)

# write.csv(df_ENRICH_leading_edge_gene_count, "leading_edge_hallmarks_all_genes_count_tabuchi.csv")

ggplot(df_ENRICH_leading_edge_gene_count, aes(x = fct_reorder(Target, count), y = count)) +
  geom_col() +
  labs(title="Leading edge analysis") + xlab("") + ylab("count") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(NA),
        panel.grid.minor = element_line(NA),
        title = element_text(colour = "black", size=19),
        axis.title.x = element_text(colour = "black", size=22),
        axis.title.y = element_text(colour = "black", size=22),
        axis.line.x = element_line(colour = 'black', size = .5),
        axis.line.y = element_line(colour = 'black', size = .5),
        axis.ticks.x = element_line(colour = "black", size = 0),
        axis.ticks.y = element_line(colour = "black", size = 0),
        axis.text.x = element_text(colour = "black", size=9),
        axis.text.y = element_text(colour = "black", size=9),
        aspect.ratio = 2) +
  coord_flip()

ggsave("leading_edge_genes_top20_tabuchi_cluster.pdf", dpi=1000)

#### LEADING EDGE ANALYSIS - OVERLAP IN GENES BETWEEN CLUSTERS ####

genecount_antitabuchi <- read.csv("leading_edge_ORA_H_cluster1.csv")
genecount_tabuchi <- read.csv("leading_edge_ORA_H_cluster2.csv")

genecount_antitabuchi_selection <- genecount_antitabuchi %>% dplyr::filter(count > 4)
genecount_tabuchi_selection <- genecount_tabuchi %>% dplyr::filter(count > 2)

df_venn <- list(antitabuchi = genecount_antitabuchi_selection$Target,
                tabuchi = genecount_tabuchi_selection$Target)
partitions <- get.venn.partitions(df_venn)

df_overlap <- unlist(partitions$..values..$`1`)

# select background

m_t2g <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>% #selecting a species and category of terms
  dplyr::select(gs_name, gene_symbol) #filtering gene name per process as input for GSEA

# ORA

df_enrich_antitabuchi <- clusterProfiler::enricher(gene = df_intersections_antitabuchi$Gene, TERM2GENE = m_t2g)
df_enrich_tabuchi <- clusterProfiler::enricher(gene = df_intersections_tabuchi$Gene, TERM2GENE = m_t2g)
df_enrich_overlap <- clusterProfiler::enricher(gene = df_overlap, TERM2GENE = m_t2g)

df_enrich_antitabuchi <- df_enrich_antitabuchi@result
df_enrich_tabuchi <- df_enrich_tabuchi@result
df_enrich_overlap <- df_enrich_overlap@result

# export

write.csv(df_enrich_antitabuchi, "ORA_SumOver3_cluster1.csv")

#### ORA - H - ALL DATASETS ####

# select DEGs

df_DEG <- df_all %>%
  group_by(Publication) %>%
  filter(Ratio > 1.5 | Ratio < 1/1.5)

# prepare gene list per dataset

df_ORA_Tabuchi2008 <- df_DEG %>% dplyr::filter(Publication %in% c("Tabuchi2008")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Tabuchi2011a <- df_DEG %>% dplyr::filter(Publication %in% c("Tabuchi2011a")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Tabuchi2011b <- df_DEG %>% dplyr::filter(Publication %in% c("Tabuchi2011b")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Tabuchi2011c <- df_DEG %>% dplyr::filter(Publication %in% c("Tabuchi2011c")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Tabuchi2011d <- df_DEG %>% dplyr::filter(Publication %in% c("Tabuchi2011d")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Tabuchi2011e <- df_DEG %>% dplyr::filter(Publication %in% c("Tabuchi2011e")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Tabuchi2011f <- df_DEG %>% dplyr::filter(Publication %in% c("Tabuchi2011f")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Tabuchi2011g <- df_DEG %>% dplyr::filter(Publication %in% c("Tabuchi2011g")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Tabuchi2011h <- df_DEG %>% dplyr::filter(Publication %in% c("Tabuchi2011h")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Amaya2014a <- df_DEG %>% dplyr::filter(Publication %in% c("Amaya2014a")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Amaya2014b <- df_DEG %>% dplyr::filter(Publication %in% c("Amaya2014b")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Amaya2014c <- df_DEG %>% dplyr::filter(Publication %in% c("Amaya2014c")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Court2017 <- df_DEG %>% dplyr::filter(Publication %in% c("Court2017")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Andocs2015a <- df_DEG %>% dplyr::filter(Publication %in% c("Andocs2015a")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Yunoki2016 <- df_DEG %>% dplyr::filter(Publication %in% c("Yunoki2016")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Scutigliani2022a <- df_DEG %>% dplyr::filter(Publication %in% c("Scutigliani2022a")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Scutigliani2022b <- df_DEG %>% dplyr::filter(Publication %in% c("Scutigliani2022b")) %>% ungroup() %>% dplyr::select(Gene)

# select background

m_t2g <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>% #selecting a species and category of terms
  dplyr::select(gs_name, gene_symbol) #filtering gene name per process as input for GSEA

# ORA

df_enrich_Tabuchi2008 <- clusterProfiler::enricher(gene = df_ORA_Tabuchi2008$Gene, TERM2GENE = m_t2g)
df_enrich_Tabuchi2011a <- clusterProfiler::enricher(gene = df_ORA_Tabuchi2011a$Gene, TERM2GENE = m_t2g)
df_enrich_Tabuchi2011b <- clusterProfiler::enricher(gene = df_ORA_Tabuchi2011b$Gene, TERM2GENE = m_t2g)
df_enrich_Tabuchi2011c <- clusterProfiler::enricher(gene = df_ORA_Tabuchi2011c$Gene, TERM2GENE = m_t2g)
df_enrich_Tabuchi2011d <- clusterProfiler::enricher(gene = df_ORA_Tabuchi2011d$Gene, TERM2GENE = m_t2g)
df_enrich_Tabuchi2011e <- clusterProfiler::enricher(gene = df_ORA_Tabuchi2011e$Gene, TERM2GENE = m_t2g)
df_enrich_Tabuchi2011f <- clusterProfiler::enricher(gene = df_ORA_Tabuchi2011f$Gene, TERM2GENE = m_t2g)
df_enrich_Tabuchi2011g <- clusterProfiler::enricher(gene = df_ORA_Tabuchi2011g$Gene, TERM2GENE = m_t2g)
df_enrich_Tabuchi2011h <- clusterProfiler::enricher(gene = df_ORA_Tabuchi2011h$Gene, TERM2GENE = m_t2g)
df_enrich_Amaya2014a <- clusterProfiler::enricher(gene = df_ORA_Amaya2014a$Gene, TERM2GENE = m_t2g)
df_enrich_Amaya2014b <- clusterProfiler::enricher(gene = df_ORA_Amaya2014b$Gene, TERM2GENE = m_t2g)
df_enrich_Amaya2014c <- clusterProfiler::enricher(gene = df_ORA_Amaya2014c$Gene, TERM2GENE = m_t2g)
df_enrich_Court2017 <- clusterProfiler::enricher(gene = df_ORA_Court2017$Gene, TERM2GENE = m_t2g)
df_enrich_Andocs2015a <- clusterProfiler::enricher(gene = df_ORA_Andocs2015a$Gene, TERM2GENE = m_t2g)
df_enrich_Yunoki2016 <- clusterProfiler::enricher(gene = df_ORA_Yunoki2016$Gene, TERM2GENE = m_t2g)
df_enrich_Scutigliani2022a <- clusterProfiler::enricher(gene = df_ORA_Scutigliani2022a$Gene, TERM2GENE = m_t2g)
df_enrich_Scutigliani2022b <- clusterProfiler::enricher(gene = df_ORA_Scutigliani2022b$Gene, TERM2GENE = m_t2g)

df_enrich_Tabuchi2008 <- df_enrich_Tabuchi2008@result
df_enrich_Tabuchi2011a <- df_enrich_Tabuchi2011a@result
df_enrich_Tabuchi2011b <- df_enrich_Tabuchi2011b@result
df_enrich_Tabuchi2011c <- df_enrich_Tabuchi2011c@result
df_enrich_Tabuchi2011d <- df_enrich_Tabuchi2011d@result
df_enrich_Tabuchi2011e <- df_enrich_Tabuchi2011e@result
df_enrich_Tabuchi2011f <- df_enrich_Tabuchi2011f@result
df_enrich_Tabuchi2011g <- df_enrich_Tabuchi2011g@result
df_enrich_Tabuchi2011h <- df_enrich_Tabuchi2011h@result
df_enrich_Amaya2014a <- df_enrich_Amaya2014a@result
df_enrich_Amaya2014b <- df_enrich_Amaya2014b@result
df_enrich_Amaya2014c <- df_enrich_Amaya2014c@result
df_enrich_Court2017 <- df_enrich_Court2017@result
df_enrich_Andocs2015a <- df_enrich_Andocs2015a@result
df_enrich_Yunoki2016 <- df_enrich_Yunoki2016@result
df_enrich_Scutigliani2022a <- df_enrich_Scutigliani2022a@result
df_enrich_Scutigliani2022b <- df_enrich_Scutigliani2022b@result

# visualization in heatmap

# filter signifcantly enriched terms

df_ENRICH_Tabuchi2008 <- dplyr::filter(df_enrich_Tabuchi2008, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2008")
df_ENRICH_Tabuchi2011a <- dplyr::filter(df_enrich_Tabuchi2011a, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011a")
df_ENRICH_Tabuchi2011b <- dplyr::filter(df_enrich_Tabuchi2011b, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011b")
df_ENRICH_Tabuchi2011c <- dplyr::filter(df_enrich_Tabuchi2011c, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011c")
df_ENRICH_Tabuchi2011d <- dplyr::filter(df_enrich_Tabuchi2011d, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011d")
df_ENRICH_Tabuchi2011e <- dplyr::filter(df_enrich_Tabuchi2011e, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011e")
df_ENRICH_Tabuchi2011f <- dplyr::filter(df_enrich_Tabuchi2011f, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011f")
df_ENRICH_Tabuchi2011g <- dplyr::filter(df_enrich_Tabuchi2011g, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011g")
df_ENRICH_Tabuchi2011h <- dplyr::filter(df_enrich_Tabuchi2011h, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust)%>%
  dplyr::mutate(condition = "Tabuchi2011h")
df_ENRICH_Amaya2014a <- dplyr::filter(df_enrich_Amaya2014a, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Amaya2014a")
df_ENRICH_Amaya2014b <- dplyr::filter(df_enrich_Amaya2014b, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Amaya2014b")
df_ENRICH_Amaya2014c <- dplyr::filter(df_enrich_Amaya2014c, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Amaya2014c")
df_ENRICH_Court2017 <- dplyr::filter(df_enrich_Court2017, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Court2017")
df_ENRICH_Andocs2015a <- dplyr::filter(df_enrich_Andocs2015a, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Andocs2015a")

df_ENRICH_Yunoki2016 <- dplyr::filter(df_enrich_Yunoki2016, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Yunoki2016")
df_ENRICH_Scutigliani2022a <- dplyr::filter(df_enrich_Scutigliani2022a, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Scutigliani2022a")
df_ENRICH_Scutigliani2022b <- dplyr::filter(df_enrich_Scutigliani2022b, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Scutigliani2022b")

df_ENRICH_all <- rbind(df_ENRICH_Tabuchi2008,
                       df_ENRICH_Tabuchi2011a,
                       df_ENRICH_Tabuchi2011b,
                       df_ENRICH_Tabuchi2011c,
                       df_ENRICH_Tabuchi2011d,
                       df_ENRICH_Tabuchi2011e,
                       df_ENRICH_Tabuchi2011f,
                       df_ENRICH_Tabuchi2011g,
                       df_ENRICH_Tabuchi2011h,
                       df_ENRICH_Amaya2014a,
                       df_ENRICH_Amaya2014b,
                       df_ENRICH_Amaya2014c,
                       df_ENRICH_Court2017,
                       df_ENRICH_Andocs2015a,
                       df_ENRICH_Yunoki2016,
                       df_ENRICH_Scutigliani2022a,
                       df_ENRICH_Scutigliani2022b)
list_all_terms <- unique(df_ENRICH_all$Description)

df_ENRICH_Tabuchi2008_sign <- dplyr::filter(df_enrich_Tabuchi2008, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2008")
df_ENRICH_Tabuchi2011a_sign <- dplyr::filter(df_enrich_Tabuchi2011a, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011a")
df_ENRICH_Tabuchi2011b_sign <- dplyr::filter(df_enrich_Tabuchi2011b, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011b")
df_ENRICH_Tabuchi2011c_sign <- dplyr::filter(df_enrich_Tabuchi2011c, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011c")
df_ENRICH_Tabuchi2011d_sign <- dplyr::filter(df_enrich_Tabuchi2011d, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011d")
df_ENRICH_Tabuchi2011e_sign <- dplyr::filter(df_enrich_Tabuchi2011e, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011e")
df_ENRICH_Tabuchi2011f_sign <- dplyr::filter(df_enrich_Tabuchi2011f, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011f")
df_ENRICH_Tabuchi2011g_sign <- dplyr::filter(df_enrich_Tabuchi2011g, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011g")
df_ENRICH_Tabuchi2011h_sign <- dplyr::filter(df_enrich_Tabuchi2011h, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust)%>%
  dplyr::mutate(condition = "Tabuchi2011h")
df_ENRICH_Amaya2014a_sign <- dplyr::filter(df_enrich_Amaya2014a, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Amaya2014a")
df_ENRICH_Amaya2014b_sign <- dplyr::filter(df_enrich_Amaya2014b, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Amaya2014b")
df_ENRICH_Amaya2014c_sign <- dplyr::filter(df_enrich_Amaya2014c, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Amaya2014c")
df_ENRICH_Court2017_sign <- dplyr::filter(df_enrich_Court2017, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Court2017")
df_ENRICH_Andocs2015a_sign <- dplyr::filter(df_enrich_Andocs2015a, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Andocs2015a")
df_ENRICH_Yunoki2016_sign <- dplyr::filter(df_enrich_Yunoki2016, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Yunoki2016")
df_ENRICH_Scutigliani2022a_sign <- dplyr::filter(df_enrich_Scutigliani2022a, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Scutigliani2022a")
df_ENRICH_Scutigliani2022b_sign <- dplyr::filter(df_enrich_Scutigliani2022b, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Scutigliani2022b")

df_ENRICH_sign <- rbind(df_ENRICH_Tabuchi2008_sign,
                        df_ENRICH_Tabuchi2011a_sign,
                        df_ENRICH_Tabuchi2011b_sign,
                        df_ENRICH_Tabuchi2011c_sign,
                        df_ENRICH_Tabuchi2011d_sign,
                        df_ENRICH_Tabuchi2011e_sign,
                        df_ENRICH_Tabuchi2011f_sign,
                        df_ENRICH_Tabuchi2011g_sign,
                        df_ENRICH_Tabuchi2011h_sign,
                        df_ENRICH_Amaya2014a_sign,
                        df_ENRICH_Amaya2014b_sign,
                        df_ENRICH_Amaya2014c_sign,
                        df_ENRICH_Court2017_sign,
                        df_ENRICH_Andocs2015a_sign,
                        df_ENRICH_Yunoki2016_sign,
                        df_ENRICH_Scutigliani2022a_sign,
                        df_ENRICH_Scutigliani2022b_sign)

write.csv(df_ENRICH_sign, "ORA_H_ALL_DATASETS.csv")

for (i in 1:nrow(df_ENRICH_sign))
  if (df_ENRICH_sign$p.adjust[i] < 0.05) df_ENRICH_sign$sign[i] = 1 else df_ENRICH_sign$sign[i] = 0

df_ENRICH_binary <- dplyr::select(df_ENRICH_sign, -p.adjust)
df_ENRICH_shaped <- pivot_wider(df_ENRICH_binary, names_from = condition, values_from = sign) %>%
  replace(is.na(.), 0)

df_heatmap_num <- df_ENRICH_shaped[2:ncol(df_ENRICH_shaped)]
rownames(df_heatmap_num) = df_ENRICH_shaped$Description
df_heatmap_input <- as.matrix(df_heatmap_num)

df_all_metrics <- df_all %>%
  dplyr::select(Publication, temperature_C, duration_m, timepoint_h, CEM43) %>%
  distinct()

publications <- df_all_metrics %>%
  dplyr::filter(Publication %in% colnames(df_heatmap_input)) %>%
  dplyr::select(Publication)
publications <- publications$Publication

temp <- df_all_metrics %>%
  dplyr::filter(Publication %in% colnames(df_heatmap_input)) %>%
  dplyr::select(temperature_C)
temp <- temp$temperature_C

time <- df_all_metrics %>%
  dplyr::filter(Publication %in% colnames(df_heatmap_input)) %>%
  dplyr::select(timepoint_h)
time <- time$timepoint_h

CEM43 <- df_all_metrics %>%
  dplyr::filter(Publication %in% colnames(df_heatmap_input)) %>%
  dplyr::select(CEM43)
CEM43 <- CEM43$CEM43

cat_df = data.frame("time" = time,
                    "CEM43" = CEM43)

rownames(cat_df) = colnames(df_heatmap_input)

ann_colors = list(
  time = c("0" = "#FEFFD9", "1" = "#EDF8BC", "3" = "#BFE8B2", "4" = "#6BC6BE", "6" = "#63C3BF", "24" = "#41B7C4"),
  CEM43 = c("1.875" = "#FFC685", "3.75" = "#F8B15B", "7.5" = "#F5A040","15" = "#b56f60","30" = "#ED7222", "120" = "#E26320", "180" = "#BA4C23", "480" = "#9E3D22"))

heatmap <- pheatmap(df_heatmap_input,
         color = c("white", "lightblue"),
         annotation_col = cat_df,
         annotation_color = ann_colors,
         cellheight = 10,
         cellwidth = 20,
         annotation_names_col = F,
         legend = F)
dev.off()

# export for figlinq

rownames <- rownames(df_heatmap_input[heatmap$tree_row[["order"]],])
colnames <- colnames(df_heatmap_input[,heatmap$tree_col[["order"]]])

df_heatmap_input.sorted <- df_heatmap_input[match(rownames, rownames(df_heatmap_input)),]
df_heatmap_input.sorted <- df_heatmap_input.sorted[,match(colnames, colnames(df_heatmap_input.sorted))]

write.csv(df_heatmap_input.sorted, "heatmap_ORA_H_alldata_figlinq.csv")



#### ORA - HALLMARKS - BARPLOTS OF NES/P ####

# filter signifcantly enriched terms

df_ENRICH_Tabuchi2008 <- dplyr::filter(df_ENRICH_Tabuchi2008, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2008")
df_ENRICH_Tabuchi2011a <- dplyr::filter(df_ENRICH_Tabuchi2011a, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011a")
df_ENRICH_Tabuchi2011b <- dplyr::filter(df_ENRICH_Tabuchi2011b, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011b")
df_ENRICH_Tabuchi2011c <- dplyr::filter(df_ENRICH_Tabuchi2011c, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011c")
df_ENRICH_Tabuchi2011d <- dplyr::filter(df_ENRICH_Tabuchi2011d, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011d")
df_ENRICH_Tabuchi2011e <- dplyr::filter(df_ENRICH_Tabuchi2011e, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011e")
df_ENRICH_Tabuchi2011f <- dplyr::filter(df_ENRICH_Tabuchi2011f, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011f")
df_ENRICH_Tabuchi2011g <- dplyr::filter(df_ENRICH_Tabuchi2011g, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011g")
df_ENRICH_Tabuchi2011h <- dplyr::filter(df_ENRICH_Tabuchi2011h, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust)%>%
  dplyr::mutate(condition = "Tabuchi2011h")
df_ENRICH_Amaya2014a <- dplyr::filter(df_ENRICH_Amaya2014a, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Amaya2014a")
df_ENRICH_Amaya2014b <- dplyr::filter(df_ENRICH_Amaya2014b, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Amaya2014b")
df_ENRICH_Amaya2014c <- dplyr::filter(df_ENRICH_Amaya2014c, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Amaya2014c")
df_ENRICH_Court2017 <- dplyr::filter(df_ENRICH_Court2017, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Court2017")
df_ENRICH_Andocs2015a <- dplyr::filter(df_ENRICH_Andocs2015a, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Andocs2015a")
df_ENRICH_Andocs2015b <- dplyr::filter(df_ENRICH_Andocs2015b, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Andocs2015b")
df_ENRICH_Yunoki2016 <- dplyr::filter(df_ENRICH_Yunoki2016, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Yunoki2016")
df_ENRICH_Scutigliani2022a <- dplyr::filter(df_ENRICH_Scutigliani2022a, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Scutigliani2022a")
df_ENRICH_Scutigliani2022b <- dplyr::filter(df_ENRICH_Scutigliani2022b, p.adjust < 0.05) %>%

  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Scutigliani2022b")

df_ENRICH_all <- rbind(df_ENRICH_Tabuchi2008,
                       df_ENRICH_Tabuchi2011a,
                       df_ENRICH_Tabuchi2011b,
                       df_ENRICH_Tabuchi2011c,
                       df_ENRICH_Tabuchi2011d,
                       df_ENRICH_Tabuchi2011e,
                       df_ENRICH_Tabuchi2011f,
                       df_ENRICH_Tabuchi2011g,
                       df_ENRICH_Tabuchi2011h,
                       df_ENRICH_Amaya2014a,
                       df_ENRICH_Amaya2014b,
                       df_ENRICH_Amaya2014c,
                       df_ENRICH_Court2017,
                       df_ENRICH_Andocs2015a,
                       df_ENRICH_Yunoki2016,
                       df_ENRICH_Scutigliani2022a,
                       df_ENRICH_Scutigliani2022b)

df_ENRICH_all <- df_ENRICH_all %>% group_by(Description)
df_ENRICH_all$Description <- as.factor(df_ENRICH_all$Description)

#write.csv(df_ENRICH_all, "ORA_H_barplot_all.csv")

# plot

textsize <- 12
theme <-   theme(plot.title = element_text(size = textsize, face = "bold", hjust = 0.5),
                 plot.subtitle = element_text(size = textsize, hjust = 0.5),
                 panel.background = element_rect(fill="white"),
                 panel.grid.major = element_line(NA),
                 panel.grid.minor = element_line(NA),
                 axis.title.x = element_text(colour = "black", size=textsize),
                 axis.title.y = element_blank(),
                 axis.line.x = element_line(colour = 'black', size = 0.5),
                 axis.line.y = element_line(colour = 'black', size = 0.5),
                 axis.ticks.x = element_line(colour = 'black', size = 0.5),
                 axis.ticks.y = element_blank(),
                 axis.text.x = element_text(colour = "black", size=textsize),
                 axis.text.y = element_text(colour = "black", size=6),
                 legend.justification = c("right", "top"),
                 legend.background = element_blank(),
                 legend.title = element_text(face = "bold", size = textsize),
                 legend.text = element_text(size = textsize, colour = "black"),
                 legend.key = element_blank(),
                 aspect.ratio = 2)

ggplot(df_ENRICH_all, aes(x = fct_reorder(Description, -p.adjust, .fun = "median"), y = log10(p.adjust))) +
  geom_boxplot(outlier.shape = NA, color = "black", fill = "white") +
  geom_point(size = 1) +
  labs(x = "Hallmark", y = "log10(p-value)") +
  geom_hline(yintercept = log10(0.05), linetype = "dashed") +
  coord_flip() + theme

ggsave("Hallmarks_sign_all_datasets_boxplot.pdf", dpi=1000)

#### LEADING EDGE OF ORA - H - ALL DATASETS - DATA PREPARATION ####

df_ENRICH_Tabuchi2008_sign <- dplyr::filter(df_enrich_Tabuchi2008, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::mutate(condition = "Tabuchi2008")
df_ENRICH_Tabuchi2011a_sign <- dplyr::filter(df_enrich_Tabuchi2011a, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::mutate(condition = "Tabuchi2011a")
df_ENRICH_Tabuchi2011b_sign <- dplyr::filter(df_enrich_Tabuchi2011b, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::mutate(condition = "Tabuchi2011b")
df_ENRICH_Tabuchi2011c_sign <- dplyr::filter(df_enrich_Tabuchi2011c, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::mutate(condition = "Tabuchi2011c")
df_ENRICH_Tabuchi2011d_sign <- dplyr::filter(df_enrich_Tabuchi2011d, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::mutate(condition = "Tabuchi2011d")
df_ENRICH_Tabuchi2011e_sign <- dplyr::filter(df_enrich_Tabuchi2011e, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::mutate(condition = "Tabuchi2011e")
df_ENRICH_Tabuchi2011f_sign <- dplyr::filter(df_enrich_Tabuchi2011f, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::mutate(condition = "Tabuchi2011f")
df_ENRICH_Tabuchi2011g_sign <- dplyr::filter(df_enrich_Tabuchi2011g, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::mutate(condition = "Tabuchi2011g")
df_ENRICH_Tabuchi2011h_sign <- dplyr::filter(df_enrich_Tabuchi2011h, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, geneID)%>%
  dplyr::mutate(condition = "Tabuchi2011h")
df_ENRICH_Amaya2014a_sign <- dplyr::filter(df_enrich_Amaya2014a, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::mutate(condition = "Amaya2014a")
df_ENRICH_Amaya2014b_sign <- dplyr::filter(df_enrich_Amaya2014b, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::mutate(condition = "Amaya2014b")
df_ENRICH_Amaya2014c_sign <- dplyr::filter(df_enrich_Amaya2014c, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::mutate(condition = "Amaya2014c")
df_ENRICH_Court2017_sign <- dplyr::filter(df_enrich_Court2017, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::mutate(condition = "Court2017")
df_ENRICH_Andocs2015a_sign <- dplyr::filter(df_enrich_Andocs2015a, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::mutate(condition = "Andocs2015a")
df_ENRICH_Yunoki2016_sign <- dplyr::filter(df_enrich_Yunoki2016, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::mutate(condition = "Yunoki2016")
df_ENRICH_Scutigliani2022a_sign <- dplyr::filter(df_enrich_Scutigliani2022a, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::mutate(condition = "Scutigliani2022a")
df_ENRICH_Scutigliani2022b_sign <- dplyr::filter(df_enrich_Scutigliani2022b, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust, geneID) %>%
  dplyr::mutate(condition = "Scutigliani2022b")

#### LEADING EDGE OF ORA - H - ALL DATASETS ####

# select data

df_ORA_leading_edge_sign <- rbind(df_ENRICH_Tabuchi2008_sign,
                        df_ENRICH_Tabuchi2011a_sign,
                        df_ENRICH_Tabuchi2011b_sign,
                        df_ENRICH_Tabuchi2011c_sign,
                        df_ENRICH_Tabuchi2011d_sign,
                        df_ENRICH_Tabuchi2011e_sign,
                        df_ENRICH_Tabuchi2011f_sign,
                        df_ENRICH_Tabuchi2011g_sign,
                        df_ENRICH_Tabuchi2011h_sign,
                        df_ENRICH_Amaya2014a_sign,
                        df_ENRICH_Amaya2014b_sign,
                        df_ENRICH_Amaya2014c_sign,
                        df_ENRICH_Court2017_sign,
                        df_ENRICH_Andocs2015a_sign,
                        df_ENRICH_Yunoki2016_sign,
                        df_ENRICH_Scutigliani2022a_sign,
                        df_ENRICH_Scutigliani2022b_sign)

# count gene occurence in gene sets

df_ENRICH_leading_edge <- df_ORA_leading_edge_sign %>%
  dplyr::select(Description, geneID) %>%
  tibble::remove_rownames() %>%
  tidyr::separate(geneID, sep = "/", paste("target", 1:500, sep = "_")) %>%
  tidyr::pivot_longer(paste("target", 1:500, sep = "_"), names_to = "Target_n", values_to = "Target") %>%
  dplyr::filter(!is.na(Target)) %>%
  dplyr::select(Description, Target) %>%
  dplyr::mutate(Presence = "1") %>%
  unique() %>%
  tidyr::pivot_wider(Target, names_from = Description, values_from = Presence) %>%
  tibble::column_to_rownames("Target") %>%
  replace(is.na(.), 0)

# save progress

write.csv(df_ENRICH_leading_edge, "leading_edge_ORA_H_all_datasets.csv")

# selection criteria for minimal amount of occurences of genes in gene sets

df_ENRICH_leading_edge[1:ncol(df_ENRICH_leading_edge)] <- as.numeric(unlist(df_ENRICH_leading_edge[1:ncol(df_ENRICH_leading_edge)]))

df_ENRICH_leading_edge <- df_ENRICH_leading_edge %>%
  dplyr::filter(rowSums(.) > 3) %>%
  dplyr::select_if(colSums(.) > 10)

# visualization in heatmap

heatmap <- pheatmap(df_ENRICH_leading_edge,
                    color = c("white", "lightblue"),
                    annotation_names_col = F,
                    show_rownames = T,
                    show_colnames = F,
                    legend = F,
                    filename = "heatmap_leading_edge_ORA_H_all_datasets.pdf")
dev.off()

# export data

rownames <- rownames(df_ENRICH_leading_edge[heatmap$tree_row[["order"]],])
colnames <- colnames(df_ENRICH_leading_edge[,heatmap$tree_col[["order"]]])

df_heatmap_input.sorted <- df_ENRICH_leading_edge[match(rownames, rownames(df_ENRICH_leading_edge)),]
df_heatmap_input.sorted <- df_heatmap_input.sorted[,match(colnames, colnames(df_heatmap_input.sorted))]

write.csv(df_heatmap_input.sorted, "heatmap_leading_edge_ORA_H_all_datasets_cutoff3-10_figlinq.csv")

# top 30 gene occurences

df_ENRICH_leading_edge_gene_count <- df_ENRICH_sign %>%
  dplyr::select(Description, geneID) %>%
  tibble::remove_rownames() %>%
  tidyr::separate(geneID, sep = "/", paste("target", 1:400, sep = "_")) %>%
  tidyr::pivot_longer(paste("target", 1:400, sep = "_"), names_to = "Target_n", values_to = "Target") %>%
  dplyr::filter(!is.na(Target)) %>%
  dplyr::select(Description, Target) %>%
  dplyr::mutate(Presence = "1") %>%
  unique() %>%
  ungroup() %>%
  group_by(Target) %>%
  summarize(count = n()) %>%
  arrange(-count)

# save data

write.csv(df_ENRICH_leading_edge_gene_count, "leading_edge_H_all_genes_count_all_datasets.csv")



#### LEADING EDGE OF ORA - H - CLUSTER 1 ####

# select data

df_ENRICH_sign <- rbind(df_ENRICH_Amaya2014a_sign,
                        df_ENRICH_Amaya2014b_sign,
                        df_ENRICH_Amaya2014c_sign,
                        df_ENRICH_Court2017_sign,
                        df_ENRICH_Andocs2015a_sign,
                        df_ENRICH_Yunoki2016_sign,
                        df_ENRICH_Scutigliani2022a_sign,
                        df_ENRICH_Scutigliani2022b_sign,
                        df_ENRICH_Tabuchi2008_sign)

# count gene occurence in gene sets

df_ENRICH_leading_edge <- df_ENRICH_sign %>%
  dplyr::select(Description, geneID) %>%
  tibble::remove_rownames() %>%
  tidyr::separate(geneID, sep = "/", paste("target", 1:400, sep = "_")) %>%
  tidyr::pivot_longer(paste("target", 1:400, sep = "_"), names_to = "Target_n", values_to = "Target") %>%
  dplyr::filter(!is.na(Target)) %>%
  dplyr::select(Description, Target) %>%
  dplyr::mutate(Presence = "1") %>%
  unique() %>%
  tidyr::pivot_wider(Target, names_from = Description, values_from = Presence) %>%
  tibble::column_to_rownames("Target") %>%
  replace(is.na(.), 0)

# save progress

write.csv(df_ENRICH_leading_edge, "leading_edge_ORA_H_cluster1.csv")

# selection criteria for minimal amount of occurences of genes in gene sets

df_ENRICH_leading_edge[1:ncol(df_ENRICH_leading_edge)] <- as.numeric(unlist(df_ENRICH_leading_edge[1:ncol(df_ENRICH_leading_edge)]))

df_ENRICH_leading_edge_test <- df_ENRICH_leading_edge %>%
  dplyr::filter(rowSums(.) > 3) %>%
  dplyr::select_if(colSums(.) > 10)

# visualization in heatmap

heatmap <- pheatmap(df_ENRICH_leading_edge_test,
         color = c("white", "lightblue"),
         cellwidth = 5,
         cellheight = 5,
         annotation_names_col = F,
         show_rownames = T,
         legend = F,
         cutree_rows = 2,
         cutree_cols = 2,
         fontsize = 5,
         filename = "heatmap_leading_edge_ORA_H_cluster1.pdf")
dev.off()

# export data

rownames <- rownames(df_ENRICH_leading_edge[heatmap$tree_row[["order"]],])
colnames <- colnames(df_ENRICH_leading_edge[,heatmap$tree_col[["order"]]])

df_heatmap_input.sorted <- df_ENRICH_leading_edge[match(rownames, rownames(df_ENRICH_leading_edge)),]
df_heatmap_input.sorted <- df_heatmap_input.sorted[,match(colnames, colnames(df_heatmap_input.sorted))]

write.csv(df_heatmap_input.sorted, "heatmap_leading_edge_ORA_H_cluster1_cutoff3-10_figlinq.csv")

# top 30 gene occurences

df_ENRICH_leading_edge_gene_count <- df_ENRICH_sign %>%
  dplyr::select(Description, geneID) %>%
  tibble::remove_rownames() %>%
  tidyr::separate(geneID, sep = "/", paste("target", 1:400, sep = "_")) %>%
  tidyr::pivot_longer(paste("target", 1:400, sep = "_"), names_to = "Target_n", values_to = "Target") %>%
  dplyr::filter(!is.na(Target)) %>%
  dplyr::select(Description, Target) %>%
  dplyr::mutate(Presence = "1") %>%
  unique() %>%
  ungroup() %>%
  group_by(Target) %>%
  summarize(count = n()) %>%
  arrange(-count)

# save data

write.csv(df_ENRICH_leading_edge_gene_count, "leading_edge_H_all_genes_count_cluster1.csv")



#### LEADING EDGE OF ORA - H - CLUSTER 2 ####

# select data

df_ENRICH_sign <- rbind(df_ENRICH_Tabuchi2011a_sign,
                        df_ENRICH_Tabuchi2011b_sign,
                        df_ENRICH_Tabuchi2011c_sign,
                        df_ENRICH_Tabuchi2011d_sign,
                        df_ENRICH_Tabuchi2011e_sign,
                        df_ENRICH_Tabuchi2011f_sign,
                        df_ENRICH_Tabuchi2011g_sign,
                        df_ENRICH_Tabuchi2011h_sign)

# count gene occurence in gene sets

df_ENRICH_leading_edge <- df_ENRICH_sign %>%
  dplyr::select(Description, geneID) %>%
  tibble::remove_rownames() %>%
  tidyr::separate(geneID, sep = "/", paste("target", 1:400, sep = "_")) %>%
  tidyr::pivot_longer(paste("target", 1:400, sep = "_"), names_to = "Target_n", values_to = "Target") %>%
  dplyr::filter(!is.na(Target)) %>%
  dplyr::select(Description, Target) %>%
  dplyr::mutate(Presence = "1") %>%
  unique() %>%
  tidyr::pivot_wider(Target, names_from = Description, values_from = Presence) %>%
  tibble::column_to_rownames("Target") %>%
  replace(is.na(.), 0)

# save progress

write.csv(df_ENRICH_leading_edge, "leading_edge_ORA_H_cluster2.csv")

# selection criteria for minimal amount of occurences of genes in gene sets

df_ENRICH_leading_edge[1:ncol(df_ENRICH_leading_edge)] <- as.numeric(unlist(df_ENRICH_leading_edge[1:ncol(df_ENRICH_leading_edge)]))

df_ENRICH_leading_edge_test <- df_ENRICH_leading_edge %>%
  dplyr::filter(rowSums(.) > 2) %>%
  dplyr::select_if(colSums(.) > 30)

# visualization in heatmap

pheatmap(df_ENRICH_leading_edge_test,
         color = c("white", "lightblue"),
         cellwidth = 5,
         cellheight = 5,
         annotation_names_col = F,
         show_rownames = T,
         legend = F,
         fontsize = 5,
         filename = "heatmap_leading_edge_H_cluster2.pdf")
dev.off()

# export data

rownames <- rownames(df_ENRICH_leading_edge[heatmap$tree_row[["order"]],])
colnames <- colnames(df_ENRICH_leading_edge[,heatmap$tree_col[["order"]]])

df_heatmap_input.sorted <- df_ENRICH_leading_edge[match(rownames, rownames(df_ENRICH_leading_edge)),]
df_heatmap_input.sorted <- df_heatmap_input.sorted[,match(colnames, colnames(df_heatmap_input.sorted))]

write.csv(df_heatmap_input.sorted, "heatmap_leading_edge_ORA_H_cluster2_cutoff2-30_figlinq.csv")

# top 30 gene occurences

df_ENRICH_leading_edge_gene_count <- df_ENRICH_sign %>%
  dplyr::select(Description, geneID) %>%
  tibble::remove_rownames() %>%
  tidyr::separate(geneID, sep = "/", paste("target", 1:400, sep = "_")) %>%
  tidyr::pivot_longer(paste("target", 1:400, sep = "_"), names_to = "Target_n", values_to = "Target") %>%
  dplyr::filter(!is.na(Target)) %>%
  dplyr::select(Description, Target) %>%
  dplyr::mutate(Presence = "1") %>%
  unique() %>%
  ungroup() %>%
  group_by(Target) %>%
  summarize(count = n()) %>%
  arrange(-count)

# save data

write.csv(df_ENRICH_leading_edge_gene_count, "leading_edge_H_all_genes_count_cluster2.csv")

#### LEADING EDGE OF ORA - H - OVERLAP GENES BETWEEN CLUSTERS ####

genecount_antitabuchi <- read.csv("leading_edge_H_all_genes_count_cluster1.csv")
genecount_tabuchi <- read.csv("leading_edge_H_all_genes_count_cluster2.csv")

genecount_antitabuchi_selection <- genecount_antitabuchi %>% dplyr::filter(count > 4)
genecount_tabuchi_selection <- genecount_tabuchi %>% dplyr::filter(count > 3)

df_venn <- list(antitabuchi = genecount_antitabuchi_selection$Target,
                tabuchi = genecount_tabuchi_selection$Target)
partitions <- get.venn.partitions(df_venn)

df_overlap <- unlist(partitions$..values..$`1`)

# select background

m_t2g <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>% #selecting a species and category of terms
  dplyr::select(gs_name, gene_symbol) #filtering gene name per process as input for GSEA

# ORA

df_enrich_overlap <- clusterProfiler::enricher(gene = df_overlap, TERM2GENE = m_t2g)

df_enrich_overlap <- df_enrich_overlap@result

write.csv(df_enrich_overlap, "ORA_core_genes_leading_edge_ORA_H.csv")


#### ORA - GOBP - ALL DATASETS ####

# select DEGs

df_DEG <- df_all %>%
  group_by(Publication) %>%
  filter(Ratio > 1.5 | Ratio < 1/1.5)

# prepare gene list per dataset

df_ORA_Tabuchi2008 <- df_DEG %>% dplyr::filter(Publication %in% c("Tabuchi2008")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Tabuchi2011a <- df_DEG %>% dplyr::filter(Publication %in% c("Tabuchi2011a")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Tabuchi2011b <- df_DEG %>% dplyr::filter(Publication %in% c("Tabuchi2011b")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Tabuchi2011c <- df_DEG %>% dplyr::filter(Publication %in% c("Tabuchi2011c")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Tabuchi2011d <- df_DEG %>% dplyr::filter(Publication %in% c("Tabuchi2011d")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Tabuchi2011e <- df_DEG %>% dplyr::filter(Publication %in% c("Tabuchi2011e")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Tabuchi2011f <- df_DEG %>% dplyr::filter(Publication %in% c("Tabuchi2011f")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Tabuchi2011g <- df_DEG %>% dplyr::filter(Publication %in% c("Tabuchi2011g")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Tabuchi2011h <- df_DEG %>% dplyr::filter(Publication %in% c("Tabuchi2011h")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Amaya2014a <- df_DEG %>% dplyr::filter(Publication %in% c("Amaya2014a")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Amaya2014b <- df_DEG %>% dplyr::filter(Publication %in% c("Amaya2014b")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Amaya2014c <- df_DEG %>% dplyr::filter(Publication %in% c("Amaya2014c")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Court2017 <- df_DEG %>% dplyr::filter(Publication %in% c("Court2017")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Andocs2015a <- df_DEG %>% dplyr::filter(Publication %in% c("Andocs2015a")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Yunoki2016 <- df_DEG %>% dplyr::filter(Publication %in% c("Yunoki2016")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Scutigliani2022a <- df_DEG %>% dplyr::filter(Publication %in% c("Scutigliani2022a")) %>% ungroup() %>% dplyr::select(Gene)
df_ORA_Scutigliani2022b <- df_DEG %>% dplyr::filter(Publication %in% c("Scutigliani2022b")) %>% ungroup() %>% dplyr::select(Gene)

# select background

m_t2g <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>% #selecting a species and category of terms
  dplyr::select(gs_name, gene_symbol) #filtering gene name per process as input for GSEA

# ORA

df_enrich_Tabuchi2008 <- clusterProfiler::enricher(gene = df_ORA_Tabuchi2008$Gene, TERM2GENE = m_t2g)
df_enrich_Tabuchi2011a <- clusterProfiler::enricher(gene = df_ORA_Tabuchi2011a$Gene, TERM2GENE = m_t2g)
df_enrich_Tabuchi2011b <- clusterProfiler::enricher(gene = df_ORA_Tabuchi2011b$Gene, TERM2GENE = m_t2g)
df_enrich_Tabuchi2011c <- clusterProfiler::enricher(gene = df_ORA_Tabuchi2011c$Gene, TERM2GENE = m_t2g)
df_enrich_Tabuchi2011d <- clusterProfiler::enricher(gene = df_ORA_Tabuchi2011d$Gene, TERM2GENE = m_t2g)
df_enrich_Tabuchi2011e <- clusterProfiler::enricher(gene = df_ORA_Tabuchi2011e$Gene, TERM2GENE = m_t2g)
df_enrich_Tabuchi2011f <- clusterProfiler::enricher(gene = df_ORA_Tabuchi2011f$Gene, TERM2GENE = m_t2g)
df_enrich_Tabuchi2011g <- clusterProfiler::enricher(gene = df_ORA_Tabuchi2011g$Gene, TERM2GENE = m_t2g)
df_enrich_Tabuchi2011h <- clusterProfiler::enricher(gene = df_ORA_Tabuchi2011h$Gene, TERM2GENE = m_t2g)
df_enrich_Amaya2014a <- clusterProfiler::enricher(gene = df_ORA_Amaya2014a$Gene, TERM2GENE = m_t2g)
df_enrich_Amaya2014b <- clusterProfiler::enricher(gene = df_ORA_Amaya2014b$Gene, TERM2GENE = m_t2g)
df_enrich_Amaya2014c <- clusterProfiler::enricher(gene = df_ORA_Amaya2014c$Gene, TERM2GENE = m_t2g)
df_enrich_Court2017 <- clusterProfiler::enricher(gene = df_ORA_Court2017$Gene, TERM2GENE = m_t2g)
df_enrich_Andocs2015a <- clusterProfiler::enricher(gene = df_ORA_Andocs2015a$Gene, TERM2GENE = m_t2g)
df_enrich_Yunoki2016 <- clusterProfiler::enricher(gene = df_ORA_Yunoki2016$Gene, TERM2GENE = m_t2g)
df_enrich_Scutigliani2022a <- clusterProfiler::enricher(gene = df_ORA_Scutigliani2022a$Gene, TERM2GENE = m_t2g)
df_enrich_Scutigliani2022b <- clusterProfiler::enricher(gene = df_ORA_Scutigliani2022b$Gene, TERM2GENE = m_t2g)

df_enrich_Tabuchi2008 <- df_enrich_Tabuchi2008@result
df_enrich_Tabuchi2011a <- df_enrich_Tabuchi2011a@result
df_enrich_Tabuchi2011b <- df_enrich_Tabuchi2011b@result
df_enrich_Tabuchi2011c <- df_enrich_Tabuchi2011c@result
df_enrich_Tabuchi2011d <- df_enrich_Tabuchi2011d@result
df_enrich_Tabuchi2011e <- df_enrich_Tabuchi2011e@result
df_enrich_Tabuchi2011f <- df_enrich_Tabuchi2011f@result
df_enrich_Tabuchi2011g <- df_enrich_Tabuchi2011g@result
df_enrich_Tabuchi2011h <- df_enrich_Tabuchi2011h@result
df_enrich_Amaya2014a <- df_enrich_Amaya2014a@result
df_enrich_Amaya2014b <- df_enrich_Amaya2014b@result
df_enrich_Amaya2014c <- df_enrich_Amaya2014c@result
df_enrich_Court2017 <- df_enrich_Court2017@result
df_enrich_Andocs2015a <- df_enrich_Andocs2015a@result
df_enrich_Yunoki2016 <- df_enrich_Yunoki2016@result
df_enrich_Scutigliani2022a <- df_enrich_Scutigliani2022a@result
df_enrich_Scutigliani2022b <- df_enrich_Scutigliani2022b@result

# visualization in heatmap

# filter signifcantly enriched terms

df_ENRICH_Tabuchi2008 <- dplyr::filter(df_enrich_Tabuchi2008, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2008")
df_ENRICH_Tabuchi2011a <- dplyr::filter(df_enrich_Tabuchi2011a, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011a")
df_ENRICH_Tabuchi2011b <- dplyr::filter(df_enrich_Tabuchi2011b, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011b")
df_ENRICH_Tabuchi2011c <- dplyr::filter(df_enrich_Tabuchi2011c, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011c")
df_ENRICH_Tabuchi2011d <- dplyr::filter(df_enrich_Tabuchi2011d, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011d")
df_ENRICH_Tabuchi2011e <- dplyr::filter(df_enrich_Tabuchi2011e, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011e")
df_ENRICH_Tabuchi2011f <- dplyr::filter(df_enrich_Tabuchi2011f, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011f")
df_ENRICH_Tabuchi2011g <- dplyr::filter(df_enrich_Tabuchi2011g, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011g")
df_ENRICH_Tabuchi2011h <- dplyr::filter(df_enrich_Tabuchi2011h, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust)%>%
  dplyr::mutate(condition = "Tabuchi2011h")
df_ENRICH_Amaya2014a <- dplyr::filter(df_enrich_Amaya2014a, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Amaya2014a")
df_ENRICH_Amaya2014b <- dplyr::filter(df_enrich_Amaya2014b, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Amaya2014b")
df_ENRICH_Amaya2014c <- dplyr::filter(df_enrich_Amaya2014c, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Amaya2014c")
df_ENRICH_Court2017 <- dplyr::filter(df_enrich_Court2017, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Court2017")
df_ENRICH_Andocs2015a <- dplyr::filter(df_enrich_Andocs2015a, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Andocs2015a")

df_ENRICH_Yunoki2016 <- dplyr::filter(df_enrich_Yunoki2016, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Yunoki2016")
df_ENRICH_Scutigliani2022a <- dplyr::filter(df_enrich_Scutigliani2022a, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Scutigliani2022a")
df_ENRICH_Scutigliani2022b <- dplyr::filter(df_enrich_Scutigliani2022b, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Scutigliani2022b")

df_ENRICH_all <- rbind(df_ENRICH_Tabuchi2008,
                       df_ENRICH_Tabuchi2011a,
                       df_ENRICH_Tabuchi2011b,
                       df_ENRICH_Tabuchi2011c,
                       df_ENRICH_Tabuchi2011d,
                       df_ENRICH_Tabuchi2011e,
                       df_ENRICH_Tabuchi2011f,
                       df_ENRICH_Tabuchi2011g,
                       df_ENRICH_Tabuchi2011h,
                       df_ENRICH_Amaya2014a,
                       df_ENRICH_Amaya2014b,
                       df_ENRICH_Amaya2014c,
                       df_ENRICH_Court2017,
                       df_ENRICH_Andocs2015a,
                       df_ENRICH_Yunoki2016,
                       df_ENRICH_Scutigliani2022a,
                       df_ENRICH_Scutigliani2022b)
list_all_terms <- unique(df_ENRICH_all$Description)

df_ENRICH_Tabuchi2008_sign <- dplyr::filter(df_enrich_Tabuchi2008, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2008")
df_ENRICH_Tabuchi2011a_sign <- dplyr::filter(df_enrich_Tabuchi2011a, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011a")
df_ENRICH_Tabuchi2011b_sign <- dplyr::filter(df_enrich_Tabuchi2011b, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011b")
df_ENRICH_Tabuchi2011c_sign <- dplyr::filter(df_enrich_Tabuchi2011c, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011c")
df_ENRICH_Tabuchi2011d_sign <- dplyr::filter(df_enrich_Tabuchi2011d, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011d")
df_ENRICH_Tabuchi2011e_sign <- dplyr::filter(df_enrich_Tabuchi2011e, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011e")
df_ENRICH_Tabuchi2011f_sign <- dplyr::filter(df_enrich_Tabuchi2011f, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011f")
df_ENRICH_Tabuchi2011g_sign <- dplyr::filter(df_enrich_Tabuchi2011g, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011g")
df_ENRICH_Tabuchi2011h_sign <- dplyr::filter(df_enrich_Tabuchi2011h, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust)%>%
  dplyr::mutate(condition = "Tabuchi2011h")
df_ENRICH_Amaya2014a_sign <- dplyr::filter(df_enrich_Amaya2014a, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Amaya2014a")
df_ENRICH_Amaya2014b_sign <- dplyr::filter(df_enrich_Amaya2014b, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Amaya2014b")
df_ENRICH_Amaya2014c_sign <- dplyr::filter(df_enrich_Amaya2014c, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Amaya2014c")
df_ENRICH_Court2017_sign <- dplyr::filter(df_enrich_Court2017, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Court2017")
df_ENRICH_Andocs2015a_sign <- dplyr::filter(df_enrich_Andocs2015a, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Andocs2015a")
df_ENRICH_Yunoki2016_sign <- dplyr::filter(df_enrich_Yunoki2016, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Yunoki2016")
df_ENRICH_Scutigliani2022a_sign <- dplyr::filter(df_enrich_Scutigliani2022a, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Scutigliani2022a")
df_ENRICH_Scutigliani2022b_sign <- dplyr::filter(df_enrich_Scutigliani2022b, p.adjust < 0.05) %>%
  dplyr::filter(Description %in% list_all_terms) %>%
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Scutigliani2022b")

df_ENRICH_sign <- rbind(df_ENRICH_Tabuchi2008_sign,
                        df_ENRICH_Tabuchi2011a_sign,
                        df_ENRICH_Tabuchi2011b_sign,
                        df_ENRICH_Tabuchi2011c_sign,
                        df_ENRICH_Tabuchi2011d_sign,
                        df_ENRICH_Tabuchi2011e_sign,
                        df_ENRICH_Tabuchi2011f_sign,
                        df_ENRICH_Tabuchi2011g_sign,
                        df_ENRICH_Tabuchi2011h_sign,
                        df_ENRICH_Amaya2014a_sign,
                        df_ENRICH_Amaya2014b_sign,
                        df_ENRICH_Amaya2014c_sign,
                        df_ENRICH_Court2017_sign,
                        df_ENRICH_Andocs2015a_sign,
                        df_ENRICH_Yunoki2016_sign,
                        df_ENRICH_Scutigliani2022a_sign,
                        df_ENRICH_Scutigliani2022b_sign)

write.csv(df_ENRICH_sign, "ORA_H_ALL_DATASETS.csv")

for (i in 1:nrow(df_ENRICH_sign))
  if (df_ENRICH_sign$p.adjust[i] < 0.05) df_ENRICH_sign$sign[i] = 1 else df_ENRICH_sign$sign[i] = 0

df_ENRICH_binary <- dplyr::select(df_ENRICH_sign, -p.adjust)
df_ENRICH_shaped <- pivot_wider(df_ENRICH_binary, names_from = condition, values_from = sign) %>%
  replace(is.na(.), 0)

df_heatmap_num <- df_ENRICH_shaped[2:ncol(df_ENRICH_shaped)]
df_heatmap_num <- as.data.frame(df_heatmap_num)
rownames(df_heatmap_num) = df_ENRICH_shaped$Description

df_heatmap_plot <- df_heatmap_num %>%
  dplyr::filter(rowSums(.) > 4) %>%
  dplyr::select_if(colSums(.) > 0)

df_heatmap_input <- as.matrix(df_heatmap_plot)

df_all_metrics <- df_all %>%
  dplyr::select(Publication, temperature_C, duration_m, timepoint_h, CEM43) %>%
  distinct()

publications <- df_all_metrics %>%
  dplyr::filter(Publication %in% colnames(df_heatmap_input)) %>%
  dplyr::select(Publication)
publications <- publications$Publication

temp <- df_all_metrics %>%
  dplyr::filter(Publication %in% colnames(df_heatmap_input)) %>%
  dplyr::select(temperature_C)
temp <- temp$temperature_C

time <- df_all_metrics %>%
  dplyr::filter(Publication %in% colnames(df_heatmap_input)) %>%
  dplyr::select(timepoint_h)
time <- time$timepoint_h

CEM43 <- df_all_metrics %>%
  dplyr::filter(Publication %in% colnames(df_heatmap_input)) %>%
  dplyr::select(CEM43)
CEM43 <- CEM43$CEM43

cat_df = data.frame("time" = time,
                    "CEM43" = CEM43)

rownames(cat_df) = colnames(df_heatmap_input)

ann_colors = list(
  time = c("0" = "#FEFFD9", "1" = "#EDF8BC", "3" = "#BFE8B2", "4" = "#6BC6BE", "6" = "#63C3BF", "24" = "#41B7C4"),
  CEM43 = c("1.875" = "#FFC685", "3.75" = "#F8B15B", "7.5" = "#F5A040","15" = "#b56f60","30" = "#ED7222", "120" = "#E26320", "180" = "#BA4C23", "480" = "#9E3D22"))

heatmap <- pheatmap(df_heatmap_input,
                    color = c("white", "lightblue"),
                    annotation_col = cat_df,
                    annotation_color = ann_colors,
                    annotation_names_col = F,
                    annotation_names_row = T,
                    legend = F)
dev.off()

# export for figlinq

rownames <- rownames(df_heatmap_input[heatmap$tree_row[["order"]],])
colnames <- colnames(df_heatmap_input[,heatmap$tree_col[["order"]]])

df_heatmap_input.sorted <- df_heatmap_input[match(rownames, rownames(df_heatmap_input)),]
df_heatmap_input.sorted <- df_heatmap_input.sorted[,match(colnames, colnames(df_heatmap_input.sorted))]

write.csv(df_heatmap_input.sorted, "heatmap_ORA_GOBP_alldata_figlinq.csv")

#### ORA - GOBP - HAND-CURATED CATEGORIZATION ####

df_cell_cycle <- df_ENRICH_sign %>% 
  dplyr::filter(stringr::str_detect(Description, "CELL_CYCLE")) %>% 
  dplyr::filter(!str_detect(Description, "ARREST")) %>% 
  dplyr::mutate(category = "cell_cycle_progression")

df_protein_folding <- df_ENRICH_sign %>% 
  dplyr::filter(stringr::str_detect(Description, "PROTEIN_FOLDING|PROTEASOMAL"))%>% 
  dplyr::mutate(category = "proteostasis")

df_stress_stimuli <- df_ENRICH_sign %>% 
  dplyr::filter(stringr::str_detect(Description, "STRESS|RESPONSE_TO"))%>% 
  dplyr::mutate(category = "response_to_stress_and_stimuli")

df_cytoskeleton <- df_ENRICH_sign %>% 
  dplyr::filter(stringr::str_detect(Description, "CYTOSKELET")) %>% 
  dplyr::filter(!str_detect(Description, "CYTOKINES")) %>% 
  dplyr::mutate(category = "cytoskeleton")

df_junction_adhesion <- df_ENRICH_sign %>% 
  dplyr::filter(stringr::str_detect(Description, "JUNCTION|CELL_CELL"))%>% 
  dplyr::mutate(category = "cell_junction_and_adhesion")

df_RAS_pathway <- df_ENRICH_sign %>% 
  dplyr::filter(stringr::str_detect(Description, "RAS_PROTEIN|ERK|ERBB2"))%>% 
  dplyr::mutate(category = "RAS_signaling")

df_mitosis <- df_ENRICH_sign %>% 
  dplyr::filter(stringr::str_detect(Description, "MITOSIS|SPINDLE|CHROMOSOME_SEGR"))%>% 
  dplyr::mutate(category = "mitosis")

df_ORA_cat <- rbind(df_cell_cycle,
                    df_protein_folding,
                    df_stress_stimuli,
                    df_cytoskeleton,
                    df_junction_adhesion,
                    df_RAS_pathway,
                    df_mitosis)

df_ORA_cat_count <- df_ORA_cat %>% ungroup() %>% group_by(condition, category) %>% summarise(sum = n())

write.csv(df_ORA_cat_count, "ORA_GOBP_categories.csv")

#### ORA - GOBP - BARPLOTS OF NES/P ####

# filter signifcantly enriched terms

df_ENRICH_Tabuchi2008 <- dplyr::filter(df_ENRICH_Tabuchi2008, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2008")
df_ENRICH_Tabuchi2011a <- dplyr::filter(df_ENRICH_Tabuchi2011a, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011a")
df_ENRICH_Tabuchi2011b <- dplyr::filter(df_ENRICH_Tabuchi2011b, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011b")
df_ENRICH_Tabuchi2011c <- dplyr::filter(df_ENRICH_Tabuchi2011c, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011c")
df_ENRICH_Tabuchi2011d <- dplyr::filter(df_ENRICH_Tabuchi2011d, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011d")
df_ENRICH_Tabuchi2011e <- dplyr::filter(df_ENRICH_Tabuchi2011e, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011e")
df_ENRICH_Tabuchi2011f <- dplyr::filter(df_ENRICH_Tabuchi2011f, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011f")
df_ENRICH_Tabuchi2011g <- dplyr::filter(df_ENRICH_Tabuchi2011g, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Tabuchi2011g")
df_ENRICH_Tabuchi2011h <- dplyr::filter(df_ENRICH_Tabuchi2011h, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust)%>%
  dplyr::mutate(condition = "Tabuchi2011h")
df_ENRICH_Amaya2014a <- dplyr::filter(df_ENRICH_Amaya2014a, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Amaya2014a")
df_ENRICH_Amaya2014b <- dplyr::filter(df_ENRICH_Amaya2014b, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Amaya2014b")
df_ENRICH_Amaya2014c <- dplyr::filter(df_ENRICH_Amaya2014c, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Amaya2014c")
df_ENRICH_Court2017 <- dplyr::filter(df_ENRICH_Court2017, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Court2017")
df_ENRICH_Andocs2015a <- dplyr::filter(df_ENRICH_Andocs2015a, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Andocs2015a")
f_ENRICH_Yunoki2016 <- dplyr::filter(df_ENRICH_Yunoki2016, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Yunoki2016")
df_ENRICH_Scutigliani2022a <- dplyr::filter(df_ENRICH_Scutigliani2022a, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Scutigliani2022a")
df_ENRICH_Scutigliani2022b <- dplyr::filter(df_ENRICH_Scutigliani2022b, p.adjust < 0.05) %>%
  
  dplyr::select(Description, p.adjust) %>%
  dplyr::mutate(condition = "Scutigliani2022b")

df_ENRICH_all <- rbind(df_ENRICH_Tabuchi2008,
                       df_ENRICH_Tabuchi2011a,
                       df_ENRICH_Tabuchi2011b,
                       df_ENRICH_Tabuchi2011c,
                       df_ENRICH_Tabuchi2011d,
                       df_ENRICH_Tabuchi2011e,
                       df_ENRICH_Tabuchi2011f,
                       df_ENRICH_Tabuchi2011g,
                       df_ENRICH_Tabuchi2011h,
                       df_ENRICH_Amaya2014a,
                       df_ENRICH_Amaya2014b,
                       df_ENRICH_Amaya2014c,
                       df_ENRICH_Court2017,
                       df_ENRICH_Andocs2015a,
                       df_ENRICH_Yunoki2016,
                       df_ENRICH_Scutigliani2022a,
                       df_ENRICH_Scutigliani2022b)

df_ENRICH_all <- df_ENRICH_all %>% group_by(Description)
df_ENRICH_all$Description <- as.factor(df_ENRICH_all$Description)

#### ORA - GOBP - VENN ANALYSIS ####

# cluster 1

df_venn <- list(df_ENRICH_Tabuchi2008 = df_ENRICH_Tabuchi2008$Description,
                df_ENRICH_Amaya2014a = df_ENRICH_Amaya2014a$Description,
                df_ENRICH_Amaya2014b = df_ENRICH_Amaya2014b$Description,
                df_ENRICH_Amaya2014c = df_ENRICH_Amaya2014c$Description,
                df_ENRICH_Court2017 = df_ENRICH_Court2017$Description,
                df_ENRICH_Andocs2015a = df_ENRICH_Andocs2015a$Description,
                df_ENRICH_Scutigliani2022a = df_ENRICH_Scutigliani2022a$Description,
                df_ENRICH_Scutigliani2022b = df_ENRICH_Scutigliani2022b$Description)
partitions <- get.venn.partitions(df_venn)



# cluster 2

df_venn <- list(df_ENRICH_Tabuchi2011a = df_ENRICH_Tabuchi2011a$Description,
                df_ENRICH_Tabuchi2011b = df_ENRICH_Tabuchi2011b$Description,
                df_ENRICH_Tabuchi2011c = df_ENRICH_Tabuchi2011c$Description,
                df_ENRICH_Tabuchi2011d = df_ENRICH_Tabuchi2011d$Description,
                df_ENRICH_Tabuchi2011e = df_ENRICH_Tabuchi2011e$Description,
                df_ENRICH_Tabuchi2011f = df_ENRICH_Tabuchi2011f$Description,
                df_ENRICH_Tabuchi2011g = df_ENRICH_Tabuchi2011g$Description,
                df_ENRICH_Tabuchi2011h = df_ENRICH_Tabuchi2011h$Description)
partitions <- get.venn.partitions(df_venn)

#### UPSET PLOT - OVER/UNDEREXPRESSED GENES - ALL ####

# dataframe preparation

df_diff <- df_heatmap_plot %>% 
  dplyr::mutate(sum = rowSums(.)) %>%
  rownames_to_column(var = "GOBP_term")

# cluster 1

df_diff_cl1 <- df_heatmap_plot %>% 
  dplyr::select(Tabuchi2008, 
                Amaya2014a,
                Amaya2014b,
                Amaya2014c,
                Court2017,
                Scutigliani2022a,
                Scutigliani2022b) %>% 
  dplyr::mutate(sum = rowSums(.)) %>%
  rownames_to_column(var = "GOBP_term")

ggplot(df_diff_cl1, aes(x=sum)) +
  geom_bar(stat="count")

# cluster 2

df_diff_cl2 <- df_heatmap_plot %>% 
  dplyr::select(Tabuchi2011a,
                Tabuchi2011b,
                Tabuchi2011c,
                Tabuchi2011d,
                Tabuchi2011e,
                Tabuchi2011f,
                Tabuchi2011g,
                Tabuchi2011h) %>% 
  dplyr::mutate(sum = rowSums(.)) %>%
  rownames_to_column(var = "GOBP_term")

ggplot(df_diff_cl2, aes(x=sum)) +
  geom_bar(stat="count")

#### EXPRESSION OF PROTEINS INVOLVED IN GO RESPONSE TO HEAT ####

# select background

protein_folding <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>% #selecting a species and category of terms
  dplyr::select(gs_name, gene_symbol) %>%  #filtering gene name per process as input for GSEA
  dplyr::filter(gs_name %in% c("GOBP_RESPONSE_TO_HEAT"))

df_protein_folding_exp <- df_all %>% 
  dplyr::filter(Gene %in% protein_folding$gene_symbol) %>% 
  dplyr::mutate(log10FC = log10(Ratio))

df_protein_folding_heatmap <- df_protein_folding_exp %>% 
  dplyr::select(Gene, log10FC, Publication) %>% 
  tidyr::pivot_wider(names_from = Gene, values_from = log10FC) 
df_protein_folding_heatmap[is.na(df_protein_folding_heatmap)] <- 0

df_heatmap_num <- df_protein_folding_heatmap[2:ncol(df_protein_folding_heatmap)]
df_heatmap_num <- as.data.frame(df_heatmap_num)
rownames(df_heatmap_num) = df_protein_folding_heatmap$Publication

df_heatmap_input <- t(as.matrix(df_heatmap_num))

df_all_metrics <- df_all %>%
  dplyr::select(Publication, temperature_C, duration_m, timepoint_h, CEM43) %>%
  distinct()

publications <- df_all_metrics %>%
  dplyr::filter(Publication %in% colnames(df_heatmap_input)) %>%
  dplyr::select(Publication)
publications <- publications$Publication

temp <- df_all_metrics %>%
  dplyr::filter(Publication %in% colnames(df_heatmap_input)) %>%
  dplyr::select(temperature_C)
temp <- temp$temperature_C

time <- df_all_metrics %>%
  dplyr::filter(Publication %in% colnames(df_heatmap_input)) %>%
  dplyr::select(timepoint_h)
time <- time$timepoint_h

CEM43 <- df_all_metrics %>%
  dplyr::filter(Publication %in% colnames(df_heatmap_input)) %>%
  dplyr::select(CEM43)
CEM43 <- CEM43$CEM43

cat_df = data.frame("time" = time,
                    "CEM43" = CEM43)

rownames(cat_df) = colnames(df_heatmap_input)

ann_colors = list(
  time = c("0" = "#FEFFD9", "1" = "#EDF8BC", "3" = "#BFE8B2", "4" = "#6BC6BE", "6" = "#63C3BF", "24" = "#41B7C4"),
  CEM43 = c("1.875" = "#FFC685", "3.75" = "#F8B15B", "7.5" = "#F5A040","15" = "#b56f60","30" = "#ED7222", "120" = "#E26320", "180" = "#BA4C23", "480" = "#9E3D22"))

heatmap <- pheatmap(df_heatmap_input,
                    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                    annotation_col = cat_df,
                    annotation_color = ann_colors,
                    annotation_names_col = T,
                    annotation_names_row = T,
                    legend = T,
                    scale = "none")
dev.off()
  
# export for figlinq

rownames <- rownames(df_heatmap_input[heatmap$tree_row[["order"]],])
colnames <- colnames(df_heatmap_input[,heatmap$tree_col[["order"]]])

df_heatmap_input.sorted <- df_heatmap_input[match(rownames, rownames(df_heatmap_input)),]
df_heatmap_input.sorted <- df_heatmap_input.sorted[,match(colnames, colnames(df_heatmap_input.sorted))]

write.csv(df_heatmap_input.sorted, "heatmap_GOBP_RESPONSE_TO_HEAT_figlinq.csv")

#### EXPRESSION OF PROTEINS INVOLVED IN PROTEIN FOLDING - PCA ####

df_all_wide <- df_protein_folding_exp %>%
  dplyr::select(-temperature_C, -duration_m, -timepoint_h, -CEM43, -log10FC) %>%
  pivot_wider(names_from = Publication, values_from = Ratio)

df_PCA <- df_all_wide %>%
  dplyr::select(-Gene) %>%
  na.omit()

df_PCA_input <- t(df_PCA)

samples.pca <- prcomp(df_PCA_input, center = TRUE, scale. = TRUE)

# plot eigenvalues

summary_PCA <- summary(samples.pca)
summary_PCA <- summary_PCA[["importance"]]
summary_PCA <- as.data.frame(summary_PCA)
summary_PCA <- summary_PCA %>% rownames_to_column(var = "factor")

summary_PCA_long <- summary_PCA %>%
  pivot_longer(!factor, names_to = "component", values_to = "value") %>%
  filter(factor %in% c("Cumulative Proportion")) %>%
  arrange(value)

summary_PCA_long$component <- factor(summary_PCA_long$component, levels = c("PC1", "PC2", "PC3", "PC4",
                                                                            "PC5", "PC6", "PC7", "PC8",
                                                                            "PC9", "PC10", "PC11", "PC12", "PC13",
                                                                            "PC14", "PC15", "PC16"))
summary_PCA_long_filter <- summary_PCA_long %>% filter(component %in% c("PC1", "PC2", "PC3", "PC4",
                                                                        "PC5", "PC6", "PC7", "PC8",
                                                                        "PC9", "PC10"))

write.csv(summary_PCA_long_filter, "eigenvalues.csv")

textsize <- 12
theme <-   theme(plot.title = element_text(size = textsize, face = "bold", hjust = 0.5),
                 plot.subtitle = element_text(size = textsize, hjust = 0.5),
                 panel.background = element_rect(fill="white"),
                 panel.grid.major = element_line(NA),
                 panel.grid.minor = element_line(NA),
                 axis.title.x = element_text(colour = "black", size=textsize),
                 axis.title.y = element_text(colour = "black", size=textsize),
                 axis.line.x = element_line(colour = 'black', size = 0.5),
                 axis.line.y = element_line(colour = 'black', size = 0.5),
                 axis.ticks.x = element_line(colour = "black", size = 0.5),
                 axis.ticks.y = element_line(colour = "black", size = 0.5),
                 axis.text.x = element_text(colour = "black", size=textsize),
                 axis.text.y = element_text(colour = "black", size=textsize),
                 legend.justification = c("right", "top"),
                 legend.background = element_blank(),
                 legend.title = element_text(face = "bold", size = textsize),
                 legend.text = element_text(size = textsize, colour = "black"),
                 legend.key = element_blank(),
                 aspect.ratio = 1)

ggplot(summary_PCA_long_filter, aes(x = component, y = value)) +
  geom_bar(stat = "identity") + theme +
  labs(title = "Eigenvalues", y = "Cumulative proportion of variance")

ggsave("eigenvalues.pdf", dpi = 1000)

# generate PCA plots

df_all_metrics <- df_protein_folding_exp %>%
  dplyr::select(Publication, temperature_C, duration_m, timepoint_h, CEM43) %>%
  distinct()

publications <- df_all_metrics %>%
  dplyr::select(Publication)
publications <- publications$Publication

temp <- df_all_metrics %>%
  dplyr::select(temperature_C)
temp <- temp$temperature_C

time <- df_all_metrics %>%
  dplyr::select(timepoint_h)
time <- time$timepoint_h

duration <- df_all_metrics %>%
  dplyr::select(duration_m)
duration <- duration$duration_m

CEM43 <- df_all_metrics %>%
  dplyr::select(CEM43)
CEM43 <- CEM43$CEM43

time <- factor(time, levels = c("0", "1", "3", "4", "6", "24"))
CEM43 <- as.numeric(CEM43)
CEM43 <- factor(CEM43, levels = c("1.875", "3.75", "7.5", "15", "30", "120", "180", "480"))

textsize <- 12
theme <-   theme(plot.title = element_text(size = textsize, face = "bold", hjust = 0.5),
                 plot.subtitle = element_text(size = textsize, hjust = 0.5),
                 panel.background = element_rect(fill="white"),
                 panel.grid.major = element_line(NA),
                 panel.grid.minor = element_line(NA),
                 axis.title.x = element_text(colour = "black", size=textsize),
                 axis.title.y = element_text(colour = "black", size=textsize),
                 axis.line.x = element_line(colour = 'black', size = 0.5),
                 axis.line.y = element_line(colour = 'black', size = 0.5),
                 axis.ticks.x = element_line(colour = "black", size = 0.5),
                 axis.ticks.y = element_line(colour = "black", size = 0.5),
                 axis.text.x = element_text(colour = "black", size=textsize),
                 axis.text.y = element_text(colour = "black", size=textsize),
                 legend.justification = c("right", "top"),
                 legend.background = element_blank(),
                 legend.title = element_text(face = "bold", size = textsize),
                 legend.text = element_text(size = textsize, colour = "black"),
                 legend.key = element_blank(),
                 aspect.ratio = 1)
scale <- scale_color_brewer(palette = "Set3")

p1 <- ggplot(samples.pca, aes(x=PC1, y=PC2)) +
  geom_point() + geom_text_repel(label=publications,
                                 color = "black",
                                 max.overlaps = 27,
                                 size = 4) + theme + scale
p2 <- ggplot(samples.pca, aes(x=PC2, y=PC3)) +
  geom_point() + geom_text_repel(label=publications,
                                 color = "black",
                                 max.overlaps = 26,
                                 size = 4) + theme + scale
cowplot::plot_grid(p1,p2,ncol=2, align = "vh")

ggsave("PCA_labels.pdf", dpi=1000)

p3 <- ggplot(samples.pca, aes(x=PC1, y=PC2, color = time)) +
  geom_point(aes(size = CEM43)) + theme + scale

p4 <- ggplot(samples.pca, aes(x=PC2, y=PC3, color = time)) +
  geom_point(aes(size = CEM43)) + theme + scale

cowplot::plot_grid(p3,p4,ncol=2, align = "vh")

ggsave("PCA_param.pdf", dpi=1000)

testp1 <- layer_data(p1) %>%
  mutate(publication = publications, CEM = CEM43, time = time, duration = duration)

write.csv(testp1, "PCA1vs2.csv")

testp2 <- layer_data(p2) %>%
  mutate(publication = publications, CEM = CEM43, time = time, duration = duration)

write.csv(testp2, "PCA2vs3.csv")