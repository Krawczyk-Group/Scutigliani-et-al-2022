
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

rm(list=ls()) 

#### DATA IMPORT/WRANGLING ####

df_amaya2014a <- read.csv("Amaya2014a.csv", sep=",") %>% 
  dplyr::rename(symbol = "Gene_Symbol", FC = "Amaya2014a") %>% 
  dplyr::mutate(publication = "Amaya2014a")

df_amaya2014b <- read.csv("Amaya2014b.csv", sep=",") %>% 
  dplyr::rename(symbol = "Gene_Symbol", FC = "Amaya2014a") %>% 
  dplyr::mutate(publication = "Amaya2014b")

df_amaya2014c <- read.csv("Amaya2014c.csv", sep=",") %>% 
  dplyr::rename(symbol = "Gene_Symbol", FC = "Amaya2014a") %>% 
  dplyr::mutate(publication = "Amaya2014c")

df_court2017 <- read.csv("Court2017.csv", sep=",") %>% 
  dplyr::rename(symbol = "Gene_Symbol", FC = "Court2017") %>% 
  dplyr::mutate(publication = "Court2017")

df_scutigliani2022a <- read.csv("Scutigliani2022a.csv", sep=",") %>% 
  dplyr::mutate(FC = 2^L2FC) %>% 
  dplyr::select(symbol, FC) %>% 
  na.omit() %>% 
  dplyr::mutate(publication = "Scutigliani2022a")

df_scutigliani2022b <- read.csv("Scutigliani2022b.csv", sep=",") %>% 
  dplyr::mutate(FC = 2^L2FC) %>% 
  dplyr::select(symbol, FC) %>% 
  na.omit() %>% 
  dplyr::mutate(publication = "Scutigliani2022b")

df_yunoki2016 <- read.csv("Yunoki2016.csv", sep=",") %>% 
  dplyr::rename(symbol = "Gene_Symbol", FC = "Yonoki2006") %>% 
  dplyr::mutate(publication = "Yunoki2016")

df_andocs2015a <- read.csv("Andocs2015a.csv", sep=",") %>% 
  dplyr::rename(symbol = "Gene_Symbol", FC = "Andocsa") %>% 
  dplyr::mutate(publication = "Andocs2015a")

df_andocs2015b <- read.csv("Andocs2015b.csv", sep=",") %>% 
  dplyr::rename(symbol = "Gene_Symbol", FC = "Andocsb") %>% 
  dplyr::mutate(publication = "Andocs2015b")

df_tabuchi2008 <- read.csv("Tabuchi2008.csv", sep=",") %>% 
  dplyr::rename(symbol = "Gene_Symbol", FC = "Tabuchi2008") %>% 
  dplyr::mutate(publication = "Tabuchi2008")

df_tabuchi2011a <- read.csv("Tabuchi2011a.csv", sep=",") %>% 
  dplyr::rename(symbol = "Gene_Symbol", FC = "Tabuchi2011a") %>% 
  dplyr::mutate(publication = "Tabuchi2011a")

df_tabuchi2011b <- read.csv("Tabuchi2011b.csv", sep=",") %>% 
  dplyr::rename(symbol = "Gene_Symbol", FC = "Tabuchi2011b") %>% 
  dplyr::mutate(publication = "Tabuchi2011b")

df_tabuchi2011c <- read.csv("Tabuchi2011c.csv", sep=",") %>% 
  dplyr::rename(symbol = "Gene_Symbol", FC = "Tabuchi2011c") %>% 
  dplyr::mutate(publication = "Tabuchi2011c")

df_tabuchi2011d <- read.csv("Tabuchi2011d.csv", sep=",") %>% 
  dplyr::rename(symbol = "Gene_Symbol", FC = "Tabuchi2011d") %>% 
  dplyr::mutate(publication = "Tabuchi2011d")

df_tabuchi2011e <- read.csv("Tabuchi2011e.csv", sep=",") %>% 
  dplyr::rename(symbol = "Gene_Symbol", FC = "Tabuchi2011e") %>% 
  dplyr::mutate(publication = "Tabuchi2011e")

df_tabuchi2011f <- read.csv("Tabuchi2011f.csv", sep=",") %>% 
  dplyr::rename(symbol = "Gene_Symbol", FC = "Tabuchi2011f") %>% 
  dplyr::mutate(publication = "Tabuchi2011f")

df_tabuchi2011g <- read.csv("Tabuchi2011g.csv", sep=",") %>% 
  dplyr::rename(symbol = "Gene_Symbol", FC = "Tabuchi2011g") %>% 
  dplyr::mutate(publication = "Tabuchi2011g")

df_tabuchi2011h <- read.csv("Tabuchi2011h.csv", sep=",") %>% 
  dplyr::rename(symbol = "Gene_Symbol", FC = "Tabuchi2011h") %>% 
  dplyr::mutate(publication = "Tabuchi2011h")

df_all <- rbind(df_amaya2014a,
                df_amaya2014b,
                df_amaya2014c,
                df_court2017,
                df_scutigliani2022a,
                df_scutigliani2022b,
                df_yunoki2016,
                df_andocs2015a,
                df_andocs2015b,
                df_tabuchi2008,
                df_tabuchi2011a,
                df_tabuchi2011b,
                df_tabuchi2011c,
                df_tabuchi2011d,
                df_tabuchi2011e,
                df_tabuchi2011f,
                df_tabuchi2011g,
                df_tabuchi2011h)

df_all <- df_all %>% dplyr::filter(FC > 0)

write.csv(df_all, "expression_all.csv")
