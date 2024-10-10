
## --------------------------------------
## Project Name: 
## Investigator(s):
## Author: Nicholas Larson PhD 
## Email:  larson.nicholas@mayo.edu 
## Script Name: 
## Purpose of script:
## Date Created: 2023-06-08
## Last Edited: 
##
## --------------------------------------

## Setup
library(tidyverse)

setwd("C:/Users/m139105/Documents/digipath/scripts/hsr_git")
dat <- read_csv('new_merged_df.csv')

dat.wide.resid <- dplyr::select(dat,Gene = Row.names,residual = epi_variance,type) %>% pivot_wider(values_from = residual,names_from = type)
dat.wide.epi  <- dat.wide.resid %>% arrange(-deep)


library(WebGestaltR)


gene.sig <- dplyr::filter(dat.wide.epi,deep>=0.30)$Gene
gene.all <- dat.wide.resid$Gene

test <- WebGestaltR(enrichMethod="ORA",
            organism="hsapiens", 
            enrichDatabase="geneontology_Biological_Process_noRedundant",
            interestGene = gene.sig, 
            interestGeneType="ensembl_gene_id",
            referenceGene = gene.all, 
            referenceGeneType="ensembl_gene_id")



dat.wide.resid <- dplyr::select(dat,Gene = Row.names,residual = residuals_variance,type) %>% pivot_wider(values_from = residual,names_from = type)
dat.wide.resid <- dat.wide.resid %>% arrange(deep)
