Author: Michael Peters

setwd("/projectnb2/bf528/users/dreadlocks/project_3/analyst")

# 5

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("limma")
# install.packages("ggplot2")

library(limma)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyverse)

# sample info dataframe with array_id and chemical columns
samples <- read.csv('/project/bf528/project_3/groups/group_6_mic_info.csv',as.is=TRUE)

# the full RMA normalized matrix of all experiments
rma <- read.table('/projectnb/bf528/project_3/samples/liver-normalization-rma.txt',
                  sep='\t',
                  as.is=TRUE,
                  header=TRUE,
                  row.names=1,
                  )


# subset the full expression matrix to just those in this comparison
rma.subset <- rma[paste0('X',samples$array_id[samples$chemical=='FLUCONAZOLE' | (samples$chemical=='Control' & samples$vehicle=='CORN_OIL_100_%')])]
rma.subset <- new("MAList",list(M=rma.subset,A=rma.subset))

# construct a design matrix modeling treatment vs control for use by limma
design <- model.matrix( ~factor( samples$chemical[samples$vehicle=='CORN_OIL_100_%'], levels=c('Control','FLUCONAZOLE') ) )
colnames(design) <- c('Intercept','FLUCONAZOLE')

# run limma
fit <- lmFit(rma.subset, design)
fit <- eBayes(fit)
t <- topTable(fit, coef=2, n=nrow(rma.subset), adjust='BH')

# write out the results to file
write.csv(t, paste0('FLUCONAZOLE', '_limma_results.csv') )

# sort and filter
t_filtered <- subset(t, adj.P.Val<0.05)
t_sorted_filtered <- t_filtered[order(t_filtered$adj.P.Val),]
dim(t_sorted_filtered)

# get top 10
top10 <- t_sorted_filtered[1:10,]
write.csv(top10, paste0('FLUCONAZOLE', '_topten.csv') )

# histogram
ggplot(t_sorted_filtered, aes(x=logFC)) + 
  geom_histogram(color="black", fill="cadetblue") +
  ggtitle('Fluconazole') +
  ylab("Count of Significant DE Genes")

ggsave( paste0('FLUCONAZOLE', '_histogram.jpeg'), plot = last_plot(), device = jpeg,
        width = 5, height = 5, units = 'in', dpi = 300, )

# scatter plot
ggplot(t_sorted_filtered, aes(x=logFC, y=P.Value)) +
  geom_point(shape=1,size = 1, color = "cadetblue") +
  ggtitle('Fluconazole') +
  ylab("Adj P-Value")

ggsave( paste0('FLUCONAZOLE', '_scatterplot.jpeg'), plot = last_plot(), device = jpeg,
        width = 5, height = 5, units = 'in', dpi = 300, )






rma.subset <- rma[paste0('X',samples$array_id[samples$chemical=='PIRINIXIC_ACID' | (samples$chemical=='Control' & samples$vehicle=='CMC_.5_%')])]
rma.subset <- new("MAList",list(M=rma.subset,A=rma.subset))

# construct a design matrix modeling treatment vs control for use by limma
design <- model.matrix( ~factor( samples$chemical[samples$vehicle=='CMC_.5_%'], levels=c('Control','PIRINIXIC_ACID') ) )
colnames(design) <- c('Intercept','PIRINIXIC_ACID')

# run limma
fit <- lmFit(rma.subset, design)
fit <- eBayes(fit)
t <- topTable(fit, coef=2, n=nrow(rma.subset), adjust='BH')

# write out the results to file
write.csv(t, paste0('PIRINIXIC_ACID', '_limma_results.csv') )

# sort and filter
t_filtered <- subset(t, adj.P.Val<0.05)
t_sorted_filtered <- t_filtered[order(t_filtered$adj.P.Val),]

# get top 10
top10 <- t_sorted_filtered[1:10,]
write.csv(top10, paste0('PIRINIXIC_ACID', '_topten.csv') )

# histogram
ggplot(t_sorted_filtered, aes(x=logFC)) + 
  geom_histogram(color="black", fill="cadetblue") +
  ggtitle('Pirinixic acid') +
  ylab("Count of Significant DE Genes")

ggsave( paste0('PIRINIXIC_ACID', '_histogram.jpeg'), plot = last_plot(), device = jpeg,
        width = 5, height = 5, units = 'in', dpi = 300, )

# scatter plot
ggplot(t_sorted_filtered, aes(x=logFC, y=P.Value)) +
  geom_point(shape=1,size = 1, color = "cadetblue") +
  ggtitle('Pirinixic acid') +
  ylab("Adj P-Value")

ggsave( paste0('PIRINIXIC_ACID', '_scatterplot.jpeg'), plot = last_plot(), device = jpeg,
        width = 5, height = 5, units = 'in', dpi = 300, )





rma.subset <- rma[paste0('X',samples$array_id[samples$chemical=='3-METHYLCHOLANTHRENE' | (samples$chemical=='Control' & samples$vehicle=='CMC_.5_%')])]
rma.subset <- new("MAList",list(M=rma.subset,A=rma.subset))

# construct a design matrix modeling treatment vs control for use by limma
design <- model.matrix( ~factor( samples$chemical[samples$vehicle=='CMC_.5_%'], levels=c('Control','3-METHYLCHOLANTHRENE') ) )
colnames(design) <- c('Intercept','3-METHYLCHOLANTHRENE')

# run limma
fit <- lmFit(rma.subset, design)
fit <- eBayes(fit)
t <- topTable(fit, coef=2, n=nrow(rma.subset), adjust='BH')

# write out the results to file
write.csv(t, paste0('3-METHYLCHOLANTHRENE', '_limma_results.csv') )

# sort and filter
t_filtered <- subset(t, adj.P.Val<0.05)
t_sorted_filtered <- t_filtered[order(t_filtered$adj.P.Val),]

# get top 10
top10 <- t_sorted_filtered[1:10,]
write.csv(top10, paste0('3-METHYLCHOLANTHRENE', '_topten.csv') )

# histogram
ggplot(t_sorted_filtered, aes(x=logFC)) + 
  geom_histogram(color="black", fill="cadetblue") +
  ggtitle('3-Methylcholanthrene') +
  ylab("Count of Significant DE Genes")

ggsave( paste0('3-METHYLCHOLANTHRENE', '_histogram.jpeg'), plot = last_plot(), device = jpeg,
        width = 5, height = 5, units = 'in', dpi = 300, )

# scatter plot
ggplot(t_sorted_filtered, aes(x=logFC, y=P.Value)) +
  geom_point(shape=1,size = 1, color = "cadetblue") +
  ggtitle('3-Methylcholanthrene') +
  ylab("Adj P-Value")

ggsave( paste0('3-METHYLCHOLANTHRENE', '_scatterplot.jpeg'), plot = last_plot(), device = jpeg,
        width = 5, height = 5, units = 'in', dpi = 300, )






# 6 ---------------------------------------------------------------------------


# Example differential expression results from both methods in 
# /project/bf528/project_3/results
# for you to use until the DESeq2 results come available for your samples.

# The input should be two sets of DE genes
# how to map Affymetrix probe IDs from the microarray analysis to refSeq identifiers used by the RNA-Seq analysis
# refSeq-to-probe id mapping in /project/bf528/project_3/refseq_affy_map


# list.files("/project/bf528/project_3/results")

# RNA Seq data
deseq_flu <- read.csv("../programmer/CAR_PXR_solo_deseq_results.csv") # CAR/PXR
deseq_pir <- read.csv("../programmer/PPARA_deseq_results.csv") # PPARA
deseq_met <- read.csv("../programmer/AhR_solo_deseq_results.csv") # AhR

deseq_flu <- subset(deseq_flu, pvalue < 0.05)
deseq_pir <- subset(deseq_pir, pvalue < 0.05)
deseq_met <- subset(deseq_met, pvalue < 0.05)

deseq_flu <- subset(deseq_flu, abs(log2FoldChange) >1.5)
deseq_pir <- subset(deseq_pir, abs(log2FoldChange) >1.5)
deseq_met <- subset(deseq_met, abs(log2FoldChange) >1.5)

# Microarray data
# limma_flu <- read.csv("/projectnb2/bf528/users/dreadlocks/project_3/analyst/FLUCONAZOLE_limma_results.csv")
limma_flu <- read.csv("/project/bf528/project_3/results/example_limma_results.csv")
limma_pir <- read.csv("/projectnb2/bf528/users/dreadlocks/project_3/analyst/PIRINIXIC_ACID_limma_results.csv")
limma_met <- read.csv("/projectnb2/bf528/users/dreadlocks/project_3/analyst/3-METHYLCHOLANTHRENE_limma_results.csv")

limma_flu <- subset(limma_flu, P.Value < 0.05)
limma_pir <- subset(limma_pir, P.Value < 0.05)
limma_met <- subset(limma_met, P.Value < 0.05)

limma_flu <- subset(limma_flu, abs(logFC) >1.5)
limma_pir <- subset(limma_pir, abs(logFC) >1.5)
limma_met <- subset(limma_met, abs(logFC) >1.5)

# Mapping
mapper <- read.csv("/project/bf528/project_3/refseq_affy_map.csv")


limma_flu_merge <- merge(limma_flu, mapper, by.x = "X", by.y = "PROBEID")
limma_flu_merge <- drop_na(limma_flu_merge)
deseq_flu_merge <- merge(deseq_flu, mapper, by.x = "X", by.y = "REFSEQ")
deseq_flu_merge <- drop_na(deseq_flu_merge)
both_merge <- merge(limma_flu_merge, deseq_flu_merge, by.x = "X", by.y = "PROBEID")
both_merge_directionality_agreement <- subset(both_merge, sign(both_merge$logFC) == sign(both_merge$log2FoldChange))
n_0 <- nrow(both_merge_directionality_agreement)
total <- sum(nrow(limma_flu_merge), nrow(deseq_flu_merge))
n_1 <- dim(deseq_flu)[1]
n_2 <- dim(limma_flu)[1]
N <- dim(mapper)[1]
n_x <- abs((n_0*N-n_1*n_2)/(n_0+N-n_1-n_2))
C_flu <- 2 * n_x / total


limma_pir_merge <- merge(limma_pir, mapper, by.x = "X", by.y = "PROBEID")
limma_pir_merge <- drop_na(limma_pir_merge)
deseq_pir_merge <- merge(deseq_pir, mapper, by.x = "X", by.y = "REFSEQ")
deseq_pir_merge <- drop_na(deseq_pir_merge)
both_merge <- merge(limma_pir_merge, deseq_pir_merge, by.x = "X", by.y = "PROBEID")
both_merge_directionality_agreement <- subset(both_merge, sign(both_merge$logFC) == sign(both_merge$log2FoldChange))
n_0 <- nrow(both_merge_directionality_agreement)
total <- sum(nrow(limma_pir_merge), nrow(deseq_pir_merge))
n_1 <- dim(deseq_pir)[1]
n_2 <- dim(limma_pir)[1]
N <- dim(mapper)[1]
n_x <- abs((n_0*N-n_1*n_2)/(n_0+N-n_1-n_2))
C_pir <- 2 * n_x / total

limma_met_merge <- merge(limma_met, mapper, by.x = "X", by.y = "PROBEID")
limma_met_merge <- drop_na(limma_met_merge)
deseq_met_merge <- merge(deseq_met, mapper, by.x = "X", by.y = "REFSEQ")
deseq_met_merge <- drop_na(deseq_met_merge)
both_merge <- merge(limma_met_merge, deseq_met_merge, by.x = "X", by.y = "PROBEID")
both_merge_directionality_agreement <- subset(both_merge, sign(both_merge$logFC) == sign(both_merge$log2FoldChange))
n_0 <- nrow(both_merge_directionality_agreement)
total <- sum(nrow(limma_met_merge), nrow(deseq_met_merge))
n_1 <- dim(deseq_met)[1]
n_2 <- dim(limma_met)[1]
N <- dim(mapper)[1]
n_x <- abs((n_0*N-n_1*n_2)/(n_0+N-n_1-n_2))
C_met <- 2 * n_x / total


# concordance vs number of DE genes - microarray
limma_scatter <- data.frame("Concordance"= c(C_flu, C_pir, C_met), 
                            "Treatment"= c(nrow(limma_flu), nrow(limma_pir), nrow(limma_met)))
p1 <- ggplot(limma_scatter, aes(x=Treatment, y=Concordance)) +
  geom_point() +
  ggtitle("Microarray Analysis") +
  ylab("Concordance of DEGs") +
  xlab("Treatment Effect") +
  geom_text(label=c("CAR", "PIR", "MET"), nudge_x = c(20, 20, 20))

# concordance vs number of DE genes - RNA-seq
deseq_scatter <- data.frame("Concordance"= c(C_flu, C_pir, C_met), 
                            "Treatment" = c(nrow(deseq_flu), nrow(deseq_pir), nrow(deseq_met)))
p2 <- ggplot(deseq_scatter, aes(x=Treatment, y=Concordance)) +
  geom_point() +
  ggtitle("RNA-Seq Analysis") +
  ylab("Concordance of DEGs") +
  xlab("Treatment Effect") +
  geom_text(label=c("CAR", "PIR", "MET"), nudge_x = c(40, 40, 40)
            )

#arrange plots into one figure
ggarrange(p1, p2, labels = c('A', 'B'), ncol = 2, nrow = 1)

ggsave( paste0('microarray_v_rnaseq_concordance_v_tx_effect', '.jpeg'), plot = last_plot(), device = jpeg,
        width = 8, height = 4, units = 'in', dpi = 300, )




# SUBDIVIDE -----------------------------------------------------

limma_flu_above <- subset(limma_flu, AveExpr > median(limma_flu$AveExpr))
limma_flu_below <- subset(limma_flu, AveExpr < median(limma_flu$AveExpr))
deseq_flu_above <- subset(deseq_flu, baseMean > median(deseq_flu$baseMean)) 
deseq_flu_below <- subset(deseq_flu, baseMean < median(deseq_flu$baseMean))

limma_pir_above <- subset(limma_pir, AveExpr > median(limma_pir$AveExpr))
limma_pir_below <- subset(limma_pir, AveExpr < median(limma_pir$AveExpr))
deseq_pir_above <- subset(deseq_pir, baseMean > median(deseq_pir$baseMean)) 
deseq_pir_below <- subset(deseq_pir, baseMean < median(deseq_pir$baseMean))

limma_met_above <- subset(limma_met, AveExpr > median(limma_met$AveExpr))
limma_met_below <- subset(limma_met, AveExpr < median(limma_met$AveExpr))
deseq_met_above <- subset(deseq_met, baseMean > median(deseq_met$baseMean)) 
deseq_met_below <- subset(deseq_met, baseMean < median(deseq_met$baseMean))


# ABOVE ----
limma_flu_merge <- merge(limma_flu_above, mapper, by.x = "X", by.y = "PROBEID")
limma_flu_merge <- drop_na(limma_flu_merge)
deseq_flu_merge <- merge(deseq_flu_above, mapper, by.x = "X", by.y = "REFSEQ")
deseq_flu_merge <- drop_na(deseq_flu_merge)
both_merge <- merge(limma_flu_merge, deseq_flu_merge, by.x = "X", by.y = "PROBEID")
both_merge_directionality_agreement <- subset(both_merge, sign(both_merge$logFC) == sign(both_merge$log2FoldChange))
n_0 <- nrow(both_merge_directionality_agreement)
total <- sum(nrow(limma_flu_merge), nrow(deseq_flu_merge))
n_1 <- dim(deseq_flu_above)[1]
n_2 <- dim(limma_flu_above)[1]
N <- dim(mapper)[1]
n_x <- abs((n_0*N-n_1*n_2)/(n_0+N-n_1-n_2))
C_flu_above <- 2 * n_x / total

limma_pir_merge <- merge(limma_pir_above, mapper, by.x = "X", by.y = "PROBEID")
limma_pir_merge <- drop_na(limma_pir_merge)
deseq_pir_merge <- merge(deseq_pir_above, mapper, by.x = "X", by.y = "REFSEQ")
deseq_pir_merge <- drop_na(deseq_pir_merge)
both_merge <- merge(limma_pir_merge, deseq_pir_merge, by.x = "X", by.y = "PROBEID")
both_merge_directionality_agreement <- subset(both_merge, sign(both_merge$logFC) == sign(both_merge$log2FoldChange))
n_0 <- nrow(both_merge_directionality_agreement)
total <- sum(nrow(limma_pir_merge), nrow(deseq_pir_merge))
n_1 <- dim(deseq_pir_above)[1]
n_2 <- dim(limma_pir_above)[1]
N <- dim(mapper)[1]
n_x <- abs((n_0*N-n_1*n_2)/(n_0+N-n_1-n_2))
C_pir_above <- 2 * n_x / total

limma_met_merge <- merge(limma_met_above, mapper, by.x = "X", by.y = "PROBEID")
limma_met_merge <- drop_na(limma_met_merge)
deseq_met_merge <- merge(deseq_met_above, mapper, by.x = "X", by.y = "REFSEQ")
deseq_met_merge <- drop_na(deseq_met_merge)
both_merge <- merge(limma_met_merge, deseq_met_merge, by.x = "X", by.y = "PROBEID")
both_merge_directionality_agreement <- subset(both_merge, sign(both_merge$logFC) == sign(both_merge$log2FoldChange))
n_0 <- nrow(both_merge_directionality_agreement)
total <- sum(nrow(limma_met_merge), nrow(deseq_met_merge))
n_1 <- dim(deseq_met_above)[1]
n_2 <- dim(limma_met_above)[1]
N <- dim(mapper)[1]
n_x <- abs((n_0*N-n_1*n_2)/(n_0+N-n_1-n_2))
C_met_above <- 2 * n_x / total

# BELOW ----
limma_flu_merge <- merge(limma_flu_below, mapper, by.x = "X", by.y = "PROBEID")
limma_flu_merge <- drop_na(limma_flu_merge)
deseq_flu_merge <- merge(deseq_flu_below, mapper, by.x = "X", by.y = "REFSEQ")
deseq_flu_merge <- drop_na(deseq_flu_merge)
both_merge <- merge(limma_flu_merge, deseq_flu_merge, by.x = "X", by.y = "PROBEID")
both_merge_directionality_agreement <- subset(both_merge, sign(both_merge$logFC) == sign(both_merge$log2FoldChange))
n_0 <- nrow(both_merge_directionality_agreement)
total <- sum(nrow(limma_flu_merge), nrow(deseq_flu_merge))
n_1 <- dim(deseq_flu_below)[1]
n_2 <- dim(limma_flu_below)[1]
N <- dim(mapper)[1]
n_x <- abs((n_0*N-n_1*n_2)/(n_0+N-n_1-n_2))
C_flu_below <- 2 * n_x / total

limma_pir_merge <- merge(limma_pir_below, mapper, by.x = "X", by.y = "PROBEID")
limma_pir_merge <- drop_na(limma_pir_merge)
deseq_pir_merge <- merge(deseq_pir_below, mapper, by.x = "X", by.y = "REFSEQ")
deseq_pir_merge <- drop_na(deseq_pir_merge)
both_merge <- merge(limma_pir_merge, deseq_pir_merge, by.x = "X", by.y = "PROBEID")
both_merge_directionality_agreement <- subset(both_merge, sign(both_merge$logFC) == sign(both_merge$log2FoldChange))
n_0 <- nrow(both_merge_directionality_agreement)
total <- sum(nrow(limma_pir_merge), nrow(deseq_pir_merge))
n_1 <- dim(deseq_pir_below)[1]
n_2 <- dim(limma_pir_below)[1]
N <- dim(mapper)[1]
n_x <- abs((n_0*N-n_1*n_2)/(n_0+N-n_1-n_2))
C_pir_below <- 2 * n_x / total

limma_met_merge <- merge(limma_met_below, mapper, by.x = "X", by.y = "PROBEID")
limma_met_merge <- drop_na(limma_met_merge)
deseq_met_merge <- merge(deseq_met_below, mapper, by.x = "X", by.y = "REFSEQ")
deseq_met_merge <- drop_na(deseq_met_merge)
both_merge <- merge(limma_met_merge, deseq_met_merge, by.x = "X", by.y = "PROBEID")
both_merge_directionality_agreement <- subset(both_merge, sign(both_merge$logFC) == sign(both_merge$log2FoldChange))
n_0 <- nrow(both_merge_directionality_agreement)
total <- sum(nrow(limma_met_merge), nrow(deseq_met_merge))
n_1 <- dim(deseq_met_below)[1]
n_2 <- dim(limma_met_below)[1]
N <- dim(mapper)[1]
n_x <- abs((n_0*N-n_1*n_2)/(n_0+N-n_1-n_2))
C_met_below <- 2 * n_x / total



#bar plot combining overall concordance measures 
df <- data.frame("Concordance"=c(C_flu, C_flu_above, C_flu_below, C_pir, C_pir_above, C_pir_below, C_met, C_met_above, C_met_below),
                       "Analysis"=rep(c("Overall", "Above", "Below"), 3),
                       "Chemical"= c(rep("Fluconazole", 3),
                                    rep("Pirinixic Acid", 3),
                                    rep("3-Methylcholanthrene", 3)))
write.csv(df, paste0('concordance', '.csv') )

# Plot data with a side-by-side bar chart
ggplot(df, aes(fill=Analysis, x=Chemical, y=Concordance)) + 
  geom_bar(stat="identity", position="dodge") +
  labs(title="Overall, Above-, and Below-Median Concordance")

ggsave( paste0('overall_above_below_concordance', '.jpeg'), plot = last_plot(), device = jpeg,
        width = 8, height = 4, units = 'in', dpi = 300, )
