##MAH Creating CPM Graphs for Manuscript

setwd("E:/DesktopDocuments/Documents/PHD/Studies/Bombus_impatiens/Bimpatiens_MaternalBrainGeneExpression/RNASeqAnalysis/BimpMB_RNAanalysis/GSEA")

#Use when not knitting whole docu
load("Results/GenesPerGO/SignificantGenesPerGO.Rdata")

#vector of go terms of interest
GOvector <- c("GO:0007631", "GO:0019563", "GO:0004371", "GO:0004970", "GO:0022008", "GO:0045202", "GO:0004890", 
              "GO:0007218", "GO:0062129", "GO:0035725", "GO:0019367", "GO:0034625", "GO:0034626", "GO:0042761", 
              "GO:0006680", "GO:0006629", "GO:0005615", "GO:0034464", "GO:0004181", "GO:0071683", "GO:0007601", 
              "GO:0007602","GO:0050911", "GO:0019230", "GO:0050954", "GO:0004984", "GO:0016499", "GO:0015179", 
              "GO:0005576", "GO:0004707", "GO:0010142")


#creatobj
GeneObj <- list(
  FLQQLQ=significant_genes_per_GO.FLQQLQ,
  FLQMLQ=significant_genes_per_GO.FLQMLQ,
  FLWFLQ=significant_genes_per_GO.FLWFLQ,
  FLWMLW=significant_genes_per_GO.FLWMLW,
  MLWMLQ=significant_genes_per_GO.MLWMLQ,
  QLWQLQ=significant_genes_per_GO.QLWQLQ
)




### Function to Isolate Genes for GOs of interest
extract_genes_by_go <- function(go_terms, objects) {
  if (!is.list(objects)) stop("`objects` must be a list of GO objects")
  
  result <- lapply(go_terms, function(go) {
    # collect genes from all objects and all slots
    all_genes <- unlist(lapply(names(objects), function(obj_name) {
      obj <- objects[[obj_name]]
      
      # loop through all slots in this object
      genes_in_obj <- unlist(lapply(names(obj), function(slot_name) {
        slot_obj <- obj[[slot_name]]
        if (go %in% names(slot_obj)) {
          return(slot_obj[[go]])
        } else {
          message("GO term ", go, " not found in ", obj_name, " (", slot_name, ")")
          return(character(0))
        }
      }))
      
      return(genes_in_obj)
    }))
    
    return(unique(all_genes))
  })
  
  # Clean names: remove ":" from GO IDs
  names(result) <- gsub(":", "", go_terms)
  return(result)
}


genes_byGO <- extract_genes_by_go(go_terms = GOvector, objects = GeneObj)









#load in the same count matrix used for WGCNA
rawCounts0 <- read.csv2("../../../RNASeqAnalysis/BimpMB_RNAanalysis/Filtering/GECountData_GenesOutliers_Filtered.csv", sep=";")
head(rawCounts0)
dim(rawCounts0)


#Read in with samples as column names and gene as column 1
#load in the sample mappings
sampleData0 <- read.delim("../../../RNASeqAnalysis/BimpMB_RNAanalysis/Filtering/sampleData_OutliersRemoved.csv", sep=";")
head(sampleData0)

#modifyformat
rawCounts <- rawCounts0
rownames(rawCounts) <- rawCounts$X
rawCounts <- rawCounts[,-1]

sampleData <- sampleData0
rownames(sampleData) <- sampleData$sample_id
sampleData <- sampleData[,1:2]


library(edgeR) #for cpm
library(dplyr)
library(tidyr)
library(tibble)
counts_pm <- cpm(rawCounts, log=TRUE)
# Make counts_pm a data frame and keep gene names
counts_df <- counts_pm %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene")






## Compare Counts (pm) of genes of interest


### GO:0007631 - FEEDING BEHAVIOR
FBGene <- genes_byGO$GO0007631

FBGeneCounts_lcpm <- counts_df %>%
  filter(gene %in% FBGene)

# Pivot longer so each row is gene × sample
FBGeneCounts_lcpm_long <- FBGeneCounts_lcpm %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "expression")

# Join metadata
FBGeneCounts_lcpm_long <- FBGeneCounts_lcpm_long %>%
  left_join(sampleData, by = c("sample_id" = "sample_id"))

FBGeneCounts_lcpm_long$tx_group <- factor(FBGeneCounts_lcpm_long$tx_group, levels= c("FLW", "FLQ", "MLW", "MLQ", "QLW", "QLQ"))

library(ggplot2)

#Make PLOT
FBfig <- ggplot(FBGeneCounts_lcpm_long, aes(x = tx_group, y = expression, fill = tx_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~gene, scales = "free_y", ncol=2) +  # one panel per gene
  theme_minimal() +
  labs(
    title = "Significant DEG(s) for GO:0007631",
    x = "Treatment Group",
    y = "Log-CPM Expression"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 10) # facet labels
  )

ggsave("results/GO_gene_Pics/GO0007631_FeedingBehavior.png", plot=FBfig, width = 4, height = 4, units = "in")










### GO:0007218 - NEUROPEPTIDE SIGNALLING
NPGene <- genes_byGO$GO0007218

NPGeneCounts_lcpm <- counts_df %>%
  filter(gene %in% NPGene)


# Pivot longer so each row is gene × sample
NPGeneCounts_lcpm_long <- NPGeneCounts_lcpm %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "expression")

# Join metadata
NPGeneCounts_lcpm_long <- NPGeneCounts_lcpm_long %>%
  left_join(sampleData, by = c("sample_id" = "sample_id"))

NPGeneCounts_lcpm_long$tx_group <- factor(NPGeneCounts_lcpm_long$tx_group, levels= c("FLW", "FLQ", "MLW", "MLQ", "QLW", "QLQ"))


#Make PLOT
NPfig <- ggplot(NPGeneCounts_lcpm_long, aes(x = tx_group, y = expression, fill = tx_group)) + #, label = sample_id)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  #geom_text(position = position_jitter(seed = 1)) +
  #facet_wrap(~Gene, scales = "free_y") + # One plot per gene
  theme_minimal() +
  labs(
    title = "Significant DEG(s) for GO:0007218",
    x = "Treatment Group",
    y = "Log-CPM Expression"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

ggsave("results/GO_gene_Pics/GO0007218_NeuropeptideSignaling_Fig.png", plot=NPfig, width = 4, height = 4, units = "in")




### GO:0022008 - NEUROGENESIS
NGGenes <- genes_byGO$GO0022008

NGGenesCounts_lcpm <- counts_df %>%
  filter(gene %in% NGGenes)


# Pivot longer so each row is gene × sample
NGGenesCounts_lcpm_long <- NGGenesCounts_lcpm %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "expression")

# Join metadata
NGGenesCounts_lcpm_long <- NGGenesCounts_lcpm_long %>%
  left_join(sampleData, by = c("sample_id" = "sample_id"))

NGGenesCounts_lcpm_long$tx_group <- factor(NGGenesCounts_lcpm_long$tx_group, levels= c("FLW", "FLQ", "MLW", "MLQ", "QLW", "QLQ"))


#Make PLOT
NGfig <- ggplot(NGGenesCounts_lcpm_long, aes(x = tx_group, y = expression, fill = tx_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~gene, scales = "free_y", ncol=2) +  # one panel per gene
  theme_minimal() +
  labs(
    title = "Significant DEG(s) for GO:0022008",
    x = "Treatment Group",
    y = "Log-CPM Expression"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 10) # facet labels
  )

ggsave("results/GO_gene_Pics/GO0022008_Neurogenesis_Fig.png", plot=NGfig, width = 9, height = 14, units = "in")




### GO:0045202 - SYNAPSE
SynGenes <- genes_byGO$GO0045202

#Group1
SynGenesCounts_lcpm <- counts_df %>%
  filter(gene %in% SynGenes)

#group1
SynGenesCounts_lcpm_Genes1to10 <- SynGenesCounts_lcpm[1:10,]

#group2
SynGenesCounts_lcpm_Genes11to12 <- SynGenesCounts_lcpm[11:12,]

# Pivot longer so each row is gene × sample
SynGenesCounts_lcpm_long_1to10 <- SynGenesCounts_lcpm_Genes1to10 %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "expression")

SynGenesCounts_lcpm_long_11to12 <- SynGenesCounts_lcpm_Genes11to12 %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "expression")

# Join metadata
SynGenesCounts_lcpm_long_1to10 <- SynGenesCounts_lcpm_long_1to10 %>%
  left_join(sampleData, by = c("sample_id" = "sample_id"))

SynGenesCounts_lcpm_long_11to12 <- SynGenesCounts_lcpm_long_11to12 %>%
  left_join(sampleData, by = c("sample_id" = "sample_id"))

SynGenesCounts_lcpm_long_1to10$tx_group <- factor(SynGenesCounts_lcpm_long_1to10$tx_group, levels= c("FLW", "FLQ", "MLW", "MLQ", "QLW", "QLQ"))

SynGenesCounts_lcpm_long_11to12$tx_group <- factor(SynGenesCounts_lcpm_long_11to12$tx_group, levels= c("FLW", "FLQ", "MLW", "MLQ", "QLW", "QLQ"))


#Make PLOT #1
Sny_fig1 <- ggplot(SynGenesCounts_lcpm_long_1to10, aes(x = tx_group, y = expression, fill = tx_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~gene, scales = "free_y", ncol=2) +  # one panel per gene
  theme_minimal() +
  labs(
    title = "Significant DEG(s) for GO:0045202",
    x = "Treatment Group",
    y = "Log-CPM Expression"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 10) # facet labels
  )

#Make PLOT #2
Sny_fig2 <- ggplot(SynGenesCounts_lcpm_long_11to12, aes(x = tx_group, y = expression, fill = tx_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~gene, scales = "free_y", ncol=2) +  # one panel per gene
  theme_minimal() +
  labs(
    title = "Continued: Significant DEG(s) GO:0045202",
    x = "Treatment Group",
    y = "Log-CPM Expression"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 10) # facet labels
  )


ggsave("results/GO_gene_Pics/GO0045202_Synapse_Fig.png", plot=Sny_fig1, width = 9, height = 14, units = "in")
ggsave("results/GO_gene_Pics/GO0045202_Synapse2_Fig.png", plot=Sny_fig2, width = 9, height = 4, units = "in")



### GO:0035725 - SODIUM ION TRANSMEMBRANE TRANSPORT
SITTGenes <- genes_byGO$GO0035725

SITTGenesCounts_lcpm <- counts_df %>%
  filter(gene %in% SITTGenes)




# Pivot longer so each row is gene × sample
SITTGenesCounts_lcpm_long <- SITTGenesCounts_lcpm %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "expression")

# Join metadata
SITTGenesCounts_lcpm_long <- SITTGenesCounts_lcpm_long %>%
  left_join(sampleData, by = c("sample_id" = "sample_id"))

SITTGenesCounts_lcpm_long$tx_group <- factor(SITTGenesCounts_lcpm_long$tx_group, levels= c("FLW", "FLQ", "MLW", "MLQ", "QLW", "QLQ"))


#Make PLOT
# Multipanel plot for all genes in this GO term
SITTfig <- ggplot(SITTGenesCounts_lcpm_long, aes(x = tx_group, y = expression, fill = tx_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~gene, scales = "free_y", ncol=2) +  # one panel per gene
  theme_minimal() +
  labs(
    title = "Significant DEG(s) for GO:0035725",
    x = "Treatment Group",
    y = "Log-CPM Expression"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 10) # facet labels
  )


ggsave("results/GO_gene_Pics/GO0035725_SodiumIonTransmembraneTransport.png", plot=SITTfig, width = 9, height = 8, units = "in")


### GO:0004970 - GLUTAMATE-GATED RECEPTOR ACTIVITY
GGRAGenes <- genes_byGO$GO0004970

GGRAGenesCounts_lcpm <- counts_df %>%
  filter(gene %in% GGRAGenes)


# Pivot longer so each row is gene × sample
GGRAGenesCounts_lcpm_long <- GGRAGenesCounts_lcpm %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "expression")

# Join metadata
GGRAGenesCounts_lcpm_long <- GGRAGenesCounts_lcpm_long %>%
  left_join(sampleData, by = c("sample_id" = "sample_id"))

GGRAGenesCounts_lcpm_long$tx_group <- factor(GGRAGenesCounts_lcpm_long$tx_group, levels= c("FLW", "FLQ", "MLW", "MLQ", "QLW", "QLQ"))


#Make PLOT
# Multipanel plot for all genes in this GO term
GGRAfig <- ggplot(GGRAGenesCounts_lcpm_long, aes(x = tx_group, y = expression, fill = tx_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~gene, scales = "free_y", ncol=2) +  # one panel per gene
  theme_minimal() +
  labs(
    title = "Significant DEG(s) for GO:0004970",
    x = "Treatment Group",
    y = "Log-CPM Expression"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 10) # facet labels
  )

ggsave("results/GO_gene_Pics/GO0004970_GlutamategatedReceptorActivity.png", plot=GGRAfig, width = 9, height = 10, units = "in")




### GO:0004890 - GABA-A RECEPTOR ACTIVITY
GARAGenes <- genes_byGO$GO0004890

GARAGenesCounts_lcpm <- counts_df %>%
  filter(gene %in% GARAGenes)

# Pivot longer so each row is gene × sample
GARAGenesCounts_lcpm_long <- GARAGenesCounts_lcpm %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "expression")

# Join metadata
GARAGenesCounts_lcpm_long <- GARAGenesCounts_lcpm_long %>%
  left_join(sampleData, by = c("sample_id" = "sample_id"))

GARAGenesCounts_lcpm_long$tx_group <- factor(GARAGenesCounts_lcpm_long$tx_group, levels= c("FLW", "FLQ", "MLW", "MLQ", "QLW", "QLQ"))


#Make PLOT
# Multipanel plot for all genes in this GO term
GARAfig <- ggplot(GARAGenesCounts_lcpm_long, aes(x = tx_group, y = expression, fill = tx_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~gene, scales = "free_y") +  # one panel per gene
  theme_minimal() +
  labs(
    title = "Significant DEG(s) for GO:0004890",
    x = "Treatment Group",
    y = "Log-CPM Expression"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 10) # facet labels
  )

ggsave("results/GO_gene_Pics/GO0004970_GabaReceptorActivity.png", plot=GARAfig, width = 9, height = 4, units = "in")





### GO:0004181 - METALLOCARBOXYPEPTIDASE ACTIVITY
MCPGenes <- genes_byGO$GO0004181

MCPGenesCounts_lcpm <- counts_df %>%
  filter(gene %in% MCPGenes)

# Pivot longer so each row is gene × sample
MCPGenesCounts_lcpm_long <- MCPGenesCounts_lcpm %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "expression")

# Join metadata
MCPGenesCounts_lcpm_long <- MCPGenesCounts_lcpm_long %>%
  left_join(sampleData, by = c("sample_id" = "sample_id"))

MCPGenesCounts_lcpm_long$tx_group <- factor(MCPGenesCounts_lcpm_long$tx_group, levels= c("FLW", "FLQ", "MLW", "MLQ", "QLW", "QLQ"))


#Make PLOT
# Multipanel plot for all genes in this GO term
MCPfig <- ggplot(MCPGenesCounts_lcpm_long, aes(x = tx_group, y = expression, fill = tx_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~gene, scales = "free_y") +  # one panel per gene
  theme_minimal() +
  labs(
    title = "Significant DEG(s) for GO:0004181",
    x = "Treatment Group",
    y = "Log-CPM Expression"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 10) # facet labels
  )

ggsave("results/GO_gene_Pics/GO0004181_MetallocarboxypeptidaseActivity.png", plot=MCPfig, width = 4, height = 4, units = "in")








### GO:0006629 - LIPID METABOLIC PROCESS

LMPGenes <- genes_byGO$GO0006629

LMPGenesCounts_lcpm <- counts_df %>%
  filter(gene %in% LMPGenes)


# Pivot longer so each row is gene × sample
LMPGenesCounts_lcpm_long <- LMPGenesCounts_lcpm %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "expression")

# Join metadata
LMPGenesCounts_lcpm_long <- LMPGenesCounts_lcpm_long %>%
  left_join(sampleData, by = c("sample_id" = "sample_id"))

LMPGenesCounts_lcpm_long$tx_group <- factor(LMPGenesCounts_lcpm_long$tx_group, levels= c("FLW", "FLQ", "MLW", "MLQ", "QLW", "QLQ"))


#Make PLOT
# Multipanel plot for all genes in this GO term
LMPfig <- ggplot(LMPGenesCounts_lcpm_long, aes(x = tx_group, y = expression, fill = tx_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~gene, scales = "free_y", ncol=2) +  # one panel per gene
  theme_minimal() +
  labs(
    title = "Significant DEG(s) for GO:0006629",
    x = "Treatment Group",
    y = "Log-CPM Expression"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 10) # facet labels
  )

ggsave("results/GO_gene_Pics/GO0006629_LipidMetabolicProcess.png", plot=LMPfig, width = 8, height = 8, units = "in")



### GO:0006680 - GLUCOSYLCERAMIDE CATABOLIC PROCESS
GluCPGenes <- genes_byGO$GO0006680

GluCPGenesCounts_lcpm <- counts_df %>%
  filter(gene %in% GluCPGenes)

#as.data.frame
#as.data.frame(FBGeneCounts_lcpm)
#colnames(FBGeneCounts_lcpm) <- FBGene[1]


# Pivot longer so each row is gene × sample
GluCPGenesCounts_lcpm_long <- GluCPGenesCounts_lcpm %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "expression")

# Join metadata
GluCPGenesCounts_lcpm_long <- GluCPGenesCounts_lcpm_long %>%
  left_join(sampleData, by = c("sample_id" = "sample_id"))

GluCPGenesCounts_lcpm_long$tx_group <- factor(GluCPGenesCounts_lcpm_long$tx_group, levels= c("FLW", "FLQ", "MLW", "MLQ", "QLW", "QLQ"))


#Make PLOT
# Multipanel plot for all genes in this GO term
GluCPfig <- ggplot(GluCPGenesCounts_lcpm_long, aes(x = tx_group, y = expression, fill = tx_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~gene, scales = "free_y") +  # one panel per gene
  theme_minimal() +
  labs(
    title = "Significant DEG(s) for GO:0006680",
    x = "Treatment Group",
    y = "Log-CPM Expression"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 10) # facet labels
  )


ggsave("results/GO_gene_Pics/GO0006680_GlucosylceramideCatabolicProcess.png", plot=GluCPfig, width = 4, height = 4, units = "in")




### GO:0019563 - GLYCEROL CATABOLIC PROCESS
GCPGenes <- genes_byGO$GO0019563

GCPGenesCounts_lcpm <- counts_df %>%
  filter(gene %in% GCPGenes)

# Pivot longer so each row is gene × sample
GCPGenesCounts_lcpm_long <- GCPGenesCounts_lcpm %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "expression")

# Join metadata
GCPGenesCounts_lcpm_long <- GCPGenesCounts_lcpm_long %>%
  left_join(sampleData, by = c("sample_id" = "sample_id"))

GCPGenesCounts_lcpm_long$tx_group <- factor(GCPGenesCounts_lcpm_long$tx_group, levels= c("FLW", "FLQ", "MLW", "MLQ", "QLW", "QLQ"))


#Make PLOT
# Multipanel plot for all genes in this GO term
GCPfig <- ggplot(GCPGenesCounts_lcpm_long, aes(x = tx_group, y = expression, fill = tx_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~gene, scales = "free_y") +  # one panel per gene
  theme_minimal() +
  labs(
    title = "Significant DEG(s) for GO:0019563",
    x = "Treatment Group",
    y = "Log-CPM Expression"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 10) # facet labels
  )

ggsave("results/GO_gene_Pics/GO0019563_GlycerolCatabolicProcess.png", plot=GCPfig, width = 4, height = 8, units = "in")











### GO:0019367, GO:0034625, GO:0034626, GO:0042761 - FATTY ACID LONGATION AND BIOSYNTHESIS (single gene: LOC100749525)
#There is a single significant gene across these 4 similar go terms. So, considering them all here.

FAEBGenes <- genes_byGO$GO0019367

FAEBGenesCounts_lcpm <- counts_df %>%
  filter(gene %in% FAEBGenes)

# Pivot longer so each row is gene × sample
FAEBGenesCounts_lcpm_long <- FAEBGenesCounts_lcpm %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "expression")

# Join metadata
FAEBGenesCounts_lcpm_long <- FAEBGenesCounts_lcpm_long %>%
  left_join(sampleData, by = c("sample_id" = "sample_id"))

FAEBGenesCounts_lcpm_long$tx_group <- factor(FAEBGenesCounts_lcpm_long$tx_group, levels= c("FLW", "FLQ", "MLW", "MLQ", "QLW", "QLQ"))


#Make PLOT
# Multipanel plot for all genes in this GO term
FAEBfig <- ggplot(FAEBGenesCounts_lcpm_long, aes(x = tx_group, y = expression, fill = tx_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~gene, scales = "free_y") +  # one panel per gene
  theme_minimal() +
  labs(
    title = "Significant DEG(s) for GO:0019367,\nGO:0034625, GO:0034626, GO:0042761",
    x = "Treatment Group",
    y = "Log-CPM Expression"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 10) # facet labels
  )


ggsave("results/GO_gene_Pics/MultipleGOs_FattyAcidElongationandBiosynthesis.png", plot=FAEBfig, width = 4, height = 4, units = "in")
















### GO:0005576 - EXTRACELLULAR REGION (cc)

ECRGenes <- genes_byGO$GO0005576

ECRGenesCounts_lcpm <- counts_df %>%
  filter(gene %in% ECRGenes)


# Pivot longer so each row is gene × sample
ECRGenesCounts_lcpm_long <- ECRGenesCounts_lcpm %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "expression")

# Join metadata
ECRGenesCounts_lcpm_long <- ECRGenesCounts_lcpm_long %>%
  left_join(sampleData, by = c("sample_id" = "sample_id"))

ECRGenesCounts_lcpm_long$tx_group <- factor(ECRGenesCounts_lcpm_long$tx_group, levels= c("FLW", "FLQ", "MLW", "MLQ", "QLW", "QLQ"))


#Make PLOT
# Multipanel plot for all genes in this GO term
ECRfig <- ggplot(ECRGenesCounts_lcpm_long, aes(x = tx_group, y = expression, fill = tx_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~gene, scales = "free_y", ncol=2) +  # one panel per gene
  theme_minimal() +
  labs(
    title = "Significant DEG(s) for GO:0005576",
    x = "Treatment Group",
    y = "Log-CPM Expression"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 10) # facet labels
  )

ggsave("results/GO_gene_Pics/GO0005576_ExtracellularSpace.png", plot=ECRfig, width = 8, height = 10, units = "in")









### GO:0010142 - farnesyl diphosphate biosynthetic process, mevalonate pathway (BP)

MevPathGenes <- genes_byGO$GO0010142

MevPathGenesCounts_lcpm <- counts_df %>%
  filter(gene %in% MevPathGenes)


# Pivot longer so each row is gene × sample
MevPathGenesCounts_lcpm_long <- MevPathGenesCounts_lcpm %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "expression")

# Join metadata
MevPathGenesCounts_lcpm_long <- MevPathGenesCounts_lcpm_long %>%
  left_join(sampleData, by = c("sample_id" = "sample_id"))

MevPathGenesCounts_lcpm_long$tx_group <- factor(MevPathGenesCounts_lcpm_long$tx_group, levels= c("FLW", "FLQ", "MLW", "MLQ", "QLW", "QLQ"))


#Make PLOT
# Multipanel plot for all genes in this GO term
MevPathfig <- ggplot(MevPathGenesCounts_lcpm_long, aes(x = tx_group, y = expression, fill = tx_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~gene, scales = "free_y", ncol=2) +  # one panel per gene
  theme_minimal() +
  labs(
    title = "Significant DEG(s) for GO:0010142",
    x = "Treatment Group",
    y = "Log-CPM Expression"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 10) # facet labels
  )

ggsave("results/GO_gene_Pics/GO0010142_MevalonatePathway.png", plot=MevPathfig, width = 4, height = 4, units = "in")


