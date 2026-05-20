###The following script will generate figures for the manuscript
###Supplemental figures are in another document


##Load in the filtered data
#load in the same count matrix used for WGCNA
rawCounts0 <- read.csv2("Filtering/GECountData_GenesOutliers_Filtered.csv", sep=";")
head(rawCounts0)
dim(rawCounts0)

rawCounts <- rawCounts0
rownames(rawCounts) <- rawCounts[,1]
rawCounts <- rawCounts[,-1]




#load in the sample mappings
sampleData0 <- read.delim("Filtering/sampleData_OutliersRemoved.csv", sep=";")
head(sampleData0)
sampleData <- sampleData0
sampleData <- sampleData[, c(1:2,5,8)]
sampleData$tx_group <- factor(sampleData$tx_group, levels=c("FLW", "FLQ", "MLW", "MLQ", "QLW", "QLQ"))




#Get PCA Ready
library(tidyverse)
library(ggfortify)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpattern)
library(grid)



##Create a PCA matrix
pca_matrix <- t(as.matrix(rawCounts))
#Look at the first five genes and samples
pca_matrix[1:5, 1:5]

pca_matrix0 <- pca_matrix

# Convert matrix to tibble for future use
as_tibble(pca_matrix, rownames = "sample")


# Perform the PCA
pca <- prcomp(pca_matrix)
pc_eigenvalues <- pca$sdev^2


#to view using ggplot2 create a tibble
# create a "tibble" manually with 
# a variable indicating the PC number
# and a variable with the variances
pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_eigenvalues

##tx_category x caste
autoplot(pca, data = sampleData, colour = "caste", shape = "tx_category", size=3.5) +
  scale_shape_manual(values=c(16,17,8), name="Treatment") +
  scale_colour_manual(values=c("queen" = "goldenrod1", "worker" = "steelblue"), name="Caste")+
  #scale_fill_manual(values=c("queen" = "goldenrod1", "worker" = "steelblue"), name="Caste") +
  #ggtitle(" ") +
  theme(text = element_text(color="black", size = 12),
      axis.title = element_text(color="black", size=12), 
      axis.text.x = element_text(color="black", size = 10), 
      axis.text.y = element_text(color="black", size = 10),
      axis.title.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0)), # Increase space on the right of the y-axis title
      axis.title.x = element_text(margin = margin(t = 8, r = 0, b = 0, l = 0)), # Increase space on the top of the x-axis title
      axis.ticks.length = unit (.2,"cm"),
      axis.ticks = element_line(color="black", size = 0.3),
      legend.background = element_rect(fill="white"),
      legend.key = element_rect(fill="white"),
      legend.text = element_text(size=8),
      legend.key.height = unit(17, "pt"),
      legend.key.width = unit(.4, "cm"),
      panel.background = element_rect(fill = "white"),
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(size = 0.3, linetype = "solid", colour = "black"),
      plot.margin= grid::unit(c(t = .03, r = 0, b = .03, l = .03), "in"),
      plot.title = element_text(hjust = 0.5))


ggsave("PCA_casteXtx.pdf",
       plot = last_plot(),
       device = "pdf",
       path = "Figures/",
       scale = 1,
       width = 6,
       height = 4,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

ggsave("PCA_casteXtx.png",
       plot = last_plot(),
       device = "png",
       path = "Figures/",
       scale = 1,
       width = 6,
       height = 4,
       units = "in",
       dpi = 500,
       limitsize = TRUE)




#######Venn Diagrams
library(ggVennDiagram)
library(VennDiagram)
library(ggvenn)
library(cowplot)
library(ggplot2)



# Read in R data
#load("DEG/BIMB_MB_SigDEGs.RData")
#load("DEG/BIMB_MB_DEGs_pLFC.RData")
load("DEG/BIMB_MB_DEGdataGeneLists.RData")
#load("DEG/DEG_data")

#All Contrasts
DEGlists_ALL <- list(
  FLWvFLQ = DEG_genelist_FLWFLQ, 
  MLWvMLQ = DEG_genelist_MLWMLQ, 
  QLWvQLQ = DEG_genelist_QLWQLQ, 
  FLWvMLW = DEG_genelist_FLWMLW, 
  FLQvMLQ = DEG_genelist_FLQMLQ, 
  FLQvQLQ = DEG_genelist_FLQQLQ
)

DEGlists_acrosscaste <- list(
  FLWvFLQ = DEG_genelist_FLWFLQ, 
  MLWvMLQ = DEG_genelist_MLWMLQ, 
  QLWvQLQ = DEG_genelist_QLWQLQ
)

DEGlists_withincastes <- list(
  FLWvMLW = DEG_genelist_FLWMLW, 
  FLQvMLQ = DEG_genelist_FLQMLQ, 
  FLQvQLQ = DEG_genelist_FLQQLQ
)

#Assign colors
#colors_ALL <- c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "green", "purple")
colors_acrosscaste <- c("lightgreen","coral", "purple1")
colors_withincaste <- c("#0073C2FF","#EFC000FF", "#CD534CFF")



# Across Caste Contrasts
#Create the Venn Diagram
acrosscaste <- venn.diagram(
  x = DEGlists_acrosscaste,
  filename = NULL,
  output = TRUE,
  height = 5.5, 
  width = 5, 
  resolution = 800,
  units = "in",
  cex = 1.75,
  fontface = "bold",
  col = colors_acrosscaste,
  fill = colors_acrosscaste,
  alpha = 0.5,
  #main = "B",
  #main.cex = 1.5,
  cat.cex = 1.25,
  cat.fontface = "plain",
  margin = 0.035,
  cat.pos = c(-27, 27, 180),
  cat.dist = c(0.06, 0.06, 0.04),
  cat.fontfamily="sans"
)


# pdf(file="venn_acrosscaste.pdf")
# grid.draw(acrosscaste)
# dev.off()




ggsave("venn_acrosscaste.png",
       plot = acrosscaste,
       device = "png",
       path = "Figures/",
       scale = 1,
       width = 4,
       height = 4,
       units = "in",
       dpi = 500,
       limitsize = TRUE)

ggsave("venn_acrosscaste.pdf",
       plot = acrosscaste,
       device = "pdf",
       path = "Figures/",
       scale = 1,
       width = 4,
       height = 4,
       units = "in",
       dpi = 500,
       limitsize = TRUE)



#Within Caste Contrasts

#Create the Venn Diagram
withincaste <- venn.diagram(
  x = DEGlists_withincastes,
  filename = NULL,
  output = TRUE,
  height = 5.5, 
  width = 5, 
  resolution = 800,
  units = "in",
  cex = 1.75,
  fontface = "bold",
  col = colors_withincaste,
  fill = colors_withincaste,
  alpha = 0.5,
  #main = "A",
  #main.cex = 1.5,
  cat.cex = 1.25,
  cat.fontface = "plain",
  margin = 0.035,
  cat.pos = c(-20, 27, 177),
  cat.dist = c(0.01, 0.05, 0.03),
  cat.fontfamily="sans"
)


ggsave("venn_withincaste.png",
       plot = withincaste,
       device = "png",
       path = "Figures/",
       scale = 1,
       width = 4,
       height = 4,
       units = "in",
       dpi = 500,
       limitsize = TRUE)

ggsave("venn_withincaste.pdf",
       plot = withincaste,
       device = "pdf",
       path = "Figures/",
       scale = 1,
       width = 4,
       height = 4,
       units = "in",
       dpi = 500,
       limitsize = TRUE)





#multipanel figure
spacer <- ggplot() + theme_void()

combined_plot <- plot_grid(
  withincaste, spacer, acrosscaste, labels = c("A"," ","B"), ncol = 3, nrow=1, rel_widths = c(4, 0.2, 4))

ggsave("venn_multipanel.pdf",
       plot = combined_plot,
       device = "pdf",
       path = "Figures/",
       scale = 1,
       width = 8,
       height = 4,
       units = "in",
       dpi = 500,
       limitsize = TRUE)


ggsave("venn_multipanel.png",
       plot = combined_plot,
       device = "png",
       path = "Figures/",
       scale = 1,
       width = 8,
       height = 4,
       units = "in",
       dpi = 500,
       limitsize = TRUE)





#Spoke with Karen, and we are going to try a stacked bar plot to show degs.


#First create a data frame:
df <- data.frame(
  treatment = c("FL v ML", "FL v ML", "FL v ML", "FL v ML",
                "FL v QL", "FL v QL", "FL v QL", "FL v QL",
                "ML v QL", "ML v QL", "ML v QL", "ML v QL"),
  caste     = c("worker", "worker", "queen", "queen",
                "worker", "worker", "queen", "queen",
                "worker", "worker", "queen", "queen"),
  deg_criteria  = c("strict", "lax", "strict", "lax",
                    "strict", "lax", "strict", "lax",
                    "strict", "lax", "strict", "lax"),
  count     = c(7,7,61,1744,0,0,28,409,0,0,0,1)
)

df$caste <- factor(df$caste, levels = c("worker", "queen"))

library(ggplot2)
library(dplyr)
library(ggbreak)
library(scales)
library(grid)
library(gtable)

# =========================
# Colors
# =========================
worker_strict <- "#0073C2FF"
queen_strict  <- "#EFC000FF"

worker_lax <- alpha(worker_strict, 0.45)
queen_lax  <- alpha(queen_strict, 0.45)

pal <- c(
  strict_worker = worker_strict,
  lax_worker    = worker_lax,
  strict_queen  = queen_strict,
  lax_queen     = queen_lax
)

# =========================
# Data prep
# =========================
df2 <- df %>%
  mutate(fill_group = paste(deg_criteria, caste, sep = "_"))

# =========================
# Build plot (WITH legend)
# =========================
p <- ggplot(df2, aes(
  x = treatment,
  y = count,
  fill = fill_group,
  group = caste
)) +
  
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.9),
    color = "black"
  ) +
  
  scale_fill_manual(
    values = pal,
    breaks = c(
      "strict_worker", "lax_worker",
      "strict_queen",  "lax_queen"
    ),
    labels = c(
      ""," Worker",
      ""," Queen"
    ),
    name = "Caste"
  ) +
  
  guides(
    fill = guide_legend(
      ncol = 2,
      byrow = TRUE
    )
  ) +
  
  scale_y_continuous(
    limits = c(-1, 1801),
    breaks = seq(0, 1800, by = 100),
    minor_breaks = seq(0, 1800, by = 50),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_break(
    breaks = c(450, 1700),
    space  = 0,
    expand = expansion(mult = c(0, 0))
  ) +
  guides(y = guide_axis(minor.ticks = TRUE)) +
  labs(x = "Treatment", y = "DEG Count", title = NULL) +
  theme(
    text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 12),
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title.x = element_text(hjust = 0.54),
    legend.position = "right",
    axis.text.y.right  = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.title.y.right = element_blank(),
    axis.ticks.x = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    #legend.key = element_rect(colour = NA, fill = "white"),
    legend.key.width  = unit(12, "pt"),
    legend.spacing.x  = unit(0, "pt"),
    legend.key.spacing.x = unit(0, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.background = element_rect(fill = "transparent", colour = NA), 
    legend.box.background = element_rect(fill = "transparent", colour = NA),
    legend.text = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) # Margin around key labels
  )



# =========================
# Extract legend
# =========================
g <- ggplotGrob(p)
legend <- gtable_filter(g, "guide-box")

# =========================
# Plot WITHOUT legend
# =========================
p_nolegend <- p + theme(legend.position = "none")

# Open device
png("Figures/stackedDEGcount.png",
    width = 4, height = 8, units = "in", res = 500)  # adjust size/res as needed

grid.newpage()
print(p_nolegend)   # plot no legend

# =========================
# Overlay legend manually
# =========================

grid.draw(
  editGrob(
    legend,
    vp = viewport(
      x = 0.75, y = 0.87,   # adjust position of the legend overlay
      width  = grobWidth(legend),
      height = grobHeight(legend),
      just = c("center", "center")
    )
  )
)


dev.off()  # save and close file
#all good except need to modify the legend to transparent







#Modified for NonStacked
df <- data.frame(
  treatment = c("FL v ML", "FL v ML", "FL v ML", "FL v ML",
                "FL v QL", "FL v QL", "FL v QL", "FL v QL",
                "ML v QL", "ML v QL", "ML v QL", "ML v QL"),
  caste     = c("worker", "worker", "queen", "queen",
                "worker", "worker", "queen", "queen",
                "worker", "worker", "queen", "queen"),
  deg_criteria  = c("strict", "lax", "strict", "lax",
                    "strict", "lax", "strict", "lax",
                    "strict", "lax", "strict", "lax"),
  count     = c(7,7,61,1744,0,0,28,409,0,0,0,1)
)

df$caste <- factor(df$caste, levels = c("worker", "queen"))


StrictOnly <- subset(df, deg_criteria=="strict")

library(ggplot2)
library(dplyr)
library(ggbreak)
library(scales)
library(grid)
library(gtable)

# =========================
# Colors
# =========================
worker <- "#0073C2FF"
queen  <- "#EFC000FF"


# =========================
# Build plot (WITH legend)
# =========================
DEGplot <- ggplot(StrictOnly, aes(
  x = treatment,
  y = count,
  fill = caste
)) +
  
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.9),
    color = "black"
  ) +
  
  scale_fill_manual(
    values = c(worker, queen),
    breaks = c(
      "worker", "queen"
    ),
    labels = c(" Worker"," Queen"),
    name = "Caste") +
  
  guides(
    fill = guide_legend(
      ncol = 1,
      byrow = TRUE
    )
  ) +
  
  scale_y_continuous(
    limits = c(0, 71),
    breaks = seq(0, 70, by = 10),
    minor_breaks = seq(0, 70, by = 5),
    expand = expansion(mult = c(0, 0))
  ) +
  guides(y = guide_axis(minor.ticks = TRUE)) +
  labs(x = "Treatment", y = "Number of DEGs", title = NULL) +
  theme(
    text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 12),
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title.x = element_text(hjust = 0.54),
    legend.position = "right",
    axis.text.y.right  = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.title.y.right = element_blank(),
    axis.ticks.x = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    #legend.key = element_rect(colour = NA, fill = "white"),
    legend.key.width  = unit(12, "pt"),
    legend.spacing.x  = unit(0, "pt"),
    legend.key.spacing.x = unit(0, "pt"),
    legend.key.spacing.y = unit(5, "pt"),
    legend.background = element_rect(fill = "transparent", colour = NA), 
    legend.box.background = element_rect(fill = "transparent", colour = NA),
    legend.text = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) # Margin around key labels
  )

ggsave("../Figures/DEGCountplot_nostack.png", plot=DEGplot, width = 4, height = 8, units = "in")
