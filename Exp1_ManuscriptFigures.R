##Figures for Bombus impatiens Maternal Brain Volume
##MA Hagadorn
##Last Modified 6/16/2022
##Most recent iteration: creating figures

#loadlibraries
library(plyr)
library(dplyr)
library(Rmisc)
library(ggplot2)
library(cowplot)
library(reshape2)
library(lemon)

#When working in RNA rProject
#setwd("../../../../../Manuscripts/Bee/Bimp_MaternalBrainVolumetric/Bimp_MomBrainVolume/")
MBdata<- read.csv2("RCODE/MBvolumetricData_analysisready.csv", sep=",", stringsAsFactors = FALSE)

#assign numeric
#assign data columns as numeric
numb_cols <- c(12:20)
MBdata[,numb_cols] <- apply(MBdata[ , numb_cols], 2, function(x) as.numeric(x))

#assign factors
MBdata$tx <- factor(MBdata$tx, levels=c("foundresslikeW", "foundresslikeQ", "maternallikeW", "maternallikeQ", "queenlikeW", "queenlikeQ"))
#MBdata$general_tx <- factor(MBdata$general_tx, levels=c("RS", "SR", "NM"))
MBdata$general_tx <- factor(MBdata$general_tx, levels=c("foundress", "maternalbrain", "queenlike"))
MBdata$caste <- factor(MBdata$caste, levels = c("worker", "queen"))
MBdata$colony <- factor(MBdata$colony)
print(str(MBdata))


##3####################################
##BOX PLOT BY GENERAL TX: Individual##
######################################


legendsize = 9.5
lettersize = 3.5
x_location <- c(.67, .91, 1.17, 1.68, 1.91, 2.17)
legend.position = c(.94,.98)
dotsize <- 2.5
seed <- 8675
generalsize <- 17
axes <- 12
fig.letterlocation <- (-0.33)

#colors
colors_caste <- c("gray60", "white")
colors_dots <- c("steelblue","darkgoldenrod1","steelblue","darkgoldenrod1","steelblue", "darkgoldenrod1")
Label <- c("Worker", "Queen")
Label_tx <- c("foundresslikeW", "foundresslikeQ", "maternallikeW", "maternallikeQ", "queenlikeW", "queenlikeQ")

#labels
labels_generaltx <- c("foundress-like\n(el-/bc-)", "maternal-like\n(el+/bc+)","queen-like\n(el+/bc-)")

#LIPS
#Changing to be W then Q: sig.code_lip <- c("ab","ab","ab","ab","b","a")
sig.code_lip <- c("ab","ab","ab","ab","a","b")
height_lip <- c(0.06736246+.0047, 0.06856684  + 0.003, 0.06823144 +.0044, 0.07036425 +.0047,  0.06398787 + 0.0028, 0.07474857 +0.0021)
x_location_lip <- c(.70, 1.08, 1.70, 2.08, 2.67, 3.05)

seed <- 8675
#Boxplot


lip_generaltx <- ggplot(MBdata, aes(x=general_tx, y=lip_CB, fill=caste)) +
  geom_boxplot(outlier.shape = NA, color="black", lwd=0.3) +
  scale_fill_manual(values = colors_caste, labels = Label) +
  #points
  #geom_point(color="black", shape=1, position=position_jitterdodge(jitter.width=.5), size=dotsize+.01) +
  geom_point(aes(color = factor(tx)), pch=16, position=position_jitterdodge(jitter.width=.95, seed=seed), size=dotsize) +
  scale_color_manual(values=c("foundresslikeQ"="darkgoldenrod1", "foundresslikeW"="steelblue", "maternallikeQ"="darkgoldenrod1", "maternallikeW"="steelblue", "queenlikeQ"="darkgoldenrod1", "queenlikeW"="steelblue"), labels=Label_tx) +
  ggtitle("A") +
  annotate("text", x=x_location_lip, label= sig.code_lip, y = height_lip, vjust = -0.5, size=lettersize) +
  scale_y_continuous(name = "Lip Volume", limits = c(0.05,.08)) +
  scale_x_discrete(name = element_blank(), labels = labels_generaltx, expand = c(0.15,0.15)) +
  theme(text = element_text(color="black",size = 15),
        axis.title = element_text(color="black"),
        axis.text.x = element_text(color="black", size = axes, margin = margin(t = 5, r =0, b = 5, l = 0)),
        axis.text.y = element_text(color="black", size = axes, margin = margin(t = 0, r =5, b = 0, l = 16)),
        axis.ticks.length = unit (.35,"cm"),
        axis.ticks = element_line(color = "black"),
        legend.position = "NONE",
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1, linetype = "solid",colour = "black"),
        plot.title = element_text(size = generalsize, margin = margin(t = 6, r = 0, b = 0, l = 0), face="bold", vjust = 3, hjust = fig.letterlocation))


ggsave("Lips.tiff",
       plot = lip_generaltx,
       device = "tiff",
       path = "RCODE/Figures/manuscript/",
       width = 5.5,
       height = 5,
       units = "in",
       dpi = 1000,
       limitsize = TRUE)



#COLLAR
sig.code_col <- c("a","cd","ac","bcd", "ab","d")
height_col <- c(0.07311177+.0024, 0.08430262 +0.0044, 0.07684370 +.0033, 0.08361708+.0013,  0.07585705 + 0.0016, 0.08956015+0.00178)
x_location_col <- c(.70, 1.05, 1.72, 2.07, 2.67, 3.08)

#Boxplot
col_generaltx <- ggplot(MBdata, aes(x=general_tx, y=col_CB, fill=caste)) +
  geom_boxplot(outlier.shape = NA, color="black", lwd=0.3) +
  scale_fill_manual(values = colors_caste, labels = Label) +
  #points
  #geom_point(color="black", shape=1, position=position_jitterdodge(jitter.width=.5), size=dotsize+.01) +
  geom_point(aes(color = factor(tx)), pch=16, position=position_jitterdodge(jitter.width=.95, seed=2022), size=dotsize) +
  scale_color_manual(values=c("foundresslikeQ"="darkgoldenrod1", "foundresslikeW"="steelblue", "maternallikeQ"="darkgoldenrod1", "maternallikeW"="steelblue", "queenlikeQ"="darkgoldenrod1", "queenlikeW"="steelblue"), labels=Label_tx) +
  ggtitle("B") +
  annotate("text", x=x_location_col, label= sig.code_col, y = height_col, vjust = -0.5, size=lettersize) +
  scale_y_continuous(name = "Collar Volume", limits = c(0.06,.10)) +
  scale_x_discrete(name = element_blank(), labels = labels_generaltx, expand = c(0.15,0.15)) +
  theme(text = element_text(color="black",size = 15),
        axis.title = element_text(color="black"),
        axis.text.x = element_text(color="black", size = axes, margin = margin(t = 5, r =0, b = 5, l = 0)),
        axis.text.y = element_text(color="black", size = axes, margin = margin(t = 0, r =5, b = 0, l = 16)),
        axis.ticks.length = unit (.35,"cm"),
        axis.ticks = element_line(color = "black"),
        legend.position = "NONE",
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1, linetype = "solid",colour = "black"),
        plot.title = element_text(size = generalsize, margin = margin(t = 6, r = 0, b = 0, l = 0), face="bold", vjust = 3, hjust = fig.letterlocation))

ggsave("Collar.tiff",
       plot = col_generaltx,
       device = "tiff",
       path = "RCODE/Figures/manuscript/",
       width = 5.5,
       height = 5,
       units = "in",
       dpi = 1000,
       limitsize = TRUE)




#CALYCES
sig.code_calyces <- c("a","ab","a","ab", "a", "b")
height_calyces <- c(0.1404742 +.008, 0.1528695 +0.006, 0.1450751 +.009, 0.1539813 +.0055, 0.1398449 + 0.0055, 0.1643087 + 0.0035)
x_location_calyces <- c(.66, 1.05, 1.66, 2.05, 2.665, 3.04)
cal_yaxis <- c(0.11,0.18)

#Boxplot
calyces_generaltx <- ggplot(MBdata, aes(x=general_tx, y=calyces_CB, fill=caste)) +
  geom_boxplot(outlier.shape = NA, color="black", lwd=0.3) +
  scale_fill_manual(values = colors_caste, labels = Label) +
  #points
  #geom_point(color="black", shape=1, position=position_jitterdodge(jitter.width=.5), size=dotsize+.01) +
  geom_point(aes(color = factor(tx)), pch=16, position=position_jitterdodge(jitter.width=.95, seed=5364), size=dotsize) +
  scale_color_manual(values=c("foundresslikeQ"="darkgoldenrod1", "foundresslikeW"="steelblue", "maternallikeQ"="darkgoldenrod1", "maternallikeW"="steelblue", "queenlikeQ"="darkgoldenrod1", "queenlikeW"="steelblue"), labels=Label_tx) +
  ggtitle("C") +
  annotate("text", x=x_location_calyces, label= sig.code_calyces, y = height_calyces, vjust = -0.5, size=lettersize) +
  scale_y_continuous(name = "Calyx Volume", limits = cal_yaxis) +
  scale_x_discrete(name = element_blank(), labels = labels_generaltx, expand = c(0.15,0.15)) +
  theme(text = element_text(color="black",size = 15),
        axis.title = element_text(color="black"),
        axis.text.x = element_text(color="black", size = axes, margin = margin(t = 5, r =0, b = 5, l = 0)),
        axis.text.y = element_text(color="black", size = axes, margin = margin(t = 0, r =5, b = 0, l = 16)),
        axis.ticks.length = unit (.35,"cm"),
        axis.ticks = element_line(color = "black"),
        legend.position = "NONE",
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1, linetype = "solid",colour = "black"),
        plot.title = element_text(size = generalsize, margin = margin(t = 6, r = 0, b = 0, l = 0), face="bold", vjust = 3, hjust = fig.letterlocation))


ggsave("Calyces.tiff",
       plot = calyces_generaltx,
       device = "tiff",
       path = "RCODE/Figures/manuscript/",
       width = 5.5,
       height = 5,
       units = "in",
       dpi = 1000,
       limitsize = TRUE)





#MUSHROOM BODY LOBES
sig.code_mbls <- c("a","a","a", "a","a","a")
height_mbls <- c(0.11001868 +.0028, 0.11294111 +0.0052, 0.10114805 +.0143, 0.11209956 +.00305, 0.10977559 + 0.0019, 0.11294111 + 0.0018)
x_location_mbls <- c(.66, 1.04, 1.66, 2.0, 2.66, 3.04)
mblobes_yaxis <- c(0.09,0.13)

#Boxplot
mblobes_generaltx <- ggplot(MBdata, aes(x=general_tx, y=mblobes_CB, fill=caste)) +
  geom_boxplot(outlier.shape = NA, color="black", lwd=0.3) +
  scale_fill_manual(values = colors_caste, labels = Label) +
  #points
  #geom_point(color="black", shape=1, position=position_jitterdodge(jitter.width=.5), size=dotsize+.01) +
  geom_point(aes(color = factor(tx)), pch=16, position=position_jitterdodge(jitter.width=.95, seed=2318), size=dotsize) +
  scale_color_manual(values=c("foundresslikeQ"="darkgoldenrod1", "foundresslikeW"="steelblue", "maternallikeQ"="darkgoldenrod1", "maternallikeW"="steelblue", "queenlikeQ"="darkgoldenrod1", "queenlikeW"="steelblue"), labels=Label_tx) +
  ggtitle("D") +
  annotate("text", x=x_location_mbls, label= sig.code_mbls, y = height_mbls, vjust = -0.5, size=lettersize) +
  scale_y_continuous(name = "Mushroom Body Lobe Volume", limits = mblobes_yaxis) +
  scale_x_discrete(name = element_blank(), labels = labels_generaltx, expand = c(0.15,0.15)) +
  theme(text = element_text(color="black",size = 15),
        axis.title = element_text(color="black"),
        axis.text.x = element_text(color="black", size = axes, margin = margin(t = 5, r =0, b = 5, l = 0)),
        axis.text.y = element_text(color="black", size = axes, margin = margin(t = 0, r =5, b = 0, l = 16)),
        axis.ticks.length = unit (.35,"cm"),
        axis.ticks = element_line(color = "black"),
        legend.position = "NONE",
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1, linetype = "solid",colour = "black"),
        plot.title = element_text(size = generalsize, margin = margin(t = 6, r = 0, b = 0, l = 0), face="bold", vjust = 3, hjust = fig.letterlocation))


ggsave("MBlobes.tiff",
       plot = mblobes_generaltx,
       device = "tiff",
       path = "RCODE/Figures/manuscript/",
       width = 5.5,
       height = 5,
       units = "in",
       dpi = 1000,
       limitsize = TRUE)



fig3 <- plot_grid(lip_generaltx, NULL, col_generaltx, calyces_generaltx, NULL, mblobes_generaltx, 
                  rel_widths = c(.55,.06,.55), nrow = 2, ncol = 3)

ggsave("Figure3.tiff",
       plot = fig3,
       device = "tiff",
       path ="RCODE/Figures/manuscript/",
       scale = 1,
       width = 9,
       height = 9,
       units = "in",
       dpi = 1000,
       limitsize = TRUE,
       bg = "white")


ggsave("Figure3.png",
       plot = fig3,
       device = "png",
       path ="RCODE/Figures/manuscript/png",
       scale = 1,
       width = 9,
       height = 9,
       units = "in",
       dpi = 1000,
       limitsize = TRUE,
       bg = "white")



# ggsave("Figure2.pdf",
#        plot = a_plot,
#        device = "pdf",
#        path = "Figures/",
#        scale = 1,
#        width = 4,
#        height = 6,
#        units = "in",
#        dpi = 1000,
#        limitsize = TRUE)
















#KENYON CELLS
sig.code_kcs <- c("ab","ab","ab", "a","b","a")
height_kcs <- c(0.09983231 +.0042, 0.09890870 +0.0012, 0.09919472 +.0028, 0.09471199 +.006, 0.10333983 + 0.0019, 0.09366606 + 0.0029)
x_location_kcs <- c(.68, 1.06, 1.68, 2.04, 2.66, 3.04)
kcs_yaxis <- c(0.08,.115)

#Boxplot
kcs_generaltx <- ggplot(MBdata, aes(x=general_tx, y=kcs_CB, fill=caste)) +
  geom_boxplot(outlier.shape = NA, color="black", lwd=0.3) +
  scale_fill_manual(values = colors_caste, labels = Label) +
  #points
  #geom_point(color="black", shape=1, position=position_jitterdodge(jitter.width=.5), size=dotsize+.01) +
  geom_point(aes(color = factor(tx)), pch=16, position=position_jitterdodge(jitter.width=.95, seed=5364), size=dotsize) +
  scale_color_manual(values=c("foundresslikeQ"="darkgoldenrod1", "foundresslikeW"="steelblue", "maternallikeQ"="darkgoldenrod1", "maternallikeW"="steelblue", "queenlikeQ"="darkgoldenrod1", "queenlikeW"="steelblue"), labels=Label_tx) +
  ggtitle("A") +
  annotate("text", x=x_location_kcs, label= sig.code_kcs, y = height_kcs, vjust = -0.5, size=lettersize) +
  scale_y_continuous(name = "Kenyon Cell Volume", limits = kcs_yaxis) +
  scale_x_discrete(name = element_blank(), labels = labels_generaltx, expand = c(0.15,0.15)) +
  theme(text = element_text(color="black",size = 15),
        axis.title = element_text(color="black"),
        axis.text.x = element_text(color="black", size = axes, margin = margin(t = 5, r =0, b = 5, l = 0)),
        axis.text.y = element_text(color="black", size = axes, margin = margin(t = 0, r =5, b = 0, l = 16)),
        axis.ticks.length = unit (.35,"cm"),
        axis.ticks = element_line(color = "black"),
        legend.position = "NONE",
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1, linetype = "solid",colour = "black"),
        plot.title = element_text(size = generalsize, margin = margin(t = 6, r = 0, b = 0, l = 0), face="bold", vjust = 3, hjust = fig.letterlocation))


ggsave("KCs.tiff",
       plot = kcs_generaltx,
       device = "tiff",
       path = "RCODE/Figures/manuscript/",
       width = 5.5,
       height = 5,
       units = "in",
       dpi = 1000,
       limitsize = TRUE)




#Neuropil
sig.code_neuropil <- c("a", "bc","ab","bc", "a", "c")
height_neuropil <- c(0.2504929 +.004, 0.2675858 +0.0045, 0.2554510 +.008, 0.2660809 +.0055, 0.2496205 + 0.0025, 0.2772498 - 0.0005)
x_location_neuropil <- c(.67, 1.06, 1.69, 2.07, 2.67, 3.05)
neuropil_yaxis <- c(0.22, 0.30)

#Boxplot
neuropil_generaltx <- ggplot(MBdata, aes(x=general_tx, y=neuropil_CB, fill=caste)) +
  geom_boxplot(outlier.shape = NA, color="black", lwd=0.3) +
  scale_fill_manual(values = colors_caste, labels = Label) +
  #points
  #geom_point(color="black", shape=1, position=position_jitterdodge(jitter.width=.5), size=dotsize+.01) +
  geom_point(aes(color = factor(tx)), pch=16, position=position_jitterdodge(jitter.width=.95, seed=2022), size=dotsize) +
  scale_color_manual(values=c("foundresslikeQ"="darkgoldenrod1", "foundresslikeW"="steelblue", "maternallikeQ"="darkgoldenrod1", "maternallikeW"="steelblue", "queenlikeQ"="darkgoldenrod1", "queenlikeW"="steelblue"), labels=Label_tx) +
  ggtitle("B") +
  annotate("text", x=x_location_neuropil, label= sig.code_neuropil, y = height_neuropil, vjust = -0.5, size=lettersize) +
  scale_y_continuous(name = "Neuropil Volume", limits = neuropil_yaxis) +
  scale_x_discrete(name = element_blank(), labels = labels_generaltx, expand = c(0.15,0.15)) +
  theme(text = element_text(color="black",size = 15),
        axis.title = element_text(color="black"),
        axis.text.x = element_text(color="black", size = axes, margin = margin(t = 5, r =0, b = 5, l = 0)),
        axis.text.y = element_text(color="black", size = axes, margin = margin(t = 0, r =5, b = 0, l = 16)),
        axis.ticks.length = unit (.35,"cm"),
        axis.ticks = element_line(color = "black"),
        legend.position = "NONE",
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1, linetype = "solid",colour = "black"),
        plot.title = element_text(size = generalsize, margin = margin(t = 6, r = 0, b = 0, l = 0), face="bold", vjust = 3, hjust = fig.letterlocation))


ggsave("TotalNeuropil.tiff",
       plot = neuropil_generaltx,
       device = "tiff",
       path = "RCODE/Figures/manuscript/",
       width = 5.5,
       height = 5,
       units = "in",
       dpi = 1000,
       limitsize = TRUE)



#NK
sig.code_nk <- c("ab","bc","ac","cd","a","d")
height_nk <- c(2.515303 +.075, 2.707984 +0.1, 2.582801 +.115, 2.819165 +.135, 2.419149 + 0.02, 2.965873 + 0.038)
x_location_nk <- c(.69, 1.065, 1.69, 2.07, 2.67, 3.05)
nk_yaxis <- c(2.0, 3.5)

#Boxplot
nk_generaltx <- ggplot(MBdata, aes(x=general_tx, y=nk, fill=caste)) +
  geom_boxplot(outlier.shape = NA, color="black", lwd=0.3) +
  scale_fill_manual(values = colors_caste, labels = Label) +
  #points
  #geom_point(color="black", shape=1, position=position_jitterdodge(jitter.width=.5), size=dotsize+.01) +
  geom_point(aes(color = factor(tx)), pch=16, position=position_jitterdodge(jitter.width=.95, seed=2022), size=dotsize) +
  scale_color_manual(values=c("foundresslikeQ"="darkgoldenrod1", "foundresslikeW"="steelblue", "maternallikeQ"="darkgoldenrod1", "maternallikeW"="steelblue", "queenlikeQ"="darkgoldenrod1", "queenlikeW"="steelblue"), labels=Label_tx) +
  ggtitle("C") +
  annotate("text", x=x_location_nk, label= sig.code_nk, y = height_nk, vjust = -0.5, size=lettersize) +
  scale_y_continuous(name = "Neuropil:Kenyon Cell Volume", limits = nk_yaxis) +
  scale_x_discrete(name = element_blank(), labels = labels_generaltx, expand = c(0.15,0.15)) +
  theme(text = element_text(color="black",size = 15),
        axis.title = element_text(color="black"),
        axis.text.x = element_text(color="black", size = axes, margin = margin(t = 5, r =0, b = 5, l = 0)),
        axis.text.y = element_text(color="black", size = axes, margin = margin(t = 0, r =5, b = 0, l = 16)),
        axis.ticks.length = unit (.35,"cm"),
        axis.ticks = element_line(color = "black"),
        legend.position = "NONE",
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1, linetype = "solid",colour = "black"),
        plot.title = element_text(size = generalsize, margin = margin(t = 6, r = 0, b = 0, l = 0), face="bold", vjust = 3, hjust = fig.letterlocation+0.04))


ggsave("NK.tiff",
       plot = nk_generaltx,
       device = "tiff",
       path = "RCODE/Figures/manuscript/",
       width = 5.5,
       height = 5,
       units = "in",
       dpi = 1000,
       limitsize = TRUE)



#ANTENNAL LOBES
sig.code_al <- c("a","ab","a","ab","a","b")
height_al <- c(0.05337774 +.0025, 0.05689339 +0.0013, 0.05261469 +.0006, 0.05677000 +.0016, 0.05419108 + 0.0015, 0.06007301 + 0.00135)
x_location_al <- c(.67, 1.05, 1.67, 2.05, 2.67, 3.05)
al_yaxis <- c(0.04,0.07)

#Boxplot
al_generaltx <- ggplot(MBdata, aes(x=general_tx, y=als_CB, fill=caste)) +
  geom_boxplot(outlier.shape = NA, color="black", lwd=0.3) +
  scale_fill_manual(values = colors_caste, labels = Label) +
  #points
  #geom_point(color="black", shape=1, position=position_jitterdodge(jitter.width=.5), size=dotsize+.01) +
  geom_point(aes(color = factor(tx)), pch=16, position=position_jitterdodge(jitter.width=.95, seed=6921), size=dotsize) +
  scale_color_manual(values=c("foundresslikeQ"="darkgoldenrod1", "foundresslikeW"="steelblue", "maternallikeQ"="darkgoldenrod1", "maternallikeW"="steelblue", "queenlikeQ"="darkgoldenrod1", "queenlikeW"="steelblue"), labels=Label_tx) +
  ggtitle("D") +
  annotate("text", x=x_location_al, label= sig.code_al, y = height_al, vjust = -0.5, size=lettersize) +
  scale_y_continuous(name = "Antennal Lobe Volume", limits = al_yaxis) +
  scale_x_discrete(name = element_blank(), labels = labels_generaltx, expand = c(0.15,0.15)) +
  theme(text = element_text(color="black",size = 15),
        axis.title = element_text(color="black"),
        axis.text.x = element_text(color="black", size = axes, margin = margin(t = 5, r =0, b = 5, l = 0)),
        axis.text.y = element_text(color="black", size = axes, margin = margin(t = 0, r =5, b = 0, l = 16)),
        axis.ticks.length = unit (.35,"cm"),
        axis.ticks = element_line(color = "black"),
        legend.position = "NONE",
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 1, linetype = "solid",colour = "black"),
        plot.title = element_text(size = generalsize, margin = margin(t = 6, r = 0, b = 0, l = 0), face="bold", vjust = 3, hjust = fig.letterlocation))


ggsave("AntennalLobes.tiff",
       plot = al_generaltx,
       device = "tiff",
       path = "RCODE/Figures/manuscript/",
       width = 5.5,
       height = 5,
       units = "in",
       dpi = 1000,
       limitsize = TRUE)




# # extract legend from AL plot
# legend <- get_legend(
#   lip_generaltx +
#     guides(color = guide_legend(nrow = 1)) +
#     theme(legend.position = "bottom",
#           legend.title = element_blank(),
#           legend.direction="horizontal",
#           legend.spacing.y = unit(0, "pt"),
#           legend.spacing.x = unit(3, "pt"),
#           legend.background = element_rect(fill="white"),
#           legend.key = element_rect(fill="white"),
#           legend.text = element_text(size=25, vjust = 1),
#           legend.key.height = unit(20, "pt"),
#           legend.key.width = unit(1.25, "cm"),)
# )
# 



fig4 <- plot_grid(kcs_generaltx, NULL, neuropil_generaltx, nk_generaltx, NULL, al_generaltx, rel_widths = c(.55,.06,.55), nrow = 2, ncol = 3)


ggsave("Figure4.tiff",
       plot = fig4,
       device = "tiff",
       path ="RCODE/Figures/manuscript/",
       scale = 1,
       width = 9,
       height = 9,
       units = "in",
       dpi = 1000,
       limitsize = TRUE,
       bg="white")


ggsave("Figure4.png",
       plot = fig4,
       device = "png",
       path ="RCODE/Figures/manuscript/png",
       scale = 1,
       width = 9,
       height = 9,
       units = "in",
       dpi = 1000,
       limitsize = TRUE,
       bg="white")



