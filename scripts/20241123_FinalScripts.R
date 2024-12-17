library(tidyverse)
library(devtools)
library(linkagemapping)
library(easysorter)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(broom)
library(ggplot2)
library(lme4)
library(multcompView)
library(ggpubr)


source("scripts/katiescripts.R")

##Figure 1
#Figure 1A
tbz_normnoregression <- rio::import("data/SFile 1A.tsv")
distributionnoregression <- tbz_normnoregression %>%
  dplyr::mutate(ben1 = case_when(cbn2 == "ben-1 variation" ~ "ben-1 variation",
                                 cbn2 == "No ben-1 varition" ~ "no variation"))%>%
  dplyr::mutate(cbn2 = case_when((cbn2 == "ben-1 variation" & cypvartype =="K267") ~ "ben-1",
                                 (cbn2 == "ben-1 variation" & cypvartype !="K267") ~ "ben-1+cyp-35D1",
                                 (cbn2 == "No ben-1 variation" & cypvartype =="K267") ~ "none",
                                 (cbn2 == "No ben-1 variation" & cypvartype !="K267") ~ "cyp-35D1"))%>%
  dplyr::mutate(mainfigure = case_when(strain =="N2" ~ "N2",
                                       strain == "CB4856" ~ "CB4856",
                                       TRUE ~ "Neither"))%>%
  ggplot()+
  aes(y = mean_norm)+
  geom_col(aes(x=reorder(strain,mean_norm),fill=mainfigure),width = 1)+
  scale_y_continuous("TBZ Resistance",expand = c(0, 0))+
  xlab("Strain")+
  scale_fill_manual(values = c("N2"="orange","CB4856"="blue"))+
  geom_point(data = dplyr::filter(tbz_normnoregression, cbn2 == "ben-1 variation") ,aes(x=strain, y= 0), shape = 24, color='red',fill="red")+
  theme_cowplot(12)+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())+
  theme(     text = element_text(size = 12),axis.text.x=element_blank(),
             legend.title = element_blank(),
             axis.ticks.x = element_blank(),
             axis.line.y = element_line(size=0.1),
             axis.line.x = element_line(size=0.1),
             legend.position = c(0.2,0.8))
##Figure 1B
total_independent_tests <- rio::import("data/total_independent_tests.txt")
independent_test_cutoff <- -log10(0.05/total_independent_tests[[1]])

# load loco mapping results
processed_mapping <- rio::import("data/SFile 1B.tsv")%>%
  dplyr::mutate(CHROM = factor(CHROM, levels = c("I","II","III","IV","V","X","MtDNA"))) %>%
  dplyr::select(-marker) %>%
  tidyr::unite("marker", CHROM, POS, sep = ":", remove = F) %>%
  dplyr::mutate(algorithm = "LOCO")

processed_mapping <- rbind(processed_mapping)

# each trait has a separate processed_mapping file now. So the plotting function and loop is removed
# but do check there is only 1 trait and if not, issue warning:
num_traits = length(unique(dplyr::select(processed_mapping,trait)))
if(num_traits > 1){
  print("WARNING: More than 1 trait in processed_mapping table. Only the first one will be plotted.")
}
# do we have mito mapping?
mito_check <- processed_mapping %>%
  na.omit()
## MANHATTAN PLOTS ##
for.plot <- processed_mapping %>%
  dplyr::mutate(CHROM = as.factor(CHROM)) %>%
  {
    if(!("MtDNA" %in% mito_check$CHROM)) dplyr::filter(., CHROM != "MtDNA") else .
  }
BF <- processed_mapping %>% 
  dplyr::group_by(trait, algorithm) %>% 
  dplyr::filter(log10p != 0) %>% 
  dplyr::distinct(marker, log10p) %>%
  dplyr::mutate(BF = -log10(0.05/sum(log10p > 0, na.rm = T))) %>%
  dplyr::ungroup() %>%
  dplyr::select(BF) %>%
  unique(.) %>%
  dplyr::slice(1) %>% # BF can be slightly different between loco and inbred... but just plot one (5.46 v 5.47...)
  as.numeric()
EIGEN <- independent_test_cutoff
BF.frame <- processed_mapping %>%
  dplyr::select(trait) %>%
  dplyr::filter(!duplicated(trait)) %>%
  dplyr::mutate(BF = BF, EIGEN  = EIGEN, user = unique(processed_mapping$BF)[1])
# if user selected a different threshold, use that, otherwise plot BF and EIGEN
if(BF.frame$user %in% c(BF.frame$BF, BF.frame$EIGEN)) {
  for.plot.ann <- for.plot %>%
    dplyr::mutate(sig = case_when(log10p > BF.frame$BF ~ "BF",
                                  log10p > BF.frame$EIGEN ~ "EIGEN",
                                  TRUE ~ "NONSIG"))
  
  sig.colors <- c("red","#EE4266", "black")
  names(sig.colors) <- c("BF","EIGEN", "NONSIG")
} else {
  for.plot.ann <- for.plot %>%
    dplyr::mutate(sig = case_when(log10p > BF.frame$user ~ "user",
                                  TRUE ~ "NONSIG"))
  
  sig.colors <- c("red", "black")
  names(sig.colors) <- c("user", "NONSIG")
}

test <- BF.frame %>%
  tidyr::pivot_longer(BF:user) %>%
  dplyr::distinct() %>%
  dplyr::filter(name %in% names(sig.colors))
# are we plotting mito or no?
if("MtDNA" %in% unique(for.plot.ann$CHROM)) {
  facet_scales <- "fixed"
} else {
  facet_scales <- "free"
}
man.plot <- ggplot() + 
  theme_bw() + 
  geom_point(data = for.plot.ann, 
             mapping = aes(x = POS/1000000, 
                           y = log10p,
                           colour = sig,
                           alpha = sig)) +
  scale_alpha_manual(values = c("BF" = 1, "EIGEN" = 1, "user" = 1, "NONSIG" = 0.25)) +
  scale_colour_manual(values = sig.colors) + 
  scale_x_continuous(expand = c(0, 0), breaks = c(5, 10, 15, 20)) +
  geom_hline(data = test, aes(yintercept = value, linetype = name)) + 
  scale_linetype_manual(values = c("BF" = 1, "EIGEN" = 3, "user" = 2)) +
  labs(x = "Genomic position (Mb)",
       y = expression(-log[10](italic(p)))) +
  facet_grid(trait ~ CHROM, scales = "free_x", space = facet_scales) +
  theme(panel.grid = element_blank(),
        axis.title.y = element_text(size=12),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12),
        legend.position = "none",
        # axis.title.x=element_blank(),
        #  axis.text.x=element_blank(),
        # axis.ticks.x=element_blank(),
        strip.text.y = element_blank()
  )

#Figure 1C
thialinkage <- rio::import("data/SFile 1C.tsv")
map1 <- thialinkage %>%
  dplyr::group_by(marker) %>%
  dplyr::filter(lod == max(lod))
tsize = 12
cis <- thialinkage %>%
  dplyr::group_by(marker) %>%
  dplyr::mutate(maxlod=max(lod))%>%
  dplyr::group_by(iteration) %>%
  dplyr::filter(!is.na(var_exp)) %>%
  dplyr::do(head(., n=1))

map1 <- linkagemapping:::cidefiner(cis, map1)
plot <- ggplot2::ggplot(map1) + 
  ggplot2::aes(x = pos/1e+06, y = lod) + 
  ggplot2::geom_ribbon(ggplot2::aes(x = pos/1e+06,  ymin = 0, ymax = ci_lod), fill = "blue", alpha = 0.5) + 
  ggplot2::geom_point(data = cis, ggplot2::aes(x = pos/1e+06, y = (1.05 * maxlod)), fill = "red", shape = 25, 
                      size = 2.2, show.legend = FALSE)

Figure_1C <- plot + ggplot2::geom_line(size = 1, alpha = 0.85) +
  ggplot2::facet_grid(.~chr, scales ="free", space = "free") +
  ggplot2::labs(x = "Genomic position (Mb)", y = "LOD") +
  ggplot2::scale_colour_discrete(name="Mapping\nIteration") +
  ggplot2::ggtitle(map1$trait[1]) +
  ylab("LOD")+ 
  ggplot2::theme_classic(12) +
  theme(strip.background = element_blank(),
        plot.title = element_blank(),
        axis.title.x = element_blank())+
  ggplot2::theme(text = element_text(size = 12),strip.background = element_blank(),
                 legend.position = "none",
                 panel.grid.minor = element_blank(),
                 panel.grid.major = element_blank(),
                 strip.text = element_blank(),
                 axis.title.y = element_text(size=12),
                 plot.margin = unit(c(0.2, 0, 0, 0), "in"),
                 axis.title.x = element_text(size = 12, face = "bold"),
                 axis.text.x = element_text(size =12),
                 axis.text.y = element_text(size=12))



figure_1 <- cowplot::plot_grid(distributionnoregression,man.plot,Figure_1C,nrow = 3, ncol = 1,align = "v",axis = "lr",labels = c("A","B","C"),label_size = 12)

ggsave("plots/Fig1.png",plot = figure_1,device = 'png',width = 7.5,height = 6.5,units = "in")


##Figure S1A
meantof <- rio::import("data/SFile S1A.tsv")

tbzdis_response <- meantof %>%
  dplyr::filter(!is.na(Thiabendazole_mean.TOF))

tbz_norm <- tbzdis_response %>%
  dplyr::mutate( mean_norm = (Thiabendazole_mean.TOF - min(Thiabendazole_mean.TOF)) / (max(Thiabendazole_mean.TOF) - min(Thiabendazole_mean.TOF)))%>%
  dplyr::mutate(cbn2 = case_when(strain == "N2" ~ "N2",
                                 strain == "CB4856" ~ "CB4856",
                                 TRUE ~ "neither"))%>%
  dplyr::mutate(ben1 = case_when(cypben1type %in% c("K267_b1","K267E_b1","K267D_b1","varnot267_b1") ~ "ben1",
                                 TRUE ~ "none"))

figure_S1A <- tbz_norm %>%
  ggplot()+
  aes(y = mean_norm)+
  geom_col(aes(x=reorder(strain,Thiabendazole_mean.TOF),fill=cbn2),width = 1)+
  theme_cowplot(12)+
  scale_y_continuous("TBZ Resistance",expand = c(0, 0))+
  xlab("Strain")+
  scale_fill_manual(values = c("neither"="grey","CB4856"="blue", "N2"="orange"))+
  geom_point(data=dplyr::filter(tbz_norm, ben1 == "ben1"),aes(x=strain, y=0), shape=24,fill="red",color="red")+
  theme_cowplot()+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())+
  theme(text = element_text(size = 12),axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(size=0.1),
        axis.line.x = element_line(size=0.1),
        legend.position = "none")

#Figure S1B
processed_mapping <- rio::import("data/SFile S1B.tsv")%>%
  dplyr::mutate(CHROM = factor(CHROM, levels = c("I","II","III","IV","V","X","MtDNA"))) %>%
  dplyr::select(-marker) %>%
  tidyr::unite("marker", CHROM, POS, sep = ":", remove = F) %>%
  dplyr::mutate(algorithm = "LOCO")

processed_mapping <- rbind(processed_mapping)

# each trait has a separate processed_mapping file now. So the plotting function and loop is removed
# but do check there is only 1 trait and if not, issue warning:
num_traits = length(unique(dplyr::select(processed_mapping,trait)))
if(num_traits > 1){
  print("WARNING: More than 1 trait in processed_mapping table. Only the first one will be plotted.")
}
# do we have mito mapping?
mito_check <- processed_mapping %>%
  na.omit()
## MANHATTAN PLOTS ##
for.plot <- processed_mapping %>%
  dplyr::mutate(CHROM = as.factor(CHROM)) %>%
  {
    if(!("MtDNA" %in% mito_check$CHROM)) dplyr::filter(., CHROM != "MtDNA") else .
  }
BF <- processed_mapping %>% 
  dplyr::group_by(trait, algorithm) %>% 
  dplyr::filter(log10p != 0) %>% 
  dplyr::distinct(marker, log10p) %>%
  dplyr::mutate(BF = -log10(0.05/sum(log10p > 0, na.rm = T))) %>%
  dplyr::ungroup() %>%
  dplyr::select(BF) %>%
  unique(.) %>%
  dplyr::slice(1) %>% # BF can be sligh
  as.numeric()
EIGEN <- independent_test_cutoff
BF.frame <- processed_mapping %>%
  dplyr::select(trait) %>%
  dplyr::filter(!duplicated(trait)) %>%
  dplyr::mutate(BF = BF, EIGEN  = EIGEN, user = unique(processed_mapping$BF)[1])
# if user selected a different threshold, use that, otherwise plot BF and EIGEN
if(BF.frame$user %in% c(BF.frame$BF, BF.frame$EIGEN)) {
  for.plot.ann <- for.plot %>%
    dplyr::mutate(sig = case_when(log10p > BF.frame$BF ~ "BF",
                                  log10p > BF.frame$EIGEN ~ "EIGEN",
                                  TRUE ~ "NONSIG"))
  
  sig.colors <- c("red","#EE4266", "black")
  names(sig.colors) <- c("BF","EIGEN", "NONSIG")
} else {
  for.plot.ann <- for.plot %>%
    dplyr::mutate(sig = case_when(log10p > BF.frame$user ~ "user",
                                  TRUE ~ "NONSIG"))
  
  sig.colors <- c("red", "black")
  names(sig.colors) <- c("user", "NONSIG")
}

test <- BF.frame %>%
  tidyr::pivot_longer(BF:user) %>%
  dplyr::distinct() %>%
  dplyr::filter(name %in% names(sig.colors))
# are we plotting mito or no?
if("MtDNA" %in% unique(for.plot.ann$CHROM)) {
  facet_scales <- "fixed"
} else {
  facet_scales <- "free"
}
man.plot1 <- ggplot() + 
  theme_bw() + 
  geom_point(data = for.plot.ann, 
             mapping = aes(x = POS/1000000, 
                           y = log10p,
                           colour = sig,
                           alpha = sig)) +
  scale_alpha_manual(values = c("BF" = 1, "EIGEN" = 1, "user" = 1, "NONSIG" = 0.25)) +
  scale_colour_manual(values = sig.colors) + 
  scale_x_continuous(expand = c(0, 0), breaks = c(5, 10, 15, 20)) +
  geom_hline(data = test, aes(yintercept = value, linetype = name)) + 
  scale_linetype_manual(values = c("BF" = 1, "EIGEN" = 3, "user" = 2)) +
  labs(x = "Genomic position (Mb)",
       y = expression(-log[10](italic(p)))) +
  facet_grid(trait ~ CHROM, scales = "free_x", space = facet_scales) +
  theme(panel.grid = element_blank(),
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        legend.position = "none",
        # axis.title.x=element_blank(),
        #  axis.text.x=element_blank(),
        # axis.ticks.x=element_blank(),
        strip.text.y = element_blank()
  )

S1 <- cowplot::plot_grid(figure_S1A, man.plot1, ncol = 1, align = "v",axis = "lr", labels = c("A","B"))
ggsave("plots/SFig1.png", plot = S1, width = 7.5, height = 6, units = "in")


##FigureS2A
thianonregressed <- rio::import("data/SFile S2.tsv")
linkagemapping::load_cross_obj("N2xCB4856cross_full")
thiaregressed <- rio::import("data/S4.tsv")

mergedpheno <- linkagemapping::mergepheno(N2xCB4856cross_full,thiaregressed, set = 2)

pxg_plot <- linkagemapping::pxgplot(mergedpheno, thialinkage)+
  ylab("Regressed animal length")+
  theme_cowplot(12)+
  theme(strip.background = element_blank(),
        axis.title.y = element_text(size = 12),
        plot.title = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

### Figure S2B
gwaspxg <- thianonregressed %>%
  dplyr::mutate(Allele = case_when(allele == -1 ~ 'REF',
                                   allele == 1 ~ 'ALT'))%>%
  dplyr::mutate(cbn2 = case_when(strain == 'N2' ~ "N2",
                                 strain == "CB4856" ~ "CB4856",
                                 TRUE ~ "neither"))
pxggwas <- gwaspxg %>%
  dplyr::filter(marker == "V_16077176")%>%
  dplyr::mutate(Allele = factor(Allele, levels = c("REF","ALT")))%>%
  ggplot()+
  aes(x=Allele, y=value)+
  geom_jitter(width=0.1,size = 0.5)+
  geom_boxplot(aes(alpha = 0.1),outlier.shape = NA)+
  theme_cowplot(12)+
  ylab("Regressed animal length")+
  theme(text = element_text(size = 12),legend.position = 'none',
        axis.title.x = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "in"),
        axis.title.y = element_text(size = 12))
pxgplots <- cowplot::plot_grid(pxggwas, pxg_plot, align = "vh", axis = "bt", labels = c("A","B"))
ggsave("plots/SFig2.png", plot = pxgplots, width = 6, height = 4, units = 'in')

##Data for Figure 2,S3,S7
cyp.metadata <- data.table::fread("data/cyp.doses.csv")
S9 <- data.table::fread("data/S9.tsv")
bg<- S9
bg1 <- data.table::fread("data/b10.txt")
genomeforallele <- bg


rd1<-rio::import(file="data/20230726_data1_processed.csv")
rd1_1 <- rd1 %>%
  dplyr::filter(strain%in% c("PHX2701","N2","PHX2883","CB4856"))
rd1_2 <- rd1 %>%
  dplyr::filter(strain%in% c("PHX2701","PHX2702","N2","PHX2883","PHX2882","CB4856"))


stats <- rd1 %>%
  dplyr::filter(concentration_um==32.5)%>%
  aov(median_wormlength_um_reg_delta ~ strain, data = .)%>%
  rstatix::tukey_hsd()

##Figure 2
#Figure 2A
genomeforallele <- bg %>%
  dplyr::mutate(sample = case_when(sample == "N2" ~ "N2",
                                   sample == "CB4856" ~ "CB4856",
                                   sample == "ECA239" ~ "PHX2883",
                                   sample == "ECA238" ~ "PHX2701"))
genomeforalleleplot <- ggplot(genomeforallele)+
  geom_segment(aes(x = start/1e6, y = factor(sample,levels = c("PHX2883","CB4856","PHX2701","N2")), xend = end/1e6, yend = sample, color = gt_name, size = 2))+
  scale_color_manual(values=c("N2"="orange","CB4856"="blue", "unknown" = "grey"))+
  # scale_color_manual(values=c("1"="#F0BA51","2"="#484DA0"))+
  facet_grid(~chrom, scales = "free",  space = "free")+
  theme_cowplot(12)+
  theme(plot.background = element_rect(fill="white"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=12, color = "black"),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks.x = element_blank())

#ggsave("test1.png", plot = genomeforalleleplot, width = 8, height = 4, units = "in")
#Figure 2B
cypforallele <- genomeforallele %>%
  dplyr::mutate(chrom = "cyp-35D1")%>%
  dplyr::mutate(gt_name =case_when(sample=="N2" ~ "N2",
                                   sample=="CB4856" ~"CB4856",
                                   sample=="PHX2701" ~"CB4856",
                                   sample=="PHX2883" ~ "N2"))
cypforalleleplot <- ggplot(cypforallele)+
  geom_segment(aes(x = start/1e6, y = factor(sample,levels = c("N2","PHX2701","CB4856","PHX2883")), xend = end/1e6, yend = sample, color = gt_name, size = 2))+
  scale_color_manual(values=c("N2"="orange","CB4856"="blue", "unknown" = "grey"))+
  # scale_color_manual(values=c("1"="#F0BA51","2"="#484DA0"))+
  facet_grid(~chrom, scales = "free",  space = "free")+
  theme_cowplot(12)+
  theme(plot.background = element_rect(fill="white"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face="italic",size=12, color = "black",hjust=0.7),
        plot.title = element_text(face="bold", hjust=0.5),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank())
#ggsave("test2.png", plot = cypforalleleplot, width = 8, height = 4, units = "in")


#Figure 2C
mainfig1 <- rd1_1
fullregressedallelesplot <- mainfig1 %>%
  dplyr::filter(concentration_um==32.5)%>%
  dplyr::mutate(strain=factor(strain, levels = c("PHX2883","CB4856","PHX2701","N2")))%>%
  ggplot+
  aes(x=strain, y=median_wormlength_um_reg_delta )+
  geom_jitter(width=0.1, size = 0.3)+
  geom_boxplot(aes(fill = strain),alpha = 0.8,outlier.shape = NA)+
  xlab("Regressed animal length")+
  ylab("Regressed animal length")+
  scale_fill_manual(values  = c("N2"="orange","PHX2701"="#E8E8E8","CB4856"="blue","PHX2883"="#696969"))+
  cowplot::theme_cowplot(12)+
  coord_flip()+
  ggpubr::geom_bracket(xmin = "N2", xmax = "PHX2701", y.position = -20,label = "****", coord.flip = TRUE)+  
  ggpubr::geom_bracket(xmin = "CB4856", xmax = "PHX2883", y.position = -20,label = "****", coord.flip = TRUE)+
  ggpubr::geom_bracket(xmin = "N2", xmax = "CB4856", y.position = 0,label = "**", coord.flip = TRUE)+
  theme(axis.text.x = element_text(size = 12),
        plot.background = element_rect(fill="white"),
        axis.title.x = element_text (size=12),
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text(face = "bold", color = "black"),
        plot.title = element_text(face="bold"),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

Assay1Fig <- cowplot::plot_grid(genomeforalleleplot,cypforalleleplot,fullregressedallelesplot,labels = c("A","","B"), nrow=1,label_size = 12,  rel_widths = c(0.35,0.2,1),align =  "h", axis = "bltr" )
ggsave("plots/Fig2.png", plot = Assay1Fig, width = 7.5, height = 4, units = "in")

#SFig3
rd12 <- rd1 %>%
  dplyr::filter(strain%in% c("ECA238","N2","ECA239","CB4856"))


#Figure A
genomeforallele <- bg %>%
  dplyr::mutate(sample = case_when(sample == "N2" ~ "N2",
                                   sample == "CB4856" ~ "CB4856",
                                   sample == "ECA239" ~ "ECA239",
                                   sample == "ECA238" ~ "ECA238"))
genomeforalleleplot <- ggplot(genomeforallele)+
  geom_segment(aes(x = start/1e6, y = factor(sample,levels = c("ECA239","CB4856","ECA238","N2")), xend = end/1e6, yend = sample, color = gt_name, size = 2))+
  scale_color_manual(values=c("N2"="orange","CB4856"="blue", "unknown" = "grey"))+
  # scale_color_manual(values=c("1"="#F0BA51","2"="#484DA0"))+
  facet_grid(~chrom, scales = "free",  space = "free")+
  theme_cowplot(12)+
  theme(plot.background = element_rect(fill="white"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text( color = "black",size="12"),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks.x = element_blank())

#ggsave("test1.png", plot = genomeforalleleplot, width = 8, height = 4, units = "in")
#Figure B
cypforallele <- genomeforallele %>%
  dplyr::mutate(chrom = "cyp-35d1")%>%
  dplyr::mutate(gt_name =case_when(sample=="N2" ~ "N2",
                                   sample=="CB4856" ~"CB4856",
                                   sample=="ECA238" ~"CB4856",
                                   sample=="ECA239" ~ "N2"))
cypforalleleplot <- ggplot(cypforallele)+
  geom_segment(aes(x = start/1e6, y = factor(sample,levels = c("N2","ECA238","CB4856","ECA239")), xend = end/1e6, yend = sample, color = gt_name, size = 2))+
  scale_color_manual(values=c("N2"="orange","CB4856"="blue", "unknown" = "grey"))+
  # scale_color_manual(values=c("1"="#F0BA51","2"="#484DA0"))+
  facet_grid(~chrom, scales = "free",  space = "free")+
  theme_cowplot(12)+
  theme(plot.background = element_rect(fill="white"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic", color = "black",hjust=0.7,size="12"),
        plot.title = element_text(face="bold", hjust=0.5),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank())
#ggsave("test2.png", plot = cypforalleleplot, width = 8, height = 4, units = "in")


#Figure C
mainfig1 <- rd12
fullregressedallelesplot <- mainfig1 %>%
  dplyr::filter(concentration_um==32.5)%>%
  dplyr::mutate(strain=factor(strain, levels = c("ECA239","CB4856","ECA238","N2")))%>%
  ggplot+
  aes(x=strain, y=median_wormlength_um_reg_delta )+
  geom_jitter(width=0.1, size = 0.3)+
  geom_boxplot(aes(fill = strain),alpha = 0.8,outlier.shape = NA)+
  xlab("Regressed animal length")+
  ylab("Regressed animal length")+
  scale_fill_manual(values  = c("N2"="orange","ECA238"="#E8E8E8","CB4856"="blue","ECA239"="#696969"))+
  cowplot::theme_cowplot(12)+
  coord_flip()+
  ggpubr::geom_bracket(xmin = "N2", xmax = "ECA238", y.position = -60,label = "****", coord.flip = TRUE)+  
  ggpubr::geom_bracket(xmin = "CB4856", xmax = "ECA239", y.position = -60,label = "****", coord.flip = TRUE)+
  theme(text=element_text(size=12),axis.text.x = element_text(size = 12),
        plot.background = element_rect(fill="white"),
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text(face = "bold", color = "black"),
        plot.title = element_text(face="bold"),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

Assay1NILFig <- cowplot::plot_grid(genomeforalleleplot,cypforalleleplot,fullregressedallelesplot,labels = c("A","","B"), nrow=1, rel_widths = c(0.35,0.2,1),align =  "h", axis = "bltr" )
ggsave("plots/Sfig3.png", plot = Assay1NILFig, width = 7.5, height = 4, units = "in")
###SFig4###
##Figure S4
tidy_genes_in_regionnob1 <- read.delim("data/SFile S4.tsv", stringsAsFactors=FALSE) %>% 
  dplyr::select(MARKER, POS, REF, ALT, WBGeneID, VARIANT_IMPACT, PEAK_MARKER, CHROM, STRAND, TRANSCRIPTION_START_POS, 
                TRANSCRIPTION_END_POS, QTL_INTERVAL_START, QTL_INTERVAL_END, VARIANT_LOG10p, STRAIN, STRAIN_GENOTYPE)
gene_infonob1 <- tidy_genes_in_regionnob1 %>%
  dplyr::filter(QTL_INTERVAL_START == 15776825 & !is.na(WBGeneID)) %>%
  dplyr::select(WBGeneID, CHROM, STRAND, TRANSCRIPTION_START_POS, TRANSCRIPTION_END_POS, QTL_INTERVAL_START, QTL_INTERVAL_END) %>% 
  dplyr::distinct()
# for each gene, keep the most significant overlapping variant log10p
# this is used as y-axix of gene position
gene_dfnob1 <- tidy_genes_in_regionnob1 %>%
  dplyr::filter(QTL_INTERVAL_START == 15776825 & !is.na(WBGeneID)) %>%
  dplyr::distinct(MARKER, WBGeneID, PEAK_MARKER, CHROM, STRAND, TRANSCRIPTION_START_POS, TRANSCRIPTION_END_POS, 
                  QTL_INTERVAL_START, QTL_INTERVAL_END, VARIANT_LOG10p) %>%
  dplyr::group_by(WBGeneID) %>% 
  dplyr::summarise(VARIANT_LOG10p = max(VARIANT_LOG10p, na.rm = T)) %>% 
  dplyr::inner_join(gene_infonob1) %>% 
  dplyr::mutate(STRAND = ifelse(STRAND=="+", "Gene_on_plus_strand", "Gene_on_minus_strand"))
peak_variantnob1 <- 16077176	

variant_dfnob1 <- tidy_genes_in_regionnob1 %>%
  dplyr::filter(QTL_INTERVAL_START == 15776825) %>%
  dplyr::distinct(CHROM, POS, VARIANT_LOG10p, VARIANT_IMPACT, REF, ALT) %>% 
  dplyr::mutate(MARKER_pos = paste0(CHROM, POS, " REF:ALT ", REF, ":", ALT))
variant_dfnob1$VARIANT_IMPACT[is.na(variant_dfnob1$VARIANT_IMPACT)] <- "Intergenic"
xsnob1 <- unique(gene_dfnob1$QTL_INTERVAL_START)
xenob1 <- unique(gene_dfnob1$QTL_INTERVAL_END)
gcnob1 <- unique(gene_dfnob1$CHROM)
variant_dfnob1high <- dplyr::filter(variant_dfnob1, VARIANT_IMPACT == "HIGH")
max_logpnob1 <- unique(max(variant_dfnob1high$VARIANT_LOG10p, na.rm = TRUE))/150
gene_plotnob1reg <- ggplot2::ggplot(gene_dfnob1) +
  ggplot2::geom_vline(aes(xintercept = peak_variantnob1/1e6),
                      linetype=3, color = "cyan")+
  ggplot2::geom_segment(aes(x = ifelse(STRAND == "Gene_on_plus_strand", TRANSCRIPTION_START_POS/1e6, TRANSCRIPTION_END_POS/1e6),
                            xend = ifelse(STRAND == "Gene_on_plus_strand", TRANSCRIPTION_END_POS/1e6, TRANSCRIPTION_START_POS/1e6),
                            y = VARIANT_LOG10p,
                            yend = VARIANT_LOG10p,
                            color = STRAND),
                        size = 1.5) +
  ggplot2::geom_segment(data = variant_dfnob1high,
                        aes(x = POS/1e6,
                            xend = POS/1e6,
                            y = VARIANT_LOG10p+max_logpnob1,
                            yend = VARIANT_LOG10p-max_logpnob1,
                            color = VARIANT_IMPACT)) +
  ggplot2::scale_color_manual(values = c("Gene_on_plus_strand" = "grey",
                                         "Gene_on_minus_strand" = "grey",
                                         "MODIFIER" = "gray50",
                                         "LOW" = "gray30",
                                         "MODERATE" = "orange",
                                         "HIGH" = "red",
                                         "Intergenic" = "tan",
                                         "Linker" = "tan"),
                              name = "Variant Impact")+
  ggplot2::labs(x = "Genomic Position (Mb)",
                y = expression(-log[10](italic(p)))) +
  ggplot2::geom_segment(x = 16.07, xend=16.07, y = 6.05, yend =5.8,arrow = arrow(length = unit(0.05, "in")))+
  ggplot2::theme_bw(12)+
  geom_text(x=16.07,y=6.2,label="cyp-35d1",fontface="italic")+
  ggplot2::theme(legend.position = "none",
                 panel.grid = element_blank())
ggsave("plots/SFig4.png", plot = pxgplots, width = 6, height = 4, units = 'in')

##SFig5
rd3<- rio::import("data/20230726_data2_processed.csv")

fig3<- rd3 %>%
  dplyr::filter(strain %in% c("N2","CB4856","ECA2859","ECA2921", "ECA2923","ECA2925"))

stats3 <- rd3 %>%
  dplyr::filter(concentration_um==32.5)%>%
  aov(median_wormlength_um_reg_delta ~ strain, data = .)%>%
  rstatix::tukey_hsd()

z1<- expression(atop(paste("WT"),paste("WT")))
z2<- expression(atop(paste(italic("cyp-35D1(ean221)")),paste("WT")))
z22<- expression(atop(paste(italic("cyp-35D1(ean222)")),paste("WT")))
z3<- expression(atop(paste("WT"),paste(italic("nhr-176(ean225)"))))
z33<- expression(atop(paste("WT"),paste(italic("nhr-176(ean226)"))))
z4<- expression(atop(paste(italic("cyp-35D1(K267E)")),paste("WT")))
z5<- expression(atop(paste(italic("cyp-35D1(ean223)")),paste("WT")))
z55<- expression(atop(paste(italic("cyp-35D1(ean224)")),paste("WT")))
z6<- expression(atop(paste(italic("cyp-35D1(K267E)")),paste(italic("nhr-176(ean227)"))))
z66<- expression(atop(paste(italic("cyp-35D1(K267E)")),paste(italic("nhr-176(ean228)"))))
#S5:
Fig5supp <- rd3 %>%
  dplyr::filter(concentration_um==32.5) %>%
  dplyr::mutate(strain=factor(strain, levels = c("N2","ECA2859","ECA2860","ECA2923","ECA2924","CB4856","ECA2921","ECA2922","ECA2925","ECA2926")))%>%
  ggplot+
  aes(x=strain, y=median_wormlength_um_reg_delta )+
  geom_jitter(width=0.1, size = 0.3)+
  geom_boxplot(aes(fill = strain),alpha = 0.8,outlier.shape = NA)+
  ylab("Regressed animal length")+
  scale_fill_manual(values  = c("N2"="orange","CB4856"="blue","ECA2859"="#E8E8E8","ECA2860"="#E8E8E8","ECA2923"="#E8E8E8","ECA2924"="#E8E8E8","ECA2921"="#696969","ECA2922"="#696969","ECA2925"="#696969","ECA2926"="#696969"))+
  scale_x_discrete(breaks=c("N2","ECA2859","ECA2860","ECA2923","ECA2924","CB4856","ECA2921","ECA2922","ECA2925","ECA2926"),labels=c(z1,z2,z2,z3,z3,z4,z5,z5,z6,z6))+
  cowplot::theme_cowplot(12)+
  ggpubr::geom_bracket(xmin = "N2", xmax = "ECA2859", y.position = -115,label = "****", coord.flip = FALSE)+  
  ggpubr::geom_bracket(xmin = "N2", xmax = "ECA2860", y.position = -85,label = "****", coord.flip = FALSE)+
  ggpubr::geom_bracket(xmin = "N2", xmax = "ECA2923", y.position = -55,label = "****", coord.flip = FALSE)+  
  ggpubr::geom_bracket(xmin = "N2", xmax = "ECA2924", y.position = -25,label = "****", coord.flip = FALSE)+
  ggpubr::geom_bracket(xmin = "CB4856", xmax = "ECA2921", y.position = -115,label = "****", coord.flip = FALSE)+
  ggpubr::geom_bracket(xmin = "CB4856", xmax = "ECA2922", y.position = -85,label = "****", coord.flip = FALSE)+
  ggpubr::geom_bracket(xmin = "CB4856", xmax = "ECA2925", y.position = -55,label = "****", coord.flip = FALSE)+
  ggpubr::geom_bracket(xmin = "CB4856", xmax = "ECA2926", y.position = -25,label = "****", coord.flip = FALSE)+
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 7.5),
        plot.background = element_rect(fill="white"),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", color = "black"),
        plot.title = element_text(face="bold"),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
blank_plot <- ggplot() + theme_void()
z0<-expression(atop(paste("CYP-35D1 Genotype:"),paste("NHR-176 Genotype:")))
Fig5supp<- cowplot::plot_grid(blank_plot,Fig5supp,labels = c("",""), nrow=1, rel_widths = c(0,1),align =  "h", axis = "bltr" )%>% 
  ggdraw() + 
  draw_label(z0, x = 0.05, y = 0.075, angle = 0,size=7.5)
Fig5supp
ggsave("plots/SFig5.png", plot = Fig5supp, width =10,height=3.5, units = "in",dpi=300)

######Double Deletion SFig6##
del<- rio::import(file="data/20230726_data4_processed.csv")
stats3 <- del %>%
  aov(median_wormlength_um_reg_delta ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::filter(group2=="N2")

z1<- expression(atop(paste("WT"),paste("WT")))
z2<- expression(atop(paste(italic("cyp-35D1(ean221)")),paste("WT")))
z3<- expression(atop(paste("WT"),paste(italic("nhr-176(ean226)"))))
z4<- expression(atop(paste(italic("cyp-35D1(ean291)")),paste(italic("nhr-176(ean226)"))))
z5<- expression(atop(paste(italic("cyp-35D1(ean292)")),paste(italic("nhr-176(ean226)"))))

#S6:
testsupp <- del %>%
  dplyr::filter(concentration_um == "32.5") %>% 
  dplyr::mutate(strain=factor(strain, levels = c("N2","ECA2859","ECA2924" ,"ECA3703", "ECA3704")))%>%
  ggplot+
  aes(x=strain, y=median_wormlength_um_reg_delta)+
  geom_jitter(width=0.1, size = 0.3)+
  geom_boxplot(aes(fill = strain),alpha = 0.8,outlier.shape = NA)+
  ylab("Regressed animal length")+
  scale_fill_manual(values  = c("N2"="orange","ECA2859"="#E8E8E8", "ECA2924"="#E8E8E8" ,"ECA3703"="#E8E8E8", "ECA3704"="#E8E8E8"))+
  scale_x_discrete(breaks=c("N2" ,"ECA2859", "ECA2924" ,"ECA3703", "ECA3704"),labels=c(z1,z2,z3,z4,z5))+
  cowplot::theme_cowplot(12)+
  ggpubr::geom_bracket(xmin = "N2", xmax = "ECA2859", y.position = -75,label = "****", coord.flip = FALSE)+
  ggpubr::geom_bracket(xmin = "N2", xmax = "ECA2924", y.position = -50,label = "****", coord.flip = FALSE)+
  ggpubr::geom_bracket(xmin = "N2", xmax = "ECA3703", y.position = -25,label = "****", coord.flip = FALSE)+
  ggpubr::geom_bracket(xmin = "N2", xmax = "ECA3704", y.position = 0,label = "****", coord.flip = FALSE)+
  ggpubr::geom_bracket(xmin = "ECA2859", xmax = "ECA3704", y.position = 25,label = "ns", coord.flip = FALSE)+
  theme(text= element_text(size=10),axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size=12),
        plot.background = element_rect(fill="white"),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", color = "black"),
        plot.margin = unit(c(0, 0, 0, 1), "cm"),
        plot.title = element_text(face="bold"),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
blank_plot <- ggplot() + theme_void()
z0<-expression(atop(paste("CYP-35D1 Genotype:"),paste("NHR-176 Genotype:")))
testsupp<- cowplot::plot_grid(blank_plot,testsupp,labels = c("",""), nrow=1, rel_widths = c(0,1),align =  "h", axis = "bltr" )%>% 
  ggdraw() + 
  draw_label(z0, x = 0.09, y = 0.058, angle = 0,size=10)
ggsave("plots/Sfig6.png", plot = testsupp, width =7.5, height = 4, units = "in")

#SFig7
bg2 <- data.table::fread("data/S10.txt")
#Figure s7A
genomeforallele <- bg2 %>%
  dplyr::mutate(sample = case_when(sample == "N2" ~ "N2",
                                   sample == "CB4856" ~ "CB4856",
                                   sample == "ECA239" ~ "PHX2883",
                                   sample == "ECA238" ~ "PHX2701", 
                                   sample == "PHX2882" ~ "PHX2882", 
                                   sample == "PHX2702" ~ "PHX2702",
  ))
genomeforalleleplot <- ggplot(genomeforallele)+
  geom_segment(aes(x = start/1e6, y = factor(sample,levels = c("PHX2883","PHX2882","CB4856","PHX2702","PHX2701","N2")), xend = end/1e6, yend = sample, color = gt_name, size = 2))+
  scale_color_manual(values=c("N2"="orange","CB4856"="blue", "unknown" = "grey"))+
  # scale_color_manual(values=c("1"="#F0BA51","2"="#484DA0"))+
  facet_grid(~chrom, scales = "free",  space = "free")+
  theme_cowplot(12)+
  theme(plot.background = element_rect(fill="white"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=12, color = "black"),
        plot.title = element_text(size=12),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks.x = element_blank())

#ggsave("test1.png", plot = genomeforalleleplot, width = 8, height = 4, units = "in")
#Figure s7B
cypforallele <- genomeforallele %>%
  dplyr::mutate(chrom = "cyp-35d1")%>%
  dplyr::mutate(gt_name =case_when(sample=="N2" ~ "N2",
                                   sample=="CB4856" ~"CB4856",
                                   sample=="PHX2701" ~"CB4856",
                                   sample=="PHX2883" ~ "N2",
                                   sample == "PHX2882" ~ "N2", 
                                   sample == "PHX2702" ~ "CB4856",))
cypforalleleplot <- ggplot(cypforallele)+
  geom_segment(aes(x = start/1e6, y = factor(sample,levels = c("PHX2883","PHX2882","CB4856","PHX2702","PHX2701","N2")), xend = end/1e6, yend = sample, color = gt_name, size = 2))+
  scale_color_manual(values=c("N2"="orange","CB4856"="blue", "unknown" = "grey"))+
  # scale_color_manual(values=c("1"="#F0BA51","2"="#484DA0"))+
  facet_grid(~chrom, scales = "free",  space = "free")+
  theme_cowplot(12)+
  theme(plot.background = element_rect(fill="white"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(face="italic",size=12, color = "black",hjust=0.7),
        plot.title = element_text(face="bold", hjust=0.5),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank())
#ggsave("test2.png", plot = cypforalleleplot, width = 8, height = 4, units = "in")


#Figure s7C
mainfig2 <- rd1_2
full2regressedallelesplot <- mainfig2 %>%
  dplyr::filter(concentration_um==32.5)%>%
  dplyr::mutate(strain=factor(strain, levels = c("PHX2883","PHX2882","CB4856","PHX2702","PHX2701","N2")))%>%
  ggplot+
  aes(x=strain, y=median_wormlength_um_reg_delta )+
  geom_jitter(width=0.1, size = 0.3)+
  geom_boxplot(aes(fill = strain),alpha = 0.8,outlier.shape = NA)+
  xlab("Regressed animal length")+
  ylab("Regressed animal length")+
  scale_fill_manual(values  = c("N2"="orange","PHX2701"="#E8E8E8","PHX2702"="#E8E8E8","CB4856"="blue","PHX2882"="#696969","PHX2883"="#696969"))+
  cowplot::theme_cowplot(12)+
  coord_flip()+
  ggpubr::geom_bracket(xmin = "N2", xmax = "PHX2701", y.position = -20,label = "****", coord.flip = TRUE)+ 
  ggpubr::geom_bracket(xmin = "N2", xmax = "PHX2702", y.position = 0,label = "****", coord.flip = TRUE)+
  ggpubr::geom_bracket(xmin = "CB4856", xmax = "PHX2882", y.position = -20,label = "*", coord.flip = TRUE)+
  ggpubr::geom_bracket(xmin = "CB4856", xmax = "PHX2883", y.position = 0,label = "****", coord.flip = TRUE)+
  theme(text=element_text(size=12),axis.text.x = element_text(size = 12),
        plot.background = element_rect(fill="white"),
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text(face = "bold", color = "black"),
        plot.title = element_text(face="bold"),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

Assay1SuppFig <- cowplot::plot_grid(genomeforalleleplot,cypforalleleplot,full2regressedallelesplot,labels = c("A","","B"),label_size = 12, nrow=1, rel_widths = c(0.35,0.2,1),align =  "h", axis = "bltr" )
ggsave("plots/SFig7.png", plot = Assay1SuppFig, width = 7.5, height = 4, units = "in")


##Fig3##
colscomp <- c("N2"="orange","CB4856"="blue","PHX2701"="orange","PHX2882"="blue")
ltype<-c("N2"="solid","CB4856"="solid","PHX2701"="longdash","PHX2882"="longdash")
''
`20230307_cypcompsummary_updated`<- rio::import("data/SFile 5A.csv") %>%
  dplyr::mutate(std=sdab/100) %>%
  dplyr::mutate(mean=meanab/100)

n2filtered<- `20230307_cypcompsummary_updated`%>%
  dplyr::filter(strain%in% c("N2","PHX2701"))

colscomp2 <- c("N2"="orange","PHX2701"="grey")


n2TBZ <- n2filtered%>%
  dplyr::filter(condition == "TBZ")%>%
  ggplot()+
  aes(x=generation,y=(meanab/100),colour = strain)+
  #geom_point(aes(color=strain))+
  ylab("Relative allele frequency")+
  xlab("Generation")+
  geom_errorbar(aes(ymin=ifelse((mean - std) < 0, 0, mean-std), ymax=(mean + std)), width=0.25,size=0.5)+
  geom_line(size=1)+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),limits = c(0,1))+
  scale_color_manual(name = "strain", labels = c("N2" = "K267","CB4856"="E267","PHX2701"="K267E","PHX2882"="E267K"), values = colscomp2)+
  scale_x_continuous(breaks = c(1,3,5,7))+
  cowplot::theme_cowplot(10)+
  theme(text= element_text(size=12),plot.background = element_rect(fill="white"),legend.position = "none")

n2DMSO <- n2filtered%>%
  dplyr::filter(condition == "DMSO")%>%
  ggplot()+
  aes(x=generation,y=(meanab/100),colour = strain)+
  #geom_point(aes(color=strain))+
  ylab("Relative allele frequency")+
  xlab("Generation")+
  geom_errorbar(aes(ymin=ifelse((mean - std) < 0, 0, mean-std), ymax=(mean + std)), width=0.25,size=0.5)+
  geom_line(size=1)+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1),limits = c(0,1))+
  scale_color_manual(name = "strain", labels = c("N2" = "K267","CB4856"="E267","PHX2701"="K267E","PHX2882"="E267K"), values = colscomp2)+
  scale_x_continuous(breaks = c(1,3,5,7))+
  cowplot::theme_cowplot(10)+
  theme(text= element_text(size=12),plot.background = element_rect(fill="white"),legend.position = "none")

#compplot2 <- cowplot::plot_grid(n2DMSO,n2TBZ,ncol = 2,nrow = 2,align = "hv", label_x = 0,axis = "lrbt",labels=c("A","B"),label_size = 12)

#ggsave(filename = "20230307_n2cypcomp.png",plot = compplot2 ,width = 5,height = 6)

##fitness##

filtfit <- rio::import("data/S File 13.csv") %>%
  dplyr::filter(fitness != "NaN")


filtd <- filtfit %>%
  dplyr::filter(condition == "DMSO")%>%
  aov(fitness ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::mutate(strain= group2)

filtt <- filtfit%>%
  dplyr::filter(condition == "TBZ")%>%
  aov(fitness ~ strain, data = .)%>%
  rstatix::tukey_hsd()%>%
  dplyr::mutate(strain = group2)

fitaplot <- filtfit %>%
  dplyr::filter(condition == "TBZ")%>%
  dplyr::mutate(strain = factor(strain, levels = c("N2","PHX2701")))%>%
  ggplot()+
  aes(x=strain, y=fitness,fill=strain)+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.15,size=.8)+
  geom_hline(yintercept = 0)+
  scale_y_continuous(breaks = c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2),limits = c(-1,0.2),sec.axis = dup_axis(name = "Thiabendazole"))+
  stat_pvalue_manual(filtt, label = "p.adj.signif",x.position="PHX2701", y.position = c(.2),remove.bracket = TRUE,size = 4)+
  scale_x_discrete(labels=c("N2" = "Lysine", "PHX2701" = "Glutamate"))+
  scale_fill_manual(name = "fancy_strain", labels = c("N2" = "Lysine","PHX2701"="Glutamate"), values = c("N2" = "orange","PHX2701"="grey"))+
  scale_color_manual(name = "Strain", labels = c("N2" = "Lysine","PHX2701"="Glutamate"), values = c("N2" = "orange","PHX2701"="grey"))+
  cowplot::theme_cowplot(10)+
  ylab("Competitive fitness")+
  theme(text= element_text(size=12),plot.background = element_rect(fill="white"),legend.position = "none",
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45,hjust = 1,vjust = 1,margin = unit(c(0,0,0,0),units = "in")),
        axis.title.x = element_blank())

fitdplot <- filtfit %>%
  dplyr::filter(condition == "DMSO")%>%
  dplyr::mutate(strain = factor(strain, levels = c("N2","PHX2701")))%>%
  ggplot()+
  aes(x=strain, y=fitness,fill=strain)+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.15,size=.8)+
  geom_hline(yintercept = 0)+
  scale_y_continuous(breaks = c(-0.8,-0.6,-0.4,-0.2,0,0.2,0.4),limits = c(-0.8,0.4),sec.axis = dup_axis(name = "DMSO"))+
  stat_pvalue_manual(filtd, label = "p.adj.signif",x.position="PHX2701", y.position = c(0.4),remove.bracket = TRUE,size = 4)+
  scale_x_discrete(labels=c("N2" = "Lysine", "PHX2701" = "Glutamate"))+
  scale_fill_manual(name = "fancy_strain", labels = c("N2" = "Lysine","PHX2701"="Glutamate"), values = c("N2" = "orange","PHX2701"="grey"))+
  scale_color_manual(name = "Strain", labels = c("N2" = "Lysine","PHX2701"="Glutamate"), values = c("N2" = "orange","PHX2701"="grey"))+
  cowplot::theme_cowplot(10)+
  ylab("Competitive fitness")+
  theme(text= element_text(size=12),plot.background = element_rect(fill="white"),legend.position = "none",
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 12, angle = 45,hjust = 1,vjust = 1,margin = unit(c(0,0,0,0),units = "in")),
        axis.title.x = element_blank())

compplot2 <- cowplot::plot_grid(n2DMSO,fitdplot,n2TBZ,fitaplot,ncol = 2,nrow = 2,align = "hv", label_x = -.03,axis = "lrbt",labels=c("A","B","C","D"),label_size = 12)

ggsave(filename = "plots/Fig3.png",plot = compplot2 ,width = 3.75,height = 5.2,units="in")

###Fig4,Sfig10,11
###########ENDOMETABOLOME####################
metab <- read.csv("data/endometab.csv")
# Convert "exposure" to numeric (if not already done)
metab$exposure <- as.numeric(metab$exposure)
metab$TBZ <- as.numeric(metab$TBZ)
metab$TBZ.OH <- as.numeric(metab$TBZ.OH)
metab$TBZ.O.glucoside <- as.numeric(metab$TBZ.O.glucoside)
metab$TBZ.O.phosphoglucoside <- as.numeric(metab$TBZ.O.phosphoglucoside)
metab$Genotype <- factor(metab$Genotype, levels = c("N2", "PHX2702", "CB4856", "PHX2883"))
metab <- metab%>%
  dplyr::filter(exposure!="0")
TBZmetab <- metab[, c(1, 2, 3)]
TBZmetab$condition <- "TBZ"
ohmetab <- metab[, c(1, 2, 4)]
ohmetab$condition <- "TBZ-OH"
ogmetab <-  metab[, c(1, 2, 5)]
ogmetab$condition <- "TBZ-O-Glu"
opmetab <-  metab[, c(1, 2, 6)]
opmetab$condition <- "TBZ-O-PGlu"
####STATS####
perform_stats_by_groups <- function(data, exposures, genotypes, columns_of_interest) {
  results <- data.frame()  # Initialize an empty dataframe to store results
  for (exposure in exposures) {
    for (genotype in genotypes) {
      # Subset the dataframe based on exposure and genotype
      subset_data <- subset(data, exposure == exposure & Genotype %in% genotype)
      # Perform statistical tests for each column of interest
      for (col in columns_of_interest) {
        result <- wilcox.test(subset_data[[col]] ~ subset_data$Genotype)
        p_value_corrected <- p.adjust(result$p.value, method = "bonferroni")
        # Determine significance level
        significance <- ifelse(p_value_corrected < 0.0001, "****",
                               ifelse(p_value_corrected < 0.001, "***",
                                      ifelse(p_value_corrected < 0.01, "**",
                                             ifelse(p_value_corrected < 0.05, "*", "ns"))))
        # Prepare a row of results
        result_row <- data.frame(
          Exposure = exposure,
          Genotype = paste(genotype, collapse = ", "),
          Column = col,
          W_statistic = result$statistic,
          P_value = result$p.value,
          P_value_corrected = p_value_corrected,
          Significance = significance
        )
        # Append the row to the results dataframe
        results <- rbind(results, result_row)
      }
    }
  }
  return(results)
}

exposures <- c("2", "6")  
genotypes <- list(c("N2", "PHX2702"), c("CB4856", "PHX2883")) 
columns_of_interest <- c("TBZ", "TBZ.OH", "TBZ.O.glucoside", "TBZ.O.phosphoglucoside")  

endo_results_df <- perform_stats_by_groups(data = metab, exposures = exposures, genotypes = genotypes, columns_of_interest = columns_of_interest)

print(endo_results_df)

###Fig4B###
cols<- c("N2"="orange","PHX2702"="#E8E8E8","CB4856"="blue","PHX2883"="#696969")
TBZmetab_filtered <- TBZmetab%>%
  dplyr::filter(exposure==6) 
exp6_annotations <- data.frame(
  exposure = 6,
  Genotype = c("N2", "CB4856"),
  xmid = c(1.5, 3.5),
  label = c("ns", "****"),
  xstart = c(1, 3),
  xend = c(2, 4),
  y = c(9.98, 9.92),
  yline = c(9.92, 9.92)
)
# Create the plot for exposure == 6
TBZ <- ggplot(TBZmetab_filtered, aes(x = as.factor(Genotype), y = TBZ, group = Genotype, fill = Genotype)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8), alpha = 0.8, size = .8,color="black") + # Move jitter to be behind the boxplot
  geom_boxplot(position = position_dodge2(width = 0.8), alpha = 0.7) +
  labs(x = "Exposure time (hrs)", y = "Norm. abundance") +
  theme_minimal() +
  scale_x_discrete( labels = c("N2" = "K267","PHX2702"="K267E","CB4856"="E267","PHX2883"="E267K"))+
  scale_fill_manual(values = cols) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(7.5, 10.05)) +  # Set y-axis limits
  cowplot::theme_cowplot(10) +
  theme(
    axis.text.x = element_text(size=12,angle = 45,hjust = 1,vjust = 1),
    axis.ticks.x = element_blank(),
    plot.background = element_rect(fill = "white"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.background = element_rect(fill = "grey80", color = "black", size = 1),
    strip.text.y = element_text(size = 12),
    strip.text = element_text(size = 12),
    axis.title.x = element_blank(),
    legend.position = "none",
    text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    plot.margin = unit(c(.1, 0.1, 0.1, 0.1), "cm")  # Reduce plot margins
  ) +
  # Adding annotations with brackets for exposure == 6
  geom_text(data = exp6_annotations, aes(x = xmid, y = y, label = label), size = 4, vjust = .1) +
  geom_segment(data = exp6_annotations, aes(x = xstart, xend = xend, y = yline, yend = yline), size = 0.5) +
  geom_segment(data = exp6_annotations, aes(x = xstart, xend = xstart, y = yline - 0.05, yend = yline), size = 0.5) +
  geom_segment(data = exp6_annotations, aes(x = xend, xend = xend, y = yline - 0.05, yend = yline), size = 0.5)

print(TBZ)



# Filter TBZmetab data for exposure == 6
TBZmetab_filtered <- ohmetab%>%
  dplyr::filter(exposure==6) 
exp6_annotations <- data.frame(
  exposure = 6,
  Genotype = c("N2", "CB4856"),
  xmid = c(1.5, 3.5),
  label = c("****", "ns"),
  xstart = c(1, 3),
  xend = c(2, 4),
  y = c(7.3, 7.34),
  yline = c(7.3, 7.3)
)

# Create the plot for exposure == 6
TBZOH <- ggplot(TBZmetab_filtered, aes(x = as.factor(Genotype), y = TBZ.OH, group = Genotype, fill = Genotype)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8), alpha = 0.8, size = .8,color="black") + # Move jitter to be behind the boxplot
  geom_boxplot(position = position_dodge2(width = 0.8), alpha = 0.7) +
  labs(x = "Exposure time (hrs)", y = "Norm. abundance") +
  theme_minimal() +
  scale_x_discrete( labels = c("N2" = "K267","PHX2702"="K267E","CB4856"="E267","PHX2883"="E267K"))+
  scale_fill_manual(values = cols) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(6, 7.6))+
  cowplot::theme_cowplot(10) +
  theme(
    axis.text.x = element_text(size=12,angle = 45,hjust = 1,vjust = 1),
    axis.ticks.x = element_blank(),
    plot.background = element_rect(fill = "white"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.background = element_rect(fill = "grey80", color = "black", size = 1),
    strip.text.y = element_text(size = 12),
    strip.text = element_text(size = 12),
    axis.title.x = element_blank(),
    legend.position = "none",
    text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_blank(),
    plot.margin = unit(c(.1, 0.1, 0.1, 0.1), "cm")  # Reduce plot margins
  ) +
  # Adding annotations with brackets for exposure == 6
  geom_text(data = exp6_annotations, aes(x = xmid, y = y, label = label), size = 4, vjust = .1) +
  geom_segment(data = exp6_annotations, aes(x = xstart, xend = xend, y = yline, yend = yline), size = 0.5) +
  geom_segment(data = exp6_annotations, aes(x = xstart, xend = xstart, y = yline - 0.05, yend = yline), size = 0.5) +
  geom_segment(data = exp6_annotations, aes(x = xend, xend = xend, y = yline - 0.05, yend = yline), size = 0.5)

print(TBZOH)





TBZmetab_filtered <- ogmetab%>%
  dplyr::filter(exposure==6) 
exp6_annotations <- data.frame(
  exposure = 6,
  Genotype = c("N2", "CB4856"),
  xmid = c(1.5, 3.5),
  label = c("*", "ns"),
  xstart = c(1, 3),
  xend = c(2, 4),
  y = c(7.34, 7.34),
  yline = c(7.3, 7.3)
)
TBZOG <- ggplot(TBZmetab_filtered, aes(x = as.factor(Genotype), y = TBZ.O.glucoside, group = Genotype, fill = Genotype)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8), alpha = 0.8, size = .8,color="black") + # Move jitter to be behind the boxplot
  geom_boxplot(position = position_dodge2(width = 0.8), alpha = 0.7) +
  labs(x = "Exposure time (hrs)", y = "Norm. abundance") +
  theme_minimal() +
  scale_x_discrete( labels = c("N2" = "K267","PHX2702"="K267E","CB4856"="E267","PHX2883"="E267K"))+
  scale_fill_manual(values = cols) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(6.1, 7.5))+
  cowplot::theme_cowplot(10) +
  theme(
    axis.text.x = element_text(size=12,angle = 45,hjust = 1,vjust = 1),
    axis.ticks.x = element_blank(),
    plot.background = element_rect(fill = "white"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.background = element_rect(fill = "grey80", color = "black", size = 1),
    strip.text.y = element_text(size = 12),
    strip.text = element_text(size = 12),
    axis.title.x = element_blank(),
    legend.position = "none",
    text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_blank(),
    plot.margin = unit(c(.1, 0.1, 0.1, 0.1), "cm")  # Reduce plot margins
  ) +
  # Adding annotations with brackets for exposure == 6
  geom_text(data = exp6_annotations, aes(x = xmid, y = y, label = label), size = 4, vjust = .1) +
  geom_segment(data = exp6_annotations, aes(x = xstart, xend = xend, y = yline, yend = yline), size = 0.5) +
  geom_segment(data = exp6_annotations, aes(x = xstart, xend = xstart, y = yline - 0.05, yend = yline), size = 0.5) +
  geom_segment(data = exp6_annotations, aes(x = xend, xend = xend, y = yline - 0.05, yend = yline), size = 0.5)

print(TBZOG)




TBZmetab_filtered <- opmetab%>%
  dplyr::filter(exposure==6)

exp6_annotations <- data.frame(
  exposure = 6,
  Genotype = c("N2", "CB4856"),
  xmid = c(1.5, 3.5),
  label = c("ns", "ns"),
  xstart = c(1, 3),
  xend = c(2, 4),
  y = c(8.02, 8.02),
  yline = c(7.98, 7.98)
)

TBZOP <- ggplot(TBZmetab_filtered, aes(x = as.factor(Genotype), y = TBZ.O.phosphoglucoside, group = Genotype, fill = Genotype)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8), alpha = 0.8, size = .8,color="black") + # Move jitter to be behind the boxplot
  geom_boxplot(position = position_dodge2(width = 0.8), alpha = 0.7) +
  labs(x = "Exposure time (hrs)", y = "Norm. abundance") +
  theme_minimal() +
  scale_x_discrete( labels = c("N2" = "K267","PHX2702"="K267E","CB4856"="E267","PHX2883"="E267K"))+
  scale_fill_manual(values = cols) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(6.8, 8.1))+
  cowplot::theme_cowplot(10) +
  theme(
    axis.text.x = element_text(size=12,angle = 45,hjust = 1,vjust = 1),
    axis.ticks.x = element_blank(),
    plot.background = element_rect(fill = "white"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.background = element_rect(fill = "grey80", color = "black", size = 1),
    strip.text.y = element_text(size = 12),
    strip.text = element_text(size = 12),
    axis.title.x = element_blank(),
    legend.position = "none",
    text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_blank(),
    plot.margin = unit(c(.1, 0.1, 0.1, 0.1), "cm")  # Reduce plot margins
  ) +
  # Adding annotations with brackets for exposure == 6
  geom_text(data = exp6_annotations, aes(x = xmid, y = y, label = label), size = 4, vjust = .1) +
  geom_segment(data = exp6_annotations, aes(x = xstart, xend = xend, y = yline, yend = yline), size = 0.5) +
  geom_segment(data = exp6_annotations, aes(x = xstart, xend = xstart, y = yline - 0.05, yend = yline), size = 0.5) +
  geom_segment(data = exp6_annotations, aes(x = xend, xend = xend, y = yline - 0.05, yend = yline), size = 0.5)

print(TBZOP)



metabplot <- plot_grid(TBZ, TBZOH, TBZOG, TBZOP, ncol = 4, align = "h",label_x=-.02,rel_widths = c(1.865,1.798,1.798,1.798), label_y = c(1, 1.1, 1.1, 1.1))
metabplot

ggsave("plots/Fig4B.png", plot = metabplot, units = "in",width = 7.36, dpi = 300)

###SupplementalPLOTS###
cols<- c("N2"="orange","PHX2702"="#E8E8E8","CB4856"="blue","PHX2883"="#696969")


# Create subsets of data for annotations
exp2_annotations <- data.frame(
  exposure = 2,
  Genotype = c("N2", "CB4856"),
  xmid = c(1.5, 3.5),
  label = c("ns", "****"),
  xstart = c(1, 3),
  xend = c(2, 4),
  y = c(9.98, 9.92),
  yline = c(9.92, 9.92)
)

exp6_annotations <- data.frame(
  exposure = 6,
  Genotype = c("N2", "CB4856"),
  xmid = c(1.5, 3.5),
  label = c("ns", "****"),
  xstart = c(1, 3),
  xend = c(2, 4),
  y = c(9.98, 9.92),
  yline = c(9.92, 9.92)
)

TBZ <- ggplot(TBZmetab, aes(x = as.factor(Genotype), y = TBZ, group = Genotype, fill = Genotype)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8), alpha = 0.8,size=.8) + # Move jitter to be behind the boxplot
  geom_boxplot(position = position_dodge2(width = 0.8, preserve = "single"), alpha = 0.7) +
  labs(x = "Exposure time (hrs)", y = "Norm. abundance") +
  theme_minimal() +
  scale_fill_manual(values = cols) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(7.5, 10.05)) +  # Set y-axis limits
  cowplot::theme_cowplot(10) +
  theme(axis.ticks.x = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(fill = "grey80", color = "black", size = 1),
        strip.text.y = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.margin = unit(c(.1, 0.1, 0.1, 0.1), "cm")  # Reduce plot margins
  )  +
  facet_grid(condition ~ exposure, scales = "free_y", labeller = labeller(exposure = c(`2` = "2 hours", `6` = "6 hours")))+
  # Adding annotations with brackets_exp2
  geom_text(data = exp2_annotations, aes(x = xmid, y = y, label = label), size = 4, vjust = .1) +
  geom_segment(data = exp2_annotations, aes(x = xstart, xend = xend, y = yline, yend = yline), size = 0.5) +
  geom_segment(data = exp2_annotations, aes(x = xstart, xend = xstart, y = yline - 0.05, yend = yline), size = 0.5) +
  geom_segment(data = exp2_annotations, aes(x = xend, xend = xend, y = yline - 0.05, yend = yline), size = 0.5) +
  # Adding annotations with brackets_exp6
  geom_text(data = exp6_annotations, aes(x = xmid, y = y, label = label), size = 4, vjust = .1) +
  geom_segment(data = exp6_annotations, aes(x = xstart, xend = xend, y = yline, yend = yline), size = 0.5) +
  geom_segment(data = exp6_annotations, aes(x = xstart, xend = xstart, y = yline - 0.05, yend = yline), size = 0.5) +
  geom_segment(data = exp6_annotations, aes(x = xend, xend = xend, y = yline - 0.05, yend = yline), size = 0.5) 

print(TBZ)

exp2_annotations <- data.frame(
  exposure = 2,
  Genotype = c("N2", "CB4856"),
  xmid = c(1.5, 3.5),
  label = c("****", "ns"),
  xstart = c(1, 3),
  xend = c(2, 4),
  y = c(7.14, 7.14),
  yline = c(7.14, 7.14)
)

# Updated exp6_annotations based on provided annotation positions
exp6_annotations <- data.frame(
  exposure = 6,
  Genotype = c("N2", "CB4856"),
  xmid = c(1.5, 3.5),
  label = c("****", "ns"),
  xstart = c(1, 3),
  xend = c(2, 4),
  y = c(7.3, 7.34),
  yline = c(7.3, 7.3)
)

TBZOH <- ggplot(ohmetab, aes(x = as.factor(Genotype), y = TBZ.OH, group = Genotype, fill = Genotype)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8), alpha = 0.8,size=.8) + # Move jitter to be behind the boxplot
  geom_boxplot(position = position_dodge2(width = 0.8, preserve = "single"), alpha = 0.7) +
  labs(x = "Exposure time (hrs)", y = "Norm. abundance") +
  theme_minimal() +
  scale_fill_manual(values = cols) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(6, 7.6))+
  cowplot::theme_cowplot(10) +
  theme(axis.ticks.x = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(fill = "grey80", color = "black", size = 1),
        strip.background.x =  element_blank(),  # Remove background of the strip
        strip.text.x= element_blank(),  # Remove Y-axis strip text
        strip.text.y = element_text( size = 12),
        strip.text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")  # Reduce plot margins
  ) +
  facet_grid(condition ~ exposure, scales = "free_y") +
  # Adding annotations with brackets_exp2
  geom_text(data = exp2_annotations, aes(x = xmid, y = y, label = label), size = 4, vjust = .1) +
  geom_segment(data = exp2_annotations, aes(x = xstart, xend = xend, y = yline, yend = yline), size = 0.5) +
  geom_segment(data = exp2_annotations, aes(x = xstart, xend = xstart, y = yline - 0.05, yend = yline), size = 0.5) +
  geom_segment(data = exp2_annotations, aes(x = xend, xend = xend, y = yline - 0.05, yend = yline), size = 0.5) +
  # Adding annotations with brackets_exp6
  geom_text(data = exp6_annotations, aes(x = xmid, y = y, label = label), size = 4, vjust = .1) +
  geom_segment(data = exp6_annotations, aes(x = xstart, xend = xend, y = yline, yend = yline), size = 0.5) +
  geom_segment(data = exp6_annotations, aes(x = xstart, xend = xstart, y = yline - 0.05, yend = yline), size = 0.5) +
  geom_segment(data = exp6_annotations, aes(x = xend, xend = xend, y = yline - 0.05, yend = yline), size = 0.5) 

print(TBZOH)

exp2_annotations <- data.frame(
  exposure = 2,
  Genotype = c("N2", "CB4856"),
  xmid = c(1.5, 3.5),
  label = c("*", "ns"),
  xstart = c(1, 3),
  xend = c(2, 4),
  y = c(7.14, 7.18),
  yline = c(7.14, 7.14)
)

# Updated exp6_annotations based on provided annotation positions
exp6_annotations <- data.frame(
  exposure = 6,
  Genotype = c("N2", "CB4856"),
  xmid = c(1.5, 3.5),
  label = c("*", "ns"),
  xstart = c(1, 3),
  xend = c(2, 4),
  y = c(7.34, 7.34),
  yline = c(7.3, 7.3)
)

TBZOG <- ggplot(ogmetab, aes(x = as.factor(Genotype), y = TBZ.O.glucoside, group = Genotype, fill = Genotype)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8), alpha = 0.8,size=.8) + # Move jitter to be behind the boxplot
  geom_boxplot(position = position_dodge2(width = 0.8, preserve = "single"), alpha = 0.7) +
  labs(x = "Exposure time (hrs)", y = "Norm. abundance") +
  theme_minimal() +
  scale_fill_manual(values = cols) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(6.1, 7.5))+
  cowplot::theme_cowplot(10) +
  theme(axis.ticks.x = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(fill = "grey80", color = "black", size = 1),
        strip.background.x =  element_blank(),  # Remove background of the strip
        strip.text.x= element_blank(),  # Remove Y-axis strip text
        strip.text.y = element_text( size = 12),
        strip.text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")  # Reduce plot margins
  )+
  facet_grid(condition ~ exposure, scales = "free_y") +
  # Adding annotations with brackets_exp2
  geom_text(data = exp2_annotations, aes(x = xmid, y = y, label = label), size = 4, vjust = .1) +
  geom_segment(data = exp2_annotations, aes(x = xstart, xend = xend, y = yline, yend = yline), size = 0.5) +
  geom_segment(data = exp2_annotations, aes(x = xstart, xend = xstart, y = yline - 0.05, yend = yline), size = 0.5) +
  geom_segment(data = exp2_annotations, aes(x = xend, xend = xend, y = yline - 0.05, yend = yline), size = 0.5) +
  # Adding annotations with brackets_exp6
  geom_text(data = exp6_annotations, aes(x = xmid, y = y, label = label), size = 4, vjust = .1) +
  geom_segment(data = exp6_annotations, aes(x = xstart, xend = xend, y = yline, yend = yline), size = 0.5) +
  geom_segment(data = exp6_annotations, aes(x = xstart, xend = xstart, y = yline - 0.05, yend = yline), size = 0.5) +
  geom_segment(data = exp6_annotations, aes(x = xend, xend = xend, y = yline - 0.05, yend = yline), size = 0.5) 

print(TBZOG)

exp2_annotations <- data.frame(
  exposure = 2,
  Genotype = c("N2", "CB4856"),
  xmid = c(1.5, 3.5),
  label = c("ns", "ns"),
  xstart = c(1, 3),
  xend = c(2, 4),
  y = c(7.48, 7.50),
  yline = c(7.46, 7.46)
)

# Updated exp6_annotations based on provided annotation positions
exp6_annotations <- data.frame(
  exposure = 6,
  Genotype = c("N2", "CB4856"),
  xmid = c(1.5, 3.5),
  label = c("ns", "ns"),
  xstart = c(1, 3),
  xend = c(2, 4),
  y = c(8.02, 8.02),
  yline = c(7.98, 7.98)
)

TBZOP <- ggplot(opmetab, aes(x = as.factor(Genotype), y = TBZ.O.phosphoglucoside, group = Genotype, fill = Genotype)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8), alpha = 0.8,size=.8) + # Move jitter to be behind the boxplot
  geom_boxplot(position = position_dodge2(width = 0.8, preserve = "single"), alpha = 0.7) +
  labs(x = "Exposure time (hrs)", y = "Norm. abundance") +
  theme_minimal() +
  scale_x_discrete( labels = c("N2" = "K267","PHX2702"="K267E","CB4856"="E267","PHX2883"="E267K"))+
  scale_fill_manual(values = cols) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(6.8, 8.1))+
  cowplot::theme_cowplot(10) +
  theme(
    plot.background = element_rect(fill = "white"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.background = element_rect(fill = "grey80", color = "black", size = 1),
    strip.background.x =  element_blank(),  # Remove background of the strip
    strip.text.x= element_blank(),  # Remove Y-axis strip text
    strip.text.y = element_text( size = 12),
    axis.text.x = element_text(size=12,angle = 45,hjust = 1,vjust = 1),
    strip.text = element_text(size = 12),
    axis.title.x = element_blank(),
    legend.position = "none",
    text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")  # Reduce plot margins
  ) +
  facet_grid(condition ~ exposure, scales = "free_y") +
  # Adding annotations with brackets_exp2
  geom_text(data = exp2_annotations, aes(x = xmid, y = y, label = label), size = 4, vjust = .1) +
  geom_segment(data = exp2_annotations, aes(x = xstart, xend = xend, y = yline, yend = yline), size = 0.5) +
  geom_segment(data = exp2_annotations, aes(x = xstart, xend = xstart, y = yline - 0.05, yend = yline), size = 0.5) +
  geom_segment(data = exp2_annotations, aes(x = xend, xend = xend, y = yline - 0.05, yend = yline), size = 0.5) +
  # Adding annotations with brackets_exp6
  geom_text(data = exp6_annotations, aes(x = xmid, y = y, label = label), size = 4, vjust = .1) +
  geom_segment(data = exp6_annotations, aes(x = xstart, xend = xend, y = yline, yend = yline), size = 0.5) +
  geom_segment(data = exp6_annotations, aes(x = xstart, xend = xstart, y = yline - 0.05, yend = yline), size = 0.5) +
  geom_segment(data = exp6_annotations, aes(x = xend, xend = xend, y = yline - 0.05, yend = yline), size = 0.5) 

print(TBZOP)

metabplot <- plot_grid(TBZ, TBZOH, TBZOG, TBZOP, ncol = 1, align = "v", labels = c("A", "B", "C", "D"), rel_heights = c(1.74, 1.64, 1.64, 1.84),label_x=-.02, label_y = c(1, 1.1, 1.1, 1.1))
metabplot
# Save the arranged plot
ggsave("plots/SFig10.png", plot = metabplot, units = "in", height = 7.67, width = 3.75, dpi = 300)







###########EXOMETABOLOME######################

xmetab <- read.csv("data/exo_metab.csv")

# Convert "exposure" to numeric (if not already done)
xmetab$exposure <- as.numeric(xmetab$exposure)
xmetab$TBZ <- as.numeric(xmetab$TBZ)
xmetab$TBZ.OH <- as.numeric(xmetab$TBZ.OH)
xmetab$TBZ.O.glucoside <- as.numeric(xmetab$TBZ.O.glucoside)
xmetab$TBZ.O.phosphoglucoside <- as.numeric(xmetab$TBZ.O.phosphoglucoside)
xmetab$Genotype <- factor(xmetab$Genotype, levels = c("N2", "PHX2702", "CB4856", "PHX2883"))



TBZmetab <- xmetab[, c(1, 2, 5)]
TBZmetab$condition <- "TBZ"

ohmetab <- xmetab[, c(1, 2, 6)]
ohmetab$condition <- "TBZ-OH"

ogmetab <-  xmetab[, c(1, 2, 8)]
ogmetab$condition <- "TBZ-O-Glu"


opmetab <-  xmetab[, c(1, 2, 9)]
opmetab$condition <- "TBZ-O-PGlu"

exposures <- c("2", "6")  
genotypes <- list(c("N2", "PHX2702"), c("CB4856", "PHX2883")) 
columns_of_interest <- c("TBZ", "TBZ.OH", "TBZ.O.glucoside", "TBZ.O.phosphoglucoside")  

exo_results_df <- perform_stats_by_groups(data = xmetab, exposures = exposures, genotypes = genotypes, columns_of_interest = columns_of_interest)

print(exo_results_df)
###PLOTS###

# Create subsets of data for annotations
exp2_annotations <- data.frame(
  exposure = 2,
  Genotype = c("N2", "CB4856"),
  xmid = c(1.5, 3.5),
  label = c("ns", "ns"),
  xstart = c(1, 3),
  xend = c(2, 4),
  y = c(11.32, 11.36),
  yline = c(11.3, 11.3)
)

exp6_annotations <- data.frame(
  exposure = 6,
  Genotype = c("N2", "CB4856"),
  xmid = c(1.5, 3.5),
  label = c("ns", "ns"),
  xstart = c(1, 3),
  xend = c(2, 4),
  y = c(11.3, 11.3),
  yline = c(11.24, 11.24)
)

TBZ <- ggplot(TBZmetab, aes(x = as.factor(Genotype), y = TBZ, group = Genotype, fill = Genotype)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8),alpha = 0.8,size=.8) + # Move jitter to be behind the boxplot
  geom_boxplot(position = position_dodge2(width = 0.8, preserve = "single"), alpha = 0.7) +
  labs(x = "Exposure time (hrs)", y = "Norm. abundance") +
  theme_minimal() +
  scale_fill_manual(values = cols) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(10.1, 11.5))+
  cowplot::theme_cowplot(10) +
  theme(axis.ticks.x = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(fill = "grey80", color = "black", size = 1),
        strip.text.y = element_text(size = 12 ),
        strip.text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.margin = unit(c(.1, 0.1, 0.1, 0.1), "cm")  # Reduce plot margins
  )  +
  facet_grid(condition ~ exposure, scales = "free_y", labeller = labeller(exposure = c(`2` = "2 hours", `6` = "6 hours")))+
  # Adding annotations with brackets_exp2
  geom_text(data = exp2_annotations, aes(x = xmid, y = y, label = label), size = 4, vjust = .1) +
  geom_segment(data = exp2_annotations, aes(x = xstart, xend = xend, y = yline, yend = yline), size = 0.5) +
  geom_segment(data = exp2_annotations, aes(x = xstart, xend = xstart, y = yline - 0.05, yend = yline), size = 0.5) +
  geom_segment(data = exp2_annotations, aes(x = xend, xend = xend, y = yline - 0.05, yend = yline), size = 0.5) +
  # Adding annotations with brackets_exp6
  geom_text(data = exp6_annotations, aes(x = xmid, y = y, label = label), size = 4, vjust = .1) +
  geom_segment(data = exp6_annotations, aes(x = xstart, xend = xend, y = yline, yend = yline), size = 0.5) +
  geom_segment(data = exp6_annotations, aes(x = xstart, xend = xstart, y = yline - 0.05, yend = yline), size = 0.5) +
  geom_segment(data = exp6_annotations, aes(x = xend, xend = xend, y = yline - 0.05, yend = yline), size = 0.5) 

print(TBZ)

exp2_annotations <- data.frame(
  exposure = 2,
  Genotype = c("N2", "CB4856"),
  xmid = c(1.5, 3.5),
  label = c("ns", "ns"),
  xstart = c(1, 3),
  xend = c(2, 4),
  y = c(8.02, 8.06),
  yline = c(8, 8)
)

# Updated exp6_annotations based on provided annotation positions
exp6_annotations <- data.frame(
  exposure = 6,
  Genotype = c("N2", "CB4856"),
  xmid = c(1.5, 3.5),
  label = c("ns", "ns"),
  xstart = c(1, 3),
  xend = c(2, 4),
  y = c(8.92, 8.96),
  yline = c(8.9, 8.9)
)

TBZOH <- ggplot(ohmetab, aes(x = as.factor(Genotype), y = TBZ.OH, group = Genotype, fill = Genotype)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8), alpha = 0.8,size=.8) + # Move jitter to be behind the boxplot
  geom_boxplot(position = position_dodge2(width = 0.8, preserve = "single"), alpha = 0.7) +
  labs(x = "Exposure time (hrs)", y = "Norm. abundance") +
  theme_minimal() +
  scale_fill_manual(values = cols) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(6.9, 9.25))+
  cowplot::theme_cowplot(10) +
  theme(axis.ticks.x = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(fill = "grey80", color = "black", size = 1),
        strip.background.x =  element_blank(),  # Remove background of the strip
        strip.text.x= element_blank(),  # Remove Y-axis strip text
        strip.text.y = element_text( size = 12),
        strip.text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")  # Reduce plot margins
  ) +
  facet_grid(condition ~ exposure, scales = "free_y") +
  # Adding annotations with brackets_exp2
  geom_text(data = exp2_annotations, aes(x = xmid, y = y, label = label), size = 4, vjust = .1) +
  geom_segment(data = exp2_annotations, aes(x = xstart, xend = xend, y = yline, yend = yline), size = 0.5) +
  geom_segment(data = exp2_annotations, aes(x = xstart, xend = xstart, y = yline - 0.05, yend = yline), size = 0.5) +
  geom_segment(data = exp2_annotations, aes(x = xend, xend = xend, y = yline - 0.05, yend = yline), size = 0.5) +
  # Adding annotations with brackets_exp6
  geom_text(data = exp6_annotations, aes(x = xmid, y = y, label = label), size = 4, vjust = .1) +
  geom_segment(data = exp6_annotations, aes(x = xstart, xend = xend, y = yline, yend = yline), size = 0.5) +
  geom_segment(data = exp6_annotations, aes(x = xstart, xend = xstart, y = yline - 0.05, yend = yline), size = 0.5) +
  geom_segment(data = exp6_annotations, aes(x = xend, xend = xend, y = yline - 0.05, yend = yline), size = 0.5) 

print(TBZOH)

exp2_annotations <- data.frame(
  exposure = 2,
  Genotype = c("N2", "CB4856"),
  xmid = c(1.5, 3.5),
  label = c("ns", "ns"),
  xstart = c(1, 3),
  xend = c(2, 4),
  y = c(8, 8.04),
  yline = c(7.98, 7.98)
)

# Updated exp6_annotations based on provided annotation positions
exp6_annotations <- data.frame(
  exposure = 6,
  Genotype = c("N2", "CB4856"),
  xmid = c(1.5, 3.5),
  label = c("ns", "ns"),
  xstart = c(1, 3),
  xend = c(2, 4),
  y = c(8.81, 8.81),
  yline = c(8.75, 8.75)
)



TBZOG <- ggplot(ogmetab, aes(x = as.factor(Genotype), y = TBZ.O.glucoside, group = Genotype, fill = Genotype)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8), alpha = 0.8,size=.8) + # Move jitter to be behind the boxplot
  geom_boxplot(position = position_dodge2(width = 0.8, preserve = "single"), alpha = 0.7) +
  labs(x = "Exposure time (hrs)", y = "Norm. abundance") +
  theme_minimal() +
  scale_fill_manual(values = cols) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(6.9, 9))+
  cowplot::theme_cowplot(10) +
  theme(axis.ticks.x = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        strip.background = element_rect(fill = "grey80", color = "black", size = 1),
        strip.background.x =  element_blank(),  # Remove background of the strip
        strip.text.x= element_blank(),  # Remove Y-axis strip text
        strip.text.y = element_text( size = 12),
        strip.text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        text = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")  # Reduce plot margins
  )+
  facet_grid(condition ~ exposure, scales = "free_y") +
  # Adding annotations with brackets_exp2
  geom_text(data = exp2_annotations, aes(x = xmid, y = y, label = label), size = 4, vjust = .1) +
  geom_segment(data = exp2_annotations, aes(x = xstart, xend = xend, y = yline, yend = yline), size = 0.5) +
  geom_segment(data = exp2_annotations, aes(x = xstart, xend = xstart, y = yline - 0.05, yend = yline), size = 0.5) +
  geom_segment(data = exp2_annotations, aes(x = xend, xend = xend, y = yline - 0.05, yend = yline), size = 0.5) +
  # Adding annotations with brackets_exp6
  geom_text(data = exp6_annotations, aes(x = xmid, y = y, label = label), size = 4, vjust = .1) +
  geom_segment(data = exp6_annotations, aes(x = xstart, xend = xend, y = yline, yend = yline), size = 0.5) +
  geom_segment(data = exp6_annotations, aes(x = xstart, xend = xstart, y = yline - 0.05, yend = yline), size = 0.5) +
  geom_segment(data = exp6_annotations, aes(x = xend, xend = xend, y = yline - 0.05, yend = yline), size = 0.5) 

print(TBZOG)


# Updated exp6_annotations based on provided annotation positions
exp6_annotations <- data.frame(
  exposure = 6,
  Genotype = c("N2", "CB4856"),
  xmid = c(1.5, 3.5),
  label = c("ns", "ns"),
  xstart = c(1, 3),
  xend = c(2, 4),
  y = c(5.65, 5.65),
  yline = c(5.59, 5.59)
)

# Assuming exposure is a factor variable with levels "2" and "6"
opmetab$exposure <- factor(opmetab$exposure, levels = c("2", "6"))

# Assuming your data frame is named 'your_data'
opmetab2 <- opmetab %>%
  mutate(TBZ.O.phosphoglucoside = ifelse(exposure == 2 & Genotype != "N2" & is.na(TBZ.O.phosphoglucoside), 0, TBZ.O.phosphoglucoside))

TBZOP <- ggplot(opmetab2, aes(x = as.factor(Genotype), y = TBZ.O.phosphoglucoside, group = Genotype, fill = Genotype)) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8),alpha = 0.8,size=.8) + # Move jitter to be behind the boxplot
  geom_boxplot(position = position_dodge2(width = 0.8, preserve = "single"), alpha = 0.7) +
  labs(x = "Exposure time (hrs)", y = "Norm. abundance") +
  theme_minimal() +
  scale_fill_manual(values = cols) +
  scale_x_discrete( labels = c("N2" = "K267","PHX2702"="K267E","CB4856"="E267","PHX2883"="E267K"))+
  scale_y_continuous(limits = c(4, 5.8), oob = scales::squish)+
  cowplot::theme_cowplot(10) +
  theme(
    plot.background = element_rect(fill = "white"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    strip.background = element_rect(fill = "grey80", color = "black", size = 1),
    strip.background.x =  element_blank(),  # Remove background of the strip
    strip.text.x= element_blank(),  # Remove Y-axis strip text
    strip.text.y = element_text( size = 12),
    axis.text.x = element_text(size=12,angle = 45,hjust = 1,vjust = 1),
    strip.text = element_text(size = 12),
    axis.title.x = element_blank(),
    legend.position = "none",
    text = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")  # Reduce plot margins
  ) +
  facet_grid(condition ~ exposure, scales = "free_y") +
  geom_text(data = exp6_annotations, aes(x = xmid, y = y, label = label), size = 4, vjust = .1) +
  geom_segment(data = exp6_annotations, aes(x = xstart, xend = xend, y = yline, yend = yline), size = 0.5) +
  geom_segment(data = exp6_annotations, aes(x = xstart, xend = xstart, y = yline - 0.05, yend = yline), size = 0.5) +
  geom_segment(data = exp6_annotations, aes(x = xend, xend = xend, y = yline - 0.05, yend = yline), size = 0.5) 

print(TBZOP)

metabplot <- plot_grid(TBZ, TBZOH, TBZOG, TBZOP, ncol = 1, align = "v", labels = c("A", "B", "C", "D"), rel_heights = c(1.74, 1.64, 1.64, 1.84),label_x=-.02, label_y = c(1, 1.1, 1.1, 1.075))
metabplot
# Save the arranged plot
ggsave("plots/Sfig11.png", plot = metabplot, units = "in", height = 7.67, width = 3.5, dpi = 300)


##Fig5 Yeast metab
test <- read.csv("data/yeast_work2.csv")
test <- test %>%
  dplyr::mutate(relative=TBZOH/TBZ)

# Control Delta
control_values <- test %>%
  dplyr::filter(strain=="K267") %>%
  dplyr::group_by(strain, trial) %>%
  dplyr::mutate(control_pheno = mean(relative)) %>% # get mean control value for each trait and strain
  dplyr::distinct(control_pheno, strain, trial) # make it neat
delta <- test %>%
  dplyr::ungroup() %>%
  dplyr::left_join(., control_values, by="trial") %>% # join control values for each trait
  dplyr::mutate(relative_delta = relative - control_pheno)

rel.reg <- lm(data = delta, formula = relative_delta ~ trial)
reg_relative <- rel.reg$residuals
delta <- delta[as.numeric(names(rel.reg$residuals)),]
coefs <- lm(data = delta, formula = relative_delta ~ trial - 1)
bleach.coefs <- data.frame(gsub(names(coefs$coefficients), pattern = "trial", replacement = ""),
                           as.numeric(coefs$coefficients))
colnames(bleach.coefs) <- c("trial","Atrial_Effect")
regressed <- cbind(delta, reg_relative)

murp <- regressed %>%
  dplyr::mutate(test=relative/control_pheno)

stats <- murp%>%
  aov(test ~ strain.x, data = .)%>%
  rstatix::tukey_hsd() 



elm <- murp%>%
  dplyr::filter(strain.x=="K267E")

mean(elm$test)



cols <- c("K267"="orange", "K267E"="grey")


mep <- murp %>%
  ggplot+
  aes(x=strain.x, y=test)+
  geom_jitter(width=0.1, size = 1)+
  geom_boxplot(aes(fill = strain.x),alpha = .6,outlier.shape = NA)+
  xlab("")+
  ylab("Activity relative to wild-type allele")+
  scale_fill_manual(name= " ", labels = c("K267" = "Lysine", "K267E" ="Glutamate"),values=cols)+
  scale_x_discrete(labels=c("K267" = "Lysine", "K267E" = "Glutamate"))+
  stat_pvalue_manual(stats, label = "p.adj.signif", y.position = c(1.3), size=6, xmax = "group2", remove.bracket = TRUE) +
  cowplot::theme_cowplot(12) +
  theme(axis.text.y = element_text(size=12),axis.text.x = element_text(size=12),axis.title.y = element_text(size=12),
        plot.background = element_rect(fill="white"),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
  )
mep
ggsave("plots/Fig5.png", plot = mep, device = "png", width = 3.75, height = 4, units = "in",dpi=300)

##Fig6
world <- map_data("world")
world <- world[world$region != "Antarctica",] # intercourse antarctica
isolation_infosev <- rio::import("data/S14.tsv")
# Reorder the levels of cypvartype so "CB4856" (blue) comes after "K267D" (green)
isolation_infosev$cypvartype <- factor(isolation_infosev$cypvartype, 
                                       levels = c("K267D", "K267", "N2", "Varnot267", "CB4856", "K267E"))

# Plot
world <- map_data("world")
world <- world[world$region != "Antarctica",] # remove Antarctica

map <- ggplot() + 
  geom_map(data = world, map = world,
           aes(x = long, y = lat, map_id = region),
           color = "white", fill = "#7f7f7f", 
           size = 0.1, alpha = 0.5) +
  geom_point(data = isolation_infosev, aes(x = long, y = lat, color = cypvartype), size = 1) +
  scale_color_manual(name = "267 Variation",
                     values = c("CB4856" = "blue", "K267D" = "cadetblue", 
                                "K267E" = "blue", "K267" = "orange", 
                                "N2" = "orange", "Varnot267" = "pink")) +
  theme_map(12) + 
  theme(text = element_text(size = 12),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position = "none")


#6B
treecyp <- ape::read.tree("data/S16.tree")%>%
  as.tibble(phytools::midpoint.root(tree))
treecyp <- treecyp %>% 
  dplyr::mutate(cypvartype = case_when(
    (label %in% dplyr::filter(tbz_norm,cypvartype == "K267")$strain) ~ "K267",
    (label %in% dplyr::filter(tbz_norm,cypvartype == "K267E")$strain) ~ "K267E",
    (label %in% dplyr::filter(tbz_norm,cypvartype == "K267D")$strain) ~ "K267D",
    TRUE ~ "K267"))
treecyp$branch.length[treecyp$branch.length < 0] <- 0  # Setting negative lengths to zero
treecypplot <- treecyp %>%
  tidytree::as.phylo(use.labels = TRUE) %>%
  ggtree::ggtree(size = 0.2, color = "gray51") +
  ggplot2::coord_flip() +
  ggplot2::scale_x_reverse() +  
  ggtree::geom_tippoint(aes(fill = treecyp$cypvartype), shape = 21, size = 1) +
  ggplot2::scale_fill_manual(name = "267 Variation", values = c("K267" = "orange", "K267E" = "blue","K267D" = "cadetblue"),labels=c("K267" = "Lysine", "K267E" = "Glutamate","K267D"="Aspartate") ) +
  ggplot2::theme(text=element_text(size=12),legend.position = c(0.1, 0.3))

evolution <- cowplot::plot_grid(map,treecypplot,ncol = 1, align = "hv", axis = "blrt", labels = c("A","B"),label_y = 1.025)
ggsave("plots/Fig6.png",plot=evolution, width = 7.5, height = 5.5, units = "in",dpi=300)

###Fig7###
if (!requireNamespace("devtools", quietly=TRUE))
  install.packages("devtools")
devtools::install_github("YuLab-SMU/ggmsa")
install.packages("ggtree")


treecaenor <- ape::read.tree("data/phylo.txt")

treecaenor<- ape::root(treecaenor, outgroup = "D.coronatus", resolve.root = TRUE)



 p<- treecaenor %>%
  tidytree::as.phylo(use.labels = T) %>%
  ggtree::ggtree(color="gray51") +
  ggplot2::theme(legend.position = "none")+
  ggtree::geom_tiplab(aes(label = label), fontface = "italic")+
  ggtree::geom_treescale(x = 0,y=-.1)+
  xlim(0,0.6)

p1<-p+theme_classic(12)+theme(text=element_text(size=12),axis.line = element_blank(),
                              axis.text = element_blank(),
                              axis.ticks = element_blank())

library(ggmsa)
ali <-ggmsa(c("data/S19.fasta"),30,50, seq_name = T, color = "Zappo_AA", by_conservation = T,position_highlight = 44)+
  theme(text=element_text(size=12),axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, face = "italic"))


alignment <- cowplot::plot_grid(p1,ali,ncol = 1, align = "v",axis = "r",labels = c("A","B","C"))
ggsave("plots/Fig7.png", plot = alignment, width = 7.5, height = 6,units = "in",dpi=300)

##SFig8
tbznob1 <- tbz_normnoregression %>%
  dplyr::mutate(cypvartype = case_when(cypvartype == "Varnot267" ~ "K267",
                                       TRUE ~ cypvartype))%>%
  ggplot()+
  aes(x=cypvartype, y= Thiabendazole_mean.TOF)+
  geom_jitter(width = 0.1)+
  ylab("Regressed animal length")+
  geom_boxplot(aes(fill=cypvartype, alpha =0.8),outlier.shape = NA)+
  scale_fill_manual(values = c("K267" = "orange", "K267D" = "cadetblue", "K267E" = "blue"))+
  scale_x_discrete(breaks = c(c("K267","K267D","K267E")),
                   labels = c("Lysine","Aspartate","Glutamate"))+cowplot::theme_cowplot(12)+
  xlab("cyp-35D1 267 allele")+
  ggtitle("Non-ben-1 regressed")+
  theme(  panel.background = element_rect(fill = "white",
                                          colour = NA_character_), # necessary to avoid drawing panel outline
          legend.position = "none",
          text = element_text(size=12),
          axis.title.x = element_blank(),
          plot.title = element_text(face = "bold",size = 12, hjust = 0.5))


statstbznoreg <- tbz_normnoregression %>%
  aov(Thiabendazole_mean.TOF ~ cypben1type, data = .)%>%
  rstatix::tukey_hsd()

tbzb1 <- tbz_norm %>%
  dplyr::mutate(cypvartype = case_when(cypvartype == "Varnot267" ~ "K267",
                                       TRUE ~ cypvartype))%>%
  ggplot()+
  aes(x=cypvartype, y= Thiabendazole_mean.TOF)+
  ylab("Regressed animal length")+
  geom_jitter(width = 0.1)+
  geom_boxplot(aes(fill=cypvartype, alpha =0.4),outlier.shape = NA)+
  scale_fill_manual(values = c("K267" = "orange", "K267D" = "cadetblue", "K267E" = "blue"))+
  scale_x_discrete(breaks = c(c("K267","K267D","K267E")),
                   labels = c("Lysine","Aspartate","Glutamate")) +cowplot::theme_cowplot(12)+
  xlab("cyp-35D1 267 allele")+
  ggtitle("ben-1 regressed")+
  theme(  panel.background = element_rect(fill = "white",
                                          colour = NA_character_), # necessary to avoid drawing panel outline
          legend.position = "none",
          text = element_text(size=12),
          axis.title.x = element_blank(),
          plot.title = element_text(face = "bold",size = 12, hjust = 0.5))
allelesplit <-cowplot::plot_grid(tbznob1,tbzb1, labels = c("A","B"))

ggsave(filename = "plots/SFig8.png",plot = allelesplit, width = 7.5, height = 4, units = "in")

##Sfig9

DL238fig <- data3 %>%
  dplyr::filter(strain %in% c("N2","PHX2701","DL238","ECA3359"))
DL238fig2 <-data3 %>%
  dplyr::filter(strain %in% c("N2","PHX2701","PHX2702","DL238","ECA3352","ECA3359"))


stats <- DL238fig %>%
  dplyr::filter(concentration_um==32.5)%>%
  aov(median_wormlength_um_reg_delta ~ strain, data = .)%>%
  rstatix::tukey_hsd()

#SupplementalFig9
#Figure A
genomeforallele <- data.table::fread("data/S12.txt")
genomeforalleleplot <- ggplot(genomeforallele)+
  geom_segment(aes(x = start/1e6, y = factor(sample,levels = c("DL238","ECA3359","ECA3352","PHX2702","PHX2701","N2")), xend = end/1e6, yend = sample, color = gt_name, size = 2))+
  scale_color_manual(values=c("N2"="orange","DL238"="cadetblue", "ECA3359" = "orange","ECA3352"="orange","PHX2702"="orange","PHX2701"="orange"))+
  # scale_color_manual(values=c("1"="#F0BA51","2"="#484DA0"))+
  facet_grid(~chrom, scales = "free",  space = "free")+
  theme_cowplot(12)+
  theme(plot.background = element_rect(fill="white"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=12,color = "black"),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks.x = element_blank())

#ggsave("test1.png", plot = genomeforalleleplot, width = 8, height = 4, units = "in")
#Figure B
cypforallele <- genomeforallele %>%
  dplyr::mutate(chrom = "cyp-35d1")
cypforalleleplot <- ggplot(cypforallele)+
  geom_segment(aes(x = start/1e6, y = factor(sample,levels = c("DL238","ECA3359","ECA3352","PHX2702","PHX2701","N2")), xend = end/1e6, yend = sample, color = gt_name, size = 2))+
  scale_color_manual(values=c("N2"="orange","DL238"="cadetblue", "ECA3359" = "cadetblue","ECA3352"="cadetblue","PHX2702"="blue","PHX2701"="blue"))+
  # scale_color_manual(values=c("1"="#F0BA51","2"="#484DA0"))+
  facet_grid(~chrom, scales = "free",  space = "free")+
  theme_cowplot(12)+
  theme(plot.background = element_rect(fill="white"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=12,face = "italic", color = "black",hjust=0.7),
        plot.title = element_text( hjust=0.5),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank())
#ggsave("test2.png", plot = cypforalleleplot, width = 8, height = 4, units = "in")


#Figure C
fullregressedallelesplot <- DL238fig2 %>%
  dplyr::mutate(strain=factor(strain, levels = c("DL238","ECA3359","ECA3352","PHX2702","PHX2701","N2")))%>%
  dplyr::filter(concentration_um==32.5)%>%
  ggplot+
  aes(x=strain, y=median_wormlength_um_reg_delta )+
  geom_jitter(width=0.1, size = 0.2)+
  geom_boxplot(aes(fill = strain),alpha = 0.4,outlier.shape = NA)+
  xlab("Regressed animal length")+
  ylab("Regressed animal length")+
  scale_fill_manual(values  = c("N2"="orange","unknown"="grey","DL238"="cadetblue"))+
  cowplot::theme_cowplot(12)+
  coord_flip()+
  ggpubr::geom_bracket(xmin = "N2", xmax = "PHX2701", y.position = -10,label = "****", coord.flip = TRUE)+
  ggpubr::geom_bracket(xmin = "N2", xmax = "PHX2702", y.position = 20,label = "**", coord.flip = TRUE)+  
  ggpubr::geom_bracket(xmin = "N2", xmax = "ECA3352", y.position = 50,label = "ns", coord.flip = TRUE)+  
  ggpubr::geom_bracket(xmin = "N2", xmax = "ECA3359", y.position = 80,label = "ns", coord.flip = TRUE)+
  ggpubr::geom_bracket(xmin = "N2", xmax = "DL238", y.position = 110,label = "****", coord.flip = TRUE)+
  theme(text=element_text(size=12),axis.text.x = element_text(size = 12),
        plot.background = element_rect(fill="white"),
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text( color = "black"),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
DL238suppFig <- cowplot::plot_grid(genomeforalleleplot,cypforalleleplot,fullregressedallelesplot,labels = c("A","","B"), nrow=1, rel_widths = c(0.35,0.2,1),align =  "h", axis = "bltr" )
ggsave("plots/SFig9.png", plot = DL238suppFig, width = 7.5, height = 4, units = "in")

##TajimasD SFig12
tdunfilt <- rio::import("data/SFile S9.csv")
tdfilt <- tdunfilt %>%
  dplyr::filter(!is.na(td)) %>%
  dplyr::filter(start > 16044238 & start < 16094238)%>%
  dplyr::mutate(start = start/1e6)%>%
  dplyr::mutate(end = end/1e6)
tajimasd <- tdfilt %>%
  ggplot()+
  aes(x=((start+end)/2), y=td)+
  geom_point()+
  ylab("Tajima's D")+
  xlab("Genomic position (Mb)")+
  geom_vline(xintercept = 16069238/1e6, colour="red")+
  geom_vline(xintercept = 16071318/1e6, colour = "red")+
  geom_rect(aes(xmin=16069238/1e6, xmax =16071318/1e6, ymax = Inf, ymin = -Inf), color = "red",fill="red", alpha = 0.01)+
  geom_vline(xintercept = 16071325/1e6, colour="blue", linetype = "dashed")+
  geom_vline(xintercept = 16072738/1e6, colour="blue")+
  geom_rect(aes(xmin=16071325/1e6, xmax =16072738/1e6, ymax = Inf, ymin = -Inf), colour = "blue",fill="blue", alpha = 0.01)+
  cowplot::theme_cowplot(12)+
  geom_rect(aes(xmin=16.09,xmax=16.095,ymin=-0.25,ymax=0),colour = "red",fill="red", alpha = 0.01)+
  geom_rect(aes(xmin=16.09,xmax=16.095,ymin=-0.5,ymax=-0.25),colour = "blue",fill="blue", alpha = 0.01)+
  geom_text(aes(x=16.098,y=-0.125,label="cyp-35d1"),fontface="italic")+
  geom_text(aes(x=16.09765,y=-0.375,label="nhr-176"),fontface="italic")+
  ylim(-2.2,0.2)+
  theme(text=element_text(size=12),legend.position = "none")
ggsave("plots/SFig12.png", plot=tajimasd, width = 7.5, height = 3, units = "in")


