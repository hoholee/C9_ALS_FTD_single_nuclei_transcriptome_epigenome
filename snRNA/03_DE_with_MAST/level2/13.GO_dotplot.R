library(tidyverse)
library(ggplot2)
library(ggdendro)
library(reshape2)
library(viridis)
library(scico)
library(scales)


# Exc up ------------------------------------------------------------------

df<-read_csv("/cndd2/jchien/project/ALS/GO/GO_postCellBender/GO_summary_full_AP_ALS.csv")%>% 
  filter(reg =='up') %>%
  filter(! celltype %in% c('Astro','OPC','Oligo','Micro')) %>%
  mutate(types=str_replace(types, "_up", "")) %>%
  mutate(types=str_replace(types, "_dn", "")) 

df$logFDR <--log10(df$FDR)
df$logFDR[!is.finite(df$logFDR)] <- 15

GOIs <- df %>% filter(FDR<0.01) %>% filter(enrichmentRatio>0) %>% filter(Affinity=='YES') %>% filter(celltype %in% c('Exc-superficial','Exc-deep')) %>%
  arrange(desc(FDR)) %>% group_by(types) %>% slice(tail(row_number(), 50)) %>% pull(description)


order_j=c(
  'protein folding',
  'heat shock protein binding',
  'unfolded protein binding',
  'Golgi-associated vesicle',
  
  'respiratory chain',
  'mitochondrial inner membrane',
  'mitochondrial protein complex',
  'nucleoside triphosphate metabolic process',
  'NADH dehydrogenase complex',
  'NADH dehydrogenase complex assembly',
  'proton transmembrane transport',
  'mitochondrial membrane organization',
  
  'protein localization to endoplasmic reticulum',
  'polysome',
  'translational elongation',
  
  'protein localization to nucleus',
  'RNA localization',
  'DNA damage response, detection of DNA damage',
  'nucleotide-excision repair',
  
  'antibiotic metabolic process',
  'nucleobase-containing compound transport',
  'ribonucleoprotein complex biogenesis',
  'process utilizing autophagic mechanism',
  'microtubule cytoskeleton organization involved in mitosis',
  'tRNA metabolic process',
  'organelle fission',
  'DNA biosynthetic process'
)


df<-df %>% filter(description %in% GOIs)

## hclust histograme
mat <- df %>% 
  distinct()%>%
  select(types, description, FDR) %>%
  pivot_wider(names_from = types, values_from = FDR,values_fill = 1) %>%
  data.frame() 
row.names(mat) <- mat$description  # put gene in `row`
is.na(mat$geneSet) <-1
clust <- hclust(dist(mat %>% as.matrix())) # hclust with distance matrix

df<-df %>%
  filter(FDR<0.25) %>%
  mutate(description = factor(description, levels=rev(order_j))) %>%
  mutate(types = factor(types, levels=c('MCX_Exc-superficial','mFCX_Exc-superficial','MCX_Exc-intermediate','mFCX_Exc-intermediate','MCX_Exc-deep','mFCX_Exc-deep',
                                        'MCX_Inh-VIP','mFCX_Inh-VIP','MCX_Inh-LAMP5','mFCX_Inh-LAMP5','MCX_Inh-PVALB','mFCX_Inh-PVALB','MCX_Inh-SST','mFCX_Inh-SST')))

df_sub<-df %>%
  filter(FDR<0.01)

ggplot(df,aes(x=types, y = description)) +
  geom_point(aes(size = enrichmentRatio, colour = logFDR)) +
  scale_color_viridis_c(name="-log10(FDR)",option="B",begin=0.1,end=0.8,limits=c(-log10(0.01),15),oob=squish,breaks=seq(2,15,length.out=5))+
  scale_size_continuous(name="Enrichment ratio", range=c(0,8),limits=c(1,12),breaks=seq(2,10,length.out=3))+
  theme_bw(base_size = 8, base_family = "Helvetica")+
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') + xlab('') +
  ggtitle("all FDR<0.01 Affinity")+
  theme(axis.ticks = element_blank()) +
  scale_x_discrete(drop = FALSE)



# Exc down --------------------------------------------------------------------

df<-read_csv("/cndd2/jchien/project/ALS/GO/GO_postCellBender/GO_summary_full_AP_ALS.csv")%>% 
  filter(reg =='dn') %>%
  filter(!grepl('Astro',celltype)) %>%
  mutate(types=str_replace(types, "_up", "")) %>%
  mutate(types=str_replace(types, "_dn", "")) 

df$logFDR <--log10(df$FDR)
df$logFDR[!is.finite(df$logFDR)] <- 15

GOIs <- df %>% filter(FDR<0.01) %>% filter(enrichmentRatio>0) %>% filter(Affinity=='YES') %>% filter(celltype %in% c('Exc-superficial','Exc-deep')) %>%
  arrange(desc(FDR)) %>% group_by(types) %>% slice(tail(row_number(), 50)) %>% pull(description)
df<-df %>% filter(description %in% GOIs)

mat <- df %>% 
  distinct()%>%
  select(types, description, FDR) %>%
  pivot_wider(names_from = types, values_from = FDR,values_fill = 1) %>%
  data.frame() 
row.names(mat) <- mat$description  
is.na(mat$geneSet) <-1
clust <- hclust(dist(mat %>% as.matrix()))
a=clust$labels[clust$order]
b=unique(df$description)
setdiff(b,a)
order_j=c(
  "synapse organization","postsynaptic specialization" ,
  "axon development",
  "regulation of neuron projection development",
  "negative regulation of nervous system development",
  "receptor complex",
  "passive transmembrane transporter activity",
  "morphogenesis of an epithelium",              
  "response to alcohol",                           
  "regulation of lipid metabolic process"
  )
df<-df %>%
  filter(FDR<0.25) %>%
  mutate(description = factor(description, levels=rev(order_j))) %>%
  mutate(types = factor(types, levels=c('MCX_Exc-superficial','mFCX_Exc-superficial','MCX_Exc-intermediate','mFCX_Exc-intermediate','MCX_Exc-deep','mFCX_Exc-deep',
                                        'MCX_Inh-VIP','mFCX_Inh-VIP','MCX_Inh-LAMP5','mFCX_Inh-LAMP5','MCX_Inh-PVALB','mFCX_Inh-PVALB','MCX_Inh-SST','mFCX_Inh-SST')))

ggplot(df,aes(x=types, y = description)) +
  geom_point(aes(size = enrichmentRatio, colour = logFDR)) +
  scale_color_viridis_c(name="-log10(FDR)",option="B",begin=0.1,end=0.8,limits=c(-log10(0.01),6),oob=squish,breaks=seq(2,6,length.out=3))+
  scale_size_continuous(name="Enrichment ratio", range=c(0,5),limits=c(1,12),breaks=seq(1,5,length.out=3))+
  theme_bw(base_size = 8, base_family = "Helvetica")+
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') + xlab('') +
  ggtitle("Downregulated, all FDR<0.01 Affinity")+
  theme(axis.ticks = element_blank()) +
  scale_x_discrete(drop = FALSE)

# Astro --------------------------------------------------------------------

df<-read_csv("/cndd2/jchien/project/ALS/GO/GO_postCellBender/GO_summary_full_AP_ALS.csv")%>% 
  filter(reg =='up') %>%
  filter(celltype %in% c('Astro','Micro','Oligo','OPC')) %>%
  mutate(types=str_replace(types, "_up", "")) %>%
  mutate(types=str_replace(types, "_dn", "")) 

df$logFDR <--log10(df$FDR)
df$logFDR[!is.finite(df$logFDR)] <- 15

GOIs <- df %>% filter(FDR<0.01) %>% filter(enrichmentRatio>0) %>% filter(Affinity=='YES')  %>% #filter(celltype=='Astro')  %>%
  arrange(desc(FDR)) %>% group_by(types) %>% slice(tail(row_number(), 50)) %>% pull(description)


order_j=c(
  "regulation of actin filament-based process",
  "actin cytoskeleton",
  
  "ameboidal-type cell migration",
  "regulation of chemotaxis",
  
  "cell-substrate adhesion",
  "cell junction organization" ,
  "cell-substrate junction", # non-astro
  
  "extracellular matrix",  
  "extracellular matrix structural constituent" , 
  
  "angiogenesis",   
  "transmembrane receptor protein serine/threonine kinase signaling pathway",
  "regulation of small GTPase mediated signal transduction" ,
  "organic hydroxy compound metabolic process" , 
  "cytokine binding",
  "regulation of body fluid levels" ,
  "extrinsic apoptotic signaling pathway" ,                                  
  "skeletal system morphogenesis", 
  "cell growth",  
  "coagulation",
  "regulation of peptide secretion",
  
  "second-messenger-mediated signaling" # non-astro
  
)


df<-df %>% filter(description %in% GOIs)

## hclust histograme
mat <- df %>% 
  distinct()%>%
  select(types, description, FDR) %>%
  pivot_wider(names_from = types, values_from = FDR,values_fill = 1) %>%
  data.frame()
row.names(mat) <- mat$description  
is.na(mat$geneSet) <-1
clust <- hclust(dist(mat %>% as.matrix())) 
a=clust$labels[clust$order]
b=unique(df$description)
setdiff(b,a)


df<-df %>%
  filter(FDR<0.25) %>%
  mutate(description = factor(description, levels=rev(order_j))) %>%
  mutate(types = factor(types, levels=c('MCX_Astro','mFCX_Astro','MCX_Micro','mFCX_Micro','MCX_Oligo','mFCX_Oligo','MCX_OPC','mFCX_OPC')))

p <- df %>% ggplot(aes(types, description))

p + geom_point(aes(color = -log10(FDR + 1e-15), size = enrichmentRatio)) +
  scale_color_viridis_c(
    name = "-log10(FDR)",
    option = "B",
    direction = 1,
    begin = 0.1,
    end = 0.9,
    breaks = seq(2, 8, length.out = 2),
    label = c("<2",  ">8"),
    limits = c(-log10(0.01), -log10(1e-8)),
    oob = squish
  ) +
  scale_size_continuous(
    name = "Enrichment ratio",
    range = c(0, 3),
    limits = c(1, 12),
    breaks = seq(1, 5, length.out = 3),
  ) +
  scale_x_discrete(drop = FALSE) +
  theme_bw(base_size = 6, base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.key.size = unit(0.5, "line")
  )
