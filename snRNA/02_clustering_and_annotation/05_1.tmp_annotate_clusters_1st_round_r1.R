
DotPlot(data_obj_sub,
        cols = c("lightgrey", "#C73D4C"),
        features = c("SNAP25", "MEF2C", "SLC17A7", "SATB2", "GAD1", "GAD2"
                     # "SOX10", "RBFOX3", "DCX", "MAP2", "TUBB3", "ENO2"
        )) +
  theme_bw(base_size = 7, base_family = "Helvetica") +
  xlab("Marker genes") +
  ylab("Cell class") +
  scale_x_discrete(position = "top") +
  theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0),
        legend.position = "bottom")

FeaturePlot(data_obj_sub,
        features = c("SNAP25", "MEF2C", "SLC17A7", "SATB2", "GAD1", "GAD2")
            )



DotPlot(data_obj_sub,
        cols = c("lightgrey", "#C73D4C"),
        
        features = c(
          "AQP4", "APOE", "GFAP", "SERPINI2", "AQP1", "PLCG1", # Astro
          "PECAM1", "VWF", "CLDN5", "NOSTRIN",  # Endo
          "C1QB", "CSF1R", # Micro
          "MOG", "MBP", "MOBP", "CNP", "ENPP6", "OPALIN", "SNAP25", "SLC17A7",  # Oligo
          "PDGFRA", "VCAN", # OPC
          "TBX18", "COLEC12", "CEMIP", "CYP1B1", "P2RY14" # VLMC
        )
) +
  theme_bw(base_size = 7, base_family = "Helvetica") +
  xlab("Marker genes") +
  ylab("Cell class") +
  scale_x_discrete(position = "top") +
  theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0),
        legend.position = "bottom")

FeaturePlot(data_obj_sub,
            features = c(
              "AQP4", "APOE", "GFAP", "SERPINI2", "AQP1", "PLCG1", # Astro
              "PECAM1", "VWF", "CLDN5", "NOSTRIN",  # Endo
              "C1QB", "CSF1R", # Micro
              "MOG", "MBP", "MOBP", "CNP", "ENPP6", "OPALIN", "SNAP25", "SLC17A7",  # Oligo
              "PDGFRA", "VCAN", # OPC
              "TBX18", "COLEC12", "CEMIP", "CYP1B1", "P2RY14" # VLMC
            )
)


DotPlot(data_obj_sub,
        cols = c("lightgrey", "#C73D4C"),
        
        features = c(
        "ACTA2", "MYH11", "TAGLN",
        "KCNMA1", "SLC4A4",
        "VWF", "PECAM1", "CLDN5", "NOSTRIN")
) +
  theme_bw(base_size = 7, base_family = "Helvetica") +
  xlab("Marker genes") +
  ylab("Cell class") +
  scale_x_discrete(position = "top") +
  theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0),
        legend.position = "bottom")


DotPlot(data_obj_sub,
        cols = c("lightgrey", "#C73D4C"),
        features = c(
          "LHX6", "SOX6",
          "PVALB", "CUX2", "MYBPC1", "PIEZO2", "SCUBE3",
          "SST", "NPY", "MAFB", "SPHKAP", "CDH12", "KLHL14",
          "ADARB2",
          "LAMP5", "SV2C", "CHST9", "CPLX3", "KIT", "CXCL14", "NDNF", "RELN", "PAX6", "FAM19A1",
          "VIP", "ABI3BP", "CLSTN2", "DACH2", "FLT1",  "ZBTB20"
        )) +
  theme_bw(base_size = 7, base_family = "Helvetica") +
  xlab("Marker genes") +
  ylab("Cell class") +
  scale_x_discrete(position = "top") +
  theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0),
        legend.position = "bottom")


DotPlot(data_obj_sub,
        cols = c("lightgrey", "#C73D4C"),
        features = c("SLC17A7",
                     "CUX2",
                     "ATP7B", "ACVR1C",
                     "LAMP5", "SERPINE2", "LINC00507",
                     "PDGFD", "CCBE1",
                     "PRKG2","GLIS3",
                     "SV2C", "SYT2",
                     "RORB",
                     "COBLL1", "PLCH1", "PRSS12", "CCDC68",
                     "OTOGL", "COL22A1", "ALDH1A1", "GABRG1", "NPFFR2",
                     "LRRK1", "TRABD2A", "RPRM",
                     "THEMIS", "SMYD1",
                     "LINC00299", "CFLAR",
                     "TLE4", "PCSK5", "MDFIC", "KCNK2",
                     "SULF1", "SEMA5A", "FEZF2", "HTR2C", "KCNIP1",
                     "CRYM", "ADRA1A", "MYO16", "NEFH", "POU3F1"
        )) +
  theme_bw(base_size = 7, base_family = "Helvetica") +
  xlab("Marker genes") +
  ylab("Cell class") +
  scale_x_discrete(position = "top") +
  theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0),
        legend.position = "bottom")


# check cluster 29 (probably NK cells? or something related to immune cells; also expressed some microglia genes, but not C1QB or CSF1R)
g1 <- WhichCells(data_obj_sub, expression = SCT_snn_res.1 == "29")
DimPlot(data_obj_sub, label=T,  cells.highlight= list(g1), cols.highlight = c( "darkred"), cols= "grey", raster = FALSE)

markers_cluster <- FindMarkers(data_obj_sub,
                               ident.1 = 29,
                               min.pct = 0.25)
cl_29_markers <- markers_cluster %>% 
  rownames_to_column("gene") %>% 
  as_tibble()

FeaturePlot(data_obj_sub,
            features = c(
             "PTPRC", "ARHGAP15", "SKAP1", "PARP8", "AOAH", "CD247", "PRKCH", "IKZF1" ,"SLFN12L"
            )
)

DotPlot(data_obj_sub,
        cols = c("lightgrey", "#C73D4C"),
          features = c(
             "PTPRC", "ARHGAP15", "SKAP1", "PARP8", "AOAH", "CD247", "PRKCH", "IKZF1" ,"SLFN12L", "C1QB", "CSF1R"
            )
        ) +
  theme_bw(base_size = 7, base_family = "Helvetica") +
  xlab("Marker genes") +
  ylab("Cell class") +
  scale_x_discrete(position = "top") +
  theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0),
        legend.position = "bottom")
write_tsv(cl_29_markers, "marker_genes_1st_round_cluster_29_r1.txt")


# check cluster 20 (exc neurons?)
g1 <- WhichCells(data_obj_sub, expression = SCT_snn_res.1 == "20")
DimPlot(data_obj_sub, label=T,  cells.highlight= list(g1), cols.highlight = c( "darkred"), cols= "grey", raster = FALSE)

markers_cluster <- FindMarkers(data_obj_sub,
                               ident.1 = 20,
                               min.pct = 0.25)
cl_20_markers <- markers_cluster %>% 
  rownames_to_column("gene") %>% 
  as_tibble()

write_tsv(cl_20_markers, "marker_genes_1st_round_cluster_20_r1.txt")


FeaturePlot(data_obj_sub,
            features = c(
             "MAP1B", "CALM1", "NRGN", "NEFL", "NEFM", "TUBA1B", "SPARCL1", "THY1",
             "VSNL1", "TUBB2A", "CCK", "RGS4", "YWHAH", "BEX1", 
             "ENC1", "CALM3", "STMN2", "UCHL1", "GAPDH"
            )
)

FeaturePlot(data_obj_sub,
            features = c(
              
              "CCDC85B","RAB3A", "NRN1", "CHGB", "MOAP1", "TPI1","TPRG1L", "PINK1", "STMN1", "NCDN", "LMO4", "ATP1A1", "CLSTN1",
              "OLFM1", "CREG2", "TSPYL1", "GABRA1", "TSPAN7", "YWHAH", "BEX1", "NDRG4", "PCSK1N", "SYP", "BEX3", "ATP6V1B2"
           
            )
)


DotPlot(data_obj_sub,
        cols = c("lightgrey", "#C73D4C"),
          features = c(
              "MAP1B", "CALM1", "NRGN", "NEFL", "NEFM", "TUBA1B", "SPARCL1", "THY1",
             "VSNL1", "TUBB2A", "CCK", "RGS4", "YWHAH", "BEX1", 
             "ENC1", "CALM3", "STMN2", "UCHL1", "GAPDH"
            )
        ) +
  theme_bw(base_size = 7, base_family = "Helvetica") +
  xlab("Marker genes") +
  ylab("Cell class") +
  scale_x_discrete(position = "top") +
  theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, hjust = 0),
        legend.position = "bottom")