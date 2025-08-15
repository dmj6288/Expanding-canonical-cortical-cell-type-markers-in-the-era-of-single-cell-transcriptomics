library(hash)

Hodge_markers     <- list( "AQP4",                                    "SLC17A7",                                            "GAD1",
                           "TYROBP",                                "OPALIN",                 "PDGFRA")
Lake_markers      <- list(c("SLC1A2", "SLC4A4"),                    c("SLC17A7", "SATB2"),                                  "GAD1",
                          "APBB1IP",                               "MOBP",                   "PCDH15")

Brenner_markers   <- list(c("GFAP", "ALDH1L1", "SLC1A2", "SLC4A4"), c("SLC17A7", "SATB2", "CAMK2A"),                      c("GAD1", "GAD2"),
                          c("P2RY12", "CSF1R", "APBB1IP"),          c("MBP", "MOBP", "PLP1"), c("VCAN", "PDGFRA", "PCDH15"))

Khrameeva_markers <- list( "GJA1",                                  c("SLC17A7", "SATB2"),                                c("GAD1", "GAD2"),
                           c("AIF1", "CX3CR1", "PTPRC", "HLA-DRA"), c("MBP", "MOBP", "MOG"),  c("PDGFRA", "CSPG4"))

Yang_markers      <- list( "SLC1A2",                                c("CBLN2", "IL1RAPL2", "TRPM3", "SEMA3E"),            c("GRIK2", "KCNC2", "CDH9", "KIT"),
                           "DOCK8",                                 "ST18",                   "VCAN")

Velmeshev_markers <- list( c("SLC1A2", "GFAP"),                     c("SYT1", "RBFOX3", "CUX2", "SATB2", "RORB", "TLE4"), c("GAD1", "GAD2", "PVALB", "SST", "VIP", "SV2C"),
                           "PTPRC",                                 "PLP1",                   "PDGFRA")

Bakken_markers    <- list( "AQP4",                                  c("SLC17A7", "SATB2"),                                c("GAD1", "GAD2"),
                           "APBB1IP",                             c("PLP1", "MOBP"),          "PDGFRA")

list_of_study_wise_markers <- list(Bakken_markers, Brenner_markers, Hodge_markers, Khrameeva_markers, Lake_markers, Velmeshev_markers, Yang_markers)

study_names <- c("Bakken", "Brenner", "Hodge", "Khrameeva", "Lake", "Velmeshev", "Yang")
cell_names  <- c("Astro", "Excite", "Inhibit", "Micro", "Oligo", "OPC")

study_name_to_marker_list_hash <- hash(study_names, list_of_study_wise_markers)

saveRDS(study_name_to_marker_list_hash, "interm/canonical_marker_hash.rds")

canonical_marker_count           <- list()
canonical_marker_superlist_count <- list()
canonical_marker_lists_test      <- list()

for (i in c(1:6)){
  
  cell_marker_length <- c()
  cell_marker_superlist <- c()
  
  for (study in names(study_name_to_marker_list_hash)){
    
    cell_marker_length <- c(cell_marker_length, length(study_name_to_marker_list_hash[[study]][[i]]))
    cell_marker_superlist <- c(cell_marker_superlist, study_name_to_marker_list_hash[[study]][[i]])
    
  }
  
  canonical_marker_count[[i]] <- round(mean(cell_marker_length))
  canonical_marker_lists_test[[i]] <- sort(unique(cell_marker_superlist))
  canonical_marker_superlist_count[[i]] <- length(unique(cell_marker_superlist))
  
}

names(canonical_marker_lists_test) <- cell_names
names(canonical_marker_superlist_count) <- cell_names

canonical_marker_count
saveRDS(canonical_marker_count,           "interm/canonical_marker_count.rds")
saveRDS(canonical_marker_superlist_count, "interm/canonical_marker_superlist_count.rds")
saveRDS(canonical_marker_lists_test,      "interm/canonical_marker_lists.rds")

