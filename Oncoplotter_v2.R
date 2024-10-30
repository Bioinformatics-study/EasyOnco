library(maftools)
library(readxl)
library(dplyr)
library(scales)
library(circlize)
#------------------------------------------------------------------------------#
args = commandArgs(trailingOnly = TRUE)
output_file <- args[1]
# output_file = 'test.maf'
Name <- sub('.maf$','',output_file)
#------------------------------------------------------------------------------#
if (file.exists('Clinical_annotation.txt')){
  #------------------------------------------------------------------------------#
  maf <- read.table(output_file, sep='\t', header=TRUE)
  maf[maf == ""] <- NA
  # filtered_maf <- maf[apply(maf, 1, function(row) all(row != "." & !is.na(row))), ]
  Clinical_annotation <- read.table('Clinical_annotation.txt', sep='\t', header=T)
  current_date <- format(Sys.Date(), "%y%m%d")
  #------------------------------------------------------------------------------#
  VariantClassification = c("Frame_Shift_Del",
                            "Frame_Shift_Ins",
                            "Frame_Shift_Indel",
                            "Splice_Site",
                            "Translation_Start_Site",
                            "Nonsense_Mutation",
                            "Nonstop_Mutation",
                            "In_Frame_Del",
                            "In_Frame_Ins",
                            "In_Frame_Indel",
                            "Missense_Mutation",
                            "Start_Codon_Del",
                            "Stop_Codon_Del",
                            "Duplication",
                            "Intron_variant",
                            "5'_UTR_variant",
                            "3'_UTR_variant",
                            "Promoter_variant",
                            "Copy_Number_Gain",
                            "Copy_Number_Loss",
                            "Translocation",
                            "Inversion",
                            "Whole_gene_deletion",
                            "Partial_gene_deletion"
)
  #------------------------------------------------------------------------------#
  Onco <-  read.maf(maf = maf, 
                    verbose=T,
                    vc_nonSyn = VariantClassification,
                    clinicalData = Clinical_annotation)
  # Filtered_Onco <-  read.maf(maf = filtered_maf,
  #                   verbose=T,
  #                   vc_nonSyn = VariantClassification,
  #                   clinicalData = Clinical_annotation)
  #------------------------------------------------------------------------------#
  output_name <- paste0(Name, "_Onco")
  write.mafSummary(maf = Onco, basename = output_name)
  #------------------------------------------------------------------------------#
  options_df <- read_excel("oncoplot_options_v2.xlsx", col_names = TRUE)
  options_list <- setNames(as.list(options_df$Value), options_df$Options)
  get_option <- function(option_name) {
    return(options_list[[option_name]])
  }
  #------------------------------------------------------------------------------#
  Oncoprint_color <- c(
    Frame_Shift_Del = get_option("Frame_Shift_Del"),
    Frame_Shift_Ins = get_option("Frame_Shift_Ins"),
    Frame_Shift_Indel = get_option("Frame_Shift_Indel"),
    Splice_Site = get_option("Splice_Site"),
    Nonsense_Mutation = get_option("Nonsense_Mutation"),
    In_Frame_Del = get_option("In_Frame_Del"),
    In_Frame_Ins = get_option("In_Frame_Ins"),
    In_Frame_Indel = get_option("In_Frame_Indel"),
    Missense_Mutation = get_option("Missense_Mutation"),
    Multi_Hit = get_option("Multi_Hit"),
    Translation_Start_Site = get_option("Translation_Start_Site"),
    Nonstop_Mutation = get_option("Nonstop_Mutation"),
    Start_Codon_Del = get_option("Start_Codon_Del"),
    Stop_Codon_Del = get_option("Stop_Codon_Del"),
    Duplication = get_option("Duplication"),
    Intron_variant = get_option("Intron_variant"),
    `5'_UTR_variant` = get_option("5'_UTR_variant"),
    `3'_UTR_variant` = get_option("3'_UTR_variant"),
    Promoter_variant = get_option("Promoter_variant"),
    Copy_Number_Gain = get_option("Copy_Number_Gain"),
    Copy_Number_Loss = get_option("Copy_Number_Loss"),
    Translocation = get_option("Translocation"),
    Inversion = get_option("Inversion"),
    Whole_gene_deletion = get_option("Whole_gene_deletion"),
    Partial_gene_deletion = get_option("Partial_gene_deletion"))
  Oncoprint_color <- Oncoprint_color[!is.na(Oncoprint_color)]
  #------------------------------------------------------------------------------#
  anno_border_col <- get_option("annotation_border")
  anno_border_col <- ifelse(anno_border_col == "yes", "white", NA)
  #------------------------------------------------------------------------------#
  top_gene <- if (get_option("Top gene") == "All") NULL else as.numeric(get_option("Top gene"))
  if (is.null(top_gene)) {
    top_genes <- getGeneSummary(Onco)$Hugo_Symbol  
  } else {
    top_genes <- getGeneSummary(Onco)$Hugo_Symbol[1:top_gene]
  }
  GeneCut <- length(top_genes)
  #------------------------------------------------------------------------------#
  plot_title <- get_option("Plot title")
  #------------------------------------------------------------------------------#
  vaf_plot <- if (get_option("VAF plot") == "yes") "yes" else FALSE
  if(vaf_plot == "yes") {
    Genes <- subsetMaf(maf = Onco,
                       genes = top_genes,
                       fields = "i_TumorVAF_WU",
                       mafObj = FALSE)[,mean(i_TumorVAF_WU, na.rm = TRUE), Hugo_Symbol]
    colnames(Genes)[2] = "VAF (%)"
  } else {
    Genes <- NULL
  }
  #------------------------------------------------------------------------------#
  if (get_option("sample_cluster_by") != "no"){
    sample_order_option <-get_option("sample_cluster_by")
  } else {
    sample_order_option <- NULL
  }
  
  if(get_option("sample_cluster_order") != "no"){
    order_option <- get_option("sample_cluster_order")
  } else {
    order_option <- NULL
  }
  

  if (!is.null(sample_order_option)){
    if (!is.null(order_option)){
      sort_order <- strsplit(order_option, ',')
      sort_order <- sort_order[[1]]
      sort_order <- gsub(' ', '', sort_order)
      
      test <- as.vector(Clinical_annotation[,sample_order_option])
      test<- factor(Clinical_annotation[,sample_order_option],levels=sort_order)
      canno_sorted <-  Clinical_annotation[order(test),]
      sample_order <- canno_sorted$Tumor_Sample_Barcode
    } else{
      sample_order <- Clinical_annotation %>%
        arrange(.[[sample_order_option]]) %>%
        pull(Tumor_Sample_Barcode)
    }
  } else{
    sample_order <- NULL
  }

  #------------------------------------------------------------------------------#
  annotation_colors <- options_df %>%
    dplyr::filter(grepl("Annotation color", options_list)) %>%
    pull(Value)
  clinical_features <- colnames(Clinical_annotation)[-1]
  categorical_features <- options_df$Value[37:41]
  
  continuous_features <- options_df$Value[42:46]
  continuous_features <- continuous_features[!continuous_features %in% c("NA", NA)]
  
  remove_na <- function(vec) {
    vec[vec != "NA"]
  }
  categorical_colors <- options_df$Value[47:96]
  categorical_groups <- split(categorical_colors, ceiling(seq_along(categorical_colors) / 10))
  categorical_groups <- lapply(categorical_groups, remove_na)
  
  valid_categorical_features <- categorical_features[!categorical_features %in% c("NA", NA)]
  
  continuous_colors <- options_df$Value[97:106]
  continuous_groups <- split(continuous_colors, ceiling(seq_along(continuous_colors) / 2))
  continuous_groups <- lapply(continuous_groups, remove_na)
  
  valid_colors <- continuous_colors[continuous_colors != "NA"]
  
  Annotation_color <- list()
  color_index <- 1
  
  for (feature in continuous_features) {
    Onco@clinical.data[[feature]] <- as.numeric(Onco@clinical.data[[feature]])
    valid_feature_values <- Clinical_annotation[[feature]]
    
    if (!is.null(valid_feature_values) && length(na.omit(valid_feature_values)) > 0) {
      color_pair <- continuous_colors[(color_index):(color_index + 1)]
      num_steps <- max(valid_feature_values, na.rm = TRUE) - min(valid_feature_values, na.rm = TRUE) + 1
      feature_gradient <- colorRampPalette(color_pair)(num_steps)
      
      feature_breaks <- seq(min(valid_feature_values, na.rm = TRUE), 
                            max(valid_feature_values, na.rm = TRUE), 
                            length.out = num_steps)
      
      Annotation_color[[feature]] <- setNames(feature_gradient, as.character(round(feature_breaks)))
      color_index <- color_index + 2
    }
  }
  
  for (i in seq_along(valid_categorical_features)) {
    feature <- valid_categorical_features[i]
    group_colors <- categorical_groups[[as.character(i)]]
    unique_values <- unique(Clinical_annotation[[feature]])
    
    if (length(group_colors) >= length(unique_values)) {
      Annotation_color[[feature]] <- setNames(group_colors[1:length(unique_values)], unique_values)
    } 
  }
  #------------------------------------------------------------------------------#
  if (tolower(options_list$`Top gene`) == 'all') {
    Oncogenes <- Onco@gene.summary$Hugo_Symbol
  } else {
    top_gene_num <- as.numeric(options_list$'Top gene')
    Oncogenes <- Onco@gene.summary$Hugo_Symbol[1:top_gene_num]
  }
  GeneCut <- length(Oncogenes)
  #------------------------------------------------------------------------------#
  Onco.titv = titv(maf = Onco, plot = FALSE, useSyn = TRUE) #missense만 그려야 경고 안뜸
  oncotitiv_filename <- paste0(current_date,"_",Name, "_onco.titv.pdf")
  pdf(oncotitiv_filename)
  plotTiTv(res = Onco.titv)
  dev.off()
  #------------------------------------------------------------------------------#
  current_date <- format(Sys.Date(), "%y%m%d")
  oncoplot_filename <- paste0(current_date,"_",Name, "_oncoplot.pdf")
  if(get_option("removeNonMutated") == "yes"){
    removeNonMutated <- TRUE
  } else{
    removeNonMutated <- FALSE
  }
  
  if(get_option("writeMatrix") == "yes"){
    writeMatrix <- TRUE
  } else{
    writeMatrix <- FALSE
  }
  
  if(get_option("showTumorSampleBarcodes") == "yes"){
    showTumorSampleBarcodes <- TRUE
  } else{
    showTumorSampleBarcodes <- FALSE
  }
  
  pdf(oncoplot_filename)
  oncoplot(maf = Onco,
           genes = top_genes,
           colors = Oncoprint_color,
           leftBarLims = c(0,100),
           fontSize = as.numeric(get_option("fontSize")),
           removeNonMutated <- removeNonMutated,
           writeMatrix = writeMatrix,
           showTumorSampleBarcodes = showTumorSampleBarcodes,
           legendFontSize = 1,
           annoBorderCol=anno_border_col,
           annotationFontSize = 1,
           anno_height = 2,
           sampleOrder = sample_order,
           titleText = plot_title,
           annotationColor = Annotation_color,
           clinicalFeatures = clinical_features,
           top = GeneCut)
  dev.off()
  #------------------------------------------------------------------------------#
  mafsummary_filename <- paste0(current_date,"_",Name, "_mafsummary.pdf")
  pdf(mafsummary_filename)
  plotmafSummary(maf = Onco,color = Oncoprint_color, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
  dev.off()
  # #------------------------------------------------------------------------------#
  maf_lolli <- read.table(output_file, sep='\t', header=T)
  gene <- subset(maf_lolli, maf_lolli$Variant_Classification!="Intron_variant"
                 & maf_lolli$Variant_Classification !="Silent_Mutation"
                 & maf_lolli$Variant_Classification != 'Splice_Site'
                 & maf_lolli$Variant_Classification != "5'_UTR_variant"
                 & maf_lolli$Variant_Classification != "3'_UTR_variant"
                 & maf_lolli$Variant_Classification != 'Promoter_variant'
                 & maf_lolli$Variant_Classification != "")
  Total_gene <- unique(gene$Hugo_Symbol)

  Onco_lolli <-  read.maf(maf = gene,
                          verbose=T,
                          vc_nonSyn = VariantClassification,
                          clinicalData = Clinical_annotation)

  maftools_path <- find.package("maftools")
  protein_db_path <- file.path(maftools_path,"extdata","protein_domains.RDs")
  proteindb <- readRDS(protein_db_path)
  DBGene <- unique(proteindb$HGNC)
  Total_gene <- intersect(DBGene, Total_gene)
  # #------------------------------------------------------------------------------#
  lollipopplot_filename <- paste0(current_date,"_",Name, "_lollipop.pdf")
  pdf(lollipopplot_filename)
  for (Gene in Total_gene) {
    lollipopPlot(maf = Onco_lolli,
                 gene = Gene,
                 AACol = NULL,
                 showMutationRate= TRUE,
                 showDomainLabel = FALSE,
                 labelPos = NULL,
                 colors = Oncoprint_color)
  }
  dev.off()
  # #------------------------------------------------------------------------------#
  vafplot_filename <- paste0(current_date,"_",Name, "_vafplot.pdf")
  pdf(vafplot_filename)
  plotVaf(maf = Onco, vafCol = 'i_TumorVAF_WU', width = 10, height = 10)
  dev.off()
  # #------------------------------------------------------------------------------#
  drugint_filename <- paste0(current_date,"_",Name, "_druginteraction.pdf")
  pdf(drugint_filename)
  drugInteractions(maf = Onco, fontSize = 0.75)
  dev.off()
  #------------------------------------------------------------------------------#
} else {
  #------------------------------------------------------------------------------#
  maf <- read.table(output_file, sep='\t', header=TRUE)
  maf[maf == ""] <- NA
  # filtered_maf <- maf[apply(maf, 1, function(row) all(row != "." & !is.na(row))), ]
  # Clinical_annotation <- read.table('Clinical_annotation.txt', sep='\t', header=T)
  current_date <- format(Sys.Date(), "%y%m%d")
  #------------------------------------------------------------------------------#
  VariantClassification = c("Frame_Shift_Del",
                            "Frame_Shift_Ins",
                            "Frame_Shift_Indel",
                            "Splice_Site",
                            "Translation_Start_Site",
                            "Nonsense_Mutation",
                            "Nonstop_Mutation",
                            "In_Frame_Del",
                            "In_Frame_Ins",
                            "In_Frame_Indel",
                            "Missense_Mutation",
                            "Start_Codon_Del",
                            "Stop_Codon_Del",
                            "Duplication",
                            "Intron_variant",
                            "5'_UTR_variant",
                            "3'_UTR_variant",
                            "Promoter_variant",
                            "Copy_Number_Gain",
                            "Copy_Number_Loss",
                            "Translocation",
                            "Inversion",
                            "Whole_gene_deletion",
                            "Partial_gene_deletion"
  )
  #------------------------------------------------------------------------------#
  Onco <-  read.maf(maf = maf, 
                    verbose=T,
                    vc_nonSyn = VariantClassification)
  # Filtered_Onco <-  read.maf(maf = filtered_maf, 
  #                   verbose=T,
  #                   vc_nonSyn = VariantClassification)
  #------------------------------------------------------------------------------#
  output_name <- paste0(Name, "_Onco")
  write.mafSummary(maf = Onco, basename = output_name)
  #------------------------------------------------------------------------------#
  options_df <- read_excel("oncoplot_options_v2.xlsx", col_names = TRUE)
  options_list <- setNames(as.list(options_df$Value), options_df$Options)
  get_option <- function(option_name) {
    return(options_list[[option_name]])
  }
  #------------------------------------------------------------------------------#
  Oncoprint_color <- c(
    Frame_Shift_Del = get_option("Frame_Shift_Del"),
    Frame_Shift_Ins = get_option("Frame_Shift_Ins"),
    Frame_Shift_Indel = get_option("Frame_Shift_Indel"),
    Splice_Site = get_option("Splice_Site"),
    Nonsense_Mutation = get_option("Nonsense_Mutation"),
    In_Frame_Del = get_option("In_Frame_Del"),
    In_Frame_Ins = get_option("In_Frame_Ins"),
    In_Frame_Indel = get_option("In_Frame_Indel"),
    Missense_Mutation = get_option("Missense_Mutation"),
    Multi_Hit = get_option("Multi_Hit"),
    Translation_Start_Site = get_option("Translation_Start_Site"),
    Nonstop_Mutation = get_option("Nonstop_Mutation"),
    Start_Codon_Del = get_option("Start_Codon_Del"),
    Stop_Codon_Del = get_option("Stop_Codon_Del"),
    Duplication = get_option("Duplication"),
    Intron_variant = get_option("Intron_variant"),
    `5'_UTR_variant` = get_option("5'_UTR_variant"),
    `3'_UTR_variant` = get_option("3'_UTR_variant"),
    Promoter_variant = get_option("Promoter_variant"),
    Copy_Number_Gain = get_option("Copy_Number_Gain"),
    Copy_Number_Loss = get_option("Copy_Number_Loss"),
    Translocation = get_option("Translocation"),
    Inversion = get_option("Inversion"),
    Whole_gene_deletion = get_option("Whole_gene_deletion"),
    Partial_gene_deletion = get_option("Partial_gene_deletion"))
  Oncoprint_color <- Oncoprint_color[!is.na(Oncoprint_color)]
  #------------------------------------------------------------------------------#
  anno_border_col <- get_option("annotation_border")
  anno_border_col <- ifelse(anno_border_col == "yes", "white", NA)
  #------------------------------------------------------------------------------#
  top_gene <- if (get_option("Top gene") == "All") NULL else as.numeric(get_option("Top gene"))
  if (is.null(top_gene)) {
    top_genes <- getGeneSummary(Onco)$Hugo_Symbol  
  } else {
    top_genes <- getGeneSummary(Onco)$Hugo_Symbol[1:top_gene]
  }
  GeneCut <- length(top_genes)
  #------------------------------------------------------------------------------#
  plot_title <- get_option("Plot title")
  #------------------------------------------------------------------------------#
  vaf_plot <- if (get_option("VAF plot") == "yes") "yes" else FALSE
  if(vaf_plot == "yes") {
    Genes <- subsetMaf(maf = Onco,
                       genes = top_genes,
                       fields = "i_TumorVAF_WU",
                       mafObj = FALSE)[,mean(i_TumorVAF_WU, na.rm = TRUE), Hugo_Symbol]
    colnames(Genes)[2] = "VAF (%)"
  } else {
    Genes <- NULL
  }
  #------------------------------------------------------------------------------#
  if (tolower(options_list$`Top gene`) == 'all') {
    Oncogenes <- Onco@gene.summary$Hugo_Symbol
  } else {
    top_gene_num <- as.numeric(options_list$'Top gene')
    Oncogenes <- Onco@gene.summary$Hugo_Symbol[1:top_gene_num]
  }
  GeneCut <- length(Oncogenes)
  #------------------------------------------------------------------------------#
  Onco.titv = titv(maf = Onco, plot = FALSE, useSyn = TRUE)
  oncotitiv_filename <- paste0(current_date,"_",Name, "_onco.titv.pdf")
  pdf(oncotitiv_filename)
  plotTiTv(res = Onco.titv)
  dev.off()
  #------------------------------------------------------------------------------#
  current_date <- format(Sys.Date(), "%y%m%d")
  oncoplot_filename <- paste0(current_date,"_",Name, "_oncoplot.pdf")
  if(get_option("removeNonMutated") == "yes"){
    removeNonMutated <- TRUE
  } else{
    removeNonMutated <- FALSE
  }
  
  if(get_option("writeMatrix") == "yes"){
    writeMatrix <- TRUE
  } else{
    writeMatrix <- FALSE
  }
  
  if(get_option("showTumorSampleBarcodes") == "yes"){
    showTumorSampleBarcodes <- TRUE
  } else{
    showTumorSampleBarcodes <- FALSE
  }
  
  pdf(oncoplot_filename)
  oncoplot(maf = Onco,
           genes = top_genes,
           colors = Oncoprint_color,
           leftBarLims = c(0,100),
           fontSize = as.numeric(get_option("fontSize")),
           removeNonMutated <- removeNonMutated,
           writeMatrix = writeMatrix,
           showTumorSampleBarcodes = showTumorSampleBarcodes,
           legendFontSize = 1,
           annoBorderCol=anno_border_col,
           annotationFontSize = 1,
           anno_height = 2,
           titleText = plot_title,
           annotationColor = Annotation_color,
           top = GeneCut)
  dev.off()
  #------------------------------------------------------------------------------#
  mafsummary_filename <- paste0(current_date,"_",Name, "_mafsummary.pdf")
  pdf(mafsummary_filename)
  plotmafSummary(maf = Onco,color = Oncoprint_color, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
  dev.off()
  # #------------------------------------------------------------------------------#
  maf_lolli <- read.table(output_file, sep='\t', header=T)
  gene <- subset(maf_lolli, maf_lolli$Variant_Classification!="Intron_variant"
                 & maf_lolli$Variant_Classification !="Silent_Mutation"
                 & maf_lolli$Variant_Classification != 'Splice_Site'
                 & maf_lolli$Variant_Classification != "5'_UTR_variant"
                 & maf_lolli$Variant_Classification != "3'_UTR_variant"
                 & maf_lolli$Variant_Classification != 'Promoter_variant'
                 & maf_lolli$Variant_Classification != "")
  Total_gene <- unique(gene$Hugo_Symbol)

  Onco_lolli <-  read.maf(maf = gene,
                          verbose=T,
                          vc_nonSyn = VariantClassification)

  maftools_path <- find.package("maftools")
  protein_db_path <- file.path(maftools_path,"extdata","protein_domains.RDs")
  proteindb <- readRDS(protein_db_path)
  DBGene <- unique(proteindb$HGNC)
  Total_gene <- intersect(DBGene, Total_gene)
  # #------------------------------------------------------------------------------#
  lollipopplot_filename <- paste0(current_date,"_",Name, "_lollipop.pdf")
  pdf(lollipopplot_filename)
  for (Gene in Total_gene) {
    lollipopPlot(maf = Onco_lolli,
                 gene = Gene,
                 AACol = NULL,
                 showMutationRate = TRUE,
                 showDomainLabel = FALSE,
                 labelPos = NULL,
                 colors = Oncoprint_color)
  }
  dev.off()
  #------------------------------------------------------------------------------#
  vafplot_filename <- paste0(current_date,"_",Name, "_vafplot.pdf")
  pdf(vafplot_filename)
  plotVaf(maf = Onco, vafCol = 'i_TumorVAF_WU', width = 10, height = 10)
  dev.off()
  #------------------------------------------------------------------------------#
  drugint_filename <- paste0(current_date,"_",Name, "_druginteraction.pdf")
  pdf(drugint_filename)
  drugInteractions(maf = Onco, fontSize = 0.75)
  dev.off()
  #------------------------------------------------------------------------------# 
}
