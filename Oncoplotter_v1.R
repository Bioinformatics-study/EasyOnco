library(maftools)
library(readxl)
library(dplyr)
library(scales)
library(circlize)
#------------------------------------------------------------------------------#
if (file.exists('Clinical_annotation.txt')){
  #------------------------------------------------------------------------------#
  maf <- read.table('EasyOnco.maf', sep='\t', header=T)
  Clinical_annotation <- read.table('Clinical_annotation.txt', sep='\t', header=T)
  current_date <- format(Sys.Date(), "%y%m%d")
  #------------------------------------------------------------------------------#
  VariantClassification = c("Frame_Shift_Del",
                            "Frame_Shift_Ins",
                            "Splice_Site",
                            "Translation_Start_Site",
                            "Nonsense_Mutation",
                            "Nonstop_Mutation",
                            "In_Frame_Del",
                            "In_Frame_Ins",
                            "Missense_Mutation",
                            "Start_Codon_Del",
                            "Stop_Codon_Del",
                            "Duplication",
                            "Intron_variant")
  #------------------------------------------------------------------------------#
  Onco <-  read.maf(maf = maf, 
                    verbose=T,
                    vc_nonSyn = VariantClassification,
                    clinicalData = Clinical_annotation)
  #------------------------------------------------------------------------------#
  write.mafSummary(maf = Onco, basename = 'Onco')
  #------------------------------------------------------------------------------#
  options_df <- read_excel("../oncoplot_options_v1.xlsx", col_names = TRUE)
  options_list <- setNames(as.list(options_df$Value), options_df$Options)
  get_option <- function(option_name) {
    return(options_list[[option_name]])
  }
  #------------------------------------------------------------------------------#
  Oncoprint_color <- c(
    Frame_Shift_Del = get_option("Frame_Shift_Del"),
    Frame_Shift_Ins = get_option("Frame_Shift_Ins"),
    Splice_Site = get_option("Splice_Site"),
    Nonsense_Mutation = get_option("Nonsense_Mutation"),
    In_Frame_Del = get_option("In_Frame_Del"),
    In_Frame_Ins = get_option("In_Frame_Ins"),
    Missense_Mutation = get_option("Missense_Mutation"),
    Multi_Hit = get_option("Multi_Hit"),
    Translation_Start_Site = get_option("Translation_Start_Site"),
    Nonstop_Mutation = get_option("Nonstop_Mutation"),
    Start_Codon_Del = get_option("Start_Codon_Del"),
    Stop_Codon_Del = get_option("Stop_Codon_Del"),
    Duplication = get_option("Duplication"),
    Intron_variant = get_option("Intron_variant"))
  Oncoprint_color <- Oncoprint_color[!is.na(Oncoprint_color)]
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
  sample_order_option <- get_option("sample_cluster_by")
  
  if (!is.null(sample_order_option) && sample_order_option %in% colnames(Clinical_annotation)) {
    sample_order <- Clinical_annotation %>%
      arrange(.[[sample_order_option]]) %>%
      pull(Tumor_Sample_Barcode)
  } else {
    sample_order <- NULL
  }
  #------------------------------------------------------------------------------#
  annotation_colors <- options_df %>%
    dplyr::filter(grepl("Annotation color", Options)) %>%
    pull(Value)
  clinical_features <- colnames(Clinical_annotation)[-1]
  categorical_features <- options_df$Value[25:29]
  
  continuous_features <- options_df$Value[30:34]
  continuous_features <- continuous_features[!continuous_features %in% c("NA", NA)]
  
  remove_na <- function(vec) {
    vec[vec != "NA"]
  }
  categorical_colors <- options_df$Value[35:84]
  categorical_groups <- split(categorical_colors, ceiling(seq_along(categorical_colors) / 10))
  categorical_groups <- lapply(categorical_groups, remove_na)
  
  valid_categorical_features <- categorical_features[!categorical_features %in% c("NA", NA)]
  
  continuous_colors <- options_df$Value[85:94]
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
    top_gene_num <- as.numeric(Option$Top_gene)
    Oncogenes <- Onco@gene.summary$Hugo_Symbol[1:top_gene_num]
  }
  GeneCut <- length(Oncogenes)
  #------------------------------------------------------------------------------#
  Onco.titv = titv(maf = Onco, plot = FALSE, useSyn = TRUE)
  oncotitiv_filename <- paste0(current_date, "_onco.titv.pdf")
  pdf(oncotitiv_filename)
  plotTiTv(res = Onco.titv)
  dev.off()
  #------------------------------------------------------------------------------#
  current_date <- format(Sys.Date(), "%y%m%d")
  oncoplot_filename <- paste0(current_date, "_oncoplot.pdf")
  pdf(oncoplot_filename)
  oncoplot(maf = Onco,
           genes = top_genes,
           colors = Oncoprint_color,
           leftBarLims = c(0,100),
           leftBarData = Genes,
           SampleNamefontSize = as.numeric(get_option("SampleNamefontSize")), 
           removeNonMutated = as.logical(get_option("removeNonMutated")),
           writeMatrix = as.logical(get_option("writeMatrix")),
           showTumorSampleBarcodes = as.logical(get_option("showTumorSampleBarcodes")),
           legendFontSize = 1,
           annotationFontSize = 1,
           anno_height = 2,
           sampleOrder = sample_order,
           titleText = plot_title,
           annotationColor = Annotation_color,
           clinicalFeatures = clinical_features,
           top = GeneCut)
  dev.off()
  #------------------------------------------------------------------------------#
  mafsummary_filename <- paste0(current_date, "_mafsummary.pdf")
  pdf(mafsummary_filename)
  plotmafSummary(maf = Onco,color = Oncoprint_color, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
  dev.off()
  #------------------------------------------------------------------------------#
  maf_lolli <- read.table('EasyOnco.maf', sep='\t', header=T)
  gene <- subset(maf_lolli, maf_lolli$Variant_Classification!="Intron_variant" & maf_lolli$Variant_Classification !="Silent_Mutation")
  Intron_gene <- unique(gene$Hugo_Symbol)
  # idx <- which(maf_lolli$Hugo_Symbol %in% Intron_gene)
  # maf_lolli <- maf_lolli[-c(idx), ]
  # lolli_gene <- unique(maf_lolli$Hugo_Symbol)
  
  Onco_lolli <-  read.maf(maf = maf_lolli, 
                          verbose=T,
                          vc_nonSyn = VariantClassification,
                          clinicalData = Clinical_annotation)
  
  #------------------------------------------------------------------------------#
  lollipopplot_filename <- paste0(current_date, "_lollipop.pdf")
  pdf(lollipopplot_filename)
  for (Gene in Intron_gene) {
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
  vafplot_filename <- paste0(current_date, "_vafplot.pdf")
  pdf(vafplot_filename)
  plotVaf(maf = Onco, vafCol = 'i_TumorVAF_WU', width = 10, height = 10)
  dev.off()
  #------------------------------------------------------------------------------#
  drugint_filename <- paste0(current_date, "_druginteraction.pdf")
  pdf(drugint_filename)
  drugInteractions(maf = Onco, fontSize = 0.75)
  dev.off()
  #------------------------------------------------------------------------------#
} else {
  #------------------------------------------------------------------------------#
  maf <- read.table('EasyOnco.maf', sep='\t', header=T)
  # Clinical_annotation <- read.table('Clinical_annotation.txt', sep='\t', header=T)
  current_date <- format(Sys.Date(), "%y%m%d")
  #------------------------------------------------------------------------------#
  VariantClassification = c("Frame_Shift_Del",
                            "Frame_Shift_Ins",
                            "Splice_Site",
                            "Translation_Start_Site",
                            "Nonsense_Mutation",
                            "Nonstop_Mutation",
                            "In_Frame_Del",
                            "In_Frame_Ins",
                            "Missense_Mutation",
                            "Start_Codon_Del",
                            "Stop_Codon_Del",
                            "Duplication",
                            "Intron_variant")
  #------------------------------------------------------------------------------#
  Onco <-  read.maf(maf = maf, 
                    verbose=T,
                    vc_nonSyn = VariantClassification)
  #------------------------------------------------------------------------------#
  write.mafSummary(maf = Onco, basename = 'Onco')
  #------------------------------------------------------------------------------#
  options_df <- read_excel("../oncoplot_options_v1.xlsx", col_names = TRUE)
  options_list <- setNames(as.list(options_df$Value), options_df$Options)
  get_option <- function(option_name) {
    return(options_list[[option_name]])
  }
  #------------------------------------------------------------------------------#
  Oncoprint_color <- c(
    Frame_Shift_Del = get_option("Frame_Shift_Del"),
    Frame_Shift_Ins = get_option("Frame_Shift_Ins"),
    Splice_Site = get_option("Splice_Site"),
    Nonsense_Mutation = get_option("Nonsense_Mutation"),
    In_Frame_Del = get_option("In_Frame_Del"),
    In_Frame_Ins = get_option("In_Frame_Ins"),
    Missense_Mutation = get_option("Missense_Mutation"),
    Multi_Hit = get_option("Multi_Hit"),
    Translation_Start_Site = get_option("Translation_Start_Site"),
    Nonstop_Mutation = get_option("Nonstop_Mutation"),
    Start_Codon_Del = get_option("Start_Codon_Del"),
    Stop_Codon_Del = get_option("Stop_Codon_Del"),
    Duplication = get_option("Duplication"),
    Intron_variant = get_option("Intron_variant"))
  Oncoprint_color <- Oncoprint_color[!is.na(Oncoprint_color)]
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
    top_gene_num <- as.numeric(Option$Top_gene)
    Oncogenes <- Onco@gene.summary$Hugo_Symbol[1:top_gene_num]
  }
  GeneCut <- length(Oncogenes)
  #------------------------------------------------------------------------------#
  Onco.titv = titv(maf = Onco, plot = FALSE, useSyn = TRUE)
  oncotitiv_filename <- paste0(current_date, "_onco.titv.pdf")
  pdf(oncotitiv_filename)
  plotTiTv(res = Onco.titv)
  dev.off()
  #------------------------------------------------------------------------------#
  current_date <- format(Sys.Date(), "%y%m%d")
  oncoplot_filename <- paste0(current_date, "_oncoplot.pdf")
  pdf(oncoplot_filename)
  oncoplot(maf = Onco, 
           genes = top_genes,
           colors = Oncoprint_color,
           leftBarLims = c(0,100),
           leftBarData = Genes,
           SampleNamefontSize = as.numeric(get_option("SampleNamefontSize")), 
           removeNonMutated = as.logical(get_option("removeNonMutated")),
           writeMatrix = as.logical(get_option("writeMatrix")),
           showTumorSampleBarcodes = as.logical(get_option("showTumorSampleBarcodes")),
           legendFontSize = 1,
           annotationFontSize = 1,
           anno_height = 2,
           titleText = plot_title,
           annotationColor = Annotation_color,
           top = GeneCut)
  dev.off()
  #------------------------------------------------------------------------------#
  mafsummary_filename <- paste0(current_date, "_mafsummary.pdf")
  pdf(mafsummary_filename)
  plotmafSummary(maf = Onco,color = Oncoprint_color, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
  dev.off()
  #------------------------------------------------------------------------------#
  maf_lolli <- read.table('EasyOnco.maf', sep='\t', header=T)
  gene <- subset(maf_lolli, maf_lolli$Variant_Classification!="Intron_variant" & maf_lolli$Variant_Classification !="Silent_Mutation")
  Intron_gene <- unique(gene$Hugo_Symbol)
  # idx <- which(maf_lolli$Hugo_Symbol %in% Intron_gene)
  # maf_lolli <- maf_lolli[-c(idx), ]
  # lolli_gene <- unique(maf_lolli$Hugo_Symbol)
  
  Onco_lolli <-  read.maf(maf = maf_lolli, 
                          verbose=T,
                          vc_nonSyn = VariantClassification)
  
  #------------------------------------------------------------------------------#
  lollipopplot_filename <- paste0(current_date, "_lollipop.pdf")
  pdf(lollipopplot_filename)
  for (Gene in Intron_gene) {
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
  vafplot_filename <- paste0(current_date, "_vafplot.pdf")
  pdf(vafplot_filename)
  plotVaf(maf = Onco, vafCol = 'i_TumorVAF_WU', width = 10, height = 10)
  dev.off()
  #------------------------------------------------------------------------------#
  drugint_filename <- paste0(current_date, "_druginteraction.pdf")
  pdf(drugint_filename)
  drugInteractions(maf = Onco, fontSize = 0.75)
  dev.off()
  #------------------------------------------------------------------------------# 
}
