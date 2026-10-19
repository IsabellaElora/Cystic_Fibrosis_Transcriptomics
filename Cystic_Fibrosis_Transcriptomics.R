# 1. Import the gene-expression count table
counts <- read.delim(
  file.choose(),
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)
# 2. Check the dimensions
dim(counts)
# 3. Look at the first six rows
head(counts)

# 4. Check the sample names
colnames(counts)

# 5. Check whether the values are numeric
str(counts)

# 6. Check for missing values
sum(is.na(counts))

# 7. Check whether gene IDs are duplicated
sum(duplicated(rownames(counts)))

# 8. Check whether all counts are non-negative integers
count_matrix <- as.matrix(counts)
is_numeric <- is.numeric(count_matrix)

if (is_numeric) {
  print(all(count_matrix >= 0))
  print(all(count_matrix == round(count_matrix)))
  
  # 9. Calculate total counts for each sample
  library_sizes <- colSums(counts)
  # Display the totals
  print(library_sizes)
  # 10. Create a data frame for plotting
  library_data <- data.frame(
    Sample = names(library_sizes),
    Total_Counts = as.numeric(library_sizes)
  )
  # 11. Plot total counts per sample
  barplot(
    library_data$Total_Counts,
    names.arg = library_data$Sample,
    las = 2,
    cex.names = 0.7,
    main = "Total Sequencing Counts per Sample",
    ylab = "Total counts",
    xlab = "Sample"
  ) 
  # 12. Create sample metadata
  sample_names <- colnames(counts)
  
  metadata <- data.frame(
    Sample = sample_names,
    Status = ifelse(grepl("_CF_", sample_names), "CF", "NCF"),
    Donor = sub(
      "^[0-9]+_(CF|NCF)_([0-9]+)_(UNC|SC)$",
      "\\2",
      sample_names
    ),
    Medium = sub(
      "^[0-9]+_(CF|NCF)_[0-9]+_(UNC|SC)$",
      "\\2",
      sample_names
    )
  )
  
  # Make the sample information easier to read
  metadata$Status <- factor(metadata$Status, levels = c("NCF", "CF"))
  metadata$Medium <- factor(metadata$Medium, levels = c("UNC", "SC"))
  
  # Display the metadata
  print(metadata)

  # 13. Select the UNC samples for the primary analysis
  unc_samples <- metadata$Sample[metadata$Medium == "UNC"]
  
  unc_counts <- counts[, unc_samples, drop = FALSE]
  unc_metadata <- metadata[metadata$Medium == "UNC", ]
  # Check the selected samples
  print(colnames(unc_counts))
  print(unc_metadata)
  
  # Count the samples in each disease group
  print(table(unc_metadata$Status))
  
  # 15. Check that counts and metadata are aligned
  stopifnot(identical(
    colnames(unc_counts),
    unc_metadata$Sample
  ))
  
  # 16. Check the dimensions of the selected data
  dim(unc_counts)
  # 17. Check the number of samples in each group
  table(unc_metadata$Status)
  
  # 18. Confirm DESeq2 is available
  library(DESeq2)
  # 19. Ensure the sample metadata matches the count columns
  stopifnot(identical(
    colnames(unc_counts),
    unc_metadata$Sample
  ))
  # 20. Make sure the counts are whole numbers
  stopifnot(all(unc_counts >= 0))
  stopifnot(all(unc_counts == round(unc_counts)))  
  # 21. Set non-CF as the reference group
  unc_metadata$Status <- relevel(
    droplevels(unc_metadata$Status),
    ref = "NCF"
  )
  
  # 22. Create the DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(unc_counts)),
    colData = as.data.frame(unc_metadata),
    design = ~ Status
  )
  
  # 23. Report the dataset dimensions
  dim(dds)  
  
  # 24. Count genes before filtering
  n_before <- nrow(dds)
  # 25. Keep genes with at least 10 counts in at least 3 samples
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep, ]
  
  # 26. Count genes after filtering
  n_after <- nrow(dds)
  
  cat("Genes before filtering:", n_before, "\n")
  cat("Genes after filtering:", n_after, "\n")
  cat("Genes removed:", n_before - n_after, "\n")  
  # 27. Estimate sequencing-depth size factors
  dds <- estimateSizeFactors(dds)
  
  # Inspect the size factors
  sizeFactors(dds)
  
  # 28. Transform counts for visualisation
  vsd <- vst(dds, blind = TRUE)
  
  # 29. Create a principal component analysis (PCA) plot
  plotPCA(vsd, intgroup = c("Status")) 
  # 30. Extract PCA coordinates and sample labels
  pca_data <- plotPCA(
    vsd,
    intgroup = "Status",
    returnData = TRUE
  )
  
  pca_data$Sample <- rownames(pca_data)
  
  # Display the sample coordinates
  print(pca_data[order(pca_data$PC1, decreasing = TRUE), ])
  
  # 31. Compare total counts for the UNC samples
  print(colSums(unc_counts))
  # 32. Inspect donor 2 in both culture media
  donor2_samples <- metadata$Sample[metadata$Donor == "2" &
                                      metadata$Status == "NCF"]
  
  print(donor2_samples)
  
  # Compare total counts for these two samples
  print(colSums(counts[, donor2_samples, drop = FALSE]))
  
  # Compare their PCA coordinates in the current UNC analysis
  print(pca_data[pca_data$Sample == "9_NCF_2_UNC", ])  
  # 33. Prepare the full dataset for visualisation
  dds_all <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(counts)),
    colData = as.data.frame(metadata),
    design = ~ Status
  )
  
  # Filter low-count genes for visualisation
  keep_all <- rowSums(counts(dds_all) >= 10) >= 3
  dds_all <- dds_all[keep_all, ]
  
  # Normalise and transform
  dds_all <- estimateSizeFactors(dds_all)
  vsd_all <- vst(dds_all, blind = TRUE)
  
  # 34. PCA coloured by disease status and shaped by medium
  pca_all <- plotPCA(
    vsd_all,
    intgroup = c("Status", "Medium"),
    returnData = TRUE
  )
  
  percentVar <- round(100 * attr(pca_all, "percentVar"))
  
  library(ggplot2)
  
  ggplot(pca_all, aes(
    x = PC1, y = PC2,
    colour = Status, shape = Medium,
    label = name
  )) +
    geom_point(size = 4) +
    geom_text(vjust = -0.8, size = 3, show.legend = FALSE) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme_minimal() +
    labs(
      title = "PCA of CF and non-CF samples",
      subtitle = "Samples shown by disease status and culture medium"
    )  
  # 35. Run the differential-expression analysis
  dds <- DESeq(dds)
  
  # 36. Compare CF against non-CF
  res <- results(
    dds,
    contrast = c("Status", "CF", "NCF"),
    alpha = 0.05
  )
  
  # 37. Summarise the results
  summary(res)
  
  # 38. Order genes by adjusted p-value
  res_ordered <- res[order(res$padj), ]
  
  # Display the 10 most significant results
  head(as.data.frame(res_ordered), 10)
  
  # 39. Count genes passing the adjusted p-value threshold
  sum(res$padj < 0.05, na.rm = TRUE)
  
  # Convert results into a standard data frame
  res_df <- as.data.frame(res_ordered)
  
  # Add gene IDs as an explicit column
  res_df$GeneID <- rownames(res_df)
  
  # Put GeneID first
  res_df <- res_df[, c("GeneID", setdiff(names(res_df), "GeneID"))]
  
  # Save the complete results table
  write.csv(
    res_df,
    "DESeq2_UNC_CF_vs_NCF_results.csv",
    row.names = FALSE
  )
  
  # Save only statistically significant genes
  significant_genes <- subset(
    res_df,
    !is.na(padj) & padj < 0.05
  )
  
  write.csv(
    significant_genes,
    "DESeq2_UNC_CF_vs_NCF_significant_genes.csv",
    row.names = FALSE
  )
  
  # Check how many significant genes were saved
  nrow(significant_genes)

  library(ggplot2)
  
  # Prepare plotting data
  volcano_data <- as.data.frame(res)
  volcano_data$GeneID <- rownames(volcano_data)
  
  # Mark genes by adjusted p-value and fold-change direction
  volcano_data$Category <- "Not significant"
  
  volcano_data$Category[
    !is.na(volcano_data$padj) &
      volcano_data$padj < 0.05 &
      volcano_data$log2FoldChange > 0
  ] <- "Higher in CF"
  
  volcano_data$Category[
    !is.na(volcano_data$padj) &
      volcano_data$padj < 0.05 &
      volcano_data$log2FoldChange < 0
  ] <- "Lower in CF"
  
  # Calculate -log10 adjusted p-value
  volcano_data$negLog10Padj <- -log10(
    pmax(volcano_data$padj, .Machine$double.xmin)
  )
  
  # Create the plot
  volcano_plot <- ggplot(
    volcano_data,
    aes(
      x = log2FoldChange,
      y = negLog10Padj,
      colour = Category
    )
  ) +
    geom_point(alpha = 0.65, size = 1.2, na.rm = TRUE) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_minimal() +
    labs(
      title = "Differential expression in CF versus non-CF",
      subtitle = "UNC culture medium; adjusted p-value < 0.05 marks significance",
      x = "Log2 fold change (CF / non-CF)",
      y = "-Log10 adjusted p-value",
      colour = "Expression category"
    )
  
  print(volcano_plot)
  
  # Save the figure
  ggsave(
    "Volcano_UNC_CF_vs_NCF.png",
    plot = volcano_plot,
    width = 9,
    height = 7,
    dpi = 300
  )  
  png(
    "MA_plot_UNC_CF_vs_NCF.png",
    width = 2400,
    height = 1800,
    res = 300
  )
  
  plotMA(
    res,
    main = "MA plot: CF versus non-CF (UNC medium)",
    ylim = c(-5, 5),
    alpha = 0.05
  )
  
  dev.off()  

  ggsave(
    "PCA_UNC_CF_vs_NCF.png",
    plot = plotPCA(vsd, intgroup = "Status"),
    width = 8,
    height = 6,
    dpi = 300
  )  
  
  
  # 44. Create a dataset excluding the unusual sample for sensitivity analysis
  # This is a diagnostic comparison only, not a decision to discard the sample
  
  sample_to_check <- "9_NCF_2_UNC"
  
  counts_sensitivity <- unc_counts[
    ,
    colnames(unc_counts) != sample_to_check,
    drop = FALSE
  ]
  
  metadata_sensitivity <- droplevels(
    unc_metadata[
      unc_metadata$Sample != sample_to_check,
      ,
      drop = FALSE
    ]
  )
  
  # Check which samples remain
  colnames(counts_sensitivity)
  metadata_sensitivity[, c("Sample", "Status")]
  
  # 45. Rebuild the DESeq2 dataset
  dds_sensitivity <- DESeqDataSetFromMatrix(
    countData = round(as.matrix(counts_sensitivity)),
    colData = as.data.frame(metadata_sensitivity),
    design = ~ Status
  )
  
  # 46. Apply the same low-count filter
  keep_sensitivity <- rowSums(
    counts(dds_sensitivity) >= 10
  ) >= 3
  
  dds_sensitivity <- dds_sensitivity[keep_sensitivity, ]
  
  # 47. Run differential expression
  dds_sensitivity <- DESeq(dds_sensitivity)
  
  res_sensitivity <- results(
    dds_sensitivity,
    contrast = c("Status", "CF", "NCF"),
    alpha = 0.05
  )
  
  # 48. Count significant genes in this sensitivity analysis
  sum(res_sensitivity$padj < 0.05, na.rm = TRUE)
  
  # Identify significant genes in each analysis
  sig_primary <- rownames(res)[
    !is.na(res$padj) & res$padj < 0.05
  ]
  
  sig_sensitivity <- rownames(res_sensitivity)[
    !is.na(res_sensitivity$padj) &
      res_sensitivity$padj < 0.05
  ]
  
  # Compare overlap
  length(sig_primary)
  length(sig_sensitivity)
  length(intersect(sig_primary, sig_sensitivity))
  
  # Jaccard similarity: overlap relative to all genes in either set
  length(intersect(sig_primary, sig_sensitivity)) /
    length(union(sig_primary, sig_sensitivity))
  
  # Compare fold changes for genes shared between analyses
  common_genes <- intersect(
    rownames(res),
    rownames(res_sensitivity)
  )
  
  comparison <- data.frame(
    GeneID = common_genes,
    log2FC_primary = res[common_genes, "log2FoldChange"],
    log2FC_sensitivity =
      res_sensitivity[common_genes, "log2FoldChange"],
    padj_primary = res[common_genes, "padj"],
    padj_sensitivity = res_sensitivity[common_genes, "padj"]
  )
  
  # Check whether fold-change directions agree
  comparison$Direction_agrees <- sign(
    comparison$log2FC_primary
  ) == sign(comparison$log2FC_sensitivity)
  
  mean(comparison$Direction_agrees, na.rm = TRUE)
  
  # Save the comparison
  write.csv(
    comparison,
    "Sensitivity_comparison_UNC.csv",
    row.names = FALSE
  )
  
  # 50. Keep genes with fold-change estimates in both analyses
  comparison_valid <- subset(
    comparison,
    is.finite(log2FC_primary) &
      is.finite(log2FC_sensitivity)
  )
  
  # Correlation between the two sets of fold-change estimates
  cor(
    comparison_valid$log2FC_primary,
    comparison_valid$log2FC_sensitivity,
    method = "spearman"
  )
  
  # Summarise the differences in fold change
  comparison_valid$FC_difference <-
    comparison_valid$log2FC_sensitivity -
    comparison_valid$log2FC_primary
  
  summary(comparison_valid$FC_difference)
  
  # Plot the estimates against one another
  library(ggplot2)
  
  fc_plot <- ggplot(
    comparison_valid,
    aes(
      x = log2FC_primary,
      y = log2FC_sensitivity
    )
  ) +
    geom_point(alpha = 0.35, size = 1) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed"
    ) +
    theme_minimal() +
    labs(
      title = "Fold-change comparison between analyses",
      subtitle = "UNC medium: all samples versus excluding one non-CF sample",
      x = "Log2 fold change: primary analysis",
      y = "Log2 fold change: sensitivity analysis"
    )
  
  print(fc_plot)
  
  ggsave(
    "Sensitivity_fold_change_comparison.png",
    plot = fc_plot,
    width = 8,
    height = 6,
    dpi = 300
  )  
  
  # Search the count matrix for the CFTR gene ID
  grep(
    "ENSG00000001626",
    rownames(counts),
    value = TRUE
  )
  
  # Check whether this gene is present in the filtered analysis
  "ENSG00000001626" %in% rownames(res)
  
  # If present, display its differential-expression result
  if ("ENSG00000001626" %in% rownames(res)) {
    print(as.data.frame(res["ENSG00000001626", ]))
  }
  # 52. Download ClinVar records for CFTR
  
  if (!requireNamespace("rentrez", quietly = TRUE)) {
    install.packages("rentrez")
  }
  
  library(rentrez)
  
  # Search ClinVar for records associated with CFTR
  clinvar_search <- entrez_search(
    db = "clinvar",
    term = '"CFTR"[Gene]',
    retmax = 10000
  )
  
  length(clinvar_search$ids)
  clinvar_xml <- entrez_fetch(
    db = "clinvar",
    id = clinvar_search$ids,
    rettype = "xml",
    retmode = "text"
  )
  
  writeLines(clinvar_xml, "ClinVar_CFTR_raw.xml")
  
  file.info("ClinVar_CFTR_raw.xml")$size  
  
  # Download ClinVar records in smaller batches
  
  library(rentrez)
  
  all_ids <- clinvar_search$ids
  batch_size <- 100
  
  # Split the record IDs into batches
  id_batches <- split(
    all_ids,
    ceiling(seq_along(all_ids) / batch_size)
  )
  
  length(id_batches)
  
  # Fetch each batch
  clinvar_batches <- vector("list", length(id_batches))
  
  for (i in seq_along(id_batches)) {
    message(
      "Downloading batch ", i, " of ", length(id_batches)
    )
    
    clinvar_batches[[i]] <- entrez_fetch(
      db = "clinvar",
      id = id_batches[[i]],
      rettype = "xml",
      retmode = "text"
    )
    
    # Brief pause between requests
    Sys.sleep(0.35)
  }
  
  # Combine and save the downloaded records
  clinvar_xml <- paste(
    clinvar_batches,
    collapse = "\n"
  )
  
  writeLines(clinvar_xml, "ClinVar_CFTR_raw.xml")
  
  # Confirm the download
  file.info("ClinVar_CFTR_raw.xml")$size
  
  # Download ClinVar records in batches of 100
  clinvar_batches <- vector("list", length(id_batches))
  
  for (i in seq_along(id_batches)) {
    message(
      "Downloading batch ", i, " of ", length(id_batches)
    )
    
    clinvar_batches[[i]] <- entrez_fetch(
      db = "clinvar",
      id = id_batches[[i]],
      rettype = "xml",
      retmode = "text"
    )
    
    Sys.sleep(0.35)
  }  # This closes the loop
  
  # Combine the downloaded responses AFTER the loop
  clinvar_xml <- paste(
    clinvar_batches,
    collapse = "\n"
  )
  
  # Save the raw responses
  writeLines(clinvar_xml, "ClinVar_CFTR_raw.xml")
  
  # Confirm the file was saved
  file.exists("ClinVar_CFTR_raw.xml")
  file.info("ClinVar_CFTR_raw.xml")$size
  
  # Download ClinVar records in batches of 100
  clinvar_batches <- vector("list", length(id_batches))
  
  for (i in seq_along(id_batches)) {
    message("Downloading batch ", i, " of ", length(id_batches))
    
    clinvar_batches[[i]] <- entrez_fetch(
      db = "clinvar",
      id = id_batches[[i]],
      rettype = "xml",
      retmode = "text"
    )
    
    Sys.sleep(0.35)
  }
  
  # The loop has now finished
  clinvar_xml <- paste(clinvar_batches, collapse = "\n")
  
  writeLines(clinvar_xml, "ClinVar_CFTR_raw.xml")
  
  file.exists("ClinVar_CFTR_raw.xml")
  file.info("ClinVar_CFTR_raw.xml")$size
  
  # Install the XML parser if needed
  if (!requireNamespace("xml2", quietly = TRUE)) {
    install.packages("xml2")
  }
  
  library(xml2)
  
  # Check that all 66 batches were downloaded
  length(clinvar_batches)
  sum(lengths(clinvar_batches))
  
  # Inspect the beginning of the first response
  cat(substr(clinvar_batches[[1]], 1, 500)) 
  
  
  # Parse the first batch as XML
  first_batch_xml <- read_xml(clinvar_batches[[1]])
  
  # Inspect the top-level structure
  xml_name(xml_root(first_batch_xml))
  
  # List the first few child elements
  xml_name(xml_children(xml_root(first_batch_xml)))[1:10]  
  
  # Retrieve a small test batch in VCF format
  test_ids <- head(all_ids, 5)
  
  test_vcf <- entrez_fetch(
    db = "clinvar",
    id = test_ids,
    rettype = "vcf",
    retmode = "text"
  )
  
  # Inspect the returned content
  cat(substr(test_vcf, 1, 1500))
  
  # Retrieve one ClinVar record in structured XML
  test_id <- all_ids[1]
  
  test_url <- paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
    "?db=clinvar&id=", test_id,
    "&rettype=clinvarset&retmode=xml"
  )
  
  response <- readLines(test_url, warn = FALSE)
  
  # Inspect the returned content
  cat(paste(head(response, 30), collapse = "\n"))
  
  # Retrieve structured summaries for five ClinVar IDs
  test_summary <- entrez_summary(
    db = "clinvar",
    id = head(all_ids, 5),
    retmode = "json"
  )
  
  # Inspect the returned structure
  str(test_summary, max.level = 2)
  
  # Display the first record
  test_summary[[1]]
  
  
  # 60. Retrieve structured ClinVar summaries in batches
  
  summary_batches <- vector("list", length(id_batches))
  
  for (i in seq_along(id_batches)) {
    message(
      "Retrieving summary batch ", i, " of ",
      length(id_batches)
    )
    
    summary_batches[[i]] <- entrez_summary(
      db = "clinvar",
      id = id_batches[[i]],
      retmode = "json"
    )
    
    Sys.sleep(0.35)
  }
  
  # Combine summaries into one list
  all_summaries <- unlist(
    summary_batches,
    recursive = FALSE
  )
  
  length(all_summaries)
  
  # Extract the fields available in each ClinVar summary
  extract_clinvar_record <- function(x) {
    
    # Extract germline clinical significance where available
    clinical_sig <- NA_character_
    
    if (!is.null(x$germline_classification$description)) {
      clinical_sig <- x$germline_classification$description
    }
    
    data.frame(
      UID = if (!is.null(x$uid)) x$uid else NA_character_,
      Accession = if (!is.null(x$accession)) x$accession else NA_character_,
      Title = if (!is.null(x$title)) x$title else NA_character_,
      Variant_type = if (!is.null(x$obj_type)) x$obj_type else NA_character_,
      Protein_change = if (!is.null(x$protein_change)) x$protein_change else NA_character_,
      Molecular_consequence = if (!is.null(x$molecular_consequence_list)) x$molecular_consequence_list else NA_character_,
      Clinical_significance = clinical_sig,
      stringsAsFactors = FALSE
    )
  }
  
  clinvar_table <- do.call(
    rbind,
    lapply(all_summaries, extract_clinvar_record)
  )
  
  # Remove duplicate records by accession
  clinvar_table <- clinvar_table[
    !duplicated(clinvar_table$Accession),
  ]
  
  # Save the extracted table
  write.csv(
    clinvar_table,
    "ClinVar_CFTR_variant_table.csv",
    row.names = FALSE
  )
  
  # Inspect the result
  dim(clinvar_table)
  head(clinvar_table, 10)
  
  # Check whether the table was created
  exists("clinvar_table")
  
  # Check its dimensions, if it exists
  if (exists("clinvar_table")) {
    print(dim(clinvar_table))
    print(head(clinvar_table, 10))
  }
  
  # Check that the downloaded summaries exist
  exists("all_summaries")
  length(all_summaries)
  
  # Inspect the structure of the first record
  str(all_summaries[[1]], max.level = 2)
  
  length(all_summaries)
2+2
exists("all_summaries")
length(all_summaries)

clinvar_table <- do.call(rbind, lapply(all_summaries, function(x) {
  data.frame(
    UID = as.character(x$uid %||% NA),
    Accession = as.character(x$accession %||% NA),
    Title = as.character(x$title %||% NA),
    Variant_type = as.character(x$obj_type %||% NA),
    Protein_change = as.character(x$protein_change %||% NA),
    Molecular_consequence = paste(x$molecular_consequence_list, collapse = "; "),
    Clinical_significance = as.character(x$germline_classification$description %||% NA),
    stringsAsFactors = FALSE
  )
}))

get_field <- function(x, field) {
  value <- x[[field]]
  if (is.null(value) || length(value) == 0) NA_character_
  else paste(as.character(value), collapse = "; ")
}

clinvar_table <- do.call(rbind, lapply(all_summaries, function(x) {
  data.frame(
    UID = get_field(x, "uid"),
    Accession = get_field(x, "accession"),
    Title = get_field(x, "title"),
    Variant_type = get_field(x, "obj_type"),
    Protein_change = get_field(x, "protein_change"),
    Molecular_consequence = get_field(x, "molecular_consequence_list"),
    Clinical_significance = get_field(x$germline_classification, "description"),
    stringsAsFactors = FALSE
  )
}))

clinvar_table <- clinvar_table[!duplicated(clinvar_table$Accession), ]
write.csv(clinvar_table, "ClinVar_CFTR_variant_table.csv", row.names = FALSE)

print(dim(clinvar_table))
View(clinvar_table)

# Inspect the first 10 variant titles
head(clinvar_table$Title, 10)

# Count the different record types
sort(table(clinvar_table$Variant_type), decreasing = TRUE)

# Search for the CFTR Phe508del variant
clinvar_table[
  grepl("Phe508del|Phe508", clinvar_table$Title, ignore.case = TRUE),
  c("Accession", "Title", "Protein_change", "Clinical_significance")
]

# Keep records whose titles explicitly identify CFTR
cftr_variants <- clinvar_table[
  grepl("\\(CFTR\\):", clinvar_table$Title, ignore.case = TRUE),
]

# Check the number of records retained
dim(cftr_variants)

# Check for missing titles and accessions
sum(is.na(cftr_variants$Title) | cftr_variants$Title == "")
sum(is.na(cftr_variants$Accession) | cftr_variants$Accession == "")

# Save the focused dataset
write.csv(
  cftr_variants,
  "ClinVar_CFTR_focused_variants.csv",
  row.names = FALSE
)

# Inspect the pathogenic records
cftr_variants[
  grepl("pathogenic", cftr_variants$Clinical_significance, ignore.case = TRUE),
  c("Accession", "Title", "Protein_change", "Clinical_significance")
] |> head(15)

# Count records by clinical significance
sort(
  table(cftr_variants$Clinical_significance, useNA = "ifany"),
  decreasing = TRUE
)

# Count variants with a protein change recorded
sum(!is.na(cftr_variants$Protein_change) &
      cftr_variants$Protein_change != "")

# Inspect all records matching the canonical F508del notation
cftr_variants[
  grepl("c\\.1522_1524del|p\\.Phe508del|F508del",
        cftr_variants$Title, ignore.case = TRUE) |
    grepl("^F508del$", cftr_variants$Protein_change, ignore.case = TRUE),
  c("Accession", "Title", "Protein_change", "Clinical_significance")
]

# Select one representative record for each variant category
examples <- rbind(
  cftr_variants[grepl("p\\.Phe508del", cftr_variants$Title), ][1, ],
  cftr_variants[grepl("frameshift", cftr_variants$Title,
                      ignore.case = TRUE) &
                  grepl("Pathogenic|Likely pathogenic",
                        cftr_variants$Clinical_significance,
                        ignore.case = TRUE), ][1, ],
  cftr_variants[grepl("Ter|\\*", cftr_variants$Title) &
                  grepl("Pathogenic|Likely pathogenic",
                        cftr_variants$Clinical_significance,
                        ignore.case = TRUE), ][1, ],
  cftr_variants[grepl("splice", cftr_variants$Molecular_consequence,
                      ignore.case = TRUE) &
                  grepl("Pathogenic|Likely pathogenic",
                        cftr_variants$Clinical_significance,
                        ignore.case = TRUE), ][1, ]
)

examples[, c("Accession", "Title", "Protein_change",
             "Molecular_consequence", "Clinical_significance")]

write.csv(examples, "CFTR_variant_examples.csv", row.names = FALSE)

# Find pathogenic or likely pathogenic frameshift variants
frameshift_candidates <- cftr_variants[
  grepl("fs", cftr_variants$Protein_change, ignore.case = TRUE) &
    grepl("Pathogenic|Likely pathogenic",
          cftr_variants$Clinical_significance,
          ignore.case = TRUE),
]

# Inspect the first five candidates
frameshift_candidates[
  1:min(5, nrow(frameshift_candidates)),
  c("Accession", "Title", "Protein_change",
    "Molecular_consequence", "Clinical_significance")
]

# Fill the missing frameshift row with p.Arg75fs
examples[2, ] <- frameshift_candidates[1, ]

# Display the four selected variants
examples[, c(
  "Accession",
  "Title",
  "Protein_change",
  "Molecular_consequence",
  "Clinical_significance"
)]

# Save the final selection
write.csv(
  examples,
  "CFTR_variant_examples.csv",
  row.names = FALSE
)

# Reset row numbers
rownames(examples) <- NULL

# Save the cleaned table
write.csv(
  examples,
  "CFTR_variant_examples.csv",
  row.names = FALSE
)

# Confirm the final table
View(examples)
variant_map <- data.frame(
  Variant = c("F508del", "R75fs", "K1292Ter", "c.2989-2A>C"),
  Protein_position = c(508, 75, 1292, NA),
  Protein_region = c(
    "NBD1",
    "N-terminal region",
    "NBD2",
    "RNA splicing; protein position not assigned"
  ),
  Variant_type = c(
    "In-frame deletion",
    "Frameshift",
    "Premature stop",
    "Splice-acceptor variant"
  )
)

print(variant_map)

write.csv(
  variant_map,
  "CFTR_variant_protein_map.csv",
  row.names = FALSE
)
files <- c(
  "ClinVar_CFTR_focused_variants.csv",
  "CFTR_variant_protein_map.csv",
  "DESeq2_UNC_CF_vs_NCF_results.csv"
)

file.exists(files)

variants <- read.csv("CFTR_variant_examples.csv")
protein_map <- read.csv("CFTR_variant_protein_map.csv")

names(variants)
names(protein_map)


combined <- combined[, c(
  "Variant",
  "Clinical_significance",
  "Molecular_consequence",
  "Protein_position",
  "Protein_region"
)]

print(combined)

write.csv(
  combined,
  "CFTR_combined_variant_analysis.csv",
  row.names = FALSE
)

variants$Variant <- c(
  "F508del",
  "R75fs",
  "K1292Ter",
  "c.2989-2A>C"
)

combined <- merge(
  variants,
  protein_map,
  by = "Variant"
)

print(combined)
final_variants <- combined[, c(
  "Variant",
  "Clinical_significance",
  "Molecular_consequence",
  "Protein_position",
  "Protein_region"
)]

write.csv(
  final_variants,
  "CFTR_combined_variant_analysis.csv",
  row.names = FALSE
)

View(final_variants)

library(ggplot2)

plot_variants <- subset(final_variants, !is.na(Protein_position))

ggplot(plot_variants, aes(x = Protein_position, y = 1)) +
  geom_segment(aes(x = 1, xend = 1480, y = 1, yend = 1),
               colour = "grey60", linewidth = 3) +
  geom_point(colour = "red", size = 5) +
  geom_text(aes(label = Variant), vjust = -1.2, size = 4) +
  scale_x_continuous(limits = c(0, 1500), breaks = seq(0, 1500, 250)) +
  labs(
    title = "Selected CFTR variants mapped to protein sequence",
    x = "Amino acid position",
    y = NULL
  ) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

domains <- data.frame(
  Domain = c("MSD1", "NBD1", "R region", "MSD2", "NBD2"),
  Start = c(81, 423, 654, 859, 1210),
  End = c(365, 646, 831, 1155, 1443)
)

print(domains)

cftr_domain_plot <- ggplot() +
  annotate("segment", x = 1, xend = 1480,
           y = 1, yend = 1,
           colour = "grey75", linewidth = 5) +
  geom_segment(
    data = domains,
    aes(x = Start, xend = End, y = 1, yend = 1,
        colour = Domain),
    linewidth = 9
  ) +
  geom_point(
    data = plot_variants,
    aes(x = Protein_position, y = 1),
    colour = "black", size = 4
  ) +
  geom_text(
    data = plot_variants,
    aes(x = Protein_position, y = 1.22, label = Variant),
    size = 4
  ) +
  scale_x_continuous(limits = c(0, 1500),
                     breaks = seq(0, 1500, 250)) +
  coord_cartesian(ylim = c(0.7, 1.5)) +
  labs(
    title = "CFTR variants and functional protein domains",
    x = "Amino acid position",
    y = NULL,
    colour = "Protein domain"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

print(cftr_domain_plot)

ggsave(
  "CFTR_variant_domain_map.png",
  plot = cftr_domain_plot,
  width = 12,
  height = 5,
  dpi = 300
)

file.exists("CFTR_variant_domain_map.png")

expression_results <- read.csv(
  "DESeq2_UNC_CF_vs_NCF_results.csv",
  row.names = 1
)

expression_results["ENSG00000001626", ]

cftr_expression <- expression_results[
  "ENSG00000001626", ,
  drop = FALSE
]

write.csv(
  cftr_expression,
  "CFTR_expression_summary.csv"
)

file.exists("CFTR_expression_summary.csv")

significant_genes <- subset(
  expression_results,
  !is.na(padj) & padj < 0.05
)

top10_genes <- head(
  significant_genes[order(significant_genes$padj), ],
  10
)

print(top10_genes)

write.csv(
  top10_genes,
  "CF_top10_differentially_expressed_genes.csv"
)


install.packages("BiocManager")

BiocManager::install(
  c("AnnotationDbi", "org.Hs.eg.db"),
  ask = FALSE,
  update = FALSE
)





library(AnnotationDbi)
library(org.Hs.eg.db)

top10_genes$Gene_name <- mapIds(
  org.Hs.eg.db,
  keys = rownames(top10_genes),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

print(top10_genes[, c(
  "Gene_name",
  "log2FoldChange",
  "padj"
)])

library(AnnotationDbi)
library(org.Hs.eg.db)

top10_genes$Gene_name <- mapIds(
  org.Hs.eg.db,
  keys = rownames(top10_genes),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

print(top10_genes[, c("Gene_name", "log2FoldChange", "padj")])

write.csv(
  top10_genes,
  "CF_top10_annotated_genes.csv",
  row.names = TRUE
)

file.exists("CF_top10_annotated_genes.csv")

significant_genes$Gene_name <- mapIds(
  org.Hs.eg.db,
  keys = rownames(significant_genes),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

print(head(significant_genes[, c(
  "Gene_name", "log2FoldChange", "padj"
)]))

genes_of_interest <- c(
  "IL6", "IL8", "CXCL8", "CXCL1", "CXCL2",
  "CXCL5", "CXCL10", "TNF", "IL1B",
  "MUC5AC", "MUC5B", "SCGB3A1",
  "S100A8", "S100A9", "DEFB1",
  "NFKB1", "IL20RB"
)

cf_related_genes <- significant_genes[
  !is.na(significant_genes$Gene_name) &
    significant_genes$Gene_name %in% genes_of_interest,
]

print(cf_related_genes[, c(
  "Gene_name", "log2FoldChange", "padj"
)])

write.csv(
  cf_related_genes,
  "CF_airway_inflammation_candidate_genes.csv",
  row.names = TRUE
)

file.exists("CF_airway_inflammation_candidate_genes.csv")

library(ggplot2)

cf_gene_plot <- ggplot(
  cf_related_genes,
  aes(x = Gene_name, y = log2FoldChange,
      fill = Gene_name)
) +
  geom_col(width = 0.6) +
  geom_text(
    aes(label = round(log2FoldChange, 2)),
    vjust = -0.5,
    size = 5
  ) +
  labs(
    title = "Airway-related gene expression in cystic fibrosis",
    subtitle = "CF versus non-CF | GSE154905, UNC medium",
    x = "Gene",
    y = "Log2 fold change (CF vs non-CF)"
  ) +
  coord_cartesian(ylim = c(0, 2.7)) +
  theme_minimal() +
  theme(legend.position = "none")

print(cf_gene_plot)

ggsave(
  "CF_airway_gene_expression.png",
  plot = cf_gene_plot,
  width = 9,
  height = 5,
  dpi = 300
)

file.exists("CF_airway_gene_expression.png")

project_summary <- data.frame(
  Analysis = c(
    "ClinVar variant records retrieved",
    "CFTR-focused ClinVar records",
    "Selected variants mapped",
    "Significant RNA-seq genes",
    "Airway-related candidate genes",
    "CFTR RNA-seq adjusted p-value"
  ),
  Result = c(
    6527,
    6234,
    4,
    256,
    3,
    0.8619395
  )
)

print(project_summary)

write.csv(
  project_summary,
  "CFTR_project_summary.csv",
  row.names = FALSE
)