library(RESET)
library(AUCell)
# library(dplyr)
# library(reticulate)

set.seed(42)
args <- commandArgs(trailingOnly = TRUE)

# gct_path="/home/people/chihs/Gene_Deconvolution/GSEApy/tests/extdata/Leukemia_hgu95av2.trim.txt"
# gmt_path="/home/people/chihs/Gene_Deconvolution/GSEApy/tests/extdata/enrichr.KEGG_2016.gmt"
gene_mat_path <- args[1]
gmt_path <- args[2]
tool_name <- args[3]
out_dir <- args[4]
base_name <- basename(gene_mat_path)


read_gep_rds <- function(rds_path){

  data <- readRDS(rds_path)

  # Check whether the indices of the matrix are unique
  # If it is not, sum the values for each gene (row-wise aggregation)
  if (length(rownames(data)) != length(unique(rownames(data)))){
    data_sum <- data %>%
    group_by(Gene) %>%
    summarise(across(where(is.numeric), \(x) sum(x, na.rm <- TRUE)))
  }
  
  return(data)
}

read_gep_csv <- function(path){

  data <- read.csv(path, row.names=1)

  # Check whether the indices of the matrix are unique
  # If it is not, sum the values for each gene (row-wise aggregation)
  if (length(rownames(data)) != length(unique(rownames(data)))){
    data_sum <- data %>%
    group_by(Gene) %>%
    summarise(across(where(is.numeric), \(x) sum(x, na.rm <- TRUE)))
  }
  
  return(data)
}


read_gct <- function(gct_path){

  # Step 1: Read the data
  data <- read.table(gct_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # Step 2: Remove the second column ("NAME")
  data <- data[, -2]

  # Step 3: Sum the values for each gene (row-wise aggregation)
  data_sum <- data %>%
    group_by(Gene) %>%
    summarise(across(where(is.numeric), \(x) sum(x, na.rm <- TRUE)))

  # Step 4: Convert the result to a data frame (to allow row names)
  data_sum <- as.data.frame(data_sum)

  # Step 5: Set "Gene" as the row names of the dataframe
  rownames(data_sum) <- data_sum$Gene
  data_sum <- data_sum[, -1]  # Remove the Gene column, as it's now the row name
  data_sum <- as.matrix(data_sum)

  return(data_sum)
}


run_gseapy <- function(gct_path, gmt_path){

  #Run in command mode
  venv_dir <- file.path(getwd(),"python_venv")
  print(venv_dir)
  cmd <- sprintf("'%s'/bin/gseapy ssgsea -d '%s' -g '%s' -o gseapy-ssgsea -p 4 --no-plot", venv_dir, gct_path, gmt_path)
  system(cmd)

  
}


run_reset <- function(X, gene_sets, out_file) {

  row_names <- rownames(X)
  gene_sets <- lapply(gene_sets, function(genes) {intersect(row_names, genes)})

  res <- reset(t(X), var.sets = gene_sets,  k = 2, random.threshold = 10)
  write.csv(t(res$S), file = out_file, row.names = TRUE)
}


run_ssgsea <- function(X, gene_sets, alpha = 0.25, scale = T, norm = F, 
                        single = T, out_file = "") {

  row_names <- rownames(X)
  num_genes <- nrow(X)
  gene_sets <- lapply(gene_sets, function(genes) {which(row_names %in% genes)})
  # Ranks for genes
  R <- matrixStats::colRanks(X, preserveShape = T, ties.method = "average")

  # Calculate enrichment score (es) for each sample (column)
  es <- apply(R, 2, function(R_col) {
    gene_ranks <- order(R_col, decreasing = TRUE)

    # Calc es for each gene set
    es_sample <- sapply(gene_sets, function(gene_set_idx) {
      # pos: match (within the gene set)
      # neg: non-match (outside the gene set)
      indicator_pos <- gene_ranks %in% gene_set_idx
      indicator_neg = !indicator_pos

      rank_alpha  <- (R_col[gene_ranks] * indicator_pos) ^ alpha

      step_cdf_pos <- cumsum(rank_alpha)    / sum(rank_alpha)
      step_cdf_neg <- cumsum(indicator_neg) / sum(indicator_neg)

      step_cdf_diff <- step_cdf_pos - step_cdf_neg

      # Normalize by gene number
      if (scale) step_cdf_diff <- step_cdf_diff / num_genes

      # Use ssGSEA or not
      if (single) {
        sum(step_cdf_diff)
      } else {
        step_cdf_diff[which.max(abs(step_cdf_diff))]
      }
    })
    #unlist(es_sample)
  })

  if (length(gene_sets) == 1) es <- matrix(es, nrow = 1)

  # Normalize by absolute diff between max and min
  if (norm) es <- es / diff(range(es))

  # Prepare output
  rownames(es) <- names(gene_sets)
  colnames(es) <- colnames(X)

  write.csv(t(es), file = out_file)
  # saveRDS(es, file = out_file)
}

run_aucell <- function(X, gene_sets, out_file = ""){
  sample_AUC <- AUCell_run(X, geneSets) # Add BPPARAM=BiocParallel::MulticoreParam(5) to run in parallel (number of cores=5)

}


#gene_mat <- read_gct(gct_path)
gene_mat <- read_gep_csv(gene_mat_path)
print("Loaded gene matrix")
gene_sets <- fgsea::gmtPathways(gmt_path)
print("Loaded gmt_path")

out_file <- sprintf("%s/%s_%s.csv", out_dir, base_name, tool_name)
print(out_file)
res <- switch(tool_name,
        "ssgsea" = run_ssgsea(gene_mat, gene_sets, scale = F, norm = T, 
                              out_file = out_file),
        "reset" = run_reset(gene_mat, gene_sets, out_file),
        "aucell" = run_aucell(gene_mat, gene_sets, out_file))