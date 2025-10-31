# ── 1. Setup ────────────────────────────────────────────────────────────────────
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("biomaRt")   # one-off install
}
library(biomaRt)
#=============================================================================
# helper (≤500 symbols / call -> avoids BioMart HTTP 413) -----------------
# Function chunk_get() is a  wrapper that breaks a long vector of query 
# values into blocks of ≤ 500 and runs biomaRt::getBM() on each block, 
# then binds the results back together:

chunk_get <- function(v, filter) {
  out <- list()
  for (i in seq(1, length(v), 500)) {                 # step through by 500
    slice <- v[i:min(i+499, length(v))]               # take that slice
    out[[length(out)+1]] <- getBM(                    # run getBM on it
      attributes = attrs,
      filters    = filter,
      values     = slice,
      mart       = ensembl
    )
  }
  do.call(rbind, out)                                 # glue chunks together
}


# ── 1. Load your gene list ──────────────────────────────────────────────────────
gene_df <- read.csv("genes.csv",  stringsAsFactors = FALSE)
#gene_df <- read.table("/Users/.../Desktop/rare_variants/data/bed_file/genes.csv", header=T) # example of full path
genes   <- unique(trimws(gene_df$gene_name))
cat(sprintf("[INFO] Input list has %d unique genes listed\n", length(genes)))  #

# ── 2. Connect to the correct Ensembl dataset ───────────────────────────────────
# Run to see all datasets:
#listDatasets(ensembl)

#    P.s.: choose the nearest mirror
ensembl <- biomaRt::useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror  = "useast"                   #  change to smth else e.g. "asia", or drop this line for default
)

# attributes: p.s. add gene_biotype so we can drop non-coding "annoyances"
attrs <- c("entrezgene_id",            # stable NCBI GeneID
           "external_gene_name",       # Ensembl-preferred symbol
           "external_synonym",
           "gene_biotype",
           "chromosome_name",
           "start_position",
           "end_position")

# ── 3. Retrieve coordinates for given symbols, primary + synonym look-up  ────────────────────────────────────
id_map   <- chunk_get(genes, "external_gene_name")
# we often see missing genes from provided gene list
missing  <- setdiff(genes, id_map$external_gene_name)
if (length(missing)) {
  cat(sprintf("[WARN] %d symbols not found by exact match\n", length(missing)))

  ## ask BioMart to return the synonym it matched
  syn_hits <- chunk_get(missing, "external_synonym")

  ## append to the master table
  id_map   <- rbind(id_map, syn_hits)

  ## ── pretty per-alias report ────────────────────────────────────────────────
  if (nrow(syn_hits)) {
    warn_df <- unique(syn_hits[ , c("external_synonym", "external_gene_name")])
    
    cat("[WARN] Aliases replaced:\n")
    invisible(
    apply(warn_df, 1, function(x)
      cat(sprintf("       %s  →  %s\n",
                  x["external_synonym"], x["external_gene_name"])))
    )
  }
  
  }
## 4. clean & drop biotype ----------------------------------------------------
keep <- id_map$gene_biotype=="protein_coding" &
  id_map$chromosome_name %in% c(as.character(1:22),"X","Y") # no MT - mitochondrial here
gene_coords <- id_map[keep, ]
gene_coords <- gene_coords[!duplicated(gene_coords$external_gene_name), ]
gene_coords$gene_biotype <- NULL
gene_coords$external_synonym <- NULL
gene_coords <- gene_coords[order(gene_coords$external_gene_name), ]
#dim(gene_coords) # just for checking
# Save for inspection
write_csv(gene_coords, "genes_with_coords.csv")

## --- make and sort the BED ---------------------------------------------
bed <- data.frame(
  chr   = paste0("chr", gene_coords$chromosome_name),
  start = gene_coords$start_position - 1,        # BED is 0-based
  end   = gene_coords$end_position,              # end is 1-based (not inclusive)
  name  = gene_coords$external_gene_name,
  score = gene_coords$entrezgene_id,
  strand= ".",                                   # placeholder
  stringsAsFactors = FALSE
)

bed <- bed[order(bed$chr, bed$start), ]          # ← coordinate sort
write.table(bed,
            file = "genes_of_interest.bed",
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)







