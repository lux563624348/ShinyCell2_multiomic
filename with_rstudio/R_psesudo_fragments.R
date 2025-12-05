## create psesudo_fragments from .h5mu
library(MuDataSeurat)
HT <- ReadH5MU("/data/HT.h5mu")

library(Matrix)
# Example: bin matrix
bin_mat <- HT@assays[["atac_cell_by_bin"]]@counts  # sparse matrix
bins <- rownames(bin_mat)
cells <- colnames(bin_mat)

# Parse bin coordinates
bin_coords <- do.call(rbind, strsplit(bins, " "))
colnames(bin_coords) <- c("chr", "bin_num", "start")
bin_coords <- as.data.frame(bin_coords)
bin_coords$start <- as.numeric(bin_coords$start)
bin_coords$end <- bin_coords$start + 499  # assume 500bp bins

# Open a file to write pseudo-fragments
frag_file <- "/data/HT_pseudo_fragments.tsv"
con <- file(frag_file, "w")

# Iterate over non-zero counts only
nz <- which(bin_mat != 0, arr.ind = TRUE)

max_frags <- 1000000L   # stop after writing this many
frag_count <- 0L

for(i in seq_len(nrow(nz))){
  r <- nz[i, 1]
  c <- nz[i, 2]
  count <- bin_mat[r, c]

  for(j in seq_len(count)){
    start <- sample(bin_coords$start[r]:bin_coords$end[r], 1)
    end <- start + 1  # minimal fragment length
    line <- paste(bin_coords$chr[r], start, end, cells[c], sep="\t")
    writeLines(line, con)

    frag_count <- frag_count + 1L
    if(frag_count >= max_frags){
      message("Reached ", max_frags, " fragments. Stopping early.")
      break
    }
  }
close(con)
  if(frag_count >= max_frags) break
}
