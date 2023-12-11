#!/usr/bin/env Rscript
# Usage
# Rscript extract_signatures.R --bs_genome <hg19/hg38> --fname_maf <FNAME> --fname_out <FNAME>

# Import packages
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
library(deconstructSigs)
library(getopt)
library(stringr)
library(tidyverse)

# Argument parsing
spec = matrix(
  c(
    "bs_genome",
    "g",
    1,
    "character",
    "fname_maf",
    "i",
    1,
    "character",
    "fname_out",
    "o",
    1,
    "character"
  ),
  byrow = TRUE,
  ncol = 4
)
opt = getopt(spec)
if (is.null(opt$bs_genome)) {
  opt$bs_genome = "hg19"
}

# Variables and constants
if (opt$bs_genome == "hg38") {
  bs_genome <- BSgenome.Hsapiens.UCSC.hg38
} else if (opt$bs_genome == "hg19") {
  bs_genome <- BSgenome.Hsapiens.UCSC.hg19
} else {
  str_interp("Invalid genome selected. Must be hg38/hg19. '{opt$bs_genome}' is not.")
}
signature_tobacco <- 4
signature_uv <- 7
signature_apobec <- c(2, 13)

# Read files
raw_maf <- read.delim(opt$fname_maf)
skip_rows <- which(!startsWith(raw_maf[, 1], "#"))[1]
mutdat <- read.delim(opt$fname_maf, skip=skip_rows)
mutdat$Tumor_Sample_Barcode = as.character(mutdat$Tumor_Sample_Barcode)
samples <- unique(mutdat["Tumor_Sample_Barcode"])
if (length(samples) > 1) {
  warning("Multiple samples found!")
}

# Extract signatures
mutsig <- mut.to.sigs.input(
  mut.ref = mutdat,
  sample.id = "Tumor_Sample_Barcode",
  chr = "Chromosome",
  pos = "Start_Position",
  ref = "Reference_Allele",
  alt = "Tumor_Seq_Allele2",
  bsg = bs_genome
)

signatures <- whichSignatures(
  tumor.ref = mutsig,
  signatures.ref = signatures.cosmic,
  sample.id = samples[1, 1],
  contexts.needed = TRUE
)
weights <- signatures[["weights"]]

# Save signatures
output <-
  matrix(
    c(
      # "fname_maf",
      # opt$fname_maf,
      # "bs_genome",
      # opt$bs_genome,
      "signature_apobec1",
      as.double(weights[signature_apobec[1]]),
      "signature_apobec2",
      as.double(weights[signature_apobec[2]]),
      "signature_tobacco",
      as.double(weights[signature_tobacco]),
      "signature_uv",
      as.double(weights[signature_uv])
    ),
    byrow = TRUE,
    ncol = 2
  )
write.table(
  output,
  file = opt$fname_out,
  sep = ":",
  row.names = FALSE,
  col.names = FALSE,
  quote=F
)
