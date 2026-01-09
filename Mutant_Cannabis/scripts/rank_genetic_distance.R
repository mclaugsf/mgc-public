#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(SNPRelate))
suppressPackageStartupMessages(library(argparse))

# ---------------------------
# Setup argument parser
# ---------------------------
parser <- ArgumentParser(description="Rank genetic distances for a query sample using SNPRelate")

parser$add_argument("-v", "--vcf", required=TRUE, help="Input multi-sample VCF file")
parser$add_argument("-q", "--query", required=TRUE, help="Query sample name")
parser$add_argument("-g", "--gds", default="tmp.gds", help="Temporary GDS file (default: tmp.gds)")
parser$add_argument("-t", "--threads", type="integer", default=1, help="Number of threads (default: 4)")

args <- parser$parse_args()

# ---------------------------
# Convert VCF → GDS
# ---------------------------
snpgdsVCF2GDS(
  vcf.fn = args$vcf,
  out.fn = args$gds,
  method = "biallelic.only",
  ignore.chr.prefix = "chr"
)

# ---------------------------
# Open GDS and compute IBS
# ---------------------------
genofile <- snpgdsOpen(args$gds)

ibs <- snpgdsIBS(
  genofile,
  num.thread = args$threads,
  autosome.only = FALSE,
  remove.monosnp = FALSE,
  missing.rate = 1          # allow SNPs with missing genotypes
)

dist_mat <- 1 - ibs$ibs
sample_ids <- ibs$sample.id
rownames(dist_mat) <- sample_ids
colnames(dist_mat) <- sample_ids

# ---------------------------
# Rank distances for query sample
# ---------------------------
if(!(args$query %in% sample_ids)){
  stop("Query sample not found in VCF")
}

ranking <- data.frame(
  sample = sample_ids,
  distance = dist_mat[args$query, ]
)

ranking <- ranking[ranking$sample != args$query, ]
ranking <- ranking[order(ranking$distance), ]

# ---------------------------
# Output to stdout
# ---------------------------
write.table(
  ranking,
  file = "",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ---------------------------
# Close GDS
# ---------------------------
snpgdsClose(genofile)
