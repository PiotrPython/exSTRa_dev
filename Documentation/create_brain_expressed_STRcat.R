# -------------------------------------------------------------------
# Overlap exome capture regions with STR catalogue
# -------------------------------------------------------------------

rpkm_expression_thresh <- 10
pad_gene_region <- 1000

exome_capture_bed <- "nexterarapidcapture_exome_targetedregions.bed"
pad_exome_capture_region <- 100

# -------------------------------------------------------------------

setwd("/home/users/allstaff/bennett.ma/work/Projects/STR/brain_expressed_STRcat")

### Get list of genes expressed in brain (GTEx)

# Load GTEx
library(readr)
GTEx_file <- "GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct.gz"
GTEx <- read_tsv(GTEx_file, skip=2)  # first two lines containing table version and dimensions

# Keep genes with median rpkm > 1 in at least one brain tissue
GTEx_brain <- GTEx[, startsWith(colnames(GTEx), "Brain")]
brain_expressed <- rowSums(GTEx_brain > rpkm_expression_thresh) > 0

# Convert Ensembl IDs to gene IDs and symbols
library(org.Hs.eg.db)
ensDB <- as.data.frame(org.Hs.egENSEMBL)
name <- GTEx$Name[brain_expressed]
name_no_dot <- sapply(name, function(x){strsplit(x, ".", fixed=TRUE)[[1]][1]})
geneid <- ensDB$gene_id[na.omit(match(name_no_dot, ensDB$ensembl_id))]


### Get coordinates of gene (chr, start, end) [+ pad regions]

# Get all gene coordinates
library(Homo.sapiens)
gene_coords <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
gene_coords_df <- data.frame(geneid=gene_coords$gene_id, chr=as.character(seqnames(gene_coords)), stringsAsFactors=FALSE)
gene_coords_df <- cbind(gene_coords_df, as.data.frame(ranges(gene_coords)))

# Write bedfile
bed <- gene_coords_df[gene_coords_df$geneid %in% geneid, c("chr", "start", "end")]
bed$start <- bed$start - pad_gene_region
bed$end <- bed$end + pad_gene_region
# Sort bedfile
bed <- bed[order(factor(bed$chr, levels=paste0("chr", c(1:22, "X", "Y"))), bed$start, bed$end), ]
write.table(bed, file=paste0("brain_expressed_genes_rpkm", rpkm_expression_thresh, "_pad", pad_gene_region, ".bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)


### Intersect brain expressed genes and STR catalogue
system(paste0("bedtools intersect -a str_reference_fix.bed -b brain_expressed_genes_rpkm", rpkm_expression_thresh, "_pad", pad_gene_region, ".bed > brain_expressed_rpkm", rpkm_expression_thresh, "_padGene", pad_gene_region, "_STRcat.bed"))


### Read exome capture region bed and write new bed file with extra padding

exome_bed <- read.delim(exome_capture_bed, sep="\t", header=FALSE)
colnames(exome_bed) <- c("chr", "start", "stop", "name")
# Pad exonic regions
exome_bed$start <- exome_bed$start - pad_exome_capture_region
exome_bed$stop <- exome_bed$stop + pad_exome_capture_region
# Merge and remove duplicate rows which now overlap after padding
exome_bed_merge <- (exome_bed$start[-1] < exome_bed$stop[-nrow(exome_bed)]) & (exome_bed$chr[-1] == exome_bed$chr[-nrow(exome_bed)])
exome_bed$start[c(FALSE, exome_bed_merge)] <- exome_bed$start[c(exome_bed_merge, FALSE)]
exome_bed$stop[c(exome_bed_merge, FALSE)] <- exome_bed$stop[c(FALSE, exome_bed_merge)]
exome_bed <- exome_bed[!duplicated(exome_bed[, 1:3]), ]
# Write padded bedfile
write.table(exome_bed, file=gsub(".bed", paste0("_pad", pad_exome_capture_region, ".bed"), exome_capture_bed, fixed=TRUE), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# Intersect brain expressed genes and STR catalogue
system(paste0("bedtools intersect -a  brain_expressed_rpkm", rpkm_expression_thresh, "_padGene", pad_gene_region, "_STRcat.bed -b ", gsub(".bed", paste0("_pad", pad_exome_capture_region, ".bed"), exome_capture_bed, fixed=TRUE), "> brain_expressed_rpkm", rpkm_expression_thresh, "_padGene", pad_gene_region, "_padExomeCapture", pad_exome_capture_region, "_STRcat.bed"))


### Convert to format required by exSTRa

# Load exome STR catalogue
STRcat <- read.delim(paste0("brain_expressed_rpkm", rpkm_expression_thresh, "_padGene", pad_gene_region, "_padExomeCapture", pad_exome_capture_region, "_STRcat.bed"), sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(STRcat) <- c("chr", "start", "stop", "period", "num_reps", "X1", "start_original", "stop_original", "TRFscore", "PercentA", "PercentC", "PercentG", "PercentT", "Entropy", "Motif", "Disease")
STRcat <- STRcat[!duplicated(STRcat), ]

dim(STRcat)
# [1] 5048     16

# Format output to create Excel sheet needed to run exSTRa
excel_STRcat <- data.frame(Disease=paste0(STRcat$chr, ":", STRcat$start_original, "-", STRcat$stop_original, "_", STRcat$Motif))
excel_STRcat$"OMIM number" <- 0
excel_STRcat$"Mode of inheritance" <- "NA"
excel_STRcat$Gene <- "NA"
excel_STRcat$"Gene location" <- "NA"
excel_STRcat$strcat_all <- paste0("http://strcat.teamerlich.org/chart/", STRcat$chr, "/", STRcat$start_original, "/", STRcat$stop_original)
excel_STRcat$"hg19 chrom" <- STRcat$chr
excel_STRcat$"hg19 start 0" <- STRcat$start_original - 1
excel_STRcat$"hg19 end" <- STRcat$stop_original
excel_STRcat$strand <- "+"
excel_STRcat$copyNum <- STRcat$num_reps
excel_STRcat$perMatch <- 100  # ?
excel_STRcat$perIndel <- 0    # ?
excel_STRcat$STR_size <- excel_STRcat$copyNum * nchar(STRcat$Motif)
excel_STRcat$read_detect_size <- excel_STRcat$STR_size
excel_STRcat$"Location of repeat within gene" <- "?"
excel_STRcat$"Repeat sequence" <- STRcat$Motif
excel_STRcat$"Stable repeat number" <- "?"
excel_STRcat$"Unstable repeat number" <- "?"
excel_STRcat$"Unpathogenic size range" <- "?"
excel_STRcat <- excel_STRcat[!duplicated(excel_STRcat), ]

write.table(excel_STRcat, file=paste0("brain_expressed_rpkm", rpkm_expression_thresh, "_padGene", pad_gene_region, "_padExomeCapture", pad_exome_capture_region, "_STRcat.xls"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


