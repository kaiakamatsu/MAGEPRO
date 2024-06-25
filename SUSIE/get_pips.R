library("optparse")
library("susieR")
library("data.table")

suppressPackageStartupMessages({
    library("Rfast")
})


arg_parser <- function() {
    option_list <- list(
        make_option(c("-g", "--genes_folder"), type = "character",
                    help = "Path to the genes folder", metavar = "GENESFOLDER"),
        make_option(c("-c", "--current_gene"), type = "character",
                    help = "Identifier of the current gene", metavar = "CURRENTGENE"),
        make_option(c("-l", "--ld_matrix"), type = "character",
                    help = "Path to the LD matrix file", metavar = "LDMATRIX"),
        make_option(c("-n", "--n_of_people"), type = "integer",
                    help = "Number of people/sample size", metavar="NSUSIE"),
        make_option(c("-o", "--output_folder", type="character", help="output folder name", metavar="OUTPUTFOLDER"))
    )

    opt_parser <- OptionParser(option_list = option_list)
    opt <- parse_args(opt_parser)

    # Check for required options
    if (is.null(opt$genes_folder)) {
        stop("Error: Genes folder (-g, --genes_folder) is required", call. = FALSE)
    }
    if (is.null(opt$current_gene)) {
        stop("Error: Current gene txt (-c, --current_gene) is required", call. = FALSE)
    }
    if (is.null(opt$ld_matrix)) {
        stop("Error: LD Matrix (-l, --ld_matrix) is required", call. = FALSE)
    }
    if (is.null(opt$n_of_people)) {
        stop("Error: Number of people/sample size (-n, --n_of_people) is required", call. = FALSE)
    }
    if (is.null(opt$output_folder)) {
        stop("Error: Output folder (-o, --output_folder) is required", call. = FALSE)
    }

    return(opt)
}
opt <- arg_parser()
GENES <- opt$g
R2 <- opt$l
CGENE = paste0(GENES, opt$c)

# read files
df <- fread(CGENE, header=TRUE, sep=" ", dec=".")       # read in current gene
colnames(df) = c("Gene", "SNP", "A1", "A2", "b", "SE", "P")
ld_stats <- fread(R2, header=TRUE, sep=" ", dec=".")    # read in ld_stats

# get unique SNPs that are both in ld matrix and gene file
df <- df[df[[2]] %in% (union(ld_stats$SNP_A, ld_stats$SNP_B))]
snp_list <- df[[2]]

# build the ld matrix
ld_matrix <- matrix(NA, nrow = length(snp_list), ncol = length(snp_list), dimnames = list(snp_list, snp_list))
index_A <- match(ld_stats$SNP_A, snp_list)
index_B <- match(ld_stats$SNP_B, snp_list)

# fill in the matrix
ld_matrix[cbind(index_A, index_B)] <- ld_stats$R2
ld_matrix[cbind(index_B, index_A)] <- ld_stats$R2
diag(ld_matrix) <- 1

# run susie
res <- susie_rss(bhat=df[[5]], shat=df[[6]], R=ld_matrix, n=opt$n, max_iter=250)
pips <- res$pip

# create credible set column
cs <- rep(-1, length(df[[5]]))
for(cname in names(res$sets$cs)) {
  number <- as.integer(gsub("L", "", cname))
  indices <- res$sets$cs[[cname]]
  cs[indices] <- number
}

# writing the table
df$pip <- pips
betas <- coef(res)[-1]  # coef.susie – extract regression
df$beta <- betas
df$cs <- cs
output = opt$o
df <- df[order(df$cs, decreasing = TRUE), ]

colnames(df) = c("Gene", "SNP", "A1", "A2", "b", "SE", "P", "pip", "beta", "cs")
fwrite(x=df, file=output, sep=" ", col.names=TRUE, row.names=FALSE, quote=FALSE)