#!/usr/bin/env Rscript

VERSION="0.0.1"


suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("circlize"))
suppressPackageStartupMessages(library("ape"))
suppressPackageStartupMessages(library("dendextend"))
suppressPackageStartupMessages(library("tidyr"))

option_list <- list(
    make_option(
        c("-f", "--faidx"),
        type="character",
        action="store",
        help="The fasta index (required)."
    ),
    make_option(
        c("-b", "--bedgraph"),
        type="character",
        action="store",
        help="The bedgraph (required)."
    ),
    make_option(
        c("-d", "--outdir"),
        type="character",
        action="store",
        help="Where to store the resulting newick tree (required)."
    ),
    make_option(
        "--version",
        type="logical",
        action="store_true",
        default=FALSE,
        help="Print version and exit.",
    )
)

parser <- OptionParser(
    usage = "%prog --bedgraph my.bg --faidx my.fasta.fai --outdir out",
    option_list = option_list
)

args <- parse_args(parser)

log_stderr <- function(...) {
  cat(sprintf(...), sep='', file=stderr())
}

quit_with_err <- function(...) {
  log_stderr(...)
  quit(save = "no", status = 1, runLast = FALSE)
}

validate_file <- function(path) {
  if (is.null(path)) {
    quit_with_err("Please provide required file")
  }
}

stripper <- function(str) {
  str <- gsub("^X", "", str, perl = TRUE)
  str
}

main <- function(args) {
  if (args$version) {
    cat(VERSION, file=stdout())
    quit(save = "no", status = 0, runLast = FALSE)
  }

  validate_file(args$faidx)
  validate_file(args$bedgraph)
  validate_file(args$outdir)

  faidx <- read.table(
    args$faidx,
    stringsAsFactors = FALSE,
    col.names = c("seqid", "length", "J", "U", "N")
  )

  faidx["start"] <- 0
  faidx["end"] <- faidx["length"]
  faidx <- faidx[c("seqid", "start", "end")]

  bg <- read.table(
    args$bedgraph,
    header = TRUE,
    stringsAsFactors = FALSE,
    comment.char = ""
  )

  sample_names <- names(bg)[4:length(names(bg))]
  names(bg)[4:length(names(bg))] <- as.character(
    sapply(
      sample_names,
      FUN = stripper
    )
  )

  names(bg)[1:3] <- c("seqid", "start", "end")


  bg_mat <- bg
  row.names(bg_mat) <- apply(
    bg_mat,
    1,
    function(x) paste0(x["seqid"], ":", x["start"], "-", x["end"])
  )

  bg_mat["seqid"] <- NULL
  bg_mat["start"] <- NULL
  bg_mat["end"] <- NULL

  bg_mat <- t(bg_mat)

  bg_dist <- dist(bg_mat, method = "euclidean")
  bg_tree <- hclust(bg_dist, method = "average")
  bg_order <- bg_tree$labels[bg_tree$order]

  # Create the output directory
  if (!dir.exists(args$outdir)) {
    dir.create(args$outdir)
  }

  write.table(
    bg_order,
    file.path(args$outdir, "isolate_order.txt"),
    row.names=FALSE,
    col.names=FALSE
  )

  pdf(file.path(args$outdir, "pav_tree.pdf"), width = 20, height = 5)
  plot(bg_tree)
  dev.off()

  write.tree(as.phylo(bg_tree), file.path(args$outdir, "pav_tree.nwk"))

  pdf(file.path(args$outdir, "pav.pdf"), width = 9, height = 9)
  circos.par(
    track.height = 0.1,
    gap.degree = 2,
    start.degree = 90,
    canvas.xlim = c(-1.1, 1.1),
    canvas.ylim = c(-1.1, 1.1)
  )

  circos.genomicInitialize(
    faidx,
    major.by = 1000000,
    track.height = 0.05,
    sector.names = faidx[,"seqid"]
  )

  col_fun = colorRamp2(c(0, 1, 2), c("#0000FF", "white", "#FF0000"))

  circos.genomicHeatmap(
    bg[, c(colnames(bg)[1:3], bg_order)],
    col = col_fun,
    connection_height = 0.0001,
    line_lwd = 0,
    heatmap_height = 0.5
  )

  circos.clear()
  dev.off()
}

main(args)
