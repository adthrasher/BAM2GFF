#!/usr/bin/env Rscript
#============================================================================
#=================================HEATMAP 3 CODE=============================
#============================================================================

heatmap.4 <- function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,"Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      xoption = c("none", "promoters", "genebody"),
                      hline = median(breaks),
                      vline = median(breaks),
                      margins = c(5,5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction=0.3,
                      cexRow = 0.2 + 1/log10(nr),
                      cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      NumColSideColors = 1,
                      NumRowSideColors = 1,
                      KeyValueName="Value",...){
  
  invalid <- function (x) {
    if (missing(x) || is.null(x) || length(x) == 0)
      return(TRUE)
    if (is.list(x))
      return(all(sapply(x, invalid)))
    else if (is.vector(x))
      return(all(is.na(x)))
    else return(FALSE)
  }
  
  x <- as.matrix(x)
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  xoption <- match.arg(xoption)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    
    if (!missing(ColSideColors)) {
      #if (!is.matrix(ColSideColors))
      #stop("'ColSideColors' must be a matrix")
      if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
        stop("'ColSideColors' must be a matrix of nrow(x) rows")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      #lhei <- c(lhei[1], 0.2, lhei[2])
      lhei=c(lhei[1], side.height.fraction*NumColSideColors, lhei[2])
    }
    
    if (!missing(RowSideColors)) {
      #if (!is.matrix(RowSideColors))
      #stop("'RowSideColors' must be a matrix")
      if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
        stop("'RowSideColors' must be a matrix of ncol(x) columns")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
      #lwid <- c(lwid[1], 0.2, lwid[2])
      lwid <- c(lwid[1], side.height.fraction*NumRowSideColors, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)){
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    } else {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc = t(RowSideColors[,rowInd, drop=F])
      rsc.colors = matrix()
      rsc.names = names(table(rsc))
      rsc.i = 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] = rsc.name
        rsc[rsc == rsc.name] = rsc.i
        rsc.i = rsc.i + 1
      }
      rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      if (length(colnames(RowSideColors)) > 0) {
        axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), colnames(RowSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  if (!missing(ColSideColors)) {
    
    if (!is.matrix(ColSideColors)){
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } else {
      par(mar = c(0.5, 0, 0, margins[2]))
      csc = ColSideColors[colInd, , drop=F]
      csc.colors = matrix()
      csc.names = names(table(csc))
      csc.i = 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
      }
      csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
      image(csc, col = as.vector(csc.colors), axes = FALSE)
      if (length(colnames(ColSideColors)) > 0) {
        axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.5)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (xoption == "promoters") {
    axis(1, at=c(1,50,100), labels=c(paste("-",distance,"kb",sep=""), "TSS", paste("+",distance,"kb",sep="")))
  }
  if (xoption == "genebody") {
    axis(1, at=c(1,50,83,116,150,198.5), labels=c(paste("-",distance,"kb",sep=""), "TSS", "33%","66%", "TES", paste("+",distance,"kb",sep="")));
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    min.raw <- min(x, na.rm = TRUE)
    max.raw <- max(x, na.rm = TRUE)
    
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    mtext(side = 1, text=NULL, line = 2)
    title("Color Key")
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
  
}

subtract_matrices <- function(sample_mx, control_mx) {
    num_of_rows = nrow(sample_mx)
    num_of_cols = ncol(sample_mx)
    final_matrix = matrix(, nrow = num_of_rows, num_of_cols)
    for(row in 1:num_of_rows) {
        for(col in 1:num_of_cols) {
            value <- sample_mx[row, col] - control_mx[row, col]
            if (value < 0) value = 0  
            final_matrix[row, col] = value
        }
    }
    return(final_matrix)
}


#============================================================================
#=================================PLOT HEATMAPS==============================
#============================================================================

suppressPackageStartupMessages(require(optparse))

option_list <- list(
    make_option(c("-f", "--folder"), type='character', default="matrix",
        help="Folder containing matrix files generated"),
    make_option(c("-c", "--control"), type='character', default=NA,
        help="Folder containing control/input matrix files generated"),
    make_option(c("-z", "--zip"), action = "store_true", default=FALSE,
        help="Folder(s) provided are in zipped format"),
    make_option(c("-n", "--name"), type='character', default=NA, 
        help="Sample Name"),
    make_option(c("-d", "--distance"), type='integer', default=2000,
        help="Distance (bp) from TSS/TES")
    )

opt = parse_args(OptionParser(option_list=option_list))
if (is.na(opt$n)) {
    stop("Sample Name must be provided. See script usage (--help)")
}

#ARGS
folder = opt$f
unzipped_folder = "UNZIPPED"
samplename = opt$n
distance = round(opt$d/1000,1)

#library(animation)
# animation package is not used because it has ImageMagick dependency which is complicated to install on ubuntu docker. 
# using pdftools

library(pdftools)

if (opt$z) {
    unzip(folder,exdir=unzipped_folder)
    sample_folder = unzipped_folder
} else {
    sample_folder = folder
}

#sample files
s_promoters <- read.table(paste(sample_folder,"/",(dir(sample_folder,pattern="*promoters.txt"))[1],sep=""), sep="\t", header=T);
s_upstream <- read.table(paste(sample_folder,"/",(dir(sample_folder,pattern="*upstream.txt"))[1],sep=""), sep="\t", header=T);
s_downstream <- read.table(paste(sample_folder,"/",(dir(sample_folder,pattern="*downstream.txt"))[1],sep=""), sep="\t", header=T);
s_genebody <- read.table(paste(sample_folder,"/",(dir(sample_folder,pattern="*genebody.txt"))[1],sep=""), sep="\t", header=T);
if (opt$z) { unlink(unzipped_folder, recursive=TRUE) }

s_promoters=s_promoters[,3:ncol(s_promoters)]
s_upstream=s_upstream[,3:ncol(s_upstream)]
s_downstream=s_downstream[,3:ncol(s_downstream)]
s_genebody=s_genebody[,3:ncol(s_genebody)]

if (!is.na(opt$c)) {
    #input/control file
    control = opt$c
    if (opt$z) {
        unzip(control,exdir=unzipped_folder)
        control_folder = unzipped_folder
    } else {
        control_folder = control
    }

    #input/control files
    c_promoters <- read.table(paste(control_folder,"/",(dir(control_folder,pattern="*promoters.txt"))[1],sep=""), sep="\t", header=T);
    c_upstream <- read.table(paste(control_folder,"/",(dir(control_folder,pattern="*upstream.txt"))[1],sep=""), sep="\t", header=T);
    c_downstream <- read.table(paste(control_folder,"/",(dir(control_folder,pattern="*downstream.txt"))[1],sep=""), sep="\t", header=T);
    c_genebody <- read.table(paste(control_folder,"/",(dir(control_folder,pattern="*genebody.txt"))[1],sep=""), sep="\t", header=T);
    if (opt$z) { unlink(unzipped_folder, recursive=TRUE) }

    c_promoters=c_promoters[,3:ncol(c_promoters)]
    c_upstream=c_upstream[,3:ncol(c_upstream)]
    c_downstream=c_downstream[,3:ncol(c_downstream)]
    c_genebody=c_genebody[,3:ncol(c_genebody)]

    #subtracting matrices
    promoters = subtract_matrices(s_promoters, c_promoters)
    upstream = subtract_matrices(s_upstream, c_upstream)
    downstream = subtract_matrices(s_downstream, c_downstream)
    genebody = subtract_matrices(s_genebody, c_genebody)

} else {
    #only sample file provided
    promoters=s_promoters
    upstream=s_upstream
    downstream=s_downstream
    genebody=s_genebody

}

#combining entire genebody
combined<-cbind(upstream, genebody, downstream);

#removing NA.
promoters<-na.omit(promoters)
combined<-na.omit(combined)

#matplot of promoters & genebody
pdf(paste(samplename, "-promoters.pdf",sep=""))
matplot(colMeans(promoters), type='l', main=paste(samplename, "Promoters",sep=" "), ylab="Average normalized mapped reads", 
    xlim=NULL, xaxt='n', xlab="Genomic Region (bp)");
axis(1, at=c(0,50,100), labels=c(paste("-",distance,"kb",sep=""), "TSS", paste("+",distance,"kb",sep="")))
dev.off();

pdf(paste(samplename, "-entiregene.pdf",sep=""));
matplot(colMeans(combined),type='l', main=paste(samplename, "MetaGenes",sep=" "), ylab="Average normalized mapped reads", 
    xlim=NULL, xaxt='n', xlab="Genomic Region (bp)");
axis(1, at=c(0,50,83,116,150,200), labels=c(paste("-",distance,"kb",sep=""), "TSS", "33%","66%", "TES", paste("+",distance,"kb",sep="")));
dev.off();

#heatmap of promoters & genebody

#extrapolate breaks & colz
remainder = round((quantile(as.vector(t(promoters)),.80)),digits=0) %% 2;
finalcount = round((quantile(as.vector(t(promoters)),.80)),digits=0) + remainder;
if (finalcount < 2) { finalcount = 2; } #adjusting for lack of variability in bam density scores
breaks=seq(0,finalcount,by=(finalcount/100));
colz=colorRampPalette(c("white", "red"))(length(breaks)-1);

pdf(paste(samplename, "-heatmap.promoters.pdf", sep=""))
heatmap.4(promoters, col=colz, breaks=breaks, dendrogram="none", Colv=NA, Rowv=NA, labRow=NA, labCol=NA, xoption="promoters", xlab="Genomic Region (bp)", main=paste(samplename, "Promoters",sep="\n"))
dev.off()

png::writePNG(pdf_render_page(paste(samplename, "-heatmap.promoters.pdf", sep=""),page=1,dpi=300), paste(samplename, "-heatmap.promoters.png2", sep=""))
jpeg::writeJPEG(pdf_render_page(paste(samplename, "-heatmap.promoters.pdf", sep=""),page=1,dpi=300), paste(samplename, "-heatmap.promoters.jpg", sep=""))
#im.convert(paste(samplename, "-heatmap.promoters.pdf", sep=""), output = paste(samplename, "-heatmap.promoters.jpg", sep=""), extra.opts="-density 300")

remainder = round((quantile(as.vector(t(combined[,3:ncol(combined)])),.80)),digits=0) %% 2
finalcount = round((quantile(as.vector(t(combined[,3:ncol(combined)])),.80) + remainder),digits=0) + remainder
if (finalcount < 2) { finalcount = 2; } #adjusting for lack of variability in bam density scores
breaks=seq(0,finalcount,by=(finalcount/100));
colz=colorRampPalette(c("white", "red"))(length(breaks)-1);

pdf(paste(samplename, "-heatmap.entiregene.pdf", sep=""))
heatmap.4(combined, col=colz, breaks=breaks, dendrogram="none", Colv=NA, Rowv=NA, labRow=NA, labCol=NA, xoption="genebody", xlab="Genomic Region (bp)", main=paste(samplename, "MetaGenes",sep="\n"))
dev.off()

png::writePNG(pdf_render_page(paste(samplename, "-heatmap.entiregene.pdf", sep=""),page=1,dpi=300), paste(samplename, "-heatmap.entiregene.png2", sep=""))
jpeg::writeJPEG(pdf_render_page(paste(samplename, "-heatmap.entiregene.pdf", sep=""),page=1,dpi=300), paste(samplename, "-heatmap.entiregene.jpg", sep=""))
#im.convert(paste(samplename, "-heatmap.entiregene.pdf", sep=""), output = paste(samplename, "-heatmap.entiregene.jpg", sep=""), extra.opts="-density 300")
