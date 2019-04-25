mod_cov.plot <- function (peak, weightCol = NULL, xlab = "Chromosome Size (bp)", 
          ylab = "", title = "ChIP Peaks over Chromosomes", chrs = NULL, 
          xlim = NULL, lower = 1, co) 
{
  isList <- FALSE
  if (is(peak, "GRanges") || length(peak) == 1) {
    tm <- getChrCov(peak = peak, weightCol = weightCol, chrs, 
                    xlim, lower = lower)
  }
  else {
    isList <- TRUE
    ltm <- lapply(peak, getChrCov, weightCol = weightCol, 
                  chrs = chrs, xlim = xlim, lower = lower)
    if (is.null(names(ltm))) {
      nn <- paste0("peak", seq_along(ltm))
      warning("input is not a named list, set the name automatically to ", 
              paste(nn, collapse = " "))
      names(ltm) <- nn
    }
    tm <- list_to_dataframe(ltm)
    chr.sorted <- sortChrName(as.character(unique(tm$chr)))
    tm$chr <- factor(tm$chr, levels = chr.sorted)
  }
  chr <- start <- end <- value <- .id <- NULL
  p <- ggplot(tm, aes(start, value))
  if (isList) {
    p <- p + geom_rect(aes(xmin = start, ymin = 0, xmax = end, 
                           ymax = value, fill = .id, color = .id), alpha=0.2)
  }
  else {
    p <- p + geom_rect(aes(xmin = start, ymin = 0, xmax = end, 
                           ymax = value), fill = co, color = co) + scale_x_continuous(expand = c(0,0))
  }
  if (length(unique(tm$chr)) > 1) {
    p <- p + facet_grid(chr ~ ., scales = "free")
  }
  p <- p + theme_classic()
  p <- p + xlab(xlab) + ylab(ylab) + ggtitle(title)
  p <- p + scale_y_continuous(expand = c(0, 0))
  p <- p + theme(strip.text.y = element_text(angle = 360))
  if (!is.null(xlim) && !is.na(xlim) && is.numeric(xlim) && 
      length(xlim) == 2) {
    p <- p + xlim(xlim)
  }
  return(p)
}






