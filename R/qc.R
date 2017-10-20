'%!in%' <- function(x,y)!('%in%'(x,y))

catverbose <- function(...){
  cat(format(Sys.time(), "%Y%m%d %H:%M:%S |"), ..., "\n")
}

plot_vaf_by_cn_state <- function(maf, sample, wgd=F){
  if ('mcn' %!in% names(maf)){
    maf[, mcn := tcn-lcn]
  }
  maf.tmp <- maf[patient == sample & !is.na(mcn) & mcn <= 6]
  phi <- unique(maf.tmp[!is.na(purity)]$purity)
  if(length(phi) == 0){
    catverbose("No FACETS annotations!")
  } else {
    gg <- ggplot(maf.tmp, aes(x=VAF)) +
      geom_histogram(col="black", fill="#41B6C4", lwd=1.5, binwidth = 0.02) +
      geom_vline(xintercept=(phi/2), linetype=2, color = "#FB6A4A") +
      facet_grid(lcn ~ mcn) +
      xlab("Variant Allele Fraction") +
      ylab("Frequency") +
      ggtitle(sample) +
      theme_bw() +
      theme(plot.title=element_text(size=25, face = "bold"),
            axis.title=element_text(size=20, face = "bold"),
            strip.text.x=element_text(size=20, face = "bold"),
            strip.text.y=element_text(size=20, face = "bold"),
            axis.text.x=element_text(size=15, angle=45, hjust=1),
            axis.text.y=element_text(size=15),
            legend.text=element_text(size=15),
            legend.title=element_text(size=15))
    if(wgd){
      gg <- gg + geom_vline(xintercept=(phi/4), linetype=2, color = "#FD8D3C")
    }
    plot(gg)
  }
}

FGA <- function(cncf){
  cncf[, total_seglen := sum(loc.end - loc.start), Tumor_Sample_Barcode]
  FGA <- cncf[, .(FGA = sum(loc.end - loc.start) / unique(total_seglen)),
              keyby = .(Tumor_Sample_Barcode, tcn = tcn, lcn = lcn)]
  FGA <- FGA[order(FGA)]
  FGA[, mcn := tcn - lcn]
  FGA[, f := 1:.N/.N, .(tcn, lcn)]
  copy(FGA)
}

qc <- function(maf, cncf){

  # cncf[, total_seglen := sum(na.rm = T, as.numeric(loc.end) - as.numeric(loc.start)),
  #      Tumor_Sample_Barcode]

  QC <- cncf[, .(
    n.dip.seg = .SD[tcn == 2 & lcn == 1, .N],
    loh = .SD[lcn == 0, sum(na.rm = T, as.numeric(loc.end) - as.numeric(loc.start))] / sum(na.rm = T, as.numeric(loc.end) - as.numeric(loc.start)),
    f_hi_mcn = .SD[(tcn - lcn) >= 2, sum(na.rm = T, as.numeric(loc.end) - as.numeric(loc.start))] / sum(na.rm = T, as.numeric(loc.end) - as.numeric(loc.start))
  ),
  keyby = Tumor_Sample_Barcode]
  QC
}
