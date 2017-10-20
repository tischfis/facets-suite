FGA <- function(cncf){
  cncf[, total_seglen := sum(loc.end - loc.start), Tumor_Sample_Barcode]
  FGA <- cncf[, .(FGA = sum(loc.end - loc.start) / unique(total_seglen)),
              keyby = .(Tumor_Sample_Barcode, tcn = tcn, lcn = lcn)]
  FGA <- FGA[order(FGA)]
  FGA[, mcn := tcn - lcn]
  FGA[, f := 1:.N/.N, .(tcn, lcn)]
  copy(FGA)
}

cna_summary <- function(maf, cncf){

  # cncf[, total_seglen := sum(na.rm = T, as.numeric(loc.end) - as.numeric(loc.start)),
  #      Tumor_Sample_Barcode]

  cna_summary <- cncf[, .(
    ploidy = unique(ploidy),
    purity = unique(purity),
    loh = .SD[lcn == 0, sum(na.rm = T, as.numeric(loc.end) - as.numeric(loc.start))] / sum(na.rm = T, as.numeric(loc.end) - as.numeric(loc.start)),
    f_hi_mcn = .SD[(tcn - lcn) >= 2, sum(na.rm = T, as.numeric(loc.end) - as.numeric(loc.start))] / sum(na.rm = T, as.numeric(loc.end) - as.numeric(loc.start))
  ),
  keyby = Tumor_Sample_Barcode]
  cna_summary
}
