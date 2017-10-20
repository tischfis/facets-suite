out_file_table <- function(filenames){
  out <- fread_rbind(
    filenames,
    fn = function(filename)
      fread(
        paste0("(sed 's/^#//' | grep -v '^$' | grep -v '^ $') < ", filename),
        sep = "=", header = F, fill = T
      )
  )
  setnames(out,c("Tumor_Sample_Barcode", "variable", "value"))
  out <- out[!variable %in% c("INPUT PARAMETERS GIVEN",
                              "LOADED MODULE INFO",
                              "FACETS OUTPUT")]
  out
}


generate_cncf_file <- function(samples){
  # samples[, cncf_exists := file.exists(cncf), cncf]

  cncf_filenames <- samples[, structure(cncf, .Names = Tumor_Sample_Barcode)]

  ## print(length(cncf_filenames))

  cncf <- rbindlist(
    idcol = "Tumor_Sample_Barcode",
    sapply(
      USE.NAMES = T,
      simplify = F,
      cncf_filenames,
      fread
    ))
  out <- rbindlist(
    idcol = "Tumor_Sample_Barcode",
    sapply(
      USE.NAMES = T,
      simplify = F,
      gsub(".cncf.txt$", ".out",  cncf_filenames),
      out_file_table
    ))
  out[, 2 := NULL]
  out_table <- dcast.data.table(
    out[variable %in% c("Purity", "Ploidy", "dipLogR", "dipt", "loglik")],
    Tumor_Sample_Barcode ~ variable,
    value.var = "value"
  )
  if("Purity" %in% names(out_table)) setnames(out_table, "Purity", "purity")
  if("Ploidy" %in% names(out_table)) setnames(out_table, "Ploidy", "ploidy")
  out_table[, purity := as.numeric(purity)]
  out_table[, ploidy := as.numeric(ploidy)]
  out_table[, dipLogR := as.numeric(dipLogR)]
  out_table[, loglik := as.numeric(loglik)]

  cncf_out <- merge(cncf, out_table, by = "Tumor_Sample_Barcode", all = T, allow.cartesian=TRUE)
}

generate_logR_seg <- function(cncf){
  cncf[, .(
    ID = Tumor_Sample_Barcode,
    chrom,
    loc.start,
    loc.end,
    num.mark,
    seg.mean = cnlr.median)]
}

generate_logR_adj_seg <- function(cncf){
  cncf[, .(
    ID = Tumor_Sample_Barcode,
    chrom,
    loc.start,
    loc.end,
    num.mark,
    seg.mean = cnlr.median - dipLogR)]
}

generate_log_tcn_seg <- function(cncf){
  cncf[, .(
    ID = Tumor_Sample_Barcode,
    chrom,
    loc.start,
    loc.end,
    num.mark,
    log_tcn = ifelse(
      is.na(purity) | cf / purity > 0.5,
      log(tcn / 2),
      0))]
}


