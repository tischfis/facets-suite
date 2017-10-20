get_arm_level_calls <- function(cncf,
                                WGD_threshold = 0.5, ### least value of frac_elev_major_cn for WGD
                                amp_threshold = 5, ### total copy number greater than this value for an amplification
                                mean_chrom_threshold = 0, ### total copy number also greater than this value multiplied by the chromosome mean for an amplification
                                fun.rename = function(filename){filename}
){

  cncf$chrom <- as.character(cncf$chrom)
  cncf[chrom == "23", chrom := "X"]
  setkey(cncf, chrom, loc.start, loc.end)

  if (!("tcn" %in% names(cncf))) {
    cncf[, c("tcn", "lcn") := list(tcn.em, lcn.em)]
  }

  ### estimate fraction of the genome with more than one copy from a parent
  ### a large fraction implies whole genome duplication
  cncf[, frac_elev_major_cn := sum(
    as.numeric((tcn - lcn) >= 2) *
      as.numeric(loc.end-loc.start), na.rm = T) /
      sum(as.numeric(loc.end-loc.start)
      ),
    by=Tumor_Sample_Barcode]

  arm_definitions <- facets.suite::arm_definitions
  ### Extract integer copy number for each probe from cncf
  fo_impact <- foverlaps(arm_definitions, cncf, nomatch=NA)
  ### Truncate segments that span two arms
  fo_impact[, loc.start := ifelse(loc.start < start, start, loc.start)]
  fo_impact[, loc.end := ifelse(loc.end > end, end, loc.end)]
  fo_impact[, arm := factor(arm, levels = levels(arm_definitions$arm))]
  #fo_impact[,Hugo_Symbol:=gsub("_.*$", "", name)]

  ### Summarize copy number for each gene
  arm_level <- fo_impact[,
                         list(frac_elev_major_cn=unique(frac_elev_major_cn),
                              Nsegments = .N,
                              length_CN = sum(as.numeric(loc.end - loc.start))),
                         by=list(Tumor_Sample_Barcode, arm, tcn=tcn, lcn=lcn)]

  setkey(arm_level, Tumor_Sample_Barcode, arm)
  ### for each CN status, calculate fraction of arm covered
  setkey(arm_definitions, arm)
  arm_level[, arm_frac := round(length_CN / arm_definitions[as.character(arm_level$arm), as.numeric(end - start)], digits = 4)]
  setkey(arm_definitions, chr, start, end)
  ### ignore CN status valuesif present in less than 10% of arm
  ## arm_level <- arm_level[arm_frac >= 0.1 ]
  arm_level <- arm_level[order(Tumor_Sample_Barcode, arm, -arm_frac)]
  ### estimate WGD from frac_elev_major_cn
  arm_level[, WGD := factor(ifelse(frac_elev_major_cn > WGD_threshold, "WGD", "no WGD"))]

  ### fix bug where lcn == NA even when tcn is 1
  arm_level[tcn == 1 & is.na(lcn), lcn := 0]


  ### annotate integer copy number
  arm_level <- annotate_integer_copy_number(arm_level, amp_threshold = amp_threshold)

  ### remove duplicate entries for partial deletions
  ### the lower value is chosen on the basis that a
  ### partial deletion is a deletion but a
  ### partial amplification is not an amplification
  #   arm_level <- arm_level[order(FACETS_CNA, -Nprobes)]
  #   arm_level <- arm_level[!duplicated(arm_level, by=c("Tumor_Sample_Barcode", "Hugo_Symbol"))]

  arm_level <- arm_level[order(Tumor_Sample_Barcode, arm, -arm_frac)]
  arm_level[,primary := as.integer(!duplicated(arm_level, by=c("Tumor_Sample_Barcode", "arm")))]

  setkey(arm_level, Tumor_Sample_Barcode, arm)
  arm_level[, tcn_summary := as.character(tcn)]
  arm_level[duplicated(arm_level) | duplicated(arm_level, fromLast = T),
            tcn_summary := paste(collapse = " ",
                                 paste(sep = ":",
                                       paste0(
                                         round(arm_frac*100, 0), "%"),
                                       tcn)),
            by = list(Tumor_Sample_Barcode, arm)]
  arm_level
}
