annotate_integer_copy_number <- function(gene_level, amp_threshold = 5){
  ### Portal/TCGA copy number calls as well as
  ### more advanced copy number calling, including CNLOH etc.

  ### lowest value of tcn for AMP
  AMP_thresh_tcn <- amp_threshold + 1

  ### assumes that WGD variable has already been set
  input_key <- key(gene_level)
  input_cols <- names(gene_level)
  gene_level[, mcn := tcn - lcn]
  setkey(gene_level, WGD, mcn, lcn)

  if("FACETS_CNA" %in% names(gene_level)) gene_level[, FACETS_CNA := NULL]
  if("FACETS_CALL" %in% names(gene_level)) gene_level[, FACETS_CALL := NULL]

  gene_level <- merge(gene_level, facets.suite::FACETS_CALL_table, sort = F, all.x = T)

  gene_level[tcn >= AMP_thresh_tcn, FACETS_CNA := 2]
  gene_level[tcn >= AMP_thresh_tcn, FACETS_CALL := "AMP"]

  output_cols <- c(input_cols, c("FACETS_CNA", "FACETS_CALL"))
  setcolorder(gene_level, output_cols)
  setkeyv(gene_level, input_key)
  gene_level
}

get_gene_level_calls <- function(cncf,
                                 gene_targets,
                                 WGD_threshold = 0.5, ### least value of frac_elev_major_cn for WGD
                                 amp_threshold = 5, ### total copy number greater than this value for an amplification
                                 mean_chrom_threshold = 0, ### total copy number also greater than this value multiplied by the chromosome mean for an amplification
                                 fun.rename = function(filename){filename}){

  cncf[, chrom := as.character(cncf$chrom)]
  cncf[chrom == "23", chrom := "X"]
  setkey(cncf, chrom, loc.start, loc.end)

  ### estimate fraction of the genome with more than one copy from a parent
  ### a large fraction implies whole genome duplication
  cncf[, frac_elev_major_cn := sum(
    as.numeric((tcn - lcn) >= 2) *
      as.numeric(loc.end-loc.start), na.rm = T) /
      sum(as.numeric(loc.end-loc.start)
      ),
    by=Tumor_Sample_Barcode]

  ### Extract integer copy number for each probe from cncf
  fo_impact <- foverlaps(gene_targets, cncf, nomatch=NA)
  fo_impact <- fo_impact[!is.na(ID)]
  fo_impact[,Hugo_Symbol:=gsub("_.*$", "", name)]

  ### Summarize copy number for each gene
  gene_level <- fo_impact[!Hugo_Symbol %in% c("Tiling", "FP", "intron"),
                          list(chr = unique(chr),
                               seg.start=unique(loc.start),
                               seg.end=unique(loc.end),
                               cnlr.median = unique(cnlr.median),
                               #                                start=unique(start),   ### with these uncommented, the per-gene summarization is broken (??)
                               #                                end=unique(end),
                               frac_elev_major_cn=unique(frac_elev_major_cn),
                               Nprobes = .N),
                          keyby=list(Tumor_Sample_Barcode, Hugo_Symbol, tcn, lcn, cf)]

  ### fix bug where lcn == NA even when tcn is 1
  gene_level[tcn == 1 & is.na(lcn), lcn := 0]

  ### apply WGD threshold
  gene_level[, WGD := factor(ifelse(frac_elev_major_cn > WGD_threshold, "WGD", "no WGD"))]
  setkey(gene_level, Tumor_Sample_Barcode, Hugo_Symbol)

  ### focality requirement
  ### get (weighted) mean total copy number for the chromosome
  #   mean_tcn_per_chr <- cncf[, list(mean_tcn =
  #                                                sum(tcn * as.numeric(loc.end - loc.start)) /
  #                                                sum(as.numeric(loc.end - loc.start))),
  #                                       keyby=list(Tumor_Sample_Barcode, chrom)]
  #   mean_tcn_per_chr_v <- with(mean_tcn_per_chr,
  #                              structure(mean_tcn,
  #                                        .Names = paste(Tumor_Sample_Barcode,
  #                                                       chrom)))
  #   gene_level[, mean_tcn_per_chr := mean_tcn_per_chr_v[paste(Tumor_Sample_Barcode, chr)]]

  ### annotate integer copy number
  gene_level <- annotate_integer_copy_number(gene_level, amp_threshold)

  ### remove duplicate entries for partial deletions
  ### the lower value is chosen on the basis that a
  ### partial deletion is a deletion but a
  ### partial amplification is not an amplification
  gene_level <- gene_level[order(FACETS_CNA, -Nprobes)]
  gene_level <- gene_level[!duplicated(gene_level, by=c("Tumor_Sample_Barcode", "Hugo_Symbol"))]

  ### Sort & set column classes
  gene_level[, FACETS_CNA := factor(FACETS_CNA, levels = c(-2:2))]
  gene_level[, chr := factor(chr, levels = c(1:22, "X", "Y"))]
  setkey(gene_level, Tumor_Sample_Barcode, Hugo_Symbol)
  copy(gene_level)
}
