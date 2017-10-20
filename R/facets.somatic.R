### write tab delimited
write.tab <- function(...){
  write.table(..., quote = F,
              col.names=T,
              row.names=F,
              sep='\t')
}

### fread a bunch of files and rbind them
fread_rbind <- function(filenames, fn = fread){
  rbindlist(idcol = "filename",
            sapply(USE.NAMES = T,
                   simplify = F,
                   filenames,
                   fn))
}

#' Title
#'
#' @param arg_line
#' @import data.table
#' @export
facets.suite <- function(arg_line = NA){

  ### process args
  if(!is.na(arg_line) | interactive()) {
    # print("reading from arg_line")
    raw_args <- unlist(stringr::str_split(arg_line, " "))
  } else {
    # print("batch")
    raw_args <- commandArgs(TRUE)
  }
  # print(raw_args)

  option_list <- list(
    optparse::make_option(c("-s", "--samplefile"),
                          type="character", default="samples.txt",
                          help="sample file"),
    optparse::make_option(c("-m", "--maf"),
                          type="character", default="samples.maf",
                          help="maf file"),
    optparse::make_option(c("-o", "--outdir"),
                          type="character", default=NA,
                          help="output directory"),
    optparse::make_option(c("-n", "--n_threads"),
                          type="integer", default=1,
                          help="number of threads"),
    optparse::make_option(c('-t', '--targetFile'),
                          type='character', default='IMPACT468',
                          help="IMPACT341/410/468, or a Picard interval list file of gene target coordinates [default IMPACT468]")
  )
  if(any(sapply(
    option_list,
    function(option){
      option@short_flag == "-g"
    }))){
    stop("cannot use short option '-g', conflicts with Rscript --gui")
  }

  opts <- optparse::parse_args(
    optparse::OptionParser(option_list=option_list),
    args = raw_args,
    positional_arguments = TRUE
  )

  options(mc.cores = opts$options$n_threads)
  
  samplefile <- opts$options$samplefile
  outdir <- opts$options$outdir
  targetFile <- opts$options$targetFile

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  maf <- suppressWarnings(fread(opts$options$maf))
  s2c <- suppressWarnings(fread(samplefile))

  original_directory <- getwd()
  setwd(dirname(samplefile))
  if(!all(file.exists(s2c[, cncf]))) {
    warning("Cannot find all paths in sample table\n")
  } else {
    s2c[, cncf := tools::file_path_as_absolute(cncf), cncf]
  }
  setwd(original_directory)

  if (targetFile == "IMPACT341") {
    gene_targets = facets.suite::IMPACT341_targets
    genes = facets.suite::msk_impact_341
  } else if (targetFile == "IMPACT410") {
    gene_targets = facets.suite::IMPACT410_targets
    genes = facets.suite::msk_impact_410
  } else if (targetFile == "IMPACT468") {
    gene_targets = facets.suite::IMPACT468_targets
    genes = facets.suite::msk_impact_468
  } else {
    # Note the target file needs to not only be in the PICARD interval list format
    # But the names must match the regex: /GENESYMBOL_.*/ (e.g. TP53_target_02)
    gene_targets <- suppressWarnings(fread(paste0('grep -v "^@" ', targetFile)))
    setnames(gene_targets, c("chr", "start", "end", "strand", "name"))
    setkey(gene_targets, chr, start, end)
  }


  #### CNCF FILES EXIST ####

  s2c[, cncf_exists := file.exists(cncf)]
  s2c <- s2c[cncf_exists == TRUE]
  write.tab(s2c, file.path(outdir, "samples.txt"))

  #### CNCF ####

  cncf <- generate_cncf_file(s2c)
  setnames(cncf, "tcn", "tcn.cncf")
  setnames(cncf, "lcn", "lcn.cncf")
  setnames(cncf, "cf", "cf.cncf")
  cncf[,`:=`(tcn = tcn.em,
             lcn = lcn.em,
             cf = cf.em)]
  write.tab(cncf, file.path(outdir, "cncf.txt"))


  #### SEG ####

  seg <- generate_logR_seg(cncf)
  write.tab(seg, file.path(outdir, "logR.seg"))
  seg <- generate_logR_adj_seg(cncf)
  write.tab(seg, file.path(outdir, "logR_adj.seg"))
  seg <- generate_log_tcn_seg(cncf)
  write.tab(seg, file.path(outdir, "log_tcn.seg"))


  #### ARM LEVEL ####

  arm_level_calls = get_arm_level_calls(cncf)
  write.tab(arm_level_calls, file.path(outdir, "armLevel.txt"))


  #### GENE LEVEL ####

  gene_level_calls = get_gene_level_calls(cncf, gene_targets = gene_targets)
  write.tab(gene_level_calls, file.path(outdir, "geneLevel.txt"))


  #### MAF ANNO ####

  maf = annotate_maf_with_facets_cf_tcn_lcn_cncf_mc(maf, cncf)

  maf[,t_alt_count := as.numeric(t_alt_count)]
  maf[,t_ref_count := as.numeric(t_ref_count)]
  maf[, paste0("ccf_Mcopies",
               c("", "_lower", "_upper", "_prob95", "_prob90")) :=
        ccf.likelihood(purity,
                       tcn,
                       t_alt_count,
                       (t_alt_count + t_ref_count),
                       copies=(tcn-lcn)), by= 1:nrow(maf)]

  maf[, paste0("ccf_1copy",
               c("", "_lower", "_upper", "_prob95", "_prob90")) :=
        ccf.likelihood(purity,
                       tcn,
                       t_alt_count,
                       (t_alt_count + t_ref_count),
                       copies=1), by= 1:nrow(maf)]
  write.tab(maf, file.path(outdir, "mafAnno.maf"))


  #### QC ####

  QC <- qc(maf, cncf)
  write.tab(QC, file.path(outdir, "QC_summary.txt"))


  #### CNA SUMMARY ####

  cna_summary <- cna_summary(maf, cncf)
  write.tab(cna_summary, file.path(outdir, "cna_summary.txt"))


  #### MUT STATUS ####

  mut_status <- maf_to_mut_status_long(maf = maf)
  gene_mut_cna <- merge(gene_level_calls, mut_status, all.x = T)
  gene_mut_cna[is.na(mutation_category), mutation_category := "wt"]
  write.tab(gene_mut_cna, file.path(outdir, "gene_mut_cna.txt"))

}
