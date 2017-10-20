estimated_af_life_history = function(purity, ns, nw, m, M, copies=1, limit=TRUE){

  if(any(is.na(c(purity, ns, nw, m, M)))){return(NA)}

  if(!copies %in% c(1, "M")) stop('copies must be 1 or "M"')
  if(copies == "M"){r = M } ## number of copies present
  else {r = 1}

  allele_fraction = ns / (ns + nw)

  if(limit){frac = min(1, allele_fraction * (purity * (M+m) + 2*(1-purity)) / purity / r)}

  else{frac = allele_fraction * (purity * (M+m) + 2*(1-purity)) / purity / r}
  as.numeric(frac)
}


################################################################################################################################
################################################################################################################################

integer_cn_table = function(out, fit, em=FALSE){

  df = out$IGV
  n.xchr <- nrow(df[df$chrom == 23,])
  if(n.xchr > 0) {
    df[df$chrom == 23,]$chrom = "X"
  }
  df$chrom = factor(df$chrom)
  if(em==TRUE){
    dt = data.table(df,
                    cf=fit$cncf$cf.em,
                    tcn=fit$cncf$tcn.em,
                    mcn=fit$cncf$tcn.em - fit$cncf$lcn.em,
                    lcn=fit$cncf$lcn.em,
                    ploidy=fit$ploidy,
                    purity=fit$purity,
                    dipLogR=out$dipLogR)
  }
  if(em==FALSE){
    dt = data.table(df,
                    cf=fit$cncf$cf,
                    tcn=fit$cncf$tcn,
                    mcn=fit$cncf$tcn - fit$cncf$lcn,
                    lcn=fit$cncf$lcn,
                    ploidy=fit$ploidy,
                    purity=fit$purity,
                    dipLogR=out$dipLogR)
  }
  setkey(dt, chrom, loc.start, loc.end)
  dt
}

################################################################################################################################
################################################################################################################################

annotate_maf_with_facets_cf_tcn_lcn_cncf = function(maf, cncf){

  maf = as.data.table(maf)
  maf <- maf[Tumor_Sample_Barcode %in% unique(cncf$Tumor_Sample_Barcode)]
  maf_cols = colnames(maf)
  maf[, Chromosome := factor(as.character(Chromosome))]
  setkey(maf,
         Tumor_Sample_Barcode,
         Chromosome,
         Start_Position,
         End_Position)
  cncf[, chrom := factor(as.character(chrom))]
  cncf[chrom == 23, chrom := "X"]

  setkey(cncf,
         Tumor_Sample_Barcode,
         chrom,
         loc.start,
         loc.end)

  #dt = integer_cn_table(out, fit, em = F)

  ### check for duplicate columns
  if(any(duplicated(names(maf)))){
    warning("duplicate columns removed from maf file")
    maf[, which(duplicated(names(maf))) := NULL, with = F]
  }

  # if(is.null(iTumor_Sample_Barcode)){
  #   maf_ann = foverlaps(maf, dt, mult="first",nomatch=NA)
  # }else{
  maf_ann = foverlaps(maf, cncf,
                      by.x = c("Tumor_Sample_Barcode",
                               "Chromosome",
                               "Start_Position",
                               "End_Position"),
                      mult="first",
                      nomatch=NA)
  # }

  setcolorder(maf_ann, c(maf_cols, setdiff(names(maf_ann), maf_cols)))
  copy(maf_ann)
  #maf_ann[,c(maf_cols, 'dipLogR', 'seg.mean', 'cf', 'tcn', 'lcn', 'purity', 'ploidy'), with=F]
}


################################################################################################################################
################################################################################################################################

annotate_maf_with_facets_cf_tcn_lcn_cncf_mc = function(maf, cncf){
  
  maf = as.data.table(maf)
  maf_cols = colnames(maf)

  cncf_remove_names <- c("ID", "seg", "num.mark", "nhet", "mafR", "segclust", 
    "cnlr.median.clust", "mafR.clust", "start", "end", "cf.em", "tcn.em", 
    "lcn.em", "cf.cncf", "tcn.cncf", "lcn.cncf", "dipt", "loglik")
  cncf <- cncf[, setdiff(names(cncf), cncf_remove_names), with = F]
  
    ### check for duplicate columns
  if(any(duplicated(names(maf)))){
    warning("duplicate columns removed from maf file")
    maf[, which(duplicated(names(maf))) := NULL, with = F]
  }
  
  ids <- sort(intersect(maf$Tumor_Sample_Barcode, cncf$Tumor_Sample_Barcode))
  maf <- maf[Tumor_Sample_Barcode %in% ids]
  cncf <- cncf[Tumor_Sample_Barcode %in% ids]
  
  maf_list <- lapply(
    ids, function(id) {
      sample_maf <- maf[Tumor_Sample_Barcode == id]
      sample_maf[, Chromosome := factor(as.character(Chromosome))]
      setkey(sample_maf,
             Chromosome,
             Start_Position,
             End_Position)
      sample_maf
    })
  cncf_list <- lapply(
    ids, 
    function(id){
      sample_cncf <- cncf[Tumor_Sample_Barcode == id]
      sample_cncf[, chrom := factor(as.character(chrom))]
      sample_cncf[chrom == 23, chrom := "X"]
      sample_cncf[, Tumor_Sample_Barcode := NULL]
      setkey(sample_cncf,
             chrom,
             loc.start,
             loc.end)
      sample_cncf
    })
  
  maf_ann = rbindlist(parallel::mcmapply(
    foverlaps,
    x = maf_list,
    y = cncf_list,
    MoreArgs = list(mult = "first", nomatch = NA),
    SIMPLIFY = F
  ))

  setcolorder(maf_ann, c(maf_cols, setdiff(names(maf_ann), maf_cols)))
  copy(maf_ann)
  #maf_ann[,c(maf_cols, 'dipLogR', 'seg.mean', 'cf', 'tcn', 'lcn', 'purity', 'ploidy'), with=F]
}


################################################################################################################################
################################################################################################################################

ccf.likelihood = function(purity, absCN, alt_allele, coverage, copies){

  #From McGranahan_and_Swanton_2015

  CCFs = seq(0.001,1,0.001)
  vac.ccf  = function(CCF, purity, absCN){purity * CCF * copies / (2*(1 - purity) + purity * absCN)}
  probs = sapply(CCFs, function(c){dbinom(alt_allele, coverage, vac.ccf(c, purity, absCN))})
  probs = probs/sum(probs)

  ccf.max = which.max(probs)
  ccf.gt.half.max = which(probs > max(probs)/2)
  ccf.lower = max(ccf.gt.half.max[1] - 1, 1) ### closest ccf value before half-max range (within 0-1 range)
  ccf.upper = min(ccf.gt.half.max[length(ccf.gt.half.max)] + 1, length(CCFs)) ### closest ccf value after half-max range (within 0-1 range)
  if(is.na(purity)){ccf.upper=NA}
  ccf.max = ccf.max/length(CCFs)
  ccf.lower = ccf.lower/length(CCFs)
  ccf.upper = ccf.upper/length(CCFs)
  prob.95 = sum(probs[950:1000])
  prob.90 = sum(probs[900:1000])
  #if(is.na(purity)){ccf.upper=NA}
  list(ccf.max,ccf.lower,ccf.upper,prob.95,prob.90)
}


################################################################################################################################
################################################################################################################################


ccf.likelihood.purity.fit = function(purity, absCN, alt_allele, coverage, copies){

  #From McGranahan_and_Swanton_2015

  CCFs = seq(0.001,1,0.001)
  vac.ccf  = function(CCF, purity, absCN){purity * CCF * copies / (2*(1 - purity) + purity * absCN)}
  probs = sapply(CCFs, function(c){dbinom(alt_allele, coverage, vac.ccf(c, purity, absCN))})
  probs = probs/sum(probs)

  ccf.max = which.max(probs)
  ccf.gt.half.max = which(probs > max(probs)/2)
  ccf.lower = max(ccf.gt.half.max[1] - 1, 1) ### closest ccf value before half-max range (within 0-1 range)
  ccf.upper = min(ccf.gt.half.max[length(ccf.gt.half.max)] + 1, length(CCFs)) ### closest ccf value after half-max range (within 0-1 range)
  if(is.na(purity)){ccf.upper=NA}
  ccf.max = ccf.max/length(CCFs)
  ccf.lower = ccf.lower/length(CCFs)
  ccf.upper = ccf.upper/length(CCFs)
  prob.95 = sum(probs[950:1000])
  prob.90 = sum(probs[900:1000])
  #if(is.na(purity)){ccf.upper=NA}
  list(ccf.max,ccf.lower,ccf.upper,prob.95,prob.90)
}


################################################################################################################################
################################################################################################################################


# mafAnno = function(maf, cncf){
#
#   maf = as.data.table(maf)
#   maf = annotate_maf_with_facets_cf_tcn_lcn_cncf(maf, cncf_out)
#
#   maf[,t_alt_count := as.numeric(t_alt_count)]
#   maf[,t_ref_count := as.numeric(t_ref_count)]
#   maf[, paste0("ccf_Mcopies",
#                c("", "_lower", "_upper", "_prob95", "_prob90")) :=
#         ccf.likelihood(purity,
#                        tcn,
#                        t_alt_count,
#                        (t_alt_count + t_ref_count),
#                        copies=(tcn-lcn)), by= 1:nrow(maf)]
#
#   maf[, paste0("ccf_1copy",
#                c("", "_lower", "_upper", "_prob95", "_prob90")) :=
#         ccf.likelihood(purity,
#                        tcn,
#                        t_alt_count,
#                        (t_alt_count + t_ref_count),
#                        copies=1), by= 1:nrow(maf)]
#   maf
# }
