mut_summary <- function(maf){
  maf[
    Consequence %like% "missense_variant",
    tm := paste(
      Hugo_Symbol,
      stringr::str_extract(
        HGVSp_Short,
        "(?<=^p.[A-Z])[0-9]+|(?<=^p.)[A-Z][0-9]+_splice"
      )
    )
    ]
  maf[, hotspot_24k := tm %in% facets.suite::hotspots_24k]
  
  maf[, mutation_category :=
            ifelse(Variant_Type %in% c("INS", "DEL") &
                     Consequence %like% paste(collapse="|", Nonsyn_Consequences), "indel (non-syn)",
                   #ifelse(Variant_Type %in% c("INS", "DEL"), "indel (non-coding)",
                   ifelse(!is.na(stringr::str_match(Consequence, paste(collapse="|", Trunc_Consequences))), "truncating",
                          ifelse(tm %in% hotspots_24k, "hotspot",
                                 ifelse(Consequence %like% "missense_variant", "missense",
                                        #ifelse(Consequence %like% "synonymous_variant", "synonymous",
                                        #ifelse(dbSNP_RS != "", "dbSNP",
                                        #       "noncoding")
                                        "silent"
                                 ))))#))
          ]
  mutation_categories <-
    c(
      "indel (non-syn)",
      # "indel (non-coding)",
      "truncating",
      "hotspot",
      "missense",
      "silent"
    )
  maf[, mutation_category := factor(mutation_category, levels = mutation_categories)]
  # maf <- maf[!(mutation_category == "indel (non-coding)" & t_depth < 500)]
  mut_N <- maf[, .(Nmut = .N), keyby = Tumor_Sample_Barcode]
  mut_summ <- dcast.data.table(
    maf,
    Tumor_Sample_Barcode ~ mutation_category,
    fun.aggregate = length,
    value.var = "HGVSp_Short"
  )
  setkey(mut_summ, Tumor_Sample_Barcode)
  mut_summ <- merge(mut_N, mut_summ)
  mut_summ
}