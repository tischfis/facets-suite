# maf_to_mut_status <- function(maf, genes = c("KRAS", "TP53", "SMAD4", "CDKN2A", "SF3B1", "U2AF1")){
#   maf[, Tumor_Sample_Barcode := factor(Tumor_Sample_Barcode)]
#   maf <- maf[Consequence %in% facets.suite::Nonsyn_Consequences & Hugo_Symbol %in% genes]
#   maf[, Hugo_Symbol := factor(Hugo_Symbol, levels = genes)]
#
#   maf[, mutation_category := ifelse(!is.na(stringr::str_match(
#     Consequence,
#     paste(collapse="|",
#           facets.suite::Trunc_Consequences))), "truncating",
#     ifelse(tm %in% facets.suite::hotspots_24k, "hotspot",
#            ifelse(!is.na(stringr::str_match(
#              Consequence,
#              paste(collapse="|",
#                    facets.suite::Inframe_Consequences))), "inframe",
#              ifelse(Consequence %like% "missense_variant", "missense",
#                     "other")
#            )))
#     ]
#
#   dc <- dcast.data.table(
#     maf[Consequence %in% facets.suite::Nonsyn_Consequences &
#           Hugo_Symbol %in% genes],
#     Tumor_Sample_Barcode ~ Hugo_Symbol,
#     value.var = "mutation_category",
#     fun.aggregate = function(x){
#       ifelse("truncating" %in% x, "truncating",
#              ifelse("hotspot" %in% x, "hotspot",
#                     ifelse("inframe" %in% x, "inframe",
#                            ifelse("missense" %in% x, "missense",
#                                   "other"))))},
#     fill = "wt", drop = F)
#   for(g in genes){
#     dc[, (g) := factor(get(g), levels = c("truncating", "hotspot", "inframe", "missense", "other", "wt"))]
#   }
#   setorderv(dc, genes)
#   copy(dc)
# }

choose_mutation_category_mutation <- function(Hugo_Symbol, HGVSp_Short, Consequence){
  tm <- paste(
    Hugo_Symbol,
    stringr::str_extract(
      HGVSp_Short,
      "(?<=^p.[A-Z])[0-9]+|(?<=^p.)[A-Z][0-9]+_splice"))
  factor(
    levels = c("truncating", "hotspot", "inframe", "missense", "other", "wt"),
    ifelse(!is.na(stringr::str_match(
      Consequence,
      paste(collapse="|",
            facets.suite::Trunc_Consequences))), "truncating",
      ifelse(tm %in% facets.suite::hotspots_24k, "hotspot",
             ifelse(!is.na(stringr::str_match(
               Consequence,
               paste(collapse="|",
                     facets.suite::Inframe_Consequences))), "inframe",
               ifelse(Consequence %like% "missense_variant", "missense",
                      "other")
             )))
  )
}

choose_mutation_category_sample <- function(mutation_category){
  sort(mutation_category)[1]
}

maf_to_mut_status_long <- function(maf, genes = NA){
  if(is.character(genes)){
    maf <- maf[Hugo_Symbol %in% genes]
  }
  maf[, mutation_category := choose_mutation_category_mutation(Hugo_Symbol, HGVSp_Short, Consequence)]
  mut_status <- maf[, .(mutation_category = choose_mutation_category_sample(mutation_category)),
                    keyby = .(Tumor_Sample_Barcode, Hugo_Symbol)]

  copy(mut_status)
}
