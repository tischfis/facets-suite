arm_gene <- function(gene_targets){
  fo_impact <- foverlaps(
    gene_targets,
    facets.suite::arm_definitions)
  fo_impact[,Hugo_Symbol:=gsub("_.*$", "", name)]
  arm_gene <- fo_impact[, .(arm = unique(arm)), Hugo_Symbol]
  arm_gene
}
