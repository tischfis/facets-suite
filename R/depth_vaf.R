#' Plot read depth against VAF
#'
#' @param maf 
#'
#' @return
#' @export
#' @import ggplot2
depth_vaf_plot <- function(
  maf, 
  min_alt_count = 10){
  
  purities <- maf[, .(purity = unique(purity)), Tumor_Sample_Barcode]
  
  maf[, t_depth := t_alt_count + t_ref_count]
  maf[, paste0("t_va r_freq", c("", "_lower", "_upper")) := as.list(
    binom::binom.confint(t_alt_count, t_depth, method = "exact",
                         conf.level = pnorm(1) - pnorm(-1))[4:6]),
    list(t_alt_count, t_depth)]
  
  ggplot(maf) +#[Consequence %like% paste(collapse = "|", autospy::Nonsyn_Consequences)]) + 
    stat_function(
      fun = function(t_depth){min_alt_count / t_depth}, 
      colour = "grey") +
    # coord_cartesian(ylim = c(0, 1)) +
    scale_colour_manual("nsyn", values = c("grey", "grey40")) +
    theme(legend.position = "bottom") +
    scale_x_log10() +
    annotation_logticks(sides = "b") +
    labs(x = "Read Depth", y = "VAF") +
    geom_errorbar(
      size = 0.2,
      alpha = 0.2,
      width = 0,
      aes(
        #col = Consequence %like% paste(collapse = "|", autospy::Nonsyn_Consequences),
        x = t_depth, 
        ymin = t_var_freq_lower,
        ymax = t_var_freq_upper
      )) +
    geom_point(
      size = 0.5,
      aes(
        #col = Consequence %like% paste(collapse = "|", autospy::Nonsyn_Consequences),
        x = t_depth, 
        y = t_var_freq
      )) + 
    geom_hline(data = purities, aes(yintercept = purity)) +
    facet_wrap(~Tumor_Sample_Barcode)
}

plot_depth_vaf_tile <- function(
  maf, 
  min_alt_count = 10){
  
  ggplot() +#[Consequence %like% paste(collapse = "|", autospy::Nonsyn_Consequences)]) + 
    theme(legend.position = "bottom") +
    scale_x_log10() +
    annotation_logticks(sides = "b") +
    labs(x = "Read Depth", y = "VAF") +
    # coord_cartesian(ylim = c(0, 1)) +
    # scale_colour_manual("nsyn", values = c("grey", "grey40")) +
    scale_fill_gradient(low = "white", high = "blue") +
    geom_bin2d(data = maf[datasim == "sim"], 
               aes(x = t_depth, 
                   y = t_var_freq
               )) +
    stat_function(data = maf,
                  fun = function(t_depth){min_alt_count / t_depth}, 
                  colour = "black") +
    geom_point(data = maf[datasim == "data"], 
      size = 0.5, 
      aes(
        x = t_depth, 
        y = t_var_freq
      ))
  
}