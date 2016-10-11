library(tidyr)
library(magrittr)
library(dplyr)
library(pracma)  ## odregress
signedln = function(x) {
  ifelse(x <= 1, 0, log(x))
}
logmean <- function (x) {
  mean <- 10000*mean(x, na.rm = T)
  signedln(mean)
}
dNSAF_input <- read.delim("201610_MS_results.txt", head = T)
outfilename_plot_all = paste0(Sys.Date(), "_MS_dNSAF_all_pairwise_plots.pdf")
dNSAF_gathered <- dNSAF_input %>% gather(sample, dNSAF, -Protein_Name)
dNSAF <- dNSAF_gathered %>% separate(sample, into = c("flag", "bait", "mutant", "technical_replicate")) 
dNSAF_mean <- with(dNSAF, aggregate(dNSAF, by = list(Protein_Name, flag, bait, mutant), mean))
names(dNSAF_mean) <- c("protein", "flag", "bait", "mutant", "dNSAF_mean")
dNSAF_logmean <- with(dNSAF, aggregate(dNSAF, by = list(Protein_Name, flag, bait, mutant), logmean))
dNSAF_mean <- cbind(dNSAF_mean, "dNSAF_logmean" = dNSAF_logmean$x)
dNSAF_mean <- dNSAF_mean %>% gather(measure, dNSAF, -protein, -flag, -bait, -mutant)
mutants <- as.character(unique(dNSAF_mean$mutant))
mutants <- mutants[mutants != "WT"]
Nterminal_WT <- dNSAF_mean %>% subset(flag == "NFLAG" & mutant == "WT")
Cterminal_WT <- dNSAF_mean %>% subset(flag == "CFLAG" & mutant == "WT")
residuals_df <- data.frame()
orthogonal_distance_regression_ms <- function(df, df_ref, measure_to_plot) {
  merged_df <- merge(df, df_ref, by = "protein") %>% filter(measure.x == measure_to_plot, measure.y == measure_to_plot)
  odregress_fit <- odregress(merged_df$dNSAF.x, merged_df$dNSAF.y)
  fit_parameters <- odregress_fit$coeff
  merged_df <- cbind(merged_df, "residuals" = odregress_fit$resid)
  if (merged_df$flag.x[1] == merged_df$flag.y[1] & measure_to_plot == "dNSAF_logmean") {
     residuals_df <- rbind(residuals_df, merged_df)
  }
  return(merged_df)
}
xy.limits <- data.frame("measure" = rep("dNSAF_mean",2), "range" = range(dNSAF_mean$dNSAF[dNSAF_mean$measure == "dNSAF_mean"]))
xy.limits <- rbind(xy.limits, data.frame("measure" = rep("dNSAF_logmean",2), "range" = range(dNSAF_mean$dNSAF[dNSAF_mean$measure == "dNSAF_logmean"])))
odregress_plot <- function (merged_df) {
  measure_being_plotted <- as.character(merged_df$measure.x[1])
  sample <- paste(as.character(merged_df$flag.x[1]), as.character(merged_df$mutant.x[1]), sep = "_")
  control <- paste(as.character(merged_df$flag.y[1]), as.character(merged_df$mutant.y[1]), sep = "_")
  xlabel <- paste(sample, measure_being_plotted, sep = "_")
  ylabel <- paste(control, measure_being_plotted, sep = "_")
  plot_title <- paste(sample, control, "orthogonal_distance_regression", measure_being_plotted, sep = "_")
  axes_limits <- xy.limits$range[xy.limits$measure == measure_being_plotted]
  plot(merged_df$dNSAF.x, merged_df$dNSAF.y, xlim = axes_limits, ylim = axes_limits, cex.main = 0.7, main = plot_title, xlab = xlabel, ylab = ylabel, col = "gray", pch = 19)
  abline(b = fit_parameters[1], a = fit_parameters[2], col = "gray")
  interesting_subset <- subset(merged_df, abs(residuals) > quantile(abs(merged_df$residuals), probs = 0.75))
  points(interesting_subset$dNSAF.x, interesting_subset$dNSAF.y, col = "darkseagreen", pch = 19)
  text(interesting_subset$dNSAF.x, interesting_subset$dNSAF.y, labels = interesting_subset$protein, pos = 3, cex = 0.4)
}
#### plot mutants versus wt in normal and log scale
opar <- par
pdf(outfilename_plot_all, width = 10, height = 10)
op <- par(mfrow = c(2,2))
for (m in mutants) {
  mutant_subset <- subset(dNSAF_mean, mutant == m)
  if (mutant_subset$flag[1] == "NFLAG") {
    orthogonal_distance_regression_ms(mutant_subset, Nterminal_WT, "dNSAF_mean") %>% odregress_plot
    orthogonal_distance_regression_ms(mutant_subset, Cterminal_WT, "dNSAF_mean") %>% odregress_plot
    orthogonal_distance_regression_ms(mutant_subset, Nterminal_WT, "dNSAF_logmean") %>% odregress_plot
    orthogonal_distance_regression_ms(mutant_subset, Cterminal_WT, "dNSAF_logmean") %>% odregress_plot
  } else {
    orthogonal_distance_regression_ms(mutant_subset, Cterminal_WT, "dNSAF_mean") %>% odregress_plot
    orthogonal_distance_regression_ms(mutant_subset, Nterminal_WT, "dNSAF_mean") %>% odregress_plot
    orthogonal_distance_regression_ms(mutant_subset, Cterminal_WT, "dNSAF_logmean") %>% odregress_plot
    orthogonal_distance_regression_ms(mutant_subset, Nterminal_WT, "dNSAF_logmean") %>% odregress_plot
  }
}
dev.off()
par(opar)






