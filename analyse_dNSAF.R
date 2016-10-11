library(tidyr)
library(magrittr)
library(dplyr)
library(pracma)  ## odregress

options(stringsAsFactors = F)

#### for taking modified log mean of dNSAF (sets log to 0 if the values are small)
signedln = function(x) {
  ifelse(x <= 1, 0, log(x))
}
logmean <- function (x) {
  mean <- 10000*mean(x, na.rm = T)
  signedln(mean)
}
#############

### load e-map correlation data
EMAP_correlations <- read.delim("correlation_network2.0.txt", head = F)
EMAP_correlations <- EMAP_correlations[!duplicated(EMAP_correlations),]
names(EMAP_correlations) <- as.character(expression(protein1, protein2, interaction, corr_weight))
################

#### load mass spec processed data for mutants and FLAG WT (dNSAF data)
dNSAF_input <- read.delim("201610_MS_results.txt", head = T)
outfilename_plot_all = paste0(Sys.Date(), "_MS_dNSAF_all_pairwise_plots.pdf")
dNSAF_gathered <- dNSAF_input %>% gather(sample, dNSAF, -Protein_Name)

##### get mean and logmean for dNSAF (combining technical replicates)
dNSAF <- dNSAF_gathered %>% separate(sample, into = c("flag", "bait", "mutant", "technical_replicate")) 
dNSAF_mean <- with(dNSAF, aggregate(dNSAF, by = list(Protein_Name, flag, bait, mutant), mean))
names(dNSAF_mean) <- as.character(expression(protein, flag, bait, mutant, dNSAF_mean))
dNSAF_logmean <- with(dNSAF, aggregate(dNSAF, by = list(Protein_Name, flag, bait, mutant), logmean))
dNSAF_mean <- cbind(dNSAF_mean, "dNSAF_logmean" = dNSAF_logmean$x)
dNSAF_mean <- dNSAF_mean %>% gather(measure, dNSAF, -protein, -flag, -bait, -mutant)

####### limiting factor is the mass spec data - so mutants are the ones for which there is mass spec data
mutants <- as.character(unique(dNSAF_mean$mutant))
mutants <- mutants[mutants != "WT"]

#### combine E-MAP correlations and mass spec data
GSP1_partners_from_ms <- as.character(unique(dNSAF_mean$protein))
EMAP_correlation_partners_from_ms <- EMAP_correlations %>% 
  filter(protein1 %in% GSP1_partners_from_ms | protein2 %in% GSP1_partners_from_ms) %>% 
  filter(protein1 %in% mutants | protein2 %in% mutants)
EMAP_correlation_partners_from_ms <- cbind(EMAP_correlation_partners_from_ms, "scaled_corr" = scale(EMAP_correlation_partners_from_ms$corr_weight))

#### fit the data for each mutant against the FLAG - WT and WT - FLAG controls
###     -> to identify interactions that are affeced by that specific mutation
Nterminal_WT <- dNSAF_mean %>% subset(flag == "NFLAG" & mutant == "WT")
Cterminal_WT <- dNSAF_mean %>% subset(flag == "CFLAG" & mutant == "WT")
orthogonal_distance_regression_ms <- function(df, df_ref, measure_to_plot) {
  merged_df <- merge(df, df_ref, by = "protein") %>% filter(measure.x == measure_to_plot, measure.y == measure_to_plot)
  odregress_fit <- odregress(merged_df$dNSAF.x, merged_df$dNSAF.y)
  fit_parameters <- odregress_fit$coeff
  merged_df <- cbind(merged_df, data.frame("residuals" = odregress_fit$resid, "b" = fit_parameters[1], "a" = fit_parameters[2]))
  return(merged_df)
}
ms_residuals_df <- data.frame()
for (m in mutants) {
  mutant_subset <- dNSAF_mean %>% filter(mutant == m)
  wt_subset <- dNSAF_mean %>% filter(mutant == "WT", flag == mutant_subset$flag[1])
  temp_residuals_df <- orthogonal_distance_regression_ms(mutant_subset, wt_subset, "dNSAF_logmean")
  temp_residuals_df <- temp_residuals_df[c(1, 2, 4, 6, 11, 12)]
  ms_residuals_df <- rbind(ms_residuals_df, temp_residuals_df)
}
ms_residuals_df <- cbind(ms_residuals_df, "scaled_ms_fit_residuals" = scale(ms_residuals_df$residuals))
names(ms_residuals_df) <- as.character(expression(protein, flag, mutant, dNSAF_mutant, dNSAF_wt, ms_fit_residual, scaled_ms_fit_residual))

### merge ms residuals and E-MAP corr data
ms_emap_merged_1 <- merge(residuals_df, EMAP_correlation_partners_from_ms, by.x = c("mutant", "protein"), by.y = c("protein1", "protein2"))
ms_emap_merged_2 <- merge(residuals_df, EMAP_correlation_partners_from_ms, by.x = c("mutant", "protein"), by.y = c("protein2", "protein1"))
ms_emap_merged <- rbind(ms_emap_merged_1, ms_emap_merged_2)
pdf(file = "ms_residuals_versus_emap_correlations.pdf")
for (m in mutants) {
  ms_emap_subset <- ms_emap_merged %>% filter(mutant == m)
  xlabel = "dNSAF logmean odregress residual"
  ylabel = "E-MAP correlation"
  plot(ms_emap_subset$ms_fit_residual, ms_emap_subset$corr_weight, main = m, xlab = xlabel, ylab = ylabel, pch = 19)
  text(ms_emap_subset$ms_fit_residual, ms_emap_subset$corr_weight, labels = ms_emap_subset$protein, pos = 4, cex = 0.5)
  legend("topright", legend = c("depleated PPI", "high E-MAP correlation"), cex = 0.5)
  legend("topleft", legend = c("increased PPI", "high E-MAP correlation"), cex = 0.5)
}
dev.off()
##### plot and output all the odregress fits
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
  abline(b = merged_df$b[1], a = merged_df$a[1], col = "gray")
  interesting_subset <- subset(merged_df, abs(residuals) > quantile(abs(merged_df$residuals), probs = 0.75))
  points(interesting_subset$dNSAF.x, interesting_subset$dNSAF.y, col = "darkseagreen", pch = 19)
  text(interesting_subset$dNSAF.x, interesting_subset$dNSAF.y, labels = interesting_subset$protein, pos = 3, cex = 0.4)
}
#### plot mutants versus wt in normal and log scale
opar <- par
pdf(outfilename_plot_all, width = 10, height = 10)
op <- par(mfrow = c(2,2))
for (m in mutants) {
  mutant_subset <- dNSAF_mean %>% filter(mutant == m)
  wt_subset <- dNSAF_mean %>% filter(mutant == "WT", flag == mutant_subset$flag[1])
  other_wt_subset <- dNSAF_mean %>% filter(mutant == "WT", flag != mutant_subset$flag[1])
  orthogonal_distance_regression_ms(mutant_subset, wt_subset, "dNSAF_mean") %>% odregress_plot
  orthogonal_distance_regression_ms(mutant_subset, other_wt_subset, "dNSAF_mean") %>% odregress_plot
  orthogonal_distance_regression_ms(mutant_subset, wt_subset, "dNSAF_logmean") %>% odregress_plot
  orthogonal_distance_regression_ms(mutant_subset, other_wt_subset, "dNSAF_logmean") %>% odregress_plot
}
dev.off()
par(opar)






