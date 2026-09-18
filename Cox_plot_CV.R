res_comp_W <- res_comp

##

res_cox <- cv100_res_full$raw_results
res_Weibull <- cv100_res_full_W$raw_results


res_cox$iBS <- res_cox$Brier
res_Weibull$iBS <- res_Weibull$Brier

summary(res_comp$iBS)

library(ggplot2)
library(tidyr)
library(dplyr)

#############################
# Transform iBS to iAUC scale
#############################


# iAUC scale
auc_min <- 0.65
auc_max <- 0.85

# iBS scale
bs_min <- 0.14
bs_max <- 0.45


transform_bs <- function(x){
  auc_min + (x - bs_min) / (bs_max - bs_min) * (auc_max - auc_min)
}


res_comp_bs <- transform_bs(res_comp$iBS)
res_comp_W_bs <- transform_bs(res_comp_W$iBS)
res_cox_bs      <- transform_bs(res_cox$iBS)
res_W_bs      <- transform_bs(res_Weibull$iBS)



#############################
# Boxplots
#############################
jpeg("CV_boxplot.jpeg", quality = 100, units = "in",  width = 14, height = 7, res = 300)
par(
  mar = c(3, 4, 3, 5)
)

boxplot(
  
  # iAUC
  res_comp$iAUC,
  res_comp_W$iAUC,
  res_cox$iAUC,
  res_Weibull$iAUC,
  
  # iBS transformed
  res_comp_bs,
  res_comp_W_bs,
  res_cox_bs,
  res_W_bs,
  
  at = c(
    1,2,3,4,
    6,7,8,9
  ),
  
  col = c(
    "red2",
          "orange2",
          "blue2",
          "skyblue",
          
          "red2",
          "orange2",
          "blue2",
          "skyblue"
  ),
  
  border = c(
    "red4",
          "darkorange4",
          "blue4",
          "skyblue4",
          
          "red4",
          "darkorange4",
          "blue4",
          "skyblue4"
  ),
  
  xaxt = "n",
  yaxt = "n",
  
  ylim = c(0.65, 0.85),
  
  ylab = "Validation performance",
  cex.lab = 1.5,
  outline = TRUE,
  lwd = 2
)



#############################
# Left y-axis: iAUC
#############################

axis(
  2,
  at = seq(0.65,0.85,0.05),
  labels = seq(0.65,0.85,0.05),
  las = 1
)



#############################
# Right y-axis: iBS
#############################

axis(
  4,
  at = seq(0.65, 0.85, length.out = 5),
  labels = seq(0.15,0.45,length.out = 5),
  las = 1
)



#############################
# x-axis metric labels
#############################

axis(
  1,
  at = c(2.5,7.5),
  labels = c("iAUC","iBS"),
  tick = FALSE,
  line = 1,
  cex.axis = 1.5
)



#############################
# separating line
#############################

abline(
  v = 5,
  lty = 2,
  lwd = 2
)



#############################
# Legend
#############################


#legend_labels <- c(
#  expression(paste("Competing Cox (", italic("\u03c1")[1], "=", italic("\u03c1")[2], "=log n)")),
#  expression(paste("Competing Cox (", italic("\u03c1")[1], "=2, ", italic("\u03c1")[2], "=log n)")),
#  "Cox (BIC)",
#  "Cox (AIC)"
#)

legend_labels <- c("Cox-type Polyhazard",
  "poly-Weibull",
  "Cox",
  "Weibull"
)


legend(
  "topright",
  legend = legend_labels,
  fill = c(
    "red2",
          "orange2",
          "blue2",
          "skyblue"
  ),
  border = c(
    "red4",
          "darkorange4",
          "blue4",
          "skyblue4"
  ),
  bty = "o",
  cex = 1.5
)



 dev.off()
 
 
 
 
 
 
 
 
 
 
 
 
 
 