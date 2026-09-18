#load("LASSO_L4_50.Rdata")


summary(final_res)
boxplot(final_res$Cox_train_iAUC,final_res$Comp_train_iAUC,final_res$Cox_test_iAUC, final_res$Comp_test_iAUC)
boxplot(final_res$Cox_train_iBS,final_res$Comp_train_iBS,final_res$Cox_test_iBS, final_res$Comp_test_iBS)






##

res_AIC <- cv100_res_AIC$raw_results
res_BIC <- cv100_res_BIC$raw_results


summary(c(res_AIC$iBS,res_BIC$iBS))

res_AIC$iBS <- res_AIC$Brier
res_BIC$iBS <- res_BIC$Brier

library(ggplot2)
library(tidyr)
library(dplyr)

#############################
# Transform iBS to iAUC scale
#############################

# iAUC scale
auc_min <- 0.7
auc_max <- 0.9

# iBS scale
bs_min <- 0.1
bs_max <- 0.2

transform_bs <- function(x){
  auc_min + (x - bs_min) / (bs_max - bs_min) * (auc_max - auc_min)
}


res_comp_BIC_bs <- transform_bs(res_comp_BIC$iBS)
res_comp_AIC_bs <- transform_bs(res_comp_AIC$iBS)
res_BIC_bs      <- transform_bs(res_BIC$iBS)
res_AIC_bs      <- transform_bs(res_AIC$iBS)



#############################
# Boxplots
#############################
jpeg("CV_boxplot.jpeg", quality = 100, units = "in",  width = 12, height = 6, res = 300)
par(
  mar = c(3, 4, 3, 5)
)

boxplot(
  
  # iAUC
  res_comp_BIC$iAUC,
  res_comp_AIC$iAUC,
  res_BIC$iAUC,
  res_AIC$iAUC,
  
  # iBS transformed
  res_comp_BIC_bs,
  res_comp_AIC_bs,
  res_BIC_bs,
  res_AIC_bs,
  
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
  
  ylim = c(auc_min, auc_max),
  
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
  at = seq(0.7,0.9,0.05),
  labels = seq(0.7,0.9,0.05),
  las = 1
)



#############################
# Right y-axis: iBS
#############################

axis(
  4,
  at = seq(0.7, 0.9, length.out = 5),
  labels = seq(0.1,0.2,length.out = 5),
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

legend_labels <- c(
  expression(paste("Cox-type Polyhazard (", italic("\u03c1")[1], "=", italic("\u03c1")[2], "= ", log(n), ")")),
  expression(paste("Cox-type Polyhazard (", italic("\u03c1")[1], "=2, ", italic("\u03c1")[2], "= ", log(n), ")")),
  "Cox (BIC)",
  "Cox (AIC)"
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
  cex = 1.2
)



 dev.off()
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 final_res4 <- final_res
 
 
 
 
 
 ## Train-test val
  jpeg("LUAD_L_train_test_comparison.jpeg",      quality = 100,      units = "in",      width = 18,
      height = 8,     res = 300)
 
 par(mar=c(3,4,3,2))
 
 
 boxplot(
   
   #############################
   # L = 1
   #############################
   
   final_res1$Cox_train_iAUC,
   final_res1$Comp_train_iAUC,
   final_res1$Cox_test_iAUC,
   final_res1$Comp_test_iAUC,
   
   final_res1$Cox_train_iBS,
   final_res1$Comp_train_iBS,
   final_res1$Cox_test_iBS,
   final_res1$Comp_test_iBS,
   
   
   #############################
   # L = 3
   #############################
   
   final_res3$Cox_train_iAUC,
   final_res3$Comp_train_iAUC,
   final_res3$Cox_test_iAUC,
   final_res3$Comp_test_iAUC,
   
   final_res3$Cox_train_iBS,
   final_res3$Comp_train_iBS,
   final_res3$Cox_test_iBS,
   final_res3$Comp_test_iBS,
   
   
   #############################
   # L = 4
   #############################
   
   final_res4$Cox_train_iAUC,
   final_res4$Comp_train_iAUC,
   final_res4$Cox_test_iAUC,
   final_res4$Comp_test_iAUC,
   
   final_res4$Cox_train_iBS,
   final_res4$Comp_train_iBS,
   final_res4$Cox_test_iBS,
   final_res4$Comp_test_iBS,
   
   
   #############################
   # positions
   #############################
   
   at=c(
     1:8,
     10:17,
     19:26
   ),
   
   
   #############################
   # colour specification
   #############################
   
   col=rep(
     c(
       "blue2",      # Cox train
              "red2",       # Comp train
              "skyblue",    # Cox test
              "orange2"     # Comp test
     ),
     6
   ),
   
   border=rep(
     c(
       "blue4",
              "red4",
              "skyblue4",
              "darkorange4"
     ),
     6
   ),
   
   
   xaxt="n",
   yaxt="n",
   
   outline=TRUE,
   
   ylim=c(0.1,0.9),
   
   ylab="Predictive performance",
   
   lwd=2
   
 )
 
 
 
 #############################
 # x-axis: L groups
 #############################
 
 axis(
   1,
   at=c(
     4.5,
     13.5,
     22.5
   ),
   labels=c(
     expression(L==1),
     expression(L==3),
     expression(L==4)
   ),
   cex.axis=1.5
 )
 
 
 
 #############################
 # y-axis
 #############################
 
 axis(
   2,
   at=seq(0.1,0.9,0.1),
   las=1
 )
 
 
 
 #############################
 # dashed separation lines
 #############################
 
 abline(
   v=9,
   lty=2,
   lwd=2
 )
 
 abline(
   v=18,
   lty=2,
   lwd=2
 )
 
 
 
 #############################
 # iAUC / iBS labels
 #############################
 
 mtext(
   c(
     "iAUC",
     "iBS",
     "iAUC",
     "iBS",
     "iAUC",
     "iBS"
   ),
   side=1,
   at=c(
     2.5,
     6.5,
     11.5,
     15.5,
     20.5,
     24.5
   ),
   line=0.7,
   cex=1
 )
 
 
 
 #############################
 # legend
 #############################
 
 legend(
   "topright",
   
   legend=c(
     "Classic Cox (Train)",
     "Cox-type Polyhazard (Train)",
     "Classic Cox (Test)",
     "Cox-type Polyhazard (Test)"
   ),
   
   fill=c(
     "blue2",
            "red2",
            "skyblue",
            "orange2"
   ),
   
   border=c(
     "blue4",
            "red4",
            "skyblue4",
            "darkorange4"
   ),
   
   bty="o",
   
   cex=1.25
 )
 
  dev.off()
  
  