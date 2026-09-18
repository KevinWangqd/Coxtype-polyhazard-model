load("COMA_Cox.Rdata")

##KM curve
library(survival)
library(timeROC)
library(survminer)

patient_strata_KM <- patient_strata

patient_strata_KM$event <- rep(1, nrow(patient_strata_KM))

fit <- survfit(
  Surv(as.numeric(patient_strata_KM$time), patient_strata_KM$event) ~ class,
  data = patient_strata_KM
)

s <- summary(fit)

df <- data.frame(
  time = s$time,
  surv = s$surv,
  lower = s$lower,
  upper = s$upper,
  strata = gsub(".*=", "", s$strata)
)

unique(df$strata)



jpeg("strata_COMA_group.jpeg", 
     quality = 100, 
     units = "in", 
     width = 10, 
     height = 8, 
     res = 300)

plot(
  NA,
  xlim = c(0, 600),
  ylim = c(0, 1),
  xlab = "Time (days)",
  ylab = "Survival",
  main = "Stratified Survival Curves"
)

cols <- c("orange", "blue", "red")

strata_order <- unique(df$strata)

for (k in strata_order) {
  
  d <- df[df$strata == k, ]
  
  col_k <- cols[which(strata_order == k)]
  
  # 95% confidence interval
  lines(
    d$time,
    d$lower,
    type = "s",
    lwd = 2,
    col = col_k,
    lty = 2
  )
  
  lines(
    d$time,
    d$upper,
    type = "s",
    lwd = 2,
    col = col_k,
    lty = 2
  )
  
  # survival curve
  lines(
    d$time,
    d$surv,
    type = "s",
    lwd = 2,
    col = col_k
  )
}

legend(
  "topright",
  legend = c(
    "Single-CF Group",
    "Double-CF Group",
    "Triple-CF Group"
  ),
  col = c("blue", "orange", "red"),
  lwd = 2
)

dev.off()





## strata COMA


patient_strata_KM <- patient_strata[which(patient_strata$class=="single"),]

patient_strata_KM$event <- rep(1,nrow(patient_strata_KM))
fit <- survfit(Surv(as.numeric(patient_strata_KM$time), patient_strata_KM$event) ~ group, data = patient_strata_KM)

s <- summary(fit)

df <- data.frame(
  time = s$time,
  surv = s$surv,
  lower = s$lower,
  upper = s$upper,
  strata = gsub(".*=", "", s$strata)
)

unique(df$strata)





jpeg("strata_COMA.jpeg", 
     quality = 100, 
     units = "in", 
     width = 10, 
     height = 8, 
     res = 300)

plot(
  NA,
  xlim = c(0, 600),
  ylim = c(0, 1),
  xlab = "Time (days)",
  ylab = "Survival",
  main = "Stratified Survival Curves"
)

cols <- c("blue", "orange",  "red")

strata_order <- unique(df$strata)

for (k in strata_order) {
  
  d <- df[df$strata == k, ]
  
  col_k <- cols[which(strata_order == k)]
  
  # 95% confidence interval
  lines(
    d$time,
    d$lower,
    type = "s",
    lwd = 2,
    col = col_k,
    lty = 2
  )
  
  lines(
    d$time,
    d$upper,
    type = "s",
    lwd = 2,
    col = col_k,
    lty = 2
  )
  
  # survival curve
  lines(
    d$time,
    d$surv,
    type = "s",
    lwd = 2,
    col = col_k
  )
}

legend(
  "topright",
  legend = c(
    "Group 1",
    "Group 2",
    "Group 3"
  ),
  col = c("blue", "orange", "red"),
  lwd = 2
)

dev.off()













# Venn diagram
library(grid)
library(eulerr)

patient_strata$id <- seq(1,nrow(patient_strata),1)
A <- patient_strata$id[which(patient_strata$group=="1"|patient_strata$group=="2,1"|patient_strata$group=="1,2"|
                         patient_strata$group=="3,1"|patient_strata$group=="1,3"|patient_strata$group=="1,2,3")]

B <- patient_strata$id[which(patient_strata$group=="2"|patient_strata$group=="2,1"|patient_strata$group=="1,2"|
                               patient_strata$group=="3,2"|patient_strata$group=="2,3"|patient_strata$group=="1,2,3")]


C <- patient_strata$id[which(patient_strata$group=="3"|patient_strata$group=="2,3"|patient_strata$group=="3,2"|
                         patient_strata$group=="3,1"|patient_strata$group=="1,3"|patient_strata$group=="1,2,3")]



library(VennDiagram)
library(grid)


venn.plot <- draw.triple.venn(
  area1 = length(A),
  area2 = length(B),
  area3 = length(C),
  n12 = length(intersect(A, B)),
  n23 = length(intersect(B, C)),
  n13 = length(intersect(A, C)),
  n123 = length(Reduce(intersect, list(A, B, C))),
  fontfamily=1,
  category = c("Group 1", "Group 2", "Group 3"),
  fill = c("#0066FF", "orange", "red"),
  alpha = 0.8,
  cat.fontfamily=1,
  cex = 2.3,
  cat.cex = 2.5,
  cat.pos = c(-35, 35, 0),
  cat.dist = c(-0.04, -0.04, -0.26)
)
jpeg("Venn_COMA.jpeg", quality = 100, units = "in", width = 10, height = 8, res = 300)

pushViewport(viewport(
  x = 0.5, y = 0.55,
  width = 0.95,
  height = 0.8   # <- reduces plot height → creates top margin
))

grid.draw(venn.plot)

# border around plotting region
grid.rect(
  gp = gpar(fill = NA, col = "black", lwd = 1)
)

# title safely inside the full page (not clipped by viewport)
popViewport()

grid.text(
  "Patient Stratification",
  y = unit(0.98, "npc"),   # safely above viewport
  gp = gpar(fontsize = 20, fontface = "bold")
)
dev.off()

