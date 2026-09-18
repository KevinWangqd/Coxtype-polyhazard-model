load("final_gene_compare_BIC.Rdata")
result$eta_matrix[transformed_data$status==1,]

eta_matrix <- result$eta_matrix
# Subset eta matrix for event patients
eta_event <- eta_matrix[transformed_data$status == 1, ]

# Number of competing groups
L <- ncol(eta_event)

# Function to compute normalized entropy for one patient
compute_H <- function(eta_i) {
  H_i <- -sum(eta_i * log(eta_i)) / log(L)
  return(H_i)
}

# Apply to all selected patients
H_values <- apply(eta_event, 1, compute_H)

# Optional: attach to data frame
result <- data.frame(
  H = H_values
)

# Quick summary
summary(H_values)

H_values

class <- ifelse(H_values < 0.3333, "single",
                       ifelse(H_values < 0.666, "double", "triple"))

max_col <- max.col(eta_matrix, ties.method = "first")
max_col

patient_strata <- as.data.frame(cbind(time=transformed_data$time[transformed_data$status==1],
                        result$eta_matrix[transformed_data$status==1,],
                        class,
                        group=max_col[transformed_data$status==1]))
patient_strata[111,"group"] <- "2,3"

table(patient_strata$group)


##KM curve
library(survival)
library(timeROC)
library(survminer)
patient_strata_KM <- patient_strata[-111,]
patient_strata_KM$event <- 1

fit <- survfit(Surv(as.numeric(patient_strata_KM$time), patient_strata_KM$event) ~ group, data = patient_strata_KM)

s <- summary(fit)

df <- data.frame(
  time = s$time,
  surv = s$surv,
  lower = s$lower,
  upper = s$upper,
  strata = gsub(".*=", "", s$strata)
)





#jpeg("strata_curve.jpeg", quality = 100, units = "in", width = 10, height = 8, res = 300)
plot(
  NA,
  xlim = c(0, 1200),
  ylim = c(0, 1),
  xlab = "Time (days)",
  ylab = "Survival", main="Stratified Survival Curves"
)

cols <- c("blue", "orange", "red")

for (k in unique(df$strata)) {
  d <- df[df$strata == k, ]
  lines(d$time, d$surv, type = "s", lwd = 2,
        col = cols[which(unique(df$strata) == k)])
  lines(
    d$time,
    d$lower,
    type = "s",
    lwd = 2,
    col = cols[which(unique(df$strata) == k)],
    lty = 2
  )
  
  lines(
    d$time,
    d$upper,
    type = "s",
    lwd = 2,
    col = cols[which(unique(df$strata) == k)],
    lty = 2
  )
}

legend("topright", legend=c("Group 1", "Group 2", "Group 3"),
       col=c("blue", "orange", "red"), lwd=2)

#dev.off()









# Venn diagram
library(grid)
library(eulerr)

patient_strata$id <- seq(1,nrow(patient_strata),1)
A <- patient_strata$id[which(patient_strata$group=="1")]
B <- patient_strata$id[which(patient_strata$group=="2"|patient_strata$group=="2,3")]
C <- patient_strata$id[which(patient_strata$group=="3"|patient_strata$group=="2,3")]



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
  
  cat.pos = c(0, 0, 0),
  cat.dist = c(-0.07, -0.07, -0.13)
)
jpeg("Venn_LUAD.jpeg", quality = 100, units = "in", width = 10, height = 8, res = 300)

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

