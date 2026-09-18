source("Cox_competing_correct.R")
source("Cox_competing_algorithm.R")
load("LUAD_process.RData")
library(riskRegression)
library(MASS)


#Function for Cox model
generate_hazard_matrix <- function(model, data, time_points, model_type = "cox") {
  
  n_samples <- nrow(data)
  hazard_matrix <- matrix(NA, nrow = n_samples, ncol = length(time_points))
  
  if (model_type == "cox") {
    # Cox PH model
    base_hazard <- basehaz(model, centered = FALSE)
    lp <- predict(model, type = "lp")
    
    for (j in 1:length(time_points)) {
      t <- time_points[j]
      idx <- max(which(base_hazard$time <= t))
      baseline_cum_hazard <- ifelse(length(idx) > 0, base_hazard$hazard[idx], 0)
      hazard_matrix[, j] <- baseline_cum_hazard * exp(lp)
    }
    
  } else if (model_type == "weibull") {
    # Weibull AFT model
    scale <- model$scale
    lp <- predict(model, newdata = data, type = "linear")
    
    for (j in 1:length(time_points)) {
      t <- time_points[j]
      z <- (log(t) - lp) / scale
      hazard_matrix[, j] <- (1/(scale * t)) * exp(z)^(1/scale) / gamma(1/scale)
    }
  }
  
  colnames(hazard_matrix) <- paste0("t_", round(time_points, 1))
  return(hazard_matrix)
}


## Standardization
#colnames(data)
X_filtered = scale(X_filtered)

R <- 20

iAUC_comp <- numeric(R)
iBS_comp  <- numeric(R)

iAUC_cox  <- numeric(R)
iBS_cox   <- numeric(R)

screen_genes <- vector("list", R)

n <- nrow(data)

for (r in 1:R) {
  
  set.seed(1234 + r)  
  
  # -----------------------
  # Screening (L = 1)
  # -----------------------
  
  
  res_high <- high_select(time = data$time_day, status = data$status, X = X_filtered, n_group = 1,
    size_group = 5, init_iauc = 0.6, in_iauc = 0.02, tar_P = 50, n_draw = 100)
  message("Iteration", r)
  X_sel <- res_high$X_selected
  screen_genes[[r]] <- colnames(X_sel)
  # =========================================================
  # 1. COMPETING COX
  # =========================================================
  
  res_forward <- forward_select(data$time_day, data$status, X_candidate = X_sel, max_add = 20, lambda1 = log(n), lambda2 = log(n))
  
  transformed_data <- data_transform(data$time_day, data$status, res_forward$X, res_forward$label)
  
  result <- Cox_com_fit(transformed_data, maxit = 1e3, tolerance = 1e-6, frailty = FALSE)
  
  iAUC_comp[r] <- compute_iAUC(transformed_data, result = result)
  iBS_comp[r]  <- compute_iBS(transformed_data, result = result)
  
  # =========================================================
  # 2. CLASSIC COX
  # =========================================================
  
  cov_names <- paste0("`", colnames(X_sel), "`")
  
  formula_str <- paste("Surv(data$time_day, data$status) ~", paste(cov_names, collapse = " + "))
  
  cox_full <- coxph(as.formula(formula_str), data = as.data.frame(X_sel))
  
  cox_aic_result <- stepAIC(cox_full,direction = "both",trace = FALSE,k = log(n))
  
  hazard_cox <- generate_hazard_matrix(cox_aic_result, data, quantile(data$time_day,probs = seq(0.01, 0.99, by = 0.05)),"cox")
  
  eval_times <- quantile(data$time_day,probs = seq(0.01, 0.99, by = 0.05))
  
  time_AUC <- Score(object = list(Cox = hazard_cox), formula = Surv(time, status) ~ 1,
    data = data.frame(time = data$time_day, status = data$status), times = eval_times, metrics = "AUC",cens.model = "km")
  
  iAUC_cox[r] <- mean(time_AUC$AUC$score$AUC)
  
  bs_res <- Score(object = list(Cox = 1 - exp(-hazard_cox)), formula = Surv(time, status) ~ 1,
    data = data.frame(time = data$time_day,status = data$status),times = eval_times,metrics = "Brier",cens.model = "km")
  
  iBS_cox[r] <- mean(bs_res$Brier$score$Brier[bs_res$Brier$score$model == "Cox"]
  )
}



screen_freq <- sort(
  table(unlist(screen_genes)),
  decreasing = TRUE
)

head(screen_freq, 20)


save(screen_genes, file = "screen_genes5.RData")

#> iAUC_comp
#[extra1] 0.7626008 0.7665306 0.8020362 0.7832854 0.7959035 0.7772209 0.7829402 0.7673342 0.7699757 0.7711624 
# 0.7871534 0.7753712 0.7854807 0.7341317 0.7505316 0.7807907 0.7868403 0.7573279 0.6944551 0.8063144
# 0.8123762 0.7112249 0.7795758 0.7660328 0.7337807 0.7454200 0.7855498 0.7927580 0.7735511 0.7607513
# 0.7486021 0.7856268 0.7952702 0.8103739 0.7680897 0.7846280 0.7689846 0.8358357 0.7640186 0.7749442 
# 0.7729393 0.8100708 0.7572898 0.7919729 0.7654921 0.7877215 0.8299627 0.7707085 0.7699497 0.7580682
# 0.7599931 0.7932133 0.7490918 0.7496140 0.7111266 0.7468507 0.7766175 0.7505552 0.8238982 0.7574449
# 0.7566997 0.8047000 0.7945701 0.7825308 0.7978571 0.7941524 0.7413790 0.7921096 0.8271432 0.7640933 
# 0.7450585 0.7657502 0.7440571 0.7919067 0.7028183 0.7848963 0.8311656 0.7703074 0.7368550 0.7831224
# 0.8045391 0.7692477 0.8139626 0.7601924 0.8352188 0.8018859 0.7841650 0.7762743 0.7441949 0.8040918
# 0.7602492 0.7825908 0.8249264 0.7782899 0.7676128 0.7955324 0.7828969 0.7453775 0.7970101 0.7862667
# 0.7777381 0.7922761 0.8019664 0.7419691 0.7782117 0.8376357 0.7827435 0.7992984 0.7900151 0.7645697

#> iAUC_cox
##[extra1] 0.7708846 0.7111242 0.7110673 0.7553296 0.7629802 0.7749567 0.7441684 0.7374366 0.7248836 0.7443386 
# 0.6622743 0.7363446 0.6779541 0.6920697 0.7287914 0.7362365 0.7396885 0.7114318 0.6789612 0.7450919
# 0.7932741 0.7311272 0.6973281 0.7196663 0.7106351 0.7539383 0.7301710 0.7567461 0.6956081 0.7109749
# 0.6853270 0.6893151 0.7839294 0.7396787 0.7024671 0.7209175 0.7295712 0.7661067 0.6922477 0.7420136 
# 0.7676890 0.7524444 0.7585549 0.7621763 0.6876532 0.7268794 0.7659084 0.7507756 0.7094337 0.7608048
# 0.7555195 0.7592119 0.7103779 0.7423176 0.7400681 0.6876211 0.7449164 0.6922470 0.7684941 0.7190827
# 0.7215384 0.7447946 0.7461056 0.7425052 0.7991067 0.7471862 0.7229494 0.6800379 0.7131762 0.7255662 
# 0.7061370 0.7434924 0.7128699 0.7212981 0.6707690 0.7354928 0.7638456 0.7317444 0.7564286 0.8032752
# 0.7700595 0.7274775 0.7468299 0.7546742 0.7990442 0.7534475 0.7221314 0.7296964 0.7284668 0.7750414
# 0.7486062 0.7372242 0.7543685 0.7421911 0.7625274 0.7264472 0.7537401 0.7809866 0.7658182 0.7495418
# 0.7065821 0.7441831 0.7410659 0.7194816 0.7085193 0.7316741 0.7561190 0.7297229 0.7596213 0.6857930

#> iBS_comp
##[extra1] 0.1335848 0.1280144 0.1186972 0.1292043 0.1198150 0.1226134 0.1269330 0.1278504 0.1286248 0.1241727 
#0.1244177 0.1151498 0.1266206 0.1372547 0.1353839 0.1320515 0.1280040 0.1315851 0.1486093 0.1200250
#0.1186388 0.1440997 0.1204120 0.1353541 0.1362521 0.1310786 0.1280749 0.1243658 0.1304423 0.1322834
#0.1370573 0.1234416 0.1267853 0.1149803 0.1322712 0.1255571 0.1183439 0.1122503 0.1139659 0.1282515 
#0.1286665 0.1204850 0.1311251 0.1265885 0.1293127 0.1238934 0.1180000 0.1311171 0.1270035 0.1355656
#0.1333301 0.1262951 0.1286059 0.1276714 0.1405389 0.1376594 0.1207596 0.1276805 0.1183763 0.1240424
#0.1307612 0.1229343 0.1223024 0.1275181 0.1143303 0.1213023 0.1363237 0.1247373 0.1148959 0.1142279
#0.1331512 0.1305451 0.1349539 0.1258180 0.1466768 0.1265256 0.1085892 0.1275653 0.1360018 0.1241237
#0.1207926 0.1343465 0.1245921 0.1344349 0.1131594 0.1215015 0.1283086 0.1294686 0.1342879 0.1138249
#0.1279519 0.1300638 0.1155681 0.1221712 0.1283835 0.1224092 0.1291673 0.1311013 0.1253130 0.1289918
#0.1251749 0.1271588 0.1151778 0.1345623 0.1283088 0.1050107 0.1300249 0.1151168 0.1198985 0.1254573

#> iBS_cox
##[extra1] 0.1298870 0.1388282 0.1361958 0.1294688 0.1313338 0.1320827 0.1359578 0.1413415 0.1375250 0.1392121 
#0.1475320 0.1269384 0.1475638 0.1423151 0.1420263 0.1409196 0.1394473 0.1362529 0.1482938 0.1351210
#0.1266493 0.1362812 0.1412623 0.1409405 0.1412446 0.1354596 0.1348508 0.1307474 0.1453854 0.1424980
#0.1429165 0.1468827 0.1283797 0.1340302 0.1443818 0.1458786 0.1371966 0.1323894 0.1390620 0.1352179 
#0.1314251 0.1272597 0.1350409 0.1357908 0.1450090 0.1399923 0.1364769 0.1354151 0.1393107 0.1365885
#0.1302688 0.1323993 0.1418188 0.1353290 0.1363902 0.1462205 0.1319018 0.1441157 0.1274703 0.1385952
#0.1422575 0.1335423 0.1354117 0.1383949 0.1231941 0.1333851 0.1395447 0.1447714 0.1413769 0.1361308
#0.1400074 0.1360103 0.1430528 0.1377829 0.1491739 0.1366658 0.1346091 0.1315699 0.1393623 0.1193683
#0.1297142 0.1407832 0.1394329 0.1410475 0.1244449 0.1324946 0.1401501 0.1416583 0.1389357 0.1280546
#0.1371165 0.1371770 0.1351826 0.1371389 0.1271058 0.1395416 0.1357056 0.1284869 0.1343253 0.1344579
# 0.1438533 0.1355569 0.1316225 0.1440610 0.1401485 0.1387609 0.1367944 0.1334086 0.1296135 0.1426438










