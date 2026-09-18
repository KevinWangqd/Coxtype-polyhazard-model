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
  
  set.seed(1234+100+ r)  
  
  # -----------------------
  # Screening (L = 1)
  # -----------------------
  
  
  res_high <- high_select(time = data$time_day, status = data$status, X = X_filtered, n_group = 3,
    size_group = 5, init_iauc = 0.64, in_iauc = 0.02, tar_P = 30, n_draw = 100)
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
  
  time_AUC <- Score(
    object = list(Cox = hazard_cox),
    formula = Surv(time, status) ~ 1,
    data = data.frame(
      time = data$time_day,
      status = data$status
    ),
    times = eval_times,
    metrics = "AUC",
    cens.model = "km"
  )
  
  iAUC_cox[r] <- mean(time_AUC$AUC$score$AUC)
  
  bs_res <- Score(
    object = list(Cox = 1 - exp(-hazard_cox)),
    formula = Surv(time, status) ~ 1,
    data = data.frame(
      time = data$time_day,
      status = data$status
    ),
    times = eval_times,
    metrics = "Brier",
    cens.model = "km"
  )
  
  iBS_cox[r] <- mean(
    bs_res$Brier$score$Brier[
      bs_res$Brier$score$model == "Cox"
    ]
  )
}



screen_freq <- sort(
  table(unlist(screen_genes)),
  decreasing = TRUE
)

head(screen_freq, 20)


save(screen_genes, file = "screen_genes_L3_4.RData")

#> iAUC_comp
#0.7633864 0.7965820 0.7530036 0.7373329 0.7554363 0.7558632 0.8010285 0.7950515 0.8280431 0.7656947
#0.7784494 0.7706692 0.7496272 0.7547910 0.8045057 0.7806005 0.7332885 0.7637678 0.8005483 0.8065243
#0.8063953 0.8080826 0.7412320 0.7600113 0.7767390 0.7903211 0.7867214 0.7493701 0.7336277 0.7665932
#0.7404717 0.7617085 0.7694875 0.8289812 0.7754568 0.7883995 0.7672980 0.7937286 0.7697246 0.7377999
#0.8081451 0.7435771 0.7279068 0.7936283 0.7961826 0.8106398 0.8074696 0.7380957 0.7251289 0.7903975
#0.7637514 0.7798154 0.7557540 0.8098647 0.8046606 0.7609453 0.7282678 0.8065028 0.7458734 0.8042045
#0.7935531 0.7690361 0.8001410 0.8027654 0.7760629 0.7612748 0.7723661 0.7774417 0.7345069 0.7404435
#0.7460830 0.7885015 0.7754270 0.8084226 0.7161814 0.7379703 0.7829153 0.8180172 0.8234694 0.7995484
#0.7186198 0.7836940 0.8274631 0.7579317 0.7455576 0.7856005 0.7643280 0.8049317 0.8049253 0.7484045
#0.7739640 0.8023055 0.7889065 0.8130569 0.7696961 0.7812724 0.7515951 0.7534380 0.7600483 0.8287395

#> iAUC_cox
#0.6906354 0.7570114 0.7201953 0.7104767 0.6897508 0.6848551 0.7304398 0.7558239 0.7669783 0.7218309
#0.7296149 0.7251585 0.6840564 0.6940562 0.7455131 0.7585277 0.6556522 0.7534995 0.7296365 0.7220607
#0.7582537 0.7159693 0.7153983 0.7136457 0.6909528 0.7577259 0.7395465 0.7140705 0.6693874 0.7415784
#0.7047361 0.7323447 0.7226541 0.7339419 0.7367800 0.7309676 0.7455060 0.7126055 0.7596915 0.6810388
#0.7420748 0.6976668 0.7211295 0.7362465 0.7448824 0.7159282 0.7274897 0.6449640 0.6851214 0.7331415
#0.7155161 0.7419939 0.7398779 0.7601283 0.7895662 0.7470932 0.7257820 0.7228685 0.7416428 0.7401236
#0.7043720 0.7132300 0.7364183 0.7469391 0.7185508 0.7030323 0.7427212 0.7345570 0.6875878 0.7492381
#0.7328520 0.7470100 0.7350165 0.7883132 0.6943755 0.7102489 0.7370263 0.7395470 0.7694011 0.7478364
#0.6977990 0.7434715 0.7247217 0.7226377 0.6688858 0.7220658 0.6955132 0.7354095 0.7824698 0.6825027
#0.7332071 0.7405073 0.7529280 0.7504806 0.7321447 0.7287386 0.7109018 0.7150961 0.7103165 0.7570188

#> iBS_comp
#0.1295265 0.1241882 0.1339111 0.1404260 0.1322039 0.1279165 0.1239184 0.1273238 0.1060112 0.1290404
#0.1238088 0.1320468 0.1359745 0.1336755 0.1121998 0.1307028 0.1368062 0.1289435 0.1262555 0.1242947
#0.1275536 0.1199001 0.1344116 0.1323402 0.1285348 0.1247405 0.1282965 0.1393414 0.1439064 0.1311128
#0.1386013 0.1319949 0.1288633 0.1127773 0.1372832 0.1295303 0.1358095 0.1230358 0.1277353 0.1424227
#0.1214668 0.1374015 0.1376489 0.1255040 0.1213892 0.1216623 0.1142248 0.1382494 0.1372896 0.1299302
#0.1315263 0.1255463 0.1278430 0.1091150 0.1215224 0.1376032 0.1377204 0.1174930 0.1376493 0.1253592
#0.1327064 0.1367901 0.1271279 0.1248439 0.1290240 0.1307271 0.1325289 0.1333559 0.1412108 0.1276312
#0.1340389 0.1286203 0.1333632 0.1232731 0.1349098 0.1325506 0.1282638 0.1183736 0.1153432 0.1265288
#0.1209117 0.1256833 0.1181728 0.1330652 0.1402201 0.1291814 0.1324603 0.1189979 0.1253500 0.1333014
#0.1303570 0.1214878 0.1296386 0.1179226 0.1229253 0.1327336 0.1360985 0.1347398 0.1377462 0.1099981

#> iBS_cox
#0.1392164 0.1337760 0.1413134 0.1401872 0.1427852 0.1420743 0.1390434 0.1379226 0.1258161 0.1357079
#0.1342837 0.1396322 0.1432173 0.1439534 0.1335594 0.1331862 0.1473997 0.1345380 0.1399275 0.1420947
#0.1320286 0.1416779 0.1435278 0.1387000 0.1459331 0.1353881 0.1442828 0.1456481 0.1474119 0.1395564
#0.1429928 0.1429915 0.1396644 0.1338900 0.1396425 0.1413353 0.1377219 0.1402032 0.1361513 0.1493875
#0.1389201 0.1457675 0.1390118 0.1343919 0.1351570 0.1392850 0.1339292 0.1498062 0.1447425 0.1389943
#0.1412545 0.1329514 0.1379489 0.1288584 0.1217648 0.1404219 0.1395549 0.1348918 0.1420842 0.1358815
#0.1453526 0.1435989 0.1391422 0.1322071 0.1390537 0.1432328 0.1348372 0.1367971 0.1465473 0.1344874
#0.1383367 0.1310806 0.1381197 0.1289415 0.1432304 0.1438449 0.1394089 0.1426813 0.1327699 0.1347225
#0.1443829 0.1370708 0.1393118 0.1420211 0.1502113 0.1410262 0.1434632 0.1345499 0.1275502 0.1411551
#0.1349918 0.1310163 0.1326426 0.1336701 0.1398054 0.1390261 0.1442873 0.1441826 0.1436514 0.1333028
