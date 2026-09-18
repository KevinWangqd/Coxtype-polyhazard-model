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

R <- 30

iAUC_comp <- numeric(R)
iBS_comp  <- numeric(R)

iAUC_cox  <- numeric(R)
iBS_cox   <- numeric(R)

screen_genes <- vector("list", R)

n <- nrow(data)

for (r in 1:R) {
  
  set.seed(1234+70+ r)  
  
  # -----------------------
  # Screening (L = 1)
  # -----------------------
  
  
  res_high <- high_select(time = data$time_day, status = data$status, X = X_filtered, n_group = 4,
    size_group = 5, init_iauc = 0.68, in_iauc = 0.02, tar_P = 30, n_draw = 100)
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


save(screen_genes, file = "screen_genes_L4_4.RData")

#> iAUC_comp
#0.8131438 0.7977907 0.8062301 0.8300397 0.7636264 0.7518021 0.7573579 0.8170304 0.7898426 0.7705475
#0.7609617 0.7806512 0.8085137 0.7672687 0.7727560 0.7947572 0.7976804 0.7828622 0.7912823 0.8260511
#0.7822839 0.7931552 0.7194220 0.7749348 0.7955915 0.8405565 0.7266478 0.7567022 0.8033978 0.8231020
#0.8516680 0.7446725 0.7564383 0.7414106 0.7861838 0.8270517 0.7559529 0.7839387 0.7753751 0.7975642
#0.7209541 0.7669597 0.7550682 0.7777258 0.8424531 0.8295877 0.7671770 0.7656993 0.7818238 0.7297188
#0.7738342 0.7739442 0.7402360 0.7776093 0.8123381 0.7790557 0.7661994 0.7661793 0.8094738 0.7732627
#0.7496809 0.7742788 0.7774752 0.8174882 0.7272768 0.7838120 0.7468866 0.7661721 0.7939808 0.8280026
#0.7867797 0.7357137 0.7042474 0.7563539 0.7609561 0.7538977 0.7585013 0.8095645 0.7494953 0.8276955
#0.8332267 0.7465034 0.8292954 0.7483633 0.8248592 0.8122345 0.7970995 0.7714669 0.7828118 0.8061732
#0.7327704 0.7691776 0.7657707 0.7682767 0.7848784 0.8012251 0.7296746 0.7743745 0.7949178 0.8175600

#> iAUC_cox
#0.7729976 0.7293842 0.7502121 0.7211444 0.7298647 0.6489102 0.7443006 0.7524350 0.7242993 0.7512588
#0.7054086 0.7450340 0.7690375 0.7315980 0.7432485 0.7286277 0.6976348 0.7872183 0.7597657 0.7991863
#0.6945518 0.6711186 0.6529653 0.6762822 0.7460450 0.7885611 0.7072615 0.7119547 0.7456446 0.7429199
#0.7791393 0.6916729 0.7063398 0.7607820 0.7378104 0.7815039 0.6982010 0.7188243 0.7068851 0.7574069
#0.6892596 0.7405161 0.7388582 0.7723887 0.7545421 0.7008414 0.6705210 0.7471817 0.7812654 0.7276981
#0.7673890 0.7324676 0.7259344 0.7325547 0.7472269 0.7338019 0.7515197 0.6731888 0.7385907 0.7497515
#0.7071249 0.7190794 0.7372087 0.8094726 0.7109287 0.7080951 0.7426021 0.7057066 0.6972642 0.7654664
#0.7658189 0.7183176 0.6913945 0.7177417 0.7276455 0.7473476 0.7333650 0.7105368 0.7336586 0.7689137
#0.8224934 0.6820588 0.7756072 0.7364896 0.7576542 0.7085811 0.7443300 0.7303713 0.7491275 0.7847340
#0.6639294 0.6980894 0.6992348 0.7319149 0.7575742 0.6659425 0.7348386 0.7109509 0.7473501 0.7176623

#> iBS_comp
#0.1198217 0.1217047 0.1224624 0.1001428 0.1272208 0.1321676 0.1386033 0.1137062 0.1231767 0.1269105
#0.1372027 0.1287543 0.1213114 0.1296481 0.1337611 0.1251134 0.1244428 0.1159633 0.1199538 0.1209714
#0.1202239 0.1259931 0.1408532 0.1312879 0.1144992 0.1033331 0.1443762 0.1318689 0.1181424 0.1152361
#0.1101712 0.1372355 0.1382077 0.1375985 0.1271017 0.1134344 0.1362874 0.1299792 0.1343694 0.1206626
#0.1427216 0.1227385 0.1300658 0.1299204 0.1099233 0.1189392 0.1269770 0.1241476 0.1296662 0.1324397
#0.1277503 0.1328128 0.1418595 0.1266978 0.1179429 0.1236720 0.1312914 0.1326119 0.1186300 0.1180918
#0.1352189 0.1295508 0.1265032 0.1125513 0.1471291 0.1333172 0.1358963 0.1305527 0.1221048 0.1148033
#0.1130083 0.1401010 0.1441883 0.1328090 0.1309068 0.1435475 0.1357357 0.1198758 0.1314294 0.1190094
#0.1082715 0.1342841 0.1132139 0.1364934 0.1200534 0.1195602 0.1265121 0.1281501 0.1346116 0.1242662
#0.1393041 0.1340408 0.1239080 0.1254901 0.1243734 0.1224992 0.1386344 0.1300990 0.1211693 0.1192542


#> iBS_cox
#0.1291999 0.1360221 0.1341446 0.1344001 0.1399864 0.1486320 0.1374915 0.1364430 0.1430061 0.1370297
#0.1445201 0.1306140 0.1312147 0.1397570 0.1353479 0.1366594 0.1407533 0.1302683 0.1328007 0.1213805
#0.1466049 0.1464162 0.1478442 0.1493018 0.1307541 0.1323311 0.1452587 0.1403482 0.1391375 0.1366493
#0.1281613 0.1476556 0.1400300 0.1344970 0.1376741 0.1227872 0.1428154 0.1398281 0.1479844 0.1324701
#0.1501693 0.1397873 0.1397516 0.1274783 0.1371706 0.1421628 0.1404526 0.1382987 0.1190842 0.1337770
#0.1350674 0.1422742 0.1455892 0.1394311 0.1344294 0.1438400 0.1351001 0.1457208 0.1374879 0.1427677
#0.1465136 0.1471648 0.1389398 0.1217141 0.1471334 0.1424033 0.1381475 0.1464947 0.1392974 0.1319798
#0.1257868 0.1432354 0.1448145 0.1415812 0.1392107 0.1390571 0.1390278 0.1460175 0.1397884 0.1335328
#0.1179091 0.1511490 0.1272508 0.1390140 0.1353193 0.1399353 0.1358513 0.1427572 0.1378696 0.1296225
#0.1453747 0.1441423 0.1416208 0.1372935 0.1281149 0.1489761 0.1395453 0.1414037 0.1325696 0.1418686
> 