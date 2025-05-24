# -------------------- CARGA LAS LIBRERÍAS --------------------
library(vars)
library(MTS)
library(tsDyn)
library(ARDL)
library(urca)
library(fUnitRoots)

# -------------------- Descarga de los precios de cierre ajustados (22-abr-2019 → 31-dic-2022) --------------------

tickers <- c("RYLD", # Global X Russell 2000 Covered Call ETF :contentReference[oaicite:0]{index=0}
             "TSLA", # Tesla Inc. :contentReference[oaicite:1]{index=1}
             "TM",   # Toyota Motor Corp. ADR :contentReference[oaicite:2]{index=2}
             "GM")   # General Motors Co. :contentReference[oaicite:3]{index=3}

start_date <- "2019-04-22"
end_date   <- "2022-12-31"

getSymbols(tickers,
           src   = "yahoo",
           from  = start_date,
           to    = end_date,
           adjust = TRUE)

# Extrae únicamente el cierre ajustado
prices_list <- lapply(tickers, function(sym) Ad(get(sym)))
names(prices_list) <- tickers

# Unifica en un solo objeto xts
prices_xts <- do.call(merge, prices_list)
colnames(prices_xts) <- tickers

# -------------------- Visualizacion de las series para chequear comportamiento y outliers --------------------


library(dygraphs)

dygraph(prices_xts, main = "Precios de cierre ajustados (2019-2022)") %>%
  dyOptions(colors = RColorBrewer::brewer.pal(4, "Set1"))


n_total  <- nrow(prices_xts)
n_test   <- ceiling(0.05 * n_total)      # ≈ 5 %
n_train  <- n_total - n_test

train_xts <- prices_xts[1:n_train, ]
test_xts  <- prices_xts[(n_train+1):n_total, ]

# --------------------  Crea los sets train / test (95 % – 5 %) --------------------


n_total  <- nrow(prices_xts)
n_test   <- ceiling(0.05 * n_total)      # ≈ 5 %
n_train  <- n_total - n_test

train_xts <- prices_xts[1:n_train, ]
test_xts  <- prices_xts[(n_train+1):n_total, ]


# --------------------  Diagnóstico de las series del train --------------------

library(urca)

adf_results <- lapply(colnames(train_xts), function(sym){
  test <- ur.df(train_xts[, sym], type = "drift", lags = 5)
  c(Stat = test@teststat[1],  # estadístico τ₁
    Pval = test@cval[1,"5pct"])   # valor crítico al 5 %
})
adf_df <- do.call(rbind, adf_results)
rownames(adf_df) <- colnames(train_xts)
adf_df

#

johansen <- ca.jo(train_xts, type = "trace", K = 2, ecdet = "const")
summary(johansen)




# ------------------------- VAR en primeras diferencias -------------------------
# ---- 2.1  Serie en diferencias ----
diff_train <- diff(na.omit(train_xts))[-1]

# ---- 2.2  Selección de orden ----
p_opt <- VARselect(diff_train, lag.max = 10)$selection["AIC(n)"]

# ---- 2.3  Ajuste y pronóstico ----
mod_var <- vars::VAR(diff_train, p = p_opt, type = "const")
var_fc  <- predict(mod_var, n.ahead = n_test)$fcst

# ---- 2.4  Reconstruir niveles ----
var_pred_mat <- sapply(names(var_fc), \(col){
  last_lvl <- as.numeric(last(train_xts[, col]))
  last_lvl + cumsum(var_fc[[col]][, 1])
})
var_pred_xts <- xts(var_pred_mat, order.by = index(test_xts))
colnames(var_pred_xts) <- tickers


# ------------------------- Pronóstico ARDL -------------------------

# ---- 3.1  Búsqueda automática ----
auto_mod <- auto_ardl(RYLD ~ TSLA + TM + GM,
                      data = data.frame(coredata(train_xts)),
                      max_order = 6)
mod_ardl <- auto_mod$best_model

# ---- 3.2  Pronóstico h = n_test ----
if (!exists("L", mode = "function"))  # wrapper para lags
  L <- function(x, k = 1) stats::lag(x, -k)

ardl_raw <- predict(mod_ardl, n.ahead = n_test)     # vector ts
ardl_vec <- as.numeric(ardl_raw)
len_fc   <- length(ardl_vec)                        # puede ser 46-47
real_vec <- coredata(test_xts$RYLD)[seq_len(len_fc)]

#-------------------------- Metricas de Error (Rmse,MAE) ---------------------

library(forecast)

###############################################################################
# 4.  Cálculo de métricas  (RMSE y MAE) ───────────────────────────────────────
###############################################################################

# ------- 4.1  VAR en diferencias (4 series) -------
rmse <- function(e) sqrt(mean(e^2,  na.rm = TRUE))
mae  <- function(e)        mean(abs(e), na.rm = TRUE)
test_mat <- coredata(test_xts)          #  matriz N × 4 con los valores reales

err_mat <- var_pred_mat - test_mat        # matriz N × 4 de errores

metrics_var <- t(apply(err_mat, 2, \(e) c(RMSE = rmse(e), MAE = mae(e))))
print(metrics_var)
#            RMSE        MAE
# RYLD   ………
# TSLA   ………
# TM     ………
# GM     ………

# ---------- 4.2  ARDL  (solo RYLD)  ----------

e_ardl <- ardl_vec - real_vec          # vector de errores (con algún NA)
metrics_ardl <- c(RMSE = rmse(e_ardl),
                  MAE  = mae(e_ardl))
print(metrics_ardl)

######## Escoger modelo Ganador -------------------------------
# ---------- 5.  Decide modelo ganador ----------
if (metrics_var["RYLD","RMSE"] < metrics_ardl["RMSE"]) {
  winner <- "VAR"
  cat("Modelo ganador: VAR (menor RMSE en RYLD)\n")
} else {
  winner <- "ARDL"
  cat("Modelo ganador: ARDL (menor RMSE en RYLD)\n")
}


#---------------------- Pronóstico 10 días del modelo ganador--------------


h <- 10
if (winner == "VAR") {
  fut_fc  <- predict(mod_var, n.ahead = h)$fcst
  fut_mat <- sapply(names(fut_fc), \(col){
    last_lvl <- as.numeric(last(prices_xts[, col]))
    last_lvl + cumsum(fut_fc[[col]][,1])
  })
  future_xts <- xts(fut_mat,
                    order.by = seq(max(index(prices_xts))+1,
                                   length.out = h, by = "days"))
} else {
  fut_vec <- predict(mod_ardl, n.ahead = h)
  future_xts <- xts(matrix(as.numeric(fut_vec), ncol = 1),
                    order.by = seq(max(index(prices_xts))+1,
                                   length.out = h, by = "days"))
  colnames(future_xts) <- "RYLD"
}
plot(future_xts, main = paste("Pronóstico 10 días –", winner))



# Guardar gráfico de función impulso-respuesta para shock en TSLA
library(vars)

png("irf_tsla.png", width = 800, height = 600)  # crea archivo PNG
plot(irf(mod_var,
         impulse = "TSLA",
         response = c("RYLD", "TM", "GM"),
         n.ahead = 10,
         boot = TRUE))
dev.off()  # cierra y guarda el archivo

getwd()


# Guardar gráfico de pronóstico a 10 días (modelo ganador)
png("forecast_10d.png", width = 800, height = 600)
plot(future_xts,
     main = paste("Pronóstico a 10 días – modelo:", winner),
     ylab = "Precio estimado",
     col = "blue")
dev.off()

# === TABLA DE MÉTRICAS COMPLETA (VAR vs. ARDL) =========================
library(dplyr) ; library(knitr)

#  a) métrica-error del VAR (RYLD)   – ya tenías metrics_var
var_err  <- metrics_var["RYLD", c("RMSE","MAE")]

#  b) métrica-error del ARDL
rmse <- function(e) sqrt(mean(e^2 , na.rm = TRUE))
mae  <- function(e)        mean(abs(e), na.rm = TRUE)
e_ardl        <- ardl_vec - real_vec
ardl_err      <- c(RMSE = rmse(e_ardl), MAE = mae(e_ardl))

#  c) fusiona en un data-frame
tabla_error <- rbind(
  VAR  = signif(var_err , 6),
  ARDL = signif(ardl_err, 6)
)
kable(tabla_error, caption = "Errores de pronóstico sobre RYLD")


# === BOUND-TEST DE COINTEGRACIÓN (ARDL) ================================
library(ARDL)
bt <- bounds_f_test(mod_ardl, case = 2)   # const sin tendencia
print(bt)

#  ——— rescata estadístico y p-valor por si quieres citarlos explícitos
Fval <- unclass(bt)$statistic
pval <- unclass(bt)$p.value


# === ESTABILIDAD DEL VAR: RAÍCES DEL POLINOMIO =========================
roots <- vars::roots(mod_var)
print(roots)                   # valores propios
all( Mod(roots) < 1 )          # TRUE = estable


