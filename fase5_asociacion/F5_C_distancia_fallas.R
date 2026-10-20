# =====================================================================
# FASE 5 - Parte C (v2). Anomalías y distancia a las fallas
#
# Frente 1 (estructura residual): la distancia a la falla de la
# microzonificación entró a la GWR, así que esta parte mide si lo que el
# modelo no explicó se concentra cerca de las fallas.
#
# Referencia: el inventario (asesor), no el territorio ni el azar espacial.
#   - Kolmogorov-Smirnov de dos muestras: distancia de las anomalías frente
#     al resto de eventos (bilateral y unilateral "más cerca").
#   - Curva: proporción de anómalos por banda de 500 m, con la proporción
#     global como referencia y una envolvente del 95 % obtenida reasignando
#     999 veces la etiqueta de anómalo entre los mismos eventos.
#   - Prueba global: chi² de homogeneidad entre bandas con p por las mismas
#     reasignaciones (¿la proporción cambia con la distancia?).
#
# Entrada: F5_eventos.csv (Parte 0 v2)
# Salidas en fase5/parteC: C_ks.csv, C_bandas.csv, C_distancia_fallas.png
#                          (la figura también se muestra en pantalla)
# =====================================================================

# ------------------------- CONFIGURACIÓN -----------------------------
BASE    <- "D:/Users/tasal/Desktop/Tesis de grado"
ENTRADA <- file.path(BASE, "fase5", "Fase5_csv", "F5_eventos.csv")
SALIDA  <- file.path(BASE, "fase5", "parteC")

F_DIST <- "D_FALLA"
BANDA  <- 500
DMAX   <- 5000          # a partir de aquí, una sola banda abierta
N_SIM  <- 999
COL_EV   <- "grey45"
COL_ANOM <- "#C0392B"
set.seed(2026)

# ------------------------- LECTURA -----------------------------------
dir.create(SALIDA, showWarnings = FALSE, recursive = TRUE)
d <- read.csv(ENTRADA, fileEncoding = "UTF-8-BOM", check.names = FALSE, stringsAsFactors = FALSE)
names(d) <- trimws(sub("^\ufeff", "", names(d)))
if (!"ANOM" %in% names(d)) d$ANOM <- as.integer(tolower(d$ANOMALIA) != "normal")
faltan <- setdiff(c("ANOM", F_DIST), names(d))
if (length(faltan)) stop("Faltan columnas: ", paste(faltan, collapse = ", "))
d <- d[!is.na(d[[F_DIST]]), ]

x  <- d[[F_DIST]]
an <- d$ANOM == 1
p0 <- mean(an)
cat("Eventos:", length(x), "| anomalías:", sum(an), "| proporción global:", round(p0, 4), "\n")

# ------------------------- KOLMOGOROV-SMIRNOV ------------------------
ks2 <- suppressWarnings(ks.test(x[an], x[!an]))
ks1 <- suppressWarnings(ks.test(x[an], x[!an], alternative = "greater"))   # anómalos más cerca
malla <- sort(unique(x))
Fa <- ecdf(x[an])(malla); Fn <- ecdf(x[!an])(malla)
i  <- which.max(abs(Fa - Fn))
ks_tab <- data.frame(
  n_anom = sum(an), n_resto = sum(!an),
  mediana_anom_m = median(x[an]), mediana_resto_m = median(x[!an]),
  D = unname(ks2$statistic), p_bilateral = ks2$p.value, p_unilateral_cerca = ks1$p.value,
  dist_max_separacion_m = malla[i], pct_anom_hasta_ahi = 100 * Fa[i], pct_resto_hasta_ahi = 100 * Fn[i]
)

# ------------------------- CURVA POR BANDAS ---------------------------
cortes <- c(seq(0, DMAX, by = BANDA), Inf)
banda  <- cut(x, cortes, right = FALSE, dig.lab = 6)
n_ev   <- as.vector(table(banda))
cuenta <- function(lab) as.vector(tapply(lab, banda, sum))
n_an   <- cuenta(an); n_an[is.na(n_an)] <- 0
prop   <- ifelse(n_ev > 0, n_an / n_ev, NA)

chi_bandas <- function(k) {                 # chi² de homogeneidad (solo bandas con eventos)
  e <- n_ev * p0; ok <- n_ev > 0
  sum((k[ok] - e[ok])^2 / (e[ok] * (1 - p0)))
}
sim <- replicate(N_SIM, { k <- cuenta(sample(an)); k[is.na(k)] <- 0; k })
prop_sim <- sweep(sim, 1, pmax(n_ev, 1), "/")
env_inf <- apply(prop_sim, 1, quantile, 0.025)
env_sup <- apply(prop_sim, 1, quantile, 0.975)
chi_obs <- chi_bandas(n_an)
chi_sim <- apply(sim, 2, chi_bandas)
p_chi   <- (1 + sum(chi_sim >= chi_obs)) / (N_SIM + 1)

bandas <- data.frame(
  Banda = levels(banda), D_INI = head(cortes, -1), N_EVENTOS = n_ev, N_ANOM = n_an,
  PROP_ANOM = prop, ENV_INF = env_inf, ENV_SUP = env_sup,
  FUERA = ifelse(!is.na(prop) & (prop < env_inf | prop > env_sup), "si", "")
)
ks_tab$chi2_bandas <- chi_obs
ks_tab$p_chi2_bandas <- p_chi
write.csv(ks_tab, file.path(SALIDA, "C_ks.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(bandas, file.path(SALIDA, "C_bandas.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# ------------------------- FIGURA ------------------------------------
dibujar <- function() {
  op  <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1), mgp = c(2.6, 0.7, 0), las = 1,
             cex.main = 1, cex.lab = 0.95, cex.axis = 0.85)
  dec <- options(OutDec = ",")
  on.exit({ par(op); options(dec) })
  
  # (a) Proporción de anómalos por banda frente al inventario
  b  <- bandas[is.finite(cortes[-1]), ]
  km <- (b$D_INI + BANDA / 2) / 1000
  ymax <- max(c(b$PROP_ANOM, b$ENV_SUP), na.rm = TRUE) * 1.45
  plot(km, b$PROP_ANOM, type = "n", ylim = c(0, ymax), xlim = c(0, DMAX / 1000),
       xlab = "Distancia a la falla (km)", ylab = "Proporción de anómalos",
       main = "(a) Proporción por banda")
  polygon(c(km, rev(km)), c(b$ENV_INF, rev(b$ENV_SUP)), col = "grey88", border = NA)
  abline(h = p0, lty = 3)
  lines(km, b$PROP_ANOM, type = "b", pch = 17, col = COL_ANOM, lwd = 1.5)
  text(km, b$PROP_ANOM, labels = paste0(b$N_ANOM, "/", b$N_EVENTOS), pos = 3,
       offset = 0.6, cex = 0.6, col = COL_ANOM)
  legend("topright", bty = "n", cex = 0.75,
         legend = c("Observada (anómalos/eventos)", "Proporción del inventario", "Envolvente 95 %"),
         col = c(COL_ANOM, "black", "grey80"), pch = c(17, NA, 15), lty = c(1, 3, NA), pt.cex = c(1, 1, 2))
  
  # (b) Distribución acumulada y separación máxima de KS
  plot(ecdf(x[!an]), do.points = FALSE, verticals = TRUE, col = COL_EV, lwd = 1.5,
       xlim = c(0, max(x)), xaxt = "n", xlab = "Distancia a la falla (km)",
       ylab = "Proporción acumulada", main = "(b) Distribución acumulada")
  ejes <- pretty(c(0, max(x)))
  axis(1, at = ejes, labels = ejes / 1000)
  plot(ecdf(x[an]), do.points = FALSE, verticals = TRUE, col = COL_ANOM, lwd = 1.5, add = TRUE)
  segments(malla[i], Fn[i], malla[i], Fa[i], lwd = 2)
  legend("bottomright", bty = "n", cex = 0.8, lwd = 1.5, col = c(COL_EV, COL_ANOM),
         legend = c(sprintf("Resto (n = %d)", sum(!an)), sprintf("Anómalos (n = %d)", sum(an))))
  legend(x = 0.45 * max(x), y = 0.42, bty = "n", cex = 0.8,
         legend = paste0("D = ", format(round(unname(ks2$statistic), 2), nsmall = 2),
                         "; p = ", format(round(ks2$p.value, 2), nsmall = 2)))
}

dibujar()                                                   # en pantalla
png(file.path(SALIDA, "C_distancia_fallas.png"), width = 2400, height = 1100, res = 300)
dibujar()
invisible(dev.off())

# ------------------------- IMPRESIÓN ----------------------------------
redondear <- function(v) {
  num <- vapply(v, is.numeric, logical(1))
  v[num] <- lapply(v[num], function(z) signif(z, 3))
  v
}
cat("\n================ Kolmogorov-Smirnov y prueba global por bandas ================\n")
print(redondear(ks_tab), row.names = FALSE)
cat("\n================ Proporción de anómalos por banda ================\n")
print(redondear(bandas), row.names = FALSE)
cat("\nArchivos guardados en:", SALIDA, "\n")
