# =============================================================================
# FASE 3 · SCRIPT 3 - PARTE A: COMPORTAMIENTO ESPERADO (OLS Y GWR MIXTA) EN R
# Proyecto de grado: comportamiento anómalo de movimientos en masa, Valle de Aburrá
# =============================================================================
#
# Lee la tabla común del script 2 (F3_S2_tabla_modelo.csv, eventos con EN_MODELO = 1)
# y, para cada dimensión del vector de comportamiento (log10 del área, de H/L y de la
# elongación), ajusta dos modelos del comportamiento esperado:
#
#   1. OLS: todos los coeficientes globales (línea base, decisión 12).
#   2. GWR mixta (gwr.mixed de GWmodel): intercepto, pendiente, distancia a la falla,
#      distancia al drenaje y densidad de lineamientos con coeficientes locales (un solo
#      ancho de banda); geología, geomorfología, cobertura, tipo de movimiento y tipo de
#      fuente con coeficientes globales.
#
# En cada dimensión se elige el modelo de menor AICc; si el otro queda a menos de 3
# unidades, se elige el de menos parámetros efectivos (decisión 12). De ese modelo
# salen los residuales estudentizados que pasan a la Parte B.
#
# La MGWR mixta no se corre en esta versión.
#
# QUÉ REPORTA, POR DIMENSIÓN
#   A. Comparación OLS / GWR mixta (AICc, parámetros efectivos, R²) y modelo elegido.
#   B. Coeficientes globales en unidades de log10 y como factor multiplicativo
#      (10^coeficiente), con su intervalo de confianza del 95 %.
#   C. Coeficientes locales: distribución, porcentaje de eventos donde son significativos
#      (positivos y negativos) y prueba Monte Carlo de variación espacial. Un mapa PNG
#      muestra dónde son significativos.
#   D. Residuales estudentizados: cuántos pasan de |2| y de |3|, y el I de Moran de los
#      residuales de la OLS y del modelo elegido.
#   E. Colinealidad local del bloque de variables continuas (número de condición y VIF
#      locales, gwr.collin.diagno de GWmodel).
#
# CÓMO CORRERLO
#   1. Instala R (4.2 o más reciente) y RStudio. Una sola vez, en la consola de R:
#        install.packages(c("GWmodel", "spdep"))
#   2. Abre este archivo en RStudio y edita las rutas del bloque CONFIGURACIÓN. En R las
#      rutas de Windows se escriben con barras normales: "D:/Users/..." (no "D:\Users\...").
#   3. Córrelo completo con el botón "Source" (o Ctrl + Shift + S).
#   Tarda del orden de 30 a 45 minutos con 999 permutaciones; la consola muestra el
#   avance. Con N_PERMUTACIONES <- 0 se omite la prueba Monte Carlo y tarda la mitad.
#
# SALIDAS (en CARPETA_SALIDA)
#   F3_S3_reporte.txt                 todo lo que se imprime en consola
#   F3_S3_resultados.csv              por evento: modelo elegido, predicción, residual,
#                                     palanca, residual estudentizado y coeficientes locales
#   F3_S3_comparacion_modelos.csv     AICc, parámetros efectivos y R² de cada modelo
#   F3_S3_coeficientes_globales.csv   coeficientes globales de los dos modelos
#   F3_S3_coeficientes_locales.csv    resumen de los coeficientes locales y prueba Monte Carlo
#   F3_S3_residuales_todos.csv        residuales estudentizados de los dos modelos
#   F3_S3_colinealidad_local.csv      número de condición y VIF locales por evento
#   F3_S3_busqueda_ancho.csv          AICc de cada ancho de banda evaluado
#   F3_S3_mapa_coeficientes_<DIM>.png mapas de los coeficientes locales
#   Para los mapas en ArcGIS Pro: XY Table To Point sobre F3_S3_resultados.csv (campos X, Y).
#
# NOTA TÉCNICA SOBRE EL KERNEL
#   gwr.mixed con kernel adaptativo recalcula el orden de las distancias dentro de cada
#   regresión local y tarda unos 8 s por ajuste, lo que llevaría la corrida a más de 10
#   horas. Por eso se le entrega una matriz de distancias en la que la columna de cada
#   evento está dividida por su ancho adaptativo (la distancia a su vecino número bw).
#   Con esa matriz, el kernel bisquare de ancho fijo 1 da exactamente los mismos pesos que
#   el adaptativo: (1 - (d/b)^2)^2 si d < b. Es un argumento documentado de gwr.mixed
#   (dMat), no una modificación del paquete. La sección de verificaciones del reporte
#   compara ambos ajustes y muestra que coinciden.
# =============================================================================

VERSION <- "F3_S3 (R) · versión 2026-10-03 (b)"

# =============================================================================
# CONFIGURACIÓN - editar solo las rutas
# =============================================================================

TABLA_MODELO   <- "D:/Users/tasal/Desktop/Tesis de grado/fase3/F3_S2_tabla_modelo.csv"
CARPETA_SALIDA <- "D:/Users/tasal/Desktop/Tesis de grado/fase3"

# Prueba Monte Carlo de variación espacial (decisión 11): 0 = no se corre.
N_PERMUTACIONES <- 999

# =============================================================================
# DE AQUÍ EN ADELANTE SON DECISIONES METODOLÓGICAS: no hace falta editarlas
# =============================================================================

RESPUESTAS <- list(                        # columna = (nombre corto, descripción)
  LOG_AREA  = c("AREA",  "log10 del área"),
  LOG_HL    = c("HL",    "log10 de H/L"),
  LOG_ELONG = c("ELONG", "log10 de la elongación")
)
# La densidad de lineamientos es la del radio DENSIDAD_RADIO_MODELO del script 2 (1 km).
# Si allí se cambia a 2 km, aquí la columna pasa a DENS_LIN_2000.
LOCALES <- list(
  PENDIENTE           = c("PEND",   "Pendiente media del entorno"),
  LOG_DIST_ESTRUCTURA = c("DFALLA", "log10(distancia a la falla + 1 m)"),
  LOG_DIST_DRENAJE    = c("DDREN",  "log10(distancia al drenaje + 1 m)"),
  DENS_LIN_1000       = c("DLIN",   "Densidad de lineamientos (km/km², 1 km)")
)
GLOBALES <- c("GEOLOGIA", "GEOMORFOLOGIA", "COBERTURA", "TIPO", "GRUPO")   # categóricas

UMBRAL_EMPATE_AICC  <- 3      # decisión 12
ALFA                <- 0.05
K_VECINOS_MORAN     <- 8      # decisión 19
PERMUTACIONES_MORAN <- 999
ANCHO_MINIMO        <- 20     # límites de búsqueda de bw.gwr en GWmodel: de 20 vecinos a n
UMBRAL_CN           <- 30     # número de condición local (Belsley, Kuh y Welsch 1980)
UMBRAL_VIF          <- 10
SEMILLA             <- 2026

DECISIONES_REGISTRADAS <- c(
  "Muestra: eventos con EN_MODELO = 1 en la tabla del script 2.",
  paste("Bloque local (coeficientes que varían en el espacio): intercepto, pendiente media del",
        "entorno, log10 de la distancia a la falla, log10 de la distancia al drenaje y densidad",
        "de lineamientos automáticos (km/km², ventana de 1 km)."),
  paste("Bloque global (un coeficiente para todo el valle): geología, geomorfología, cobertura, tipo",
        "de movimiento y tipo de fuente, como indicadoras; la referencia de cada una es su clase más",
        "frecuente. Las categóricas no pueden tener coeficiente local: dentro de un vecindario casi",
        "todos los eventos caen en la misma clase y el coeficiente no se puede estimar."),
  paste("Las respuestas y las variables continuas se estandarizan (media 0, desviación 1) para que",
        "los coeficientes locales sean comparables; los coeficientes globales se reportan en unidades",
        "de log10 y como factor multiplicativo (10^coeficiente)."),
  paste("GWR mixta con gwr.mixed de GWmodel (estimador de Speckman). Kernel bisquare adaptativo;",
        "ancho de banda (número de vecinos) por AICc, con la función gold de GWmodel entre 20",
        "vecinos y n, como bw.gwr. El AICc de cada ancho es el que calcula gwr.mixed."),
  paste("Modelo (decisión 12, modificada por indicación del asesor): se usa siempre la GWR mixta; la",
        "OLS se ajusta solo como comparación (AICc con la misma fórmula y 16 parámetros, R² e I de Moran)."),
  paste("Palanca, errores estándar y residuales estudentizados: como el ajuste es lineal en la",
        "respuesta, se obtienen aplicando gwr.mixed a respuestas indicadoras (1 en un evento, 0 en",
        "los demás), lo que da la matriz sombrero completa."),
  paste("Significancia de los coeficientes locales: t con el alfa ajustado de da Silva y",
        "Fotheringham (2016): alfa · k / ENP del bloque local, con k = 5 términos locales."),
  paste("Prueba Monte Carlo (decisión 11; Brunsdon, Fotheringham y Charlton 1998): se permutan las",
        "ubicaciones, se reajusta el modelo con el mismo ancho (como gwr.montecarlo de GWmodel) y",
        "se compara la desviación estándar de cada coeficiente local con la observada.",
        "p = (simulaciones con desviación >= observada + 1) / (permutaciones + 1). Las pruebas de",
        "las tres dimensiones se ajustan juntas por Holm."),
  paste("I de Moran de los residuales estudentizados: k = 8 vecinos más cercanos, pesos",
        "estandarizados por filas, 999 permutaciones (moran.mc de spdep)."),
  paste("Colinealidad local del bloque continuo: gwr.collin.diagno de GWmodel con el ancho de la",
        "GWR mixta de cada dimensión. Umbrales: número de condición > 30, VIF > 10.")
)

# =============================================================================
# PAQUETES
# =============================================================================

for (paquete in c("GWmodel", "spdep")) {
  if (!requireNamespace(paquete, quietly = TRUE)) {
    stop("Falta el paquete ", paquete, ". Instálalo una vez con:\n",
         "  install.packages(c(\"GWmodel\", \"spdep\"))", call. = FALSE)
  }
}
suppressPackageStartupMessages({
  library(GWmodel)   # carga también sp
  library(spdep)
})
if (getRversion() < "4.2.0") {
  warning("Tu R es anterior a 4.2: las tildes del reporte pueden verse mal en Windows.")
}
if (.Platform$OS.type == "unix") options(bitmapType = "cairo")
options(width = 200)

# =============================================================================
# REPORTE (consola + archivo de texto)
# =============================================================================

REPORTE <- character(0)

rep_ <- function(...) {
  texto <- paste0(...)
  cat(texto, "\n", sep = "")
  REPORTE <<- c(REPORTE, texto)
}
titulo <- function(texto, nivel = 1) {
  caracter <- if (nivel == 1) "=" else "-"
  rep_("")
  rep_(strrep(caracter, 78))
  rep_(texto)
  rep_(strrep(caracter, 78))
}
tabla <- function(df) {
  for (linea in capture.output(print(df, row.names = FALSE, right = TRUE))) rep_(linea)
}
aviso <- function(texto) rep_("  AVISO: ", texto)
guardar_reporte <- function() {
  if (!dir.exists(CARPETA_SALIDA)) return(invisible(NULL))
  con <- file(file.path(CARPETA_SALIDA, "F3_S3_reporte.txt"), open = "w", encoding = "UTF-8")
  writeLines(REPORTE, con)
  close(con)
}
escribir_csv <- function(df, nombre) {
  # UTF-8 con BOM, igual que los CSV de Python (pandas, encoding "utf-8-sig")
  con <- file(file.path(CARPETA_SALIDA, nombre), open = "w", encoding = "UTF-8")
  on.exit(close(con))
  cat("\ufeff", file = con)
  write.csv(df, con, row.names = FALSE, na = "")
}
f2 <- function(x, d = 2) formatC(x, format = "f", digits = d)

# =============================================================================
# PREPARACIÓN DE DATOS
# =============================================================================

estandarizar <- function(x) {
  media <- mean(x)
  de <- sqrt(mean((x - media)^2))       # desviación poblacional, como numpy
  list(z = (x - media) / de, media = media, de = de)
}

indicadoras <- function(df, columnas) {
  # Indicadoras 0/1 con la clase más frecuente como referencia.
  X <- NULL; etiquetas <- character(0); refs <- character(0)
  for (col in columnas) {
    v <- as.character(df[[col]])
    conteo <- sort(table(v), decreasing = TRUE)
    ref <- names(conteo)[1]
    refs[col] <- ref
    for (nivel in sort(setdiff(unique(v), ref), method = "radix")) {
      X <- cbind(X, as.numeric(v == nivel))
      etiquetas <- c(etiquetas, paste0(col, "=", nivel))
    }
  }
  colnames(X) <- paste0("G", seq_len(ncol(X)))
  list(X = X, etiquetas = etiquetas, refs = refs)
}

# =============================================================================
# GWR MIXTA CON gwr.mixed
# =============================================================================

aicc_formula <- function(rss, enp, n) {
  # La misma fórmula de gwr.mixed (AICc gaussiano)
  if (!is.finite(rss) || n - 2 - enp <= 0) return(Inf)
  n * log(rss / n) + n * log(2 * pi) + n * (n + enp) / (n - 2 - enp)
}

distancias_escaladas <- function(DM_ord, DM, bw) {
  # Columna i dividida por la distancia al vecino número bw del evento i (contándose a
  # sí mismo, como gw.weight de GWmodel). Con kernel bisquare de ancho fijo 1 da los
  # mismos pesos que el kernel adaptativo de ancho bw.
  sweep(DM, 2, DM_ord[bw, ], "/")
}

ajustar_mixta <- function(spd, DMs, diagnostico = FALSE) {
  # gwr.mixed imprime los coeficientes globales en consola: se silencia.
  invisible(capture.output(
    m <- gwr.mixed(FORMULA, data = spd, fixed.vars = NOMBRES_G, bw = 1,
                   diagnostic = diagnostico, kernel = "bisquare", adaptive = FALSE,
                   dMat = DMs)
  ))
  m
}

ajustados_mixta <- function(m, XL, XG) {
  as.vector(rowSums(XL * m$local) + XG %*% t(m$global))
}

# Los parámetros efectivos dependen solo del ancho (no de la respuesta): se guardan para
# no recalcularlos cuando otra dimensión evalúa el mismo ancho.
CACHE_ENP <- new.env()
CACHE_OPERADOR <- new.env()
VERIF <- new.env()

buscar_ancho <- function(spd, y, XL, XG, DM, DM_ord, dimension) {
  n <- length(y)
  registro <- list()
  evaluar <- function(bw) {
    bw <- as.integer(round(bw))
    clave <- as.character(bw)
    previo <- Filter(function(r) r$bw == bw, registro)
    if (length(previo) > 0) return(previo[[1]]$AICc)
    t0 <- Sys.time()
    DMs <- distancias_escaladas(DM_ord, DM, bw)
    res <- tryCatch({
      if (exists(clave, envir = CACHE_ENP, inherits = FALSE)) {
        m <- ajustar_mixta(spd, DMs, FALSE)
        enp <- get(clave, envir = CACHE_ENP)
        rss <- sum((y - ajustados_mixta(m, XL, XG))^2)
        list(aicc = aicc_formula(rss, enp, n), enp = enp)
      } else {
        m <- ajustar_mixta(spd, DMs, TRUE)
        assign(clave, m$df.used, envir = CACHE_ENP)
        if (is.null(VERIF$aicc)) {
          VERIF$aicc <- abs(aicc_formula(m$r.ss, m$df.used, n) - m$aic)
        }
        a <- if (n - 2 - m$df.used <= 0) Inf else m$aic
        list(aicc = a, enp = m$df.used)
      }
    }, error = function(e) list(aicc = Inf, enp = NA))
    seg <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("    ancho %3d vecinos: AICc %s | ENP %s | %.0f s\n", bw,
                if (is.finite(res$aicc)) f2(res$aicc) else "no estimable",
                if (is.na(res$enp)) "-" else f2(res$enp), seg))
    registro[[length(registro) + 1]] <<- list(bw = bw, AICc = res$aicc, ENP = res$enp)
    # gold necesita valores finitos para comparar
    if (is.finite(res$aicc)) res$aicc else 1e12
  }
  bw <- as.integer(round(gold(evaluar, ANCHO_MINIMO, n, adapt.bw = TRUE)))
  tabla_bw <- do.call(rbind, lapply(registro, as.data.frame))
  tabla_bw <- tabla_bw[order(tabla_bw$bw), ]
  tabla_bw$AICc[!is.finite(tabla_bw$AICc)] <- NA
  list(bw = bw, tabla = cbind(Dimension = dimension, tabla_bw))
}

operador_mixta <- function(spd, XL, XG, DMs, bw) {
  # Matriz sombrero H y operadores de los coeficientes: como gwr.mixed es lineal en la
  # respuesta, aplicarlo a la respuesta indicadora del evento k da la columna k de cada
  # operador. C[i, j, k]: peso de y_k en el coeficiente local j del evento i.
  clave <- as.character(bw)
  if (exists(clave, envir = CACHE_OPERADOR, inherits = FALSE)) return(get(clave, envir = CACHE_OPERADOR))
  n <- nrow(XL); p <- ncol(XL); g <- ncol(XG)
  C <- array(0, dim = c(n, p, n))
  A <- matrix(0, g, n)
  for (k in seq_len(n)) {
    e <- numeric(n); e[k] <- 1
    spd$Y_RESP <- e
    m <- ajustar_mixta(spd, DMs, FALSE)
    C[, , k] <- m$local
    A[, k] <- m$global
    if (k %% 100 == 0) cat(sprintf("    palanca: %d de %d eventos\n", k, n))
  }
  H <- XG %*% A
  ENP_j <- numeric(p)
  for (j in seq_len(p)) {
    H <- H + XL[, j] * C[, j, ]
    ENP_j[j] <- sum(XL[, j] * diag(C[, j, ]))
  }
  op <- list(C = C, A = A, H = H, ENP_j = ENP_j, ENP_global = sum(diag(XG %*% A)))
  assign(clave, op, envir = CACHE_OPERADOR)
  op
}

resultado_mixta <- function(op, y, XL, XG, bw) {
  n <- length(y); p <- ncol(XL)
  pred <- as.vector(op$H %*% y)
  e <- y - pred
  rss <- sum(e^2)
  tr <- sum(diag(op$H))
  s2 <- rss / (n - tr)
  h <- diag(op$H)
  B <- sapply(seq_len(p), function(j) as.vector(op$C[, j, ] %*% y))
  SE <- sapply(seq_len(p), function(j) sqrt(s2 * rowSums(op$C[, j, ]^2)))
  list(modelo = "GWR mixta", pred = pred, res = e, rss = rss, tr = tr, s2 = s2, h = h,
       aicc = aicc_formula(rss, tr, n), beta_local = B, se_local = SE,
       beta_global = as.vector(op$A %*% y), se_global = sqrt(s2 * rowSums(op$A^2)),
       ENP_j = op$ENP_j, bw = bw, re = e / sqrt(s2 * (1 - h)))
}

ajustar_ols <- function(y, XL, XG) {
  X <- cbind(XL, XG)
  m <- lm(y ~ X - 1)
  n <- length(y); p <- ncol(X)
  rss <- sum(residuals(m)^2)
  cf <- summary(m)$coefficients
  list(modelo = "OLS", pred = as.vector(fitted(m)), res = as.vector(residuals(m)), rss = rss,
       tr = p, s2 = rss / (n - p), h = as.vector(hatvalues(m)), aicc = aicc_formula(rss, p, n),
       beta = cf[, 1], se = cf[, 2], re = as.vector(rstandard(m)))
}

monte_carlo <- function(spd, DMs, sd_obs, n_perm) {
  # Brunsdon, Fotheringham y Charlton (1998), como gwr.montecarlo de GWmodel pero sobre
  # gwr.mixed: se reubican los eventos al azar y se conserva el ancho de banda.
  n <- nrow(DMs)
  mayores <- numeric(length(sd_obs))
  t0 <- Sys.time()
  for (k in seq_len(n_perm)) {
    o <- sample.int(n)
    m <- ajustar_mixta(spd, DMs[o, o], FALSE)
    mayores <- mayores + (apply(m$local, 2, sd) >= sd_obs)
    if (k %% 100 == 0) {
      seg <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      cat(sprintf("    Monte Carlo: %d de %d permutaciones (faltan unos %.0f min)\n",
                  k, n_perm, seg / k * (n_perm - k) / 60))
    }
  }
  (mayores + 1) / (n_perm + 1)
}

moran_residuales <- function(r, lw) {
  m <- moran.mc(r, lw, nsim = PERMUTACIONES_MORAN)
  c(I = unname(m$statistic), p = m$p.value)
}

t_critico_local <- function(ENP_j, n) {
  # da Silva y Fotheringham (2016), como en mgwr: alfa · k / ENP, gl = n - 1
  k <- length(ENP_j)
  qt(1 - (ALFA * k / sum(ENP_j)) / 2, n - 1)
}

tabla_globales <- function(beta, se, etiquetas, sd_y, gl) {
  coef <- beta * sd_y
  ee <- se * sd_y
  t <- beta / se
  tc <- qt(0.975, gl)
  data.frame("Término" = etiquetas,
             "Coef. (log10)" = round(coef, 4), "EE" = round(ee, 4),
             "t" = round(t, 2), "p" = round(2 * pt(-abs(t), gl), 4),
             "Factor 10^coef" = round(10^coef, 3),
             "IC 95 % inf." = round(10^(coef - tc * ee), 3),
             "IC 95 % sup." = round(10^(coef + tc * ee), 3),
             check.names = FALSE)
}

mapa_coeficientes <- function(ruta, coords, B, SIG, nombres, titulo_dim, p_mc) {
  png(ruta, width = 2100, height = 1400, res = 150)
  on.exit(dev.off())
  par(mfrow = c(2, 3), mar = c(1, 1, 3.2, 1), oma = c(0, 0, 2, 0))
  paleta <- hcl.colors(11, "Blue-Red 3")
  for (j in seq_len(ncol(B))) {
    b <- B[, j]
    lim <- max(abs(b))
    color <- paleta[findInterval(b, seq(-lim, lim, length.out = 12), all.inside = TRUE)]
    sig <- SIG[, j] == 1
    sub <- if (is.null(p_mc)) "" else sprintf("p Monte Carlo = %.3f · ", p_mc[j])
    plot(coords, type = "n", asp = 1, axes = FALSE, xlab = "", ylab = "",
         main = sprintf("%s\n%ssignif. en %.0f %% · rango %s a %s", nombres[j], sub,
                        100 * mean(sig), f2(min(b)), f2(max(b))),
         cex.main = 0.95)
    points(coords[!sig, , drop = FALSE], pch = 16, col = "grey75", cex = 0.55)
    points(coords[sig, , drop = FALSE], pch = 21, bg = color[sig], col = "grey15", cex = 1.05)
    box(col = "grey60")
  }
  plot.new()
  legend("center", bty = "n", cex = 1,
         legend = c("Significativo, positivo", "Significativo, negativo", "No significativo"),
         pch = c(21, 21, 16), pt.bg = c(paleta[10], paleta[2], NA),
         col = c("grey15", "grey15", "grey75"), pt.cex = c(1.4, 1.4, 0.9))
  text(0.5, 0.18, "Coeficientes estandarizados.\nEl color se satura hacia el valor\nextremo de cada panel.", cex = 0.85)
  mtext(titulo_dim, outer = TRUE, cex = 1.1, font = 2)
}

# =============================================================================
# PROGRAMA PRINCIPAL
# =============================================================================

main <- function() {
  inicio <- Sys.time()
  if (!dir.exists(CARPETA_SALIDA)) dir.create(CARPETA_SALIDA, recursive = TRUE)
  set.seed(SEMILLA)

  titulo("FASE 3 · SCRIPT 3 - PARTE A: COMPORTAMIENTO ESPERADO (OLS Y GWR MIXTA, R)")
  rep_("Ejecución: ", format(inicio, "%Y-%m-%d %H:%M"), " | ", VERSION)
  rep_(R.version.string, " | GWmodel ", as.character(packageVersion("GWmodel")),
       " | spdep ", as.character(packageVersion("spdep")))
  rep_("")
  rep_("Decisiones registradas:")
  for (k in seq_along(DECISIONES_REGISTRADAS)) rep_("  ", k, ". ", DECISIONES_REGISTRADAS[k])

  # ---------------------------------------------------------------- datos
  if (!file.exists(TABLA_MODELO)) stop("No encuentro la tabla del script 2: ", TABLA_MODELO)
  tabla_in <- read.csv(TABLA_MODELO, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE,
                       colClasses = c(ID = "character"))
  # Los nombres de las clases llevan tildes: se marcan como UTF-8 para que el orden y las
  # comparaciones no dependan de la configuración regional del computador.
  for (cc in names(tabla_in)) {
    if (is.character(tabla_in[[cc]])) tabla_in[[cc]] <- enc2utf8(tabla_in[[cc]])
  }
  necesarias <- c("ID", "X", "Y", "EN_MODELO", names(RESPUESTAS), names(LOCALES), GLOBALES)
  faltan <- setdiff(necesarias, names(tabla_in))
  if (length(faltan) > 0) stop("A la tabla del script 2 le faltan las columnas: ",
                               paste(faltan, collapse = ", "))
  df <- tabla_in[tabla_in$EN_MODELO == 1, ]
  rownames(df) <- NULL
  vacios <- !complete.cases(df[, c(names(RESPUESTAS), names(LOCALES))]) |
    apply(df[, GLOBALES], 1, function(v) any(is.na(v) | v == ""))
  if (any(vacios)) stop("Hay ", sum(vacios), " eventos del modelo con datos vacíos: ",
                        paste(df$ID[vacios], collapse = ", "))
  n <- nrow(df)
  coords <- as.matrix(df[, c("X", "Y")])

  cortos_loc <- c("INT", sapply(LOCALES, `[`, 1))
  nombres_loc <- c("Intercepto", sapply(LOCALES, `[`, 2))
  est <- lapply(names(LOCALES), function(col) estandarizar(df[[col]]))
  XL <- cbind(1, sapply(est, `[[`, "z"))
  colnames(XL) <- c("Intercept", cortos_loc[-1])
  ind <- indicadoras(df, GLOBALES)
  XG <- ind$X
  assign("NOMBRES_G", colnames(XG), envir = globalenv())
  assign("FORMULA", as.formula(paste("Y_RESP ~", paste(c(cortos_loc[-1], colnames(XG)),
                                                         collapse = " + "))), envir = globalenv())

  titulo("DATOS")
  grupos <- table(df$GRUPO)
  rep_("Eventos en el modelo: ", n, " (", paste(names(grupos), grupos, collapse = ", "), ")")
  rep_("Bloque local: ", ncol(XL), " términos | bloque global: ", ncol(XG),
       " indicadoras | parámetros nominales: ", ncol(XL) + ncol(XG))
  rep_("Referencias de las categóricas: ",
       paste(names(ind$refs), ind$refs, sep = ": ", collapse = "; "))
  rep_("")
  rep_("Estandarización de las variables continuas locales:")
  tabla(data.frame(Variable = names(LOCALES),
                   Media = round(sapply(est, `[[`, "media"), 4),
                   "Desv. est." = round(sapply(est, `[[`, "de"), 4), check.names = FALSE))

  DM <- gw.dist(dp.locat = coords)
  DM_ord <- apply(DM, 2, sort)
  datos <- data.frame(Y_RESP = 0, XL[, -1, drop = FALSE], XG)
  spd <- SpatialPointsDataFrame(coords, datos, match.ID = FALSE)
  lw <- nb2listw(knn2nb(knearneigh(coords, k = K_VECINOS_MORAN)), style = "W")

  resultados <- df[, c("ID", "X", "Y", "GRUPO", "TIPO")]
  residuales_todos <- df[, "ID", drop = FALSE]
  colinealidad <- df[, "ID", drop = FALSE]
  filas_comp <- list(); filas_glob <- list(); filas_loc <- list(); tablas_bw <- list()
  verif_kernel <- NULL; verif_oper <- list()

  # ---------------------------------------------------------------- una dimensión a la vez
  for (col in names(RESPUESTAS)) {
    corto <- RESPUESTAS[[col]][1]
    desc <- RESPUESTAS[[col]][2]
    titulo(sprintf("DIMENSIÓN: %s (%s)", desc, col))
    ey <- estandarizar(df[[col]])
    y <- ey$z
    sd_y <- ey$de
    spd$Y_RESP <- y
    rep_(sprintf("Media %.4f | desviación %.4f (la respuesta se estandariza para ajustar)",
                 ey$media, sd_y))

    # OLS
    ols <- ajustar_ols(y, XL, XG)

    # GWR mixta: ancho por AICc
    rep_("")
    rep_("Búsqueda del ancho de banda por AICc (cada ancho nuevo tarda unos 30 s):")
    busq <- buscar_ancho(spd, y, XL, XG, DM, DM_ord, col)
    bw <- busq$bw
    tablas_bw[[col]] <- busq$tabla
    DMs <- distancias_escaladas(DM_ord, DM, bw)
    rep_(sprintf("Ancho elegido: %d vecinos (%d anchos evaluados).", bw, nrow(busq$tabla)))
    if (bw >= 0.9 * n) {
      aviso("el ancho quedó cerca del total de eventos: la GWR mixta es casi un modelo global.")
    } else if (bw <= ANCHO_MINIMO + 2) {
      aviso("el ancho quedó en el límite inferior de la búsqueda: revisa F3_S3_busqueda_ancho.csv.")
    }

    rep_("Palanca y errores estándar (gwr.mixed aplicado a respuestas indicadoras)...")
    op <- operador_mixta(spd, XL, XG, DMs, bw)
    spd$Y_RESP <- y
    gwr <- resultado_mixta(op, y, XL, XG, bw)

    # Verificaciones internas
    m_directo <- ajustar_mixta(spd, DMs, FALSE)
    verif_oper[[col]] <- c(
      pred = max(abs(ajustados_mixta(m_directo, XL, XG) - gwr$pred)),
      enp = abs(get(as.character(bw), envir = CACHE_ENP) - gwr$tr))
    if (is.null(verif_kernel)) {
      cat("    verificando el kernel adaptativo original (unos segundos)...\n")
      invisible(capture.output(
        m_adapt <- gwr.mixed(FORMULA, data = spd, fixed.vars = NOMBRES_G, bw = bw,
                             diagnostic = FALSE, kernel = "bisquare", adaptive = TRUE, dMat = DM)))
      verif_kernel <- max(abs(m_adapt$local - m_directo$local), abs(m_adapt$global - m_directo$global))
    }

    # A. Comparación
    modelos <- list(OLS = ols, "GWR mixta" = gwr)
    tss <- sum((y - mean(y))^2)
    comp <- do.call(rbind, lapply(names(modelos), function(nom) {
      r <- modelos[[nom]]
      r2 <- 1 - r$rss / tss
      data.frame(Modelo = nom, AICc = round(r$aicc, 2), ENP = round(r$tr, 2),
                 "R²" = round(r2, 3), "R² ajust." = round(1 - (1 - r2) * (n - 1) / (n - r$tr - 1), 3),
                 "Sigma (estand.)" = round(sqrt(r$s2), 3),
                 "Ancho (vecinos)" = if (nom == "OLS") "-" else as.character(bw),
                 check.names = FALSE)
    }))
    comp <- data.frame(comp[, 1:2], "Dif. AICc" = round(comp$AICc - min(comp$AICc), 2), comp[, -(1:2)],
                  check.names = FALSE)
    # Por indicación del asesor se usa siempre la GWR mixta; la OLS queda solo como comparación
    elegido <- "GWR mixta"
    dif_aicc <- comp$AICc[comp$Modelo == "OLS"] - comp$AICc[comp$Modelo == "GWR mixta"]
    motivo <- sprintf("por indicación del asesor se usa siempre la GWR mixta; el AICc de la OLS es %.2f unidades %s",
                      abs(dif_aicc), if (dif_aicc >= 0) "mayor" else "menor")
    titulo("A. Comparación de modelos", 2)
    tabla(comp)
    rep_(sprintf("Modelo elegido: %s (%s).", elegido, motivo))
    filas_comp[[col]] <- data.frame("Dimensión" = col, comp, Elegido = comp$Modelo == elegido,
                               check.names = FALSE)
    res <- modelos[[elegido]]
    resultados[[paste0("MODELO_", corto)]] <- elegido
    resultados[[paste0("PRED_", corto)]] <- res$pred * sd_y + ey$media
    resultados[[paste0("RES_", corto)]] <- res$res * sd_y
    resultados[[paste0("H_", corto)]] <- res$h
    resultados[[paste0("RE_", corto)]] <- res$re

    # B. Coeficientes globales
    k_loc <- ncol(XL)
    tg_ols <- tabla_globales(ols$beta[-(1:k_loc)], ols$se[-(1:k_loc)], ind$etiquetas, sd_y, n - ols$tr)
    tg_gwr <- tabla_globales(gwr$beta_global, gwr$se_global, ind$etiquetas, sd_y, n - gwr$tr)
    filas_glob[[col]] <- rbind(data.frame("Dimensión" = col, Modelo = "OLS", tg_ols, check.names = FALSE),
                               data.frame("Dimensión" = col, Modelo = "GWR mixta", tg_gwr, check.names = FALSE))
    titulo(sprintf("B. Coeficientes globales (%s)", elegido), 2)
    tg <- if (elegido == "OLS") tg_ols else tg_gwr
    tabla(tg)
    rep_("Coeficiente en unidades de log10 de la respuesta. El factor 10^coef es multiplicativo")
    rep_("respecto a la clase de referencia, con lo demás constante: 1,5 = 50 % mayor; 0,8 = 20 % menor.")
    sig_g <- tg[tg$p < ALFA, ]
    if (nrow(sig_g) > 0) {
      for (i in seq_len(nrow(sig_g))) {
        rep_(sprintf("  · %s: factor %s (IC 95 %%: %s a %s), p = %s", sig_g[i, 1],
                     f2(sig_g[i, "Factor 10^coef"]), f2(sig_g[i, "IC 95 % inf."]),
                     f2(sig_g[i, "IC 95 % sup."]), f2(sig_g$p[i], 4)))
      }
    } else {
      rep_("  Ningún coeficiente global es significativo con p < 0,05.")
    }

    # C. Coeficientes locales
    p_mc <- NULL
    if (elegido == "OLS") {
      titulo("C. Coeficientes de las variables continuas (OLS, estandarizados)", 2)
      cf <- ols$beta[1:k_loc]; se <- ols$se[1:k_loc]
      tabla(data.frame("Término" = nombres_loc, "Coef." = round(cf, 4), EE = round(se, 4),
                       t = round(cf / se, 2), p = round(2 * pt(-abs(cf / se), n - ols$tr), 4),
                       check.names = FALSE))
      rep_("La OLS fue elegida: los efectos de las variables continuas se leen como constantes en el valle.")
    } else {
      tcrit <- t_critico_local(gwr$ENP_j, n)
      tvals <- gwr$beta_local / gwr$se_local
      SIG <- (abs(tvals) > tcrit) * 1
      if (N_PERMUTACIONES > 0) {
        rep_(sprintf("Prueba Monte Carlo con %d permutaciones (la consola muestra el avance)...",
                     N_PERMUTACIONES))
        p_mc <- monte_carlo(spd, DMs, apply(gwr$beta_local, 2, sd), N_PERMUTACIONES)
      }
      titulo(sprintf("C. Coeficientes locales (%s, estandarizados)", elegido), 2)
      b <- gwr$beta_local
      tl <- data.frame("Término" = nombres_loc,
                       "Mín" = round(apply(b, 2, min), 3),
                       "P25" = round(apply(b, 2, quantile, 0.25), 3),
                       "Mediana" = round(apply(b, 2, median), 3),
                       "P75" = round(apply(b, 2, quantile, 0.75), 3),
                       "Máx" = round(apply(b, 2, max), 3),
                       "ENP_j" = round(gwr$ENP_j, 2),
                       "% signif. (+)" = round(100 * colMeans(SIG * (b > 0)), 1),
                       "% signif. (-)" = round(100 * colMeans(SIG * (b < 0)), 1),
                       check.names = FALSE)
      if (!is.null(p_mc)) tl[["p Monte Carlo"]] <- round(p_mc, 3)
      tabla(tl)
      rep_(sprintf("Ancho: %d vecinos. Significancia con t > %.2f (alfa ajustado %.4f; da Silva y Fotheringham 2016).",
                   bw, tcrit, ALFA * k_loc / sum(gwr$ENP_j)))
      for (j in seq_len(k_loc)) {
        txt_mc <- if (is.null(p_mc)) "" else if (p_mc[j] < ALFA)
          sprintf("; su variación en el espacio supera la esperada por azar (p = %.3f)", p_mc[j]) else
          sprintf("; su variación en el espacio no se distingue del azar (p = %.3f)", p_mc[j])
        rep_(sprintf("  · %s: significativo en %.0f %% de los eventos (%.0f %% positivo, %.0f %% negativo)%s.",
                     nombres_loc[j], 100 * mean(SIG[, j]), tl[j, "% signif. (+)"],
                     tl[j, "% signif. (-)"], txt_mc))
      }
      rep_("La decisión final sobre la variación espacial se toma con la p ajustada por Holm (RESUMEN).")
      png_nombre <- sprintf("F3_S3_mapa_coeficientes_%s.png", corto)
      mapa_coeficientes(file.path(CARPETA_SALIDA, png_nombre), coords, b, SIG, nombres_loc,
                        sprintf("%s · GWR mixta, %d vecinos", desc, bw), p_mc)
      rep_("Mapa de dónde son significativos: ", png_nombre)
      filas_loc[[col]] <- data.frame("Dimensión" = col, tl, "Ancho (vecinos)" = bw, check.names = FALSE)
      for (j in seq_len(k_loc)) {
        resultados[[sprintf("B_%s_%s", corto, cortos_loc[j])]] <- b[, j]
        resultados[[sprintf("T_%s_%s", corto, cortos_loc[j])]] <- tvals[, j]
        resultados[[sprintf("SIG_%s_%s", corto, cortos_loc[j])]] <- SIG[, j]
      }
    }

    # D. Residuales estudentizados y Moran
    titulo(sprintf("D. Residuales estudentizados (%s)", elegido), 2)
    conteo <- function(r) data.frame(
      "|r| > 2" = sum(abs(r) > 2), "% > 2" = round(100 * mean(abs(r) > 2), 1),
      "|r| > 3" = sum(abs(r) > 3), "% > 3" = round(100 * mean(abs(r) > 3), 1),
      "Máx |r|" = round(max(abs(r)), 2), check.names = FALSE)
    mor_ols <- moran_residuales(ols$re, lw)
    mor_gwr <- moran_residuales(gwr$re, lw)
    td <- rbind(data.frame(Modelo = "OLS", conteo(ols$re), "Palanca máx." = round(max(ols$h), 3),
                      "I de Moran" = round(mor_ols["I"], 3), "p Moran" = round(mor_ols["p"], 3),
                      check.names = FALSE),
                data.frame(Modelo = "GWR mixta", conteo(gwr$re), "Palanca máx." = round(max(gwr$h), 3),
                      "I de Moran" = round(mor_gwr["I"], 3), "p Moran" = round(mor_gwr["p"], 3),
                      check.names = FALSE))
    td$Modelo[td$Modelo == elegido] <- paste0(elegido, " (elegido)")
    tabla(td)
    rep_(sprintf("Esperado si los residuales fueran normales: %.1f eventos con |r| > 2 y %.1f con |r| > 3.",
                 n * 2 * pnorm(-2), n * 2 * pnorm(-3)))
    extremos <- order(-abs(res$re))
    extremos <- extremos[abs(res$re[extremos]) > 3]
    if (length(extremos) > 0) {
      rep_("Eventos con |r| > 3 en el modelo elegido:")
      tabla(data.frame(ID = df$ID[extremos], GRUPO = df$GRUPO[extremos], TIPO = df$TIPO[extremos],
                       r = round(res$re[extremos], 2), Palanca = round(res$h[extremos], 3)))
    }
    rep_(sprintf("I de Moran (k = %d, %d permutaciones): un valor positivo y significativo en la OLS",
                 K_VECINOS_MORAN, PERMUTACIONES_MORAN))
    rep_("indica que el comportamiento esperado cambia en el espacio; si baja en la GWR mixta, el")
    rep_("modelo local recoge esa variación y los residuales quedan listos para la Parte B.")

    # E. Colinealidad local del bloque continuo
    titulo(sprintf("E. Colinealidad local del bloque continuo (ancho de la GWR mixta: %d vecinos)", bw), 2)
    cd <- gwr.collin.diagno(as.formula(paste("Y_RESP ~", paste(cortos_loc[-1], collapse = " + "))),
                            data = spd, bw = 1, kernel = "bisquare", adaptive = FALSE, dMat = DMs)
    cdf <- as.data.frame(cd$SDF)
    cn <- cdf$local_CN
    vifs <- cdf[, paste0(cortos_loc[-1], "_VIF"), drop = FALSE]
    te <- rbind(
      data.frame(Medida = "Número de condición local", Mediana = round(median(cn), 2),
                 P90 = round(quantile(cn, 0.9), 2), "Máx" = round(max(cn), 2),
                 "Eventos > umbral" = sum(cn > UMBRAL_CN), Umbral = UMBRAL_CN, check.names = FALSE),
      do.call(rbind, lapply(seq_along(vifs), function(j) data.frame(
        Medida = paste("VIF local:", nombres_loc[j + 1]), Mediana = round(median(vifs[[j]]), 2),
        P90 = round(quantile(vifs[[j]], 0.9), 2), "Máx" = round(max(vifs[[j]]), 2),
        "Eventos > umbral" = sum(vifs[[j]] > UMBRAL_VIF), Umbral = UMBRAL_VIF, check.names = FALSE))))
    tabla(te)
    if (sum(cn > UMBRAL_CN) == 0 && all(sapply(vifs, max) <= UMBRAL_VIF)) {
      rep_("Sin colinealidad local preocupante: los coeficientes locales se pueden interpretar.")
    } else {
      aviso(sprintf("%d eventos superan el número de condición %d. Sus coeficientes locales son inestables;",
                    sum(cn > UMBRAL_CN), UMBRAL_CN))
      rep_("  conviene marcarlos en los mapas y no interpretar allí los coeficientes uno por uno.")
    }
    colinealidad[[paste0("CN_", corto)]] <- cn
    for (j in seq_along(vifs)) colinealidad[[sprintf("VIF_%s_%s", corto, cortos_loc[j + 1])]] <- vifs[[j]]

    residuales_todos[[paste0("RE_", corto, "_OLS")]] <- ols$re
    residuales_todos[[paste0("RE_", corto, "_GWR")]] <- gwr$re
  }

  # ---------------------------------------------------------------- Monte Carlo conjunto
  loc_todos <- do.call(rbind, filas_loc)
  titulo("RESUMEN")
  if (!is.null(loc_todos) && "p Monte Carlo" %in% names(loc_todos)) {
    loc_todos[["p Holm"]] <- round(p.adjust(loc_todos[["p Monte Carlo"]], method = "holm"), 3)
    loc_todos[["Varía en el espacio"]] <- ifelse(loc_todos[["p Holm"]] < ALFA, "sí", "no")
    rep_(sprintf("Prueba Monte Carlo de las %d pruebas, ajustada por Holm:", nrow(loc_todos)))
    tabla(loc_todos[, c("Dimensión", "Término", "% signif. (+)", "% signif. (-)",
                        "p Monte Carlo", "p Holm", "Varía en el espacio")])
    rep_(sprintf("Con %d permutaciones, la p más pequeña posible es %.4f.", N_PERMUTACIONES,
                 1 / (N_PERMUTACIONES + 1)))
    rep_("")
  }
  resumen <- data.frame("Dimensión" = names(RESPUESTAS),
                        "Modelo elegido" = sapply(names(RESPUESTAS), function(c)
                          resultados[[paste0("MODELO_", RESPUESTAS[[c]][1])]][1]),
                        "Columna para la Parte B" = paste0("RE_", sapply(RESPUESTAS, `[`, 1)),
                        check.names = FALSE)
  tabla(resumen)

  titulo("VERIFICACIONES INTERNAS (deben ser del orden de 1e-10 o menores)", 2)
  rep_(sprintf("AICc: fórmula del script frente al de gwr.mixed: %.1e", VERIF$aicc))
  rep_(sprintf("Kernel: distancias escaladas frente al kernel adaptativo de gwr.mixed: %.1e", verif_kernel))
  for (col in names(verif_oper)) {
    rep_(sprintf("%s: predicción por matriz sombrero frente a gwr.mixed: %.1e | traza frente a gwr.mixed: %.1e",
                 col, verif_oper[[col]]["pred"], verif_oper[[col]]["enp"]))
  }

  # ---------------------------------------------------------------- salidas
  escribir_csv(resultados, "F3_S3_resultados.csv")
  escribir_csv(do.call(rbind, filas_comp), "F3_S3_comparacion_modelos.csv")
  escribir_csv(do.call(rbind, filas_glob), "F3_S3_coeficientes_globales.csv")
  if (!is.null(loc_todos)) escribir_csv(loc_todos, "F3_S3_coeficientes_locales.csv")
  escribir_csv(residuales_todos, "F3_S3_residuales_todos.csv")
  escribir_csv(colinealidad, "F3_S3_colinealidad_local.csv")
  escribir_csv(do.call(rbind, tablas_bw), "F3_S3_busqueda_ancho.csv")
  titulo("SALIDAS", 2)
  rep_(sprintf("Resultados por evento: F3_S3_resultados.csv (%d eventos)", n))
  rep_("Otras: F3_S3_comparacion_modelos.csv, F3_S3_coeficientes_globales.csv, F3_S3_coeficientes_locales.csv,")
  rep_("       F3_S3_residuales_todos.csv, F3_S3_colinealidad_local.csv, F3_S3_busqueda_ancho.csv, mapas PNG")
  rep_("Todo en: ", CARPETA_SALIDA)
  rep_(sprintf("Tiempo total: %.1f min", as.numeric(difftime(Sys.time(), inicio, units = "mins"))))
}

tryCatch(main(),
         error = function(e) {
           rep_("")
           rep_("ERROR: ", conditionMessage(e))
           stop(e)
         },
         finally = guardar_reporte())
