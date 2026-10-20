# =============================================================================
# FASE 3 · SCRIPT 5 - PARTE C: ESTRUCTURA LOCAL DE LAS ANOMALIAS (LISA) E INTEGRACION
# Proyecto de grado: comportamiento anomalo de movimientos en masa, Valle de Aburra
# =============================================================================
#
# Lee la salida de la Parte B (F3_S4_resultados.csv) y pregunta DONDE se concentran las
# anomalias: el I de Moran local (LISA, Anselin 1995) de la distancia robusta de cada
# evento frente a la de sus 8 vecinos mas cercanos.
#
#   Alto-Alto: evento anomalo rodeado de eventos anomalos -> zona que se aparta de lo esperado
#   Alto-Bajo: evento anomalo rodeado de eventos normales -> anomalia aislada
#   Bajo-Alto: evento normal rodeado de anomalos
#   Bajo-Bajo: evento normal rodeado de normales -> zona de comportamiento tipico
#
# Al final integra las tres partes en una tabla por evento (y una capa para ArcGIS Pro):
# anomalia contextual, perfil, rareza en el inventario y contexto local.
#
# COMO CORRERLO
#   En RStudio: editar las rutas de CONFIGURACION y correr con Source (Ctrl + Shift + S).
#   Usa spdep y sf, ya instalados con GWmodel. Tarda menos de un minuto.
#
# SALIDAS (en CARPETA_SALIDA)
#   F3_S5_reporte.txt              todo lo que se imprime en consola
#   F3_S5_resultados.csv           LISA por evento
#   F3_fase3_integrado.csv         tabla final de la Fase 3 (Partes A, B y C)
#   F3_fase3_integrado.gpkg        la misma tabla como capa de puntos para ArcGIS Pro
#   F3_S5_mapa_lisa.png, F3_S5_dispersion_moran.png
# =============================================================================

VERSION <- "F3_S5 (R) - version 2026-10-03 (c)"

# =============================================================================
# CONFIGURACION - editar solo las rutas
# =============================================================================

RESULTADOS_B   <- "D:/Users/tasal/Desktop/Tesis de grado/fase3/F3_S4_resultados.csv"
CARPETA_SALIDA <- "D:/Users/tasal/Desktop/Tesis de grado/fase3"

# Sistema de coordenadas de X e Y (el de las coronas): MAGNA-SIRGAS 2018 / Origen-Nacional
CRS_EPSG <- 9377

# Reactivaciones (decision 21): campo de la capa de poligonos que las marca. Si se deja
# CAMPO_REACTIVACION <- "", la sensibilidad sin reactivaciones no se corre.
# El nombre puede estar recortado en la geodatabase (p. ej. "Reactivaci"): se busca por su
# comienzo. Si el campo tiene un dominio codificado, el script traduce cada codigo a su
# descripcion (la que se ve en ArcGIS) antes de compararla con VALOR_REACTIVACION.
POLIGONOS_GDB      <- "C:/Users/villa/Documents/ArcGIS/Projects/GDB_proyecto de grado/GDB_proyecto de grado.gdb"
POLIGONOS_CAPA     <- "Inventario_integrado"
CAMPO_ID_POLIGONOS <- "ID_evento"
CAMPO_REACTIVACION <- "Reactivación"
VALOR_REACTIVACION <- "si"          # descripcion que marca una reactivacion (sin importar tildes)

# =============================================================================
# DECISIONES METODOLOGICAS (no hace falta editarlas)
# =============================================================================

K_VECINOS      <- 8          # decision 19
PERMUTACIONES  <- 9999       # con FDR sobre 430 pruebas, 999 no alcanzan (p minima 0,002)
ALFA           <- 0.05
SEMILLA        <- 2026

DECISIONES_REGISTRADAS <- c(
  "Variable (decision 18): la distancia robusta contextual de la Parte B, sin elevar al cuadrado para reducir la asimetria.",
  "Vecindad (decision 19): los 8 vecinos mas cercanos de cada corona, con pesos estandarizados por filas.",
  "Inferencia: permutaciones condicionales (Anselin 1995), 9999 por evento, con p bilateral. Con FDR sobre 430 pruebas, la p mas pequena debe bajar de 0,05/430 = 0,00012: con 999 permutaciones la p minima es 0,002 y nada podria resultar significativo; con 9999 es 0,0002.",
  "Comparaciones multiples (decision 20): FDR de Benjamini y Hochberg (1995) sobre todos los eventos.",
  "Cuadrantes por la media de la distancia y de su promedio en los vecinos (Alto-Alto, Alto-Bajo, Bajo-Alto, Bajo-Bajo).",
  "Sensibilidades: sin los eventos con errores conocidos (su distancia extrema pesa en el LISA de sus vecinos) y, si hay campo, sin reactivaciones (decision 21)."
)

# =============================================================================
# PAQUETES Y REPORTE
# =============================================================================

for (paq in c("spdep", "sf")) {
  if (!requireNamespace(paq, quietly = TRUE)) stop("Falta el paquete ", paq, ". Instalalo con install.packages(\"", paq, "\")", call. = FALSE)
}
suppressPackageStartupMessages({ library(spdep); library(sf) })
if (.Platform$OS.type == "unix") options(bitmapType = "cairo")
options(width = 200)

REPORTE <- character(0)
rep_ <- function(...) { texto <- paste0(...); cat(texto, "\n", sep = ""); REPORTE <<- c(REPORTE, texto) }
titulo <- function(texto, nivel = 1) {
  car <- if (nivel == 1) "=" else "-"
  rep_(""); rep_(strrep(car, 78)); rep_(texto); rep_(strrep(car, 78))
}
tabla <- function(df) for (l in capture.output(print(df, row.names = FALSE, right = TRUE))) rep_(l)
aviso <- function(texto) rep_("  AVISO: ", texto)
guardar_reporte <- function() {
  if (!dir.exists(CARPETA_SALIDA)) return(invisible(NULL))
  con <- file(file.path(CARPETA_SALIDA, "F3_S5_reporte.txt"), open = "w", encoding = "UTF-8")
  writeLines(REPORTE, con); close(con)
}
escribir_csv <- function(df, nombre) {
  con <- file(file.path(CARPETA_SALIDA, nombre), open = "w", encoding = "UTF-8")
  on.exit(close(con)); cat("\ufeff", file = con)
  write.csv(df, con, row.names = FALSE, na = "")
}
f <- function(x, d = 3) formatC(x, format = "f", digits = d)

CATS <- c("Alto-Alto", "Alto-Bajo", "Bajo-Alto", "Bajo-Bajo", "No significativo")
COLS <- c("Alto-Alto" = "#B04A2E", "Alto-Bajo" = "#E0A15E", "Bajo-Alto" = "#8DB3D6",
          "Bajo-Bajo" = "#2D5F8B", "No significativo" = "#CFCBC2")
C_OSC <- "#1F2A33"

# =============================================================================
# FUNCIONES
# =============================================================================

pesos <- function(xy) nb2listw(knn2nb(knearneigh(xy, k = K_VECINOS)), style = "W")

lisa <- function(x, lw) {
  set.seed(SEMILLA)
  r <- localmoran_perm(x, lw, nsim = PERMUTACIONES, alternative = "two.sided", iseed = SEMILLA)
  col_p <- grep("Sim$", colnames(r), value = TRUE)
  col_p <- col_p[!grepl("folded", col_p)][1]           # p bilateral por rangos de la permutacion
  p <- as.numeric(r[, col_p])
  q <- as.character(attr(r, "quadr")$mean)
  q <- c("High-High" = "Alto-Alto", "High-Low" = "Alto-Bajo", "Low-High" = "Bajo-Alto",
         "Low-Low" = "Bajo-Bajo")[q]
  p_fdr <- p.adjust(p, method = "BH")
  data.frame(Ii = as.numeric(r[, "Ii"]), Z = as.numeric(r[, "Z.Ii"]), P = p, P_FDR = p_fdr,
             CUADRANTE = unname(q),
             LISA = ifelse(p_fdr < ALFA, unname(q), "No significativo"),
             LISA_P05 = ifelse(p < ALFA, unname(q), "No significativo"), stringsAsFactors = FALSE)
}

conteo <- function(v) as.integer(table(factor(v, levels = CATS)))

sin_tildes <- function(x) chartr("\u00e1\u00e9\u00ed\u00f3\u00fa\u00fc", "aeiouu", tolower(trimws(x)))

leer_reactivaciones <- function(ids) {
  # Lee el campo de reactivaciones de la capa de poligonos. GDAL (-resolveDomains) agrega una
  # columna <campo>_resolved con la descripcion del dominio de cada codigo; si el campo no tiene
  # dominio, se usa el valor guardado.
  if (!dir.exists(POLIGONOS_GDB)) stop("No encuentro la geodatabase: ", POLIGONOS_GDB)
  tmp <- tempfile(fileext = ".csv")
  traducido <- tryCatch({
    gdal_utils("vectortranslate", POLIGONOS_GDB, tmp, options = c("-f", "CSV", "-resolveDomains", POLIGONOS_CAPA))
    TRUE
  }, error = function(e) FALSE)
  tab <- if (traducido && file.exists(tmp)) {
    suppressWarnings(st_read(tmp, quiet = TRUE, stringsAsFactors = FALSE))
  } else {
    st_drop_geometry(st_read(POLIGONOS_GDB, layer = POLIGONOS_CAPA, quiet = TRUE))
  }
  nombres <- names(tab); low <- tolower(nombres); obj <- tolower(CAMPO_REACTIVACION)
  k <- which(low == obj)
  if (length(k) == 0) k <- which(startsWith(obj, low) & nchar(low) >= 6)      # nombre recortado
  if (length(k) == 0) k <- which(startsWith(low, obj) & !grepl("_resolved$", low))
  if (length(k) == 0) stop("La capa ", POLIGONOS_CAPA, " no tiene un campo que empiece como ",
                           CAMPO_REACTIVACION, ". Campos: ", paste(nombres[!grepl("_resolved$", low)], collapse = ", "))
  campo <- nombres[k[1]]
  kid <- which(low == tolower(CAMPO_ID_POLIGONOS))
  if (length(kid) == 0) stop("La capa ", POLIGONOS_CAPA, " no tiene el campo ", CAMPO_ID_POLIGONOS)
  id_pol <- trimws(as.character(tab[[nombres[kid[1]]]]))
  id_pol <- sub("\\.0+$", "", id_pol)
  # GDAL entrega el texto en UTF-8: se declara asi (convertirlo lo codificaria dos veces)
  utf8 <- function(x) { x <- as.character(x); Encoding(x) <- "UTF-8"; trimws(x) }
  codigo <- utf8(tab[[campo]])
  res <- paste0(campo, "_resolved")
  desc <- if (res %in% nombres) utf8(tab[[res]]) else codigo
  faltan <- setdiff(ids, id_pol)
  if (length(faltan) > 0) stop(length(faltan), " eventos no estan en la capa de poligonos: ",
                               paste(head(faltan, 20), collapse = ", "))
  m <- match(ids, id_pol)
  vals <- table(paste0(codigo[m], " = ", desc[m]))
  rep_("Campo de reactivaciones: ", campo, if (res %in% nombres) " (con dominio: se usa la descripcion de cada codigo)" else " (sin dominio)")
  rep_("Valores en los eventos del modelo (codigo = descripcion): ",
       paste0(names(vals), " (", as.integer(vals), ")", collapse = "; "))
  es <- sin_tildes(desc[m]) == sin_tildes(VALOR_REACTIVACION)
  es[is.na(es)] <- FALSE
  if (!any(es)) aviso(paste0("ningun evento tiene la descripcion '", VALOR_REACTIVACION,
                             "'; revisa VALOR_REACTIVACION con los valores de arriba."))
  es
}

comparar <- function(Lp, Ls, ids_sens, ids, etiqueta) {
  # Compara la corrida principal con una sensibilidad, con FDR y sin corregir (exploratorio)
  m <- match(ids_sens, ids)
  for (col in c("LISA", "LISA_P05")) {
    pm <- Lp[[col]][m]; ps <- Ls[[col]]
    sig_p <- pm != "No significativo"; sig_s <- ps != "No significativo"
    rep_(sprintf("%s, %s: significativos %d con la sensibilidad y %d en la principal; coinciden en categoria %d.",
                 etiqueta, if (col == "LISA") "con FDR" else "con p < 0,05 sin corregir",
                 sum(sig_s), sum(sig_p), sum(sig_p & ps == pm)))
    cambian <- which(ps != pm)
    if (length(cambian) > 0) tabla(data.frame(ID = ids_sens[cambian], Principal = pm[cambian], Sensibilidad = ps[cambian]))
  }
}

# =============================================================================
# PROGRAMA PRINCIPAL
# =============================================================================

main <- function() {
  inicio <- Sys.time()
  if (!dir.exists(CARPETA_SALIDA)) dir.create(CARPETA_SALIDA, recursive = TRUE)
  titulo("FASE 3 · SCRIPT 5 - PARTE C: ESTRUCTURA LOCAL DE LAS ANOMALIAS (LISA)")
  rep_("Ejecucion: ", format(inicio, "%Y-%m-%d %H:%M"), " | ", VERSION)
  rep_(R.version.string, " | spdep ", as.character(packageVersion("spdep")), " | sf ",
       as.character(packageVersion("sf")))
  rep_(""); rep_("Decisiones registradas:")
  for (k in seq_along(DECISIONES_REGISTRADAS)) rep_("  ", k, ". ", DECISIONES_REGISTRADAS[k])

  # ------------------------------------------------------------------ datos
  if (!file.exists(RESULTADOS_B)) stop("No encuentro la salida de la Parte B: ", RESULTADOS_B)
  d <- read.csv(RESULTADOS_B, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE,
                colClasses = c(ID = "character"), na.strings = c("NA"))
  for (cc in names(d)) if (is.character(d[[cc]])) d[[cc]] <- enc2utf8(d[[cc]])
  d$NOTA[is.na(d$NOTA)] <- ""; d$PERFIL[is.na(d$PERFIL)] <- ""
  necesarias <- c("ID", "X", "Y", "RD", "NIVEL", "ANOMALIA", "PERFIL", "NIVEL_ORIGINAL", "NOTA")
  faltan <- setdiff(necesarias, names(d))
  if (length(faltan) > 0) stop("A F3_S4_resultados.csv le faltan columnas: ", paste(faltan, collapse = ", "))
  n <- nrow(d)
  xy <- as.matrix(d[, c("X", "Y")])
  lw <- pesos(xy)
  comp <- n.comp.nb(lw$neighbours)

  titulo("1. DATOS Y VECINDAD")
  rep_("Eventos: ", n, " | variable: distancia robusta de la Parte B (mediana ", f(median(d$RD), 2),
       ", maximo ", f(max(d$RD), 2), " en el evento ", d$ID[which.max(d$RD)], ")")
  tam <- sort(table(comp$comp.id), decreasing = TRUE)
  rep_("Red de ", K_VECINOS, " vecinos: ", comp$nc, " grupos conectados (tamanos: ", paste(tam, collapse = ", "), ").")
  rep_("Cada evento solo se compara con eventos de su propio grupo; un grupo pequeno da resultados")
  rep_("locales menos estables.")

  # ------------------------------------------------------------------ global
  titulo("2. AUTOCORRELACION GLOBAL DE LA DISTANCIA ROBUSTA")
  set.seed(SEMILLA)
  mg <- moran.mc(d$RD, lw, nsim = PERMUTACIONES)
  rep_(sprintf("I de Moran = %s (p = %s, %d permutaciones). Sin patron se esperaria %s.",
               f(mg$statistic), f(mg$p.value, 4), PERMUTACIONES, f(-1 / (n - 1), 4)))
  rep_("Positivo y significativo: los eventos anomalos tienden a estar cerca de otros anomalos.")

  # ------------------------------------------------------------------ LISA
  titulo("3. LISA DE LA DISTANCIA ROBUSTA")
  L <- lisa(d$RD, lw)
  tabla(data.frame(Categoria = CATS, "Con FDR (principal)" = conteo(L$LISA),
                   "Con p < 0,05 sin corregir" = conteo(L$LISA_P05), check.names = FALSE))
  rep_("Con p < 0,05 sin corregir se esperarian unos ", round(n * ALFA), " falsos positivos; por eso la")
  rep_("clasificacion principal usa el FDR. Eventos con FDR menor que 0,10: ", sum(L$P_FDR < 0.10),
       "; menor que 0,20: ", sum(L$P_FDR < 0.20), ". Con p < 0,01 sin corregir: ", sum(L$P < 0.01),
       " (por azar se esperarian ", f(n * 0.01, 1), ").")
  sig <- which(L$LISA != "No significativo")
  if (length(sig) > 0) {
    o <- sig[order(L$P_FDR[sig], -d$RD[sig])]
    rep_("")
    rep_("Eventos con LISA significativo (FDR):")
    tabla(data.frame(ID = d$ID[o], GRUPO = d$GRUPO[o], TIPO = d$TIPO[o], RD = round(d$RD[o], 2),
                     "Parte B" = d$ANOMALIA[o], LISA = L$LISA[o], p = signif(L$P[o], 2),
                     "p FDR" = signif(L$P_FDR[o], 2), Nota = d$NOTA[o], check.names = FALSE))
  } else {
    rep_("Ningun evento tiene LISA significativo con FDR.")
  }
  rep_("")
  rep_("LISA (FDR) frente a la clasificacion de la Parte B:")
  ct <- table(factor(d$ANOMALIA, levels = c("normal", "moderada", "fuerte")), factor(L$LISA, levels = CATS))
  tabla(cbind("Parte B" = rownames(ct), as.data.frame.matrix(ct)))

  # ------------------------------------------------------------------ sensibilidades
  titulo("4. SENSIBILIDADES")
  err <- d$NOTA != ""
  if (any(err)) {
    keep <- !err
    nb_all <- lw$neighbours
    for (i in which(err)) {
      vec <- which(sapply(nb_all, function(v) i %in% v))
      rep_(sprintf("El %s esta entre los %d vecinos de: %s. Su clasificacion sin corregir: %s.", d$ID[i], K_VECINOS,
                   paste(d$ID[vec], collapse = ", "), paste(paste0(d$ID[vec], " ", L$LISA_P05[vec]), collapse = "; ")))
    }
    lw_s <- pesos(xy[keep, , drop = FALSE])
    set.seed(SEMILLA); mgs <- moran.mc(d$RD[keep], lw_s, nsim = PERMUTACIONES)
    rep_(sprintf("Sin errores conocidos: I de Moran global = %s (p = %s); con ellos, %s.",
                 f(mgs$statistic), f(mgs$p.value, 4), f(mg$statistic)))
    Ls <- lisa(d$RD[keep], lw_s)
    comparar(L, Ls, d$ID[keep], d$ID, paste0("Sin errores conocidos (", paste(d$ID[err], collapse = ", "), ")"))
  } else {
    rep_("No hay eventos con errores conocidos marcados en la Parte B.")
  }
  rep_("")
  lisa_sin_react <- rep(NA_character_, n)
  if (nzchar(CAMPO_REACTIVACION)) {
    react <- leer_reactivaciones(d$ID)
    rep_("Reactivaciones marcadas (campo ", CAMPO_REACTIVACION, "): ", sum(react))
    if (sum(react) > 0) {
      keep <- !react
      Lr <- lisa(d$RD[keep], pesos(xy[keep, , drop = FALSE]))
      lisa_sin_react[keep] <- Lr$LISA
      comparar(L, Lr, d$ID[keep], d$ID, "Sin reactivaciones")
    }
  } else {
    aviso("CAMPO_REACTIVACION esta vacio: la sensibilidad sin reactivaciones (decision 21) no se corrio.")
  }

  # ------------------------------------------------------------------ integracion
  titulo("5. INTEGRACION DE LA FASE 3")
  contexto <- c("Alto-Alto" = "en un agrupamiento de anomalias", "Alto-Bajo" = "aislada entre eventos normales",
                "Bajo-Alto" = "normal rodeado de anomalias", "Bajo-Bajo" = "en zona de comportamiento tipico",
                "No significativo" = "sin patron local significativo")
  sintesis <- ifelse(d$NIVEL > 0,
                     paste0("Anomalia ", d$ANOMALIA, "; ", d$PERFIL, "; ", contexto[L$LISA]),
                     ifelse(L$LISA == "No significativo", "Normal",
                            paste0("Normal, ", contexto[L$LISA])))
  final <- data.frame(ID = d$ID, X = d$X, Y = d$Y, GRUPO = d$GRUPO, TIPO = d$TIPO,
                      RE_AREA = d$RE_AREA, RE_HL = d$RE_HL, RE_ELONG = d$RE_ELONG,
                      RD = d$RD, ANOMALIA = d$ANOMALIA, PERFIL = d$PERFIL,
                      ADAPTATIVO = if ("ADAPTATIVO" %in% names(d)) d$ADAPTATIVO else NA,
                      RARO_INVENTARIO = as.integer(d$NIVEL_ORIGINAL > 0),
                      LISA = L$LISA, LISA_P05 = L$LISA_P05, LISA_P = L$P, LISA_P_FDR = L$P_FDR, CUADRANTE = L$CUADRANTE,
                      LISA_SIN_REACT = lisa_sin_react,
                      BAJO_UMBRAL = if ("BAJO_UMBRAL" %in% names(d)) d$BAJO_UMBRAL else NA,
                      FRONTERA = if ("FRONTERA" %in% names(d)) d$FRONTERA else NA,
                      NOTA = d$NOTA, SINTESIS = sintesis, stringsAsFactors = FALSE)
  rs <- sort(table(sintesis), decreasing = TRUE)
  tabla(data.frame(Sintesis = names(rs), Eventos = as.integer(rs)))

  # ------------------------------------------------------------------ graficos
  cat_ <- factor(L$LISA, levels = CATS)
  png(file.path(CARPETA_SALIDA, "F3_S5_mapa_lisa.png"), width = 1500, height = 1500, res = 160)
  par(mar = c(1, 1, 4, 1))
  cat05 <- factor(L$LISA_P05, levels = CATS)
  plot(xy, asp = 1, axes = FALSE, xlab = "", ylab = "", type = "n",
       main = "LISA de la distancia robusta\nrelleno: significativo con FDR | solo contorno: p < 0,05 sin corregir (exploratorio)",
       cex.main = 0.9)
  ns <- cat05 == "No significativo"
  points(xy[ns, , drop = FALSE], pch = 21, bg = COLS["No significativo"], col = "#8A8F94", cex = 0.55)
  ex <- !ns & cat_ == "No significativo"
  points(xy[ex, , drop = FALSE], pch = 21, bg = NA, col = COLS[as.character(cat05[ex])], lwd = 2.2, cex = 1.3)
  fd <- cat_ != "No significativo"
  points(xy[fd, , drop = FALSE], pch = 21, bg = COLS[as.character(cat_[fd])], col = C_OSC, cex = 1.5)
  s <- which(!ns & cat05 %in% c("Alto-Alto", "Alto-Bajo"))
  if (length(s) > 0) text(xy[s, 1], xy[s, 2], d$ID[s], pos = 4, cex = 0.65, col = C_OSC)
  legend("bottomright", bty = "n", pch = 21, pt.bg = COLS, col = C_OSC,
         legend = paste0(CATS, " (FDR: ", conteo(L$LISA), "; sin corregir: ", conteo(L$LISA_P05), ")"), cex = 0.8)
  box(col = "#CFCBC2")
  dev.off()

  png(file.path(CARPETA_SALIDA, "F3_S5_dispersion_moran.png"), width = 1500, height = 1300, res = 160)
  par(mar = c(5, 5, 3, 1))
  z <- (d$RD - mean(d$RD)) / sd(d$RD); wz <- lag.listw(lw, z)
  plot(z, wz, pch = 21, bg = COLS[as.character(cat_)], col = C_OSC, cex = ifelse(cat_ == "No significativo", 0.7, 1.3),
       xlab = "Distancia robusta (estandarizada)", ylab = "Promedio en los 8 vecinos",
       main = sprintf("Diagrama de Moran de la distancia robusta (I = %s)", f(mg$statistic)))
  abline(h = 0, v = 0, col = "#8A8F94", lty = 2)
  abline(lm(wz ~ z), col = "#B04A2E", lwd = 2)
  dev.off()

  # ------------------------------------------------------------------ salidas
  escribir_csv(data.frame(ID = d$ID, L), "F3_S5_resultados.csv")
  escribir_csv(final, "F3_fase3_integrado.csv")
  gpkg <- file.path(CARPETA_SALIDA, "F3_fase3_integrado.gpkg")
  ok_gpkg <- tryCatch({
    puntos <- st_as_sf(final, coords = c("X", "Y"), crs = CRS_EPSG, remove = FALSE)
    if (file.exists(gpkg)) file.remove(gpkg)
    st_write(puntos, gpkg, layer = "fase3_eventos", quiet = TRUE)
    TRUE
  }, error = function(e) { aviso(paste("no se pudo escribir la capa GPKG:", conditionMessage(e))); FALSE })
  titulo("SALIDAS", 2)
  rep_("F3_S5_resultados.csv (LISA por evento) y F3_fase3_integrado.csv (tabla final de la Fase 3)")
  if (ok_gpkg) rep_("F3_fase3_integrado.gpkg: capa de puntos fase3_eventos (EPSG ", CRS_EPSG, "); se abre en ArcGIS Pro con Agregar datos.")
  rep_("Graficos: F3_S5_mapa_lisa.png, F3_S5_dispersion_moran.png | carpeta: ", CARPETA_SALIDA)
  rep_(sprintf("Tiempo: %.1f s", as.numeric(difftime(Sys.time(), inicio, units = "secs"))))
}

tryCatch(main(),
         error = function(e) { rep_(""); rep_("ERROR: ", conditionMessage(e)); stop(e) },
         finally = guardar_reporte())
