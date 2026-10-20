# =====================================================================
# FASE 5 - Parte B (v2). Pesos de evidencia sobre la clase anómala
#
# Anómalos frente a no anómalos dentro de los mismos eventos (asesor).
# Para cada clase B de un factor y D = evento anómalo:
#   W+ = ln[ P(B | D) / P(B | no D) ]      W- = ln[ P(no B | D) / P(no B | no D) ]
#   C  = W+ - W-     s(C) = sqrt(1/N11 + 1/N12 + 1/N21 + 1/N22)
#   Se interpreta la clase si |C / s(C)| > 2. Con una celda en cero se suma
#   0,5 a las cuatro (Haldane) y queda marcado.
# Prueba global por factor: Fisher exacto (clase x anómalo), Holm dentro de
# cada frente. Los pesos no se suman en un índice, así que no se requiere
# independencia condicional entre factores.
#
# Frentes (asesor):
#   1. Estructura residual: todo lo que entró a la GWR (terciles de pendiente,
#      distancia a falla, distancia a drenaje y densidad de lineamientos;
#      grupos de geología, geomorfología y cobertura; tipo y fuente).
#   2. Independiente: unidades sin agrupar (>= MIN_EV eventos).
#
# Chequeos de la pendiente: perfil de las anomalías por tercil y magnitud
# de los residuales de todos los eventos por tercil.
#
# Entrada: F5_eventos.csv (Parte 0 v2)
# Salidas en fase5/parteB: B_pesos_evidencia.csv, B_fisher.csv,
#   B_perfil_por_pendiente.csv, B_residuales_por_pendiente.csv
# =====================================================================

# ------------------------- CONFIGURACIÓN -----------------------------
BASE    <- "D:/Users/tasal/Desktop/Tesis de grado"
ENTRADA <- file.path(BASE, "fase5", "Fase5_csv", "F5_eventos.csv")
SALIDA  <- file.path(BASE, "fase5", "parteB")

CONTINUOS <- c(Pendiente = "PENDIENTE", Distancia_falla = "D_FALLA",
               Distancia_drenaje = "DIST_DRENAJE", Densidad_lineamientos = "DENS_LIN_1000")
CATEG_GWR <- c(Geologia = "GEOLOGIA", Geomorfologia = "GEOMORFOLOGIA",
               Cobertura = "COBERTURA", Tipo = "TIPO", Fuente = "GRUPO")
SIN_AGRUPAR <- c(Geologia_unidades = "GEO_UNI", Geomorfologia_unidades = "GEOM_UNI",
                 Cobertura_unidades = "COB_UNI")
EXCLUIR  <- c("Nube")
MIN_EV   <- 10
N_CLASES <- 3
UMBRAL_C <- 2
set.seed(2026)

# ------------------------- LECTURA Y LIMPIEZA ------------------------
dir.create(SALIDA, showWarnings = FALSE, recursive = TRUE)
d <- read.csv(ENTRADA, fileEncoding = "UTF-8-BOM", check.names = FALSE,
              stringsAsFactors = FALSE)
names(d) <- trimws(sub("^\ufeff", "", names(d)))
if (!"ANOM" %in% names(d)) d$ANOM <- as.integer(tolower(d$ANOMALIA) != "normal")

faltan <- setdiff(c("ANOM", CONTINUOS, CATEG_GWR, SIN_AGRUPAR), names(d))
if (length(faltan)) {
  stop("Faltan columnas: ", paste(faltan, collapse = ", "),
       "\nDisponibles: ", paste(names(d), collapse = ", "))
}

limpiar <- function(x) {            # misma limpieza de la Parte A
  x <- trimws(gsub("\\s+", " ", as.character(x)))
  x[x %in% c("", "NA", "NaN", "None")] <- NA
  x <- gsub("incisadas", "incisados", x, fixed = TRUE)
  x <- gsub("suelo residuales", "suelo residual", x, fixed = TRUE)
  frec <- table(x)
  nom  <- names(frec)
  if (!length(nom)) return(x)
  etiq <- vapply(split(nom, tolower(nom)), function(v) v[which.max(frec[v])], "")
  x <- unname(etiq[tolower(x)])
  x[x %in% EXCLUIR] <- NA
  x
}
for (cc in SIN_AGRUPAR) {
  x <- limpiar(d[[cc]])
  tb <- table(x)
  x[x %in% names(tb)[tb < MIN_EV]] <- NA          # corte de eventos por unidad
  d[[cc]] <- x
}
cat("Eventos:", nrow(d), "| anomalías:", sum(d$ANOM == 1), "\n")

# Terciles de las continuas (sobre todos los eventos)
for (nm in names(CONTINUOS)) {
  x <- d[[CONTINUOS[[nm]]]]
  cortes <- quantile(x, probs = seq(0, 1, length.out = N_CLASES + 1), na.rm = TRUE)
  etiquetas <- sprintf("%s-%s", formatC(head(cortes, -1), format = "f", digits = 1),
                       formatC(cortes[-1], format = "f", digits = 1))
  d[[paste0("CL_", nm)]] <- cut(x, breaks = cortes, labels = etiquetas, include.lowest = TRUE)
}

factores <- rbind(
  data.frame(Frente = "1. Estructura residual", Factor = names(CONTINUOS),
             col = paste0("CL_", names(CONTINUOS))),
  data.frame(Frente = "1. Estructura residual", Factor = names(CATEG_GWR), col = unname(CATEG_GWR)),
  data.frame(Frente = "2. Independiente", Factor = names(SIN_AGRUPAR), col = unname(SIN_AGRUPAR))
)

# ------------------------- PESOS DE UNA CLASE ------------------------
woe_clase <- function(en_clase, anom) {
  n11 <- sum(en_clase & anom == 1); n12 <- sum(!en_clase & anom == 1)
  n21 <- sum(en_clase & anom == 0); n22 <- sum(!en_clase & anom == 0)
  corr <- any(c(n11, n12, n21, n22) == 0)
  k <- if (corr) 0.5 else 0
  a <- n11 + k; b <- n12 + k; c <- n21 + k; e <- n22 + k
  wp <- log((a / (a + b)) / (c / (c + e)))
  wm <- log((b / (a + b)) / (e / (c + e)))
  s  <- sqrt(1 / a + 1 / b + 1 / c + 1 / e)
  data.frame(N_EV = n11 + n21, N_ANOM = n11, PCT_ANOM = 100 * n11 / (n11 + n21),
             W_mas = wp, W_menos = wm, C = wp - wm, s_C = s, C_std = (wp - wm) / s,
             Correccion_0_5 = corr)
}

# ------------------------- CÁLCULO ------------------------------------
woe <- list(); fis <- list()
for (i in seq_len(nrow(factores))) {
  f  <- factores[i, ]
  g  <- d[[f$col]]
  ok <- !is.na(g) & as.character(g) != ""
  g  <- droplevels(factor(g[ok]))
  an <- d$ANOM[ok]
  for (cl in levels(g)) {
    woe[[paste(f$Factor, cl)]] <- cbind(Frente = f$Frente, Factor = f$Factor, Clase = cl,
                                        woe_clase(g == cl, an))
  }
  tab <- table(g, an)
  ft  <- tryCatch(fisher.test(tab, workspace = 2e8),
                  error = function(e) fisher.test(tab, simulate.p.value = TRUE, B = 1e5))
  fis[[f$Factor]] <- data.frame(Frente = f$Frente, Factor = f$Factor, Clases = nlevels(g),
                                n = length(g), N_ANOM = sum(an), p_Fisher = ft$p.value)
}
woe <- do.call(rbind, woe); rownames(woe) <- NULL
woe$Supera_umbral <- abs(woe$C_std) > UMBRAL_C
fis <- do.call(rbind, fis); rownames(fis) <- NULL
fis$p_Holm_frente <- ave(fis$p_Fisher, fis$Frente, FUN = function(x) p.adjust(x, "holm"))

write.csv(woe, file.path(SALIDA, "B_pesos_evidencia.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(fis, file.path(SALIDA, "B_fisher.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# ------------------------- CHEQUEOS DE LA PENDIENTE ------------------
an <- d$ANOM == 1
if ("PERFIL" %in% names(d)) {
  pp <- addmargins(table(Perfil = d$PERFIL[an], Tercil_pendiente = d$CL_Pendiente[an]))
  write.csv(as.data.frame.matrix(pp), file.path(SALIDA, "B_perfil_por_pendiente.csv"),
            fileEncoding = "UTF-8")
}
res_cols <- intersect(c("RD", "RE_AREA", "RE_HL", "RE_ELONG"), names(d))
if (length(res_cols)) {
  mags <- data.frame(row.names = seq_len(nrow(d)))
  for (cc in res_cols) {
    if (cc == "RD") mags$RD <- d$RD else mags[[paste0("abs_", cc)]] <- abs(d[[cc]])
  }
  het <- aggregate(mags, by = list(Tercil_pendiente = d$CL_Pendiente), FUN = median)
  het$Tercil_pendiente <- as.character(het$Tercil_pendiente)
  het_p <- sapply(mags, function(v) kruskal.test(v ~ d$CL_Pendiente)$p.value)
  het <- rbind(het, cbind(Tercil_pendiente = "p Kruskal-Wallis", as.data.frame(t(het_p))))
  write.csv(het, file.path(SALIDA, "B_residuales_por_pendiente.csv"), row.names = FALSE,
            fileEncoding = "UTF-8")
}

# ------------------------- IMPRESIÓN ----------------------------------
redondear <- function(x) {
  num <- vapply(x, is.numeric, logical(1))
  x[num] <- lapply(x[num], function(v) signif(v, 3))
  x
}
cat("\n================ Fisher por factor (Holm dentro de cada frente) ================\n")
print(redondear(fis), row.names = FALSE)
cat("\n================ Pesos de evidencia por clase ================\n")
print(redondear(woe[, c("Factor", "Clase", "N_EV", "N_ANOM", "PCT_ANOM", "C", "s_C",
                        "C_std", "Supera_umbral", "Correccion_0_5")]), row.names = FALSE)
if (exists("pp")) {
  cat("\n================ Perfil de las anomalías por tercil de pendiente ================\n")
  print(pp)
}
if (exists("het")) {
  cat("\n================ Residuales de todos los eventos por tercil de pendiente (medianas) ================\n")
  print(redondear(het), row.names = FALSE)
}
cat("\nArchivos guardados en:", SALIDA, "\n")
