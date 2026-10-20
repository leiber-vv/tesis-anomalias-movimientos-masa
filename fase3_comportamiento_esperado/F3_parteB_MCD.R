# =============================================================================
# FASE 3 · SCRIPT 4 - PARTE B: ANOMALIA CONTEXTUAL CON MCD (R)
# Proyecto de grado: comportamiento anomalo de movimientos en masa, Valle de Aburra
# =============================================================================
#
# Lee la salida de la Parte A (F3_S3_resultados.csv) y mide, para cada evento, que tan
# raro es su comportamiento PARA SUS CONDICIONES: la distancia robusta de sus tres
# residuales estudentizados (area, H/L y elongacion) al centro de los demas eventos,
# estimado con el MCD (minimo determinante de la covarianza).
#
# Por que el MCD: la distancia de Mahalanobis clasica usa la media y la covarianza de
# todos los eventos, y las propias anomalias las deforman y se esconden (enmascaramiento).
# El MCD estima el centro y la covarianza con el 75 % de eventos mas concentrado, asi que
# las anomalias no las mueven (Rousseeuw y Van Driessen 1999).
#
# QUE ENTREGA
#   1. Distancia robusta (RD) y clasica (MD) de cada evento; nivel de anomalia en dos
#      umbrales: moderada (chi-cuadrado 0,975) y fuerte (0,99).
#   2. Perfil de cada anomalia: la dimension con el residual robusto mas extremo y su
#      signo (por ejemplo, "H/L menor de lo esperado (mas movil)").
#   3. Grafico distancia-distancia (Rousseeuw y Van Zomeren 1990) y grafico de cuantiles.
#   4. Referencia: MCD sobre el vector original (sin modelo), para separar lo raro en el
#      inventario de lo raro para sus condiciones.
#   5. Concentracion de anomalias por fuente, tipo y eventos bajo el umbral de
#      completitud (prueba exacta de Fisher).
#   6. Sensibilidad: el MCD sin los eventos bajo el umbral de completitud, y estabilidad
#      del algoritmo (FAST-MCD frente al MCD determinista).
#
# COMO CORRERLO
#   En RStudio: editar las rutas de CONFIGURACION y correr con Source (Ctrl + Shift + S).
#   Usa robustbase y sf, que ya quedaron instalados con GWmodel. Tarda unos segundos.
#   sf lee la marca de completitud directamente de la geodatabase de ArcGIS (solo lectura;
#   puede estar abierta en ArcGIS Pro mientras corre).
#
# SALIDAS (en CARPETA_SALIDA)
#   F3_S4_reporte.txt            todo lo que se imprime en consola
#   F3_S4_resultados.csv         por evento: distancias, nivel, perfil (entrada de la Parte C)
#   F3_S4_anomalias.csv          solo las anomalias, de la mas fuerte a la mas debil
#   F3_S4_distancia_distancia.png, F3_S4_cuantiles.png, F3_S4_pares.png,
#   F3_S4_mapa.png, F3_S4_contexto_vs_global.png
# =============================================================================

VERSION <- "F3_S4 (R) - version 2026-10-03 (d)"

# =============================================================================
# CONFIGURACION - editar solo las rutas
# =============================================================================

RESULTADOS_A   <- "D:/Users/tasal/Desktop/Tesis de grado/fase3/F3_S3_resultados.csv"
CARPETA_SALIDA <- "D:/Users/tasal/Desktop/Tesis de grado/fase3"

# Marca de completitud de la Fase 0: campo de la capa de poligonos del inventario integrado.
# Rutas con barras normales. La capa es la misma que POLIGONOS del script 2.
POLIGONOS_GDB      <- "C:/Users/villa/Documents/ArcGIS/Projects/GDB_proyecto de grado/GDB_proyecto de grado.gdb"
POLIGONOS_CAPA     <- "Inventario_integrado"     # nombre de la capa de poligonos en la geodatabase
CAMPO_ID_POLIGONOS <- "ID_evento"         # el mismo campo de ID que usa el script 2
CAMPO_COMPLETITUD  <- "completitud"
VALOR_BAJO_UMBRAL  <- "bajo_umbral"

# Eventos con errores conocidos: NO se excluyen (decision del autor); se marcan en las
# tablas y graficos para leerlos en la discusion.
ERRORES_CONOCIDOS <- c("73" = "Error de digitalizacion (revisado por el autor)")

# =============================================================================
# DECISIONES METODOLOGICAS (no hace falta editarlas)
# =============================================================================

ALFA_MCD  <- 0.75                  # h = 0,75 n (decision 15)
NIVELES   <- c(0.975, 0.99)        # anomalia moderada y fuerte (decision 16)
SEMILLAS  <- 1:10                  # corridas del FAST-MCD; se elige la de menor determinante
N_SUBMUESTRAS <- 5000              # subconjuntos iniciales de cada corrida
# Umbrales de completitud de la Fase 0 (m2), por valor del campo GRUPO. Solo se usan para
# comparar con la marca oficial de la Fase 0 (o como respaldo si POLIGONOS_GDB es "").
UMBRAL_COMPLETITUD <- c(Mapeo = 421, Institucional = 123)

DIMENSIONES <- list(
  AREA  = c("Area",       "area mayor de lo esperado",            "area menor de lo esperado"),
  HL    = c("H/L",        "H/L mayor de lo esperado (menos movil)", "H/L menor de lo esperado (mas movil)"),
  ELONG = c("Elongacion", "mas alargado de lo esperado",          "menos alargado de lo esperado")
)

DECISIONES_REGISTRADAS <- c(
  "Anomalia contextual (decision 14): distancia robusta de los tres residuales estudentizados del modelo elegido en la Parte A (GWR mixta en las tres dimensiones). El MCD del vector original (sin modelo) es solo una referencia.",
  "MCD reponderado con h = 0,75 n (decisiones 13 y 15): resiste hasta 25 % de eventos anomalos y es mas eficiente que h = 0,5 n. FAST-MCD de robustbase (covMcd), con correccion de consistencia y de muestra pequena.",
  "El MCD es el subconjunto de h eventos con el menor determinante de la covarianza; los algoritmos lo aproximan. Se hacen 10 corridas del FAST-MCD con semillas distintas y una del MCD determinista, y se usa la de menor determinante. Los eventos cuyo nivel cambia entre corridas se marcan como de frontera.",
  "Umbrales en dos niveles (decision 16): anomalia moderada si la distancia robusta al cuadrado supera el cuantil 0,975 de una chi-cuadrado con 3 grados de libertad (9,35); fuerte si supera el 0,99 (11,34). Rousseeuw y Van Zomeren (1990).",
  "Comprobacion de los umbrales (reserva de la decision 16): corte adaptativo de Filzmoser, Garrett y Reimann (2005), que compara la cola de las distancias con la chi-cuadrado y marca solo el exceso. Se reporta junto a los dos niveles, que siguen siendo la clasificacion principal.",
  "Perfil (decision 17): la dimension con el residual robusto estandarizado de mayor valor absoluto, con su signo.",
  "Concentracion por fuente, tipo y eventos bajo el umbral de completitud: prueba exacta de Fisher (decision 22). La marca de completitud es la de la Fase 0 (campo completitud de la capa de poligonos).",
  "Sensibilidad minima (decision 23): MCD sin los eventos bajo el umbral de completitud de la Fase 0. Los residuales son los de la Parte A con todos los eventos; solo cambia la estimacion del centro y la covarianza.",
  "Los eventos con errores conocidos no se excluyen (decision del autor): se marcan para la discusion. El MCD es robusto, asi que no deforman el centro ni la covarianza."
)

# =============================================================================
# PAQUETES Y REPORTE
# =============================================================================

if (!requireNamespace("robustbase", quietly = TRUE)) {
  stop("Falta el paquete robustbase. Instalalo una vez con: install.packages(\"robustbase\")", call. = FALSE)
}
suppressPackageStartupMessages(library(robustbase))
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
  con <- file(file.path(CARPETA_SALIDA, "F3_S4_reporte.txt"), open = "w", encoding = "UTF-8")
  writeLines(REPORTE, con); close(con)
}
escribir_csv <- function(df, nombre) {
  con <- file(file.path(CARPETA_SALIDA, nombre), open = "w", encoding = "UTF-8")
  on.exit(close(con)); cat("\ufeff", file = con)
  write.csv(df, con, row.names = FALSE, na = "")
}
f <- function(x, d = 2) formatC(x, format = "f", digits = d)

# Colores (los mismos de los mapas de la Parte A)
C_OSC <- "#1F2A33"; C_GRIS <- "#B9B5AB"; C_MOD <- "#E0A15E"; C_FUE <- "#B04A2E"; C_AZUL <- "#2D5F8B"

# =============================================================================
# FUNCIONES
# =============================================================================

mcd <- function(X, alfa = ALFA_MCD) {
  # Varias corridas; se queda con la de menor determinante (el criterio que define el MCD)
  corridas <- lapply(SEMILLAS, function(s) { set.seed(s); covMcd(X, alpha = alfa, nsamp = N_SUBMUESTRAS) })
  corridas[[length(corridas) + 1]] <- covMcd(X, alpha = alfa, nsamp = "deterministic")
  crit <- sapply(corridas, `[[`, "crit")
  mejor <- corridas[[which.min(crit)]]
  mejor$corridas <- corridas
  mejor$n_optimo <- sum(abs(crit - min(crit)) < 1e-9)
  mejor
}

distancia <- function(X, centro, cov) sqrt(mahalanobis(X, centro, cov))

nivel_de <- function(rd, cortes) ifelse(rd > cortes[2], 2L, ifelse(rd > cortes[1], 1L, 0L))

perfil <- function(Z) {
  # Dimension con el residual robusto estandarizado de mayor valor absoluto, con su signo
  j <- apply(abs(Z), 1, which.max)
  vapply(seq_len(nrow(Z)), function(i) {
    d <- DIMENSIONES[[colnames(Z)[j[i]]]]
    if (Z[i, j[i]] >= 0) d[2] else d[3]
  }, character(1))
}

elipse <- function(centro, S, q, k = 200) {
  # Puntos de la elipse de tolerancia {x : (x - c)' S^-1 (x - c) = q}
  e <- eigen(S); t <- seq(0, 2 * pi, length.out = k)
  u <- cbind(cos(t), sin(t)) %*% diag(sqrt(e$values * q)) %*% t(e$vectors)
  sweep(u, 2, centro, "+")
}

leer_completitud <- function(ids) {
  # Marca oficial de la Fase 0, leida de la geodatabase con sf. Devuelve TRUE si esta bajo umbral.
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Falta el paquete sf para leer la geodatabase. Instalalo con: install.packages(\"sf\")")
  }
  if (!dir.exists(POLIGONOS_GDB)) stop("No encuentro la geodatabase: ", POLIGONOS_GDB)
  capas <- sf::st_layers(POLIGONOS_GDB)$name
  if (!POLIGONOS_CAPA %in% capas) {
    stop("La geodatabase no tiene la capa ", POLIGONOS_CAPA, ". Capas disponibles: ", paste(capas, collapse = ", "))
  }
  tab <- sf::st_drop_geometry(sf::st_read(POLIGONOS_GDB, layer = POLIGONOS_CAPA, quiet = TRUE))
  campo <- function(nombre) {
    k <- which(tolower(names(tab)) == tolower(nombre))
    if (length(k) == 0) stop("La capa ", POLIGONOS_CAPA, " no tiene el campo ", nombre,
                             ". Campos: ", paste(names(tab), collapse = ", "))
    names(tab)[k[1]]
  }
  id_pol <- tab[[campo(CAMPO_ID_POLIGONOS)]]
  if (is.numeric(id_pol)) id_pol <- format(id_pol, scientific = FALSE, trim = TRUE, drop0trailing = TRUE)
  id_pol <- trimws(as.character(id_pol))
  val <- tolower(trimws(enc2utf8(as.character(tab[[campo(CAMPO_COMPLETITUD)]]))))
  faltan <- setdiff(ids, id_pol)
  if (length(faltan) > 0) stop(length(faltan), " eventos del modelo no estan en la capa de poligonos: ",
                               paste(head(faltan, 20), collapse = ", "))
  if (anyDuplicated(id_pol[id_pol %in% ids])) stop("Hay IDs repetidos en la capa de poligonos.")
  valores <- unique(val[id_pol %in% ids])
  raros <- setdiff(valores, c(tolower(VALOR_BAJO_UMBRAL), "completo"))
  if (length(raros) > 0) aviso(paste0("valores inesperados en ", CAMPO_COMPLETITUD, ": ",
                                      paste(raros, collapse = ", "), " (se tratan como completos)"))
  val[match(ids, id_pol)] == tolower(VALOR_BAJO_UMBRAL)
}

corte_adaptativo <- function(d2, p) {
  # Filzmoser, Garrett y Reimann (2005), con la formulacion de arw (paquete mvoutlier)
  n <- length(d2); delta <- qchisq(0.975, p)
  pcrit <- if (p <= 10) (0.24 - 0.003 * p) / sqrt(n) else (0.252 - 0.0018 * p) / sqrt(n)
  d2ord <- sort(d2)
  dif <- pchisq(d2ord, p) - (0.5:n) / n
  i <- (d2ord >= delta) & (dif > 0)
  alfan <- if (sum(i) == 0) 0 else max(dif[i])
  if (alfan < pcrit) alfan <- 0
  list(corte = if (alfan > 0) max(d2ord[n - ceiling(n * alfan)], delta) else Inf, alfan = alfan, pcrit = pcrit)
}

fisher_tabla <- function(grupo, anomalo, nombre) {
  tab <- table(grupo, factor(anomalo, levels = c(FALSE, TRUE)))
  pr <- tryCatch(fisher.test(tab, workspace = 2e7)$p.value,
                 error = function(e) { set.seed(1); fisher.test(tab, simulate.p.value = TRUE, B = 1e5)$p.value })
  df <- data.frame(Clase = rownames(tab), Eventos = as.vector(rowSums(tab)),
                   Anomalias = as.vector(tab[, "TRUE"]),
                   "% anomalias" = round(100 * tab[, "TRUE"] / rowSums(tab), 1), check.names = FALSE)
  list(tabla = df, p = pr, nombre = nombre)
}

# =============================================================================
# PROGRAMA PRINCIPAL
# =============================================================================

main <- function() {
  inicio <- Sys.time()
  if (!dir.exists(CARPETA_SALIDA)) dir.create(CARPETA_SALIDA, recursive = TRUE)
  titulo("FASE 3 · SCRIPT 4 - PARTE B: ANOMALIA CONTEXTUAL CON MCD (R)")
  rep_("Ejecucion: ", format(inicio, "%Y-%m-%d %H:%M"), " | ", VERSION)
  rep_(R.version.string, " | robustbase ", as.character(packageVersion("robustbase")))
  rep_(""); rep_("Decisiones registradas:")
  for (k in seq_along(DECISIONES_REGISTRADAS)) rep_("  ", k, ". ", DECISIONES_REGISTRADAS[k])

  # ------------------------------------------------------------------ datos
  if (!file.exists(RESULTADOS_A)) stop("No encuentro la salida de la Parte A: ", RESULTADOS_A)
  d <- read.csv(RESULTADOS_A, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE,
                colClasses = c(ID = "character"))
  for (cc in names(d)) if (is.character(d[[cc]])) d[[cc]] <- enc2utf8(d[[cc]])
  dims <- names(DIMENSIONES)
  faltan <- setdiff(c("ID", "X", "Y", "GRUPO", "TIPO", paste0("RE_", dims), paste0("PRED_", dims),
                      paste0("RES_", dims)), names(d))
  if (length(faltan) > 0) stop("A F3_S3_resultados.csv le faltan columnas: ", paste(faltan, collapse = ", "))
  n <- nrow(d); p <- length(dims)
  X <- as.matrix(d[, paste0("RE_", dims)]); colnames(X) <- dims
  if (any(!is.finite(X))) stop("Hay residuales vacios o no finitos.")
  Xo <- as.matrix(d[, paste0("PRED_", dims)] + d[, paste0("RES_", dims)]); colnames(Xo) <- dims
  area_m2 <- 10^Xo[, "AREA"]
  bajo_calc <- area_m2 < UMBRAL_COMPLETITUD[d$GRUPO]
  if (any(is.na(bajo_calc))) stop("Valores de GRUPO sin umbral de completitud: ",
                                  paste(unique(d$GRUPO[is.na(bajo_calc)]), collapse = ", "))
  oficial <- nzchar(POLIGONOS_GDB)
  bajo <- if (oficial) leer_completitud(d$ID) else bajo_calc
  nota <- ifelse(d$ID %in% names(ERRORES_CONOCIDOS), ERRORES_CONOCIDOS[d$ID], "")

  titulo("1. DATOS")
  rep_("Eventos: ", n, " | dimensiones: ", paste(sapply(DIMENSIONES, `[`, 1), collapse = ", "),
       " (residuales estudentizados de la Parte A)")
  if (oficial) {
    rep_("Eventos bajo el umbral de completitud: ", sum(bajo), " (campo ", CAMPO_COMPLETITUD, " de ",
         POLIGONOS_CAPA, ", marca de la Fase 0)")
    rep_("Comparacion con el calculo desde el area observada (", paste(names(UMBRAL_COMPLETITUD),
         UMBRAL_COMPLETITUD, "m2", collapse = "; "), "):")
    tabla(data.frame(" " = c("Fase 0: bajo umbral", "Fase 0: completo"),
                     "Calculo: bajo umbral" = c(sum(bajo & bajo_calc), sum(!bajo & bajo_calc)),
                     "Calculo: completo" = c(sum(bajo & !bajo_calc), sum(!bajo & !bajo_calc)), check.names = FALSE))
    if (any(bajo != bajo_calc)) {
      rep_("Las diferencias pueden venir de como se midio el area en la Fase 0; se usa la marca de la Fase 0.")
    }
  } else {
    rep_("Eventos bajo el umbral de completitud: ", sum(bajo), " (calculado con el area observada y los")
    rep_("umbrales de la Fase 0, porque POLIGONOS_GDB esta vacio)")
    aviso("sin la marca oficial de la Fase 0; solo afecta a la seccion 6 y a la sensibilidad.")
  }
  rep_("Errores conocidos marcados (no se excluyen): ", paste(names(ERRORES_CONOCIDOS), collapse = ", "))

  # ------------------------------------------------------------------ MCD
  m <- mcd(X)
  cortes <- sqrt(qchisq(NIVELES, df = p))
  RD <- distancia(X, m$center, m$cov)
  MD <- distancia(X, colMeans(X), cov(X))
  nivel <- nivel_de(RD, cortes)
  Z <- sweep(sweep(X, 2, m$center, "-"), 2, sqrt(diag(m$cov)), "/")
  ca <- corte_adaptativo(RD^2, p)
  adapt <- RD^2 > ca$corte
  perf <- ifelse(nivel > 0, perfil(Z), "")
  niveles_corridas <- sapply(m$corridas, function(mc) nivel_de(distancia(X, mc$center, mc$cov), cortes))
  frontera <- apply(niveles_corridas, 1, function(v) length(unique(v)) > 1)

  titulo("2. ESTIMACION ROBUSTA (MCD)")
  rep_(sprintf("h = %d de %d eventos (alfa = %.2f). Eventos con peso 0 en la reponderacion: %d.",
               m$quan, n, ALFA_MCD, sum(m$mcd.wt == 0)))
  rep_(sprintf("Corridas: %d (10 FAST-MCD y 1 determinista); %d alcanzan el menor determinante (log det = %s).",
               length(m$corridas), m$n_optimo, f(m$crit, 6)))
  rep_("")
  rep_("Centro y dispersion de los residuales: clasico frente a robusto")
  tabla(data.frame(Dimension = sapply(DIMENSIONES, `[`, 1),
                   "Media clasica" = round(colMeans(X), 3), "Centro MCD" = round(m$center, 3),
                   "DE clasica" = round(apply(X, 2, sd), 3), "DE robusta" = round(sqrt(diag(m$cov)), 3),
                   check.names = FALSE))
  rep_("")
  rep_("Correlaciones entre los residuales (clasica / robusta):")
  cc <- cor(X); cr <- cov2cor(m$cov)
  for (i in 1:(p - 1)) for (j in (i + 1):p) {
    rep_(sprintf("  %s - %s: %s / %s", DIMENSIONES[[i]][1], DIMENSIONES[[j]][1], f(cc[i, j], 3), f(cr[i, j], 3)))
  }
  rep_("Una correlacion robusta alta indicaria que dos dimensiones se desvian juntas; el MCD la tiene en cuenta.")

  # ------------------------------------------------------------------ anomalias
  titulo("3. ANOMALIAS CONTEXTUALES")
  rep_(sprintf("Umbrales de la distancia robusta: moderada > %s (chi-cuadrado 0,975 con 3 gl) | fuerte > %s (0,99)",
               f(cortes[1], 3), f(cortes[2], 3)))
  esp <- n * (1 - NIVELES)
  tabla(data.frame(Nivel = c("Moderada o fuerte (> 0,975)", "Fuerte (> 0,99)"),
                   Eventos = c(sum(nivel >= 1), sum(nivel == 2)),
                   "%" = round(100 * c(sum(nivel >= 1), sum(nivel == 2)) / n, 1),
                   "Esperados por azar" = round(esp, 1), check.names = FALSE))
  rep_("Esperados por azar: los eventos que pasarian el umbral aunque no hubiera anomalias, si los")
  rep_("residuales fueran normales. Si los observados los superan con claridad, hay anomalias reales.")
  rep_("")
  qs <- c(0.25, 0.5, 0.75, 0.9)
  rep_("Calibracion (distancia robusta al cuadrado frente a la chi-cuadrado con 3 gl, por cuantiles):")
  tabla(data.frame(Cuantil = qs, Observado = round(quantile(RD^2, qs), 2), Teorico = round(qchisq(qs, p), 2),
                   Cociente = round(quantile(RD^2, qs) / qchisq(qs, p), 2), row.names = NULL))
  rep_("Cociente cercano a 1 en el centro: el MCD esta bien calibrado. Si crece hacia la cola, la cola es")
  rep_("mas pesada que la normal y parte de las anomalias moderadas puede venir de esa forma de la cola.")
  rep_("")
  rep_(sprintf("Corte adaptativo (Filzmoser, Garrett y Reimann 2005): exceso de la cola %.3f (critico %.3f);",
               ca$alfan, ca$pcrit))
  if (is.finite(ca$corte)) {
    rep_(sprintf("corte en la distancia %s -> %d anomalias: %d de las %d fuertes y %d moderadas (%s).",
                 f(sqrt(ca$corte), 3), sum(adapt), sum(adapt & nivel == 2), sum(nivel == 2),
                 sum(adapt & nivel == 1), paste(d$ID[adapt & nivel == 1], collapse = ", ")))
  } else {
    rep_("la cola no se aparta de la chi-cuadrado: el metodo no marca anomalias.")
  }
  rep_("")
  o <- order(-RD); o <- o[nivel[o] > 0]
  lista <- data.frame(ID = d$ID[o], GRUPO = d$GRUPO[o], TIPO = d$TIPO[o],
                      RD = round(RD[o], 2), MD = round(MD[o], 2),
                      Nivel = c("", "moderada", "fuerte")[nivel[o] + 1],
                      "z area" = round(Z[o, "AREA"], 2), "z H/L" = round(Z[o, "HL"], 2),
                      "z elong." = round(Z[o, "ELONG"], 2), Perfil = perf[o],
                      Adaptativo = ifelse(adapt[o], "si", ""),
                      "Bajo umbral" = ifelse(bajo[o], "si", ""), Frontera = ifelse(frontera[o], "si", ""),
                      Nota = nota[o], check.names = FALSE)
  rep_("Anomalias, de la mas fuerte a la mas debil (z = residual robusto estandarizado):")
  tabla(lista)
  rep_("")
  rep_("Perfiles de las anomalias:")
  tp <- sort(table(perf[nivel > 0]), decreasing = TRUE)
  for (k in names(tp)) rep_("  ", tp[[k]], "  ", k)

  # ------------------------------------------------------------------ distancia-distancia
  titulo("4. GRAFICO DISTANCIA-DISTANCIA (clasica frente a robusta)")
  c1 <- cortes[1]
  q_ambas <- sum(RD > c1 & MD > c1); q_rob <- sum(RD > c1 & MD <= c1); q_cla <- sum(RD <= c1 & MD > c1)
  rep_(sprintf("Anomalos con las dos distancias: %d | solo con la robusta (enmascarados para la clasica): %d | solo con la clasica: %d",
               q_ambas, q_rob, q_cla))
  rep_("Los enmascarados son eventos que la distancia clasica no detecta porque las propias anomalias")
  rep_("inflan la covarianza clasica: es la razon para usar el MCD.")

  # ------------------------------------------------------------------ referencia global
  titulo("5. REFERENCIA: MCD SOBRE EL VECTOR ORIGINAL (sin modelo)")
  mo <- mcd(Xo)
  RDo <- distancia(Xo, mo$center, mo$cov)
  nivel_o <- nivel_de(RDo, cortes)
  ctx <- nivel >= 1; glo <- nivel_o >= 1
  tabla(data.frame(" " = c("Anomalo para sus condiciones", "Normal para sus condiciones"),
                   "Raro en el inventario" = c(sum(ctx & glo), sum(!ctx & glo)),
                   "Comun en el inventario" = c(sum(ctx & !glo), sum(!ctx & !glo)), check.names = FALSE))
  rep_("Anomalo para sus condiciones pero comun en el inventario: su comportamiento es normal en el valle,")
  rep_("pero no donde ocurrio. Raro en el inventario pero normal para sus condiciones: es extremo, pero")
  rep_("sus condiciones lo explican. La tesis usa el primer criterio (decision 14).")

  # ------------------------------------------------------------------ concentracion
  titulo("6. CONCENTRACION DE LAS ANOMALIAS (prueba exacta de Fisher)")
  for (ft in list(fisher_tabla(d$GRUPO, ctx, "Fuente"), fisher_tabla(d$TIPO, ctx, "Tipo de movimiento"),
                  fisher_tabla(ifelse(bajo, "Bajo el umbral", "Sobre el umbral"), ctx, "Completitud"))) {
    rep_(""); rep_(ft$nombre, " (p de Fisher = ", f(ft$p, 3), ")")
    tabla(ft$tabla)
  }
  rep_("")
  rep_("p < 0,05: las anomalias no se reparten al azar entre las clases.")

  # ------------------------------------------------------------------ sensibilidad
  titulo("7. SENSIBILIDAD Y ESTABILIDAD")
  ms <- mcd(X[!bajo, , drop = FALSE])
  RDs <- distancia(X, ms$center, ms$cov)
  ctx_s <- RDs > c1
  rep_(sprintf("Sin los %d eventos bajo umbral en la estimacion: %d anomalias entre los %d restantes (antes, %d de ellos); coinciden %d.",
               sum(bajo), sum(ctx_s[!bajo]), sum(!bajo), sum(ctx[!bajo]), sum(ctx_s[!bajo] & ctx[!bajo])))
  rep_(sprintf("Estabilidad del algoritmo: %d de %d corridas llegan al mismo minimo. Eventos de frontera (su nivel cambia en alguna corrida): %s.",
               m$n_optimo, length(m$corridas), if (any(frontera)) paste(d$ID[frontera], collapse = ", ") else "ninguno"))
  rep_("Los eventos de frontera estan justo en el limite de un umbral: su clasificacion es menos firme.")

  # ------------------------------------------------------------------ graficos
  png(file.path(CARPETA_SALIDA, "F3_S4_distancia_distancia.png"), width = 1600, height = 1300, res = 160)
  par(mar = c(5, 5, 3, 1))
  col <- c(C_GRIS, C_MOD, C_FUE)[nivel + 1]
  plot(MD, RD, pch = 21, bg = col, col = C_OSC, cex = ifelse(nivel > 0, 1.3, 0.8),
       xlab = "Distancia clasica (Mahalanobis)", ylab = "Distancia robusta (MCD)",
       main = "Grafico distancia-distancia de los residuales")
  abline(0, 1, col = C_OSC, lty = 3)
  abline(h = cortes, v = cortes[1], col = c(C_MOD, C_FUE, C_GRIS), lty = 2, lwd = 1.5)
  top <- which(nivel == 2)
  text(MD[top], RD[top], d$ID[top], pos = 4, cex = 0.75, col = C_OSC)
  legend("topleft", bty = "n", pch = 21, pt.bg = c(C_GRIS, C_MOD, C_FUE), col = C_OSC,
         legend = c("Normal", "Anomalia moderada", "Anomalia fuerte (con su ID)"))
  dev.off()

  png(file.path(CARPETA_SALIDA, "F3_S4_cuantiles.png"), width = 1400, height = 1300, res = 160)
  par(mar = c(5, 5, 3, 1))
  teo <- qchisq(ppoints(n), df = p)
  plot(teo, sort(RD^2), pch = 21, bg = C_GRIS, col = C_OSC, xlab = "Cuantiles de chi-cuadrado (3 gl)",
       ylab = "Distancia robusta al cuadrado (ordenada)", main = "Grafico de cuantiles de la distancia robusta")
  abline(0, 1, col = C_FUE, lwd = 2)
  abline(h = cortes^2, col = c(C_MOD, C_FUE), lty = 2)
  if (is.finite(ca$corte)) abline(h = ca$corte, col = C_AZUL, lty = 3, lwd = 2)
  legend("topleft", bty = "n", lty = c(2, 2, 3), col = c(C_MOD, C_FUE, C_AZUL), lwd = c(1, 1, 2),
         legend = c("Umbral 0,975", "Umbral 0,99", "Corte adaptativo"))
  dev.off()

  png(file.path(CARPETA_SALIDA, "F3_S4_pares.png"), width = 2400, height = 850, res = 160)
  par(mfrow = c(1, 3), mar = c(5, 5, 3, 1), oma = c(2.5, 0, 0, 0))
  for (par_ in list(c(1, 2), c(1, 3), c(2, 3))) {
    i <- par_[1]; j <- par_[2]
    plot(X[, i], X[, j], pch = 21, bg = col, col = C_OSC, cex = ifelse(nivel > 0, 1.2, 0.7),
         xlab = paste("Residual", DIMENSIONES[[i]][1]), ylab = paste("Residual", DIMENSIONES[[j]][1]))
    lines(elipse(m$center[c(i, j)], m$cov[c(i, j), c(i, j)], qchisq(0.975, 2)), col = C_FUE, lwd = 2)
    lines(elipse(colMeans(X)[c(i, j)], cov(X)[c(i, j), c(i, j)], qchisq(0.975, 2)), col = C_AZUL, lwd = 2, lty = 2)
  }
  mtext("Elipses del 97,5 %: robusta (continua) y clasica (punteada). Color: nivel en las tres dimensiones.",
        side = 1, line = 1, outer = TRUE, cex = 0.8)
  dev.off()

  png(file.path(CARPETA_SALIDA, "F3_S4_mapa.png"), width = 1500, height = 1500, res = 160)
  par(mar = c(1, 1, 3, 1))
  o2 <- order(nivel)
  plot(d$X[o2], d$Y[o2], asp = 1, axes = FALSE, xlab = "", ylab = "", pch = 21, bg = col[o2], col = C_OSC,
       cex = c(0.6, 1.3, 1.7)[nivel[o2] + 1], main = "Anomalias contextuales (Parte B)")
  text(d$X[top], d$Y[top], d$ID[top], pos = 4, cex = 0.7, col = C_OSC)
  legend("bottomright", bty = "n", pch = 21, pt.bg = c(C_GRIS, C_MOD, C_FUE), col = C_OSC,
         legend = c("Normal", "Moderada", "Fuerte"))
  box(col = C_GRIS)
  dev.off()

  png(file.path(CARPETA_SALIDA, "F3_S4_contexto_vs_global.png"), width = 1500, height = 1300, res = 160)
  par(mar = c(5, 5, 3, 1))
  plot(RDo, RD, pch = 21, bg = col, col = C_OSC, xlab = "Distancia robusta del vector original (sin modelo)",
       ylab = "Distancia robusta de los residuales", main = "Raro para sus condiciones frente a raro en el inventario")
  abline(h = cortes[1], v = cortes[1], col = C_OSC, lty = 2)
  dev.off()

  # ------------------------------------------------------------------ salidas
  salida <- data.frame(ID = d$ID, X = d$X, Y = d$Y, GRUPO = d$GRUPO, TIPO = d$TIPO,
                       RE_AREA = d$RE_AREA, RE_HL = d$RE_HL, RE_ELONG = d$RE_ELONG,
                       Z_AREA = Z[, "AREA"], Z_HL = Z[, "HL"], Z_ELONG = Z[, "ELONG"],
                       MD = MD, RD = RD, NIVEL = nivel,
                       ANOMALIA = c("normal", "moderada", "fuerte")[nivel + 1], PERFIL = perf,
                       RD_ORIGINAL = RDo, NIVEL_ORIGINAL = nivel_o, BAJO_UMBRAL = as.integer(bajo),
                       ADAPTATIVO = as.integer(adapt), PESO_MCD = m$mcd.wt, FRONTERA = as.integer(frontera), NOTA = nota)
  escribir_csv(salida, "F3_S4_resultados.csv")
  escribir_csv(lista, "F3_S4_anomalias.csv")
  titulo("SALIDAS", 2)
  rep_("F3_S4_resultados.csv (", n, " eventos; la columna RD es la entrada de la Parte C)")
  rep_("F3_S4_anomalias.csv (", nrow(lista), " anomalias) y graficos PNG en: ", CARPETA_SALIDA)
  rep_(sprintf("Tiempo: %.1f s", as.numeric(difftime(Sys.time(), inicio, units = "secs"))))
}

tryCatch(main(),
         error = function(e) { rep_(""); rep_("ERROR: ", conditionMessage(e)); stop(e) },
         finally = guardar_reporte())
