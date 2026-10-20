# =====================================================================
# FASE 5 - Parte A (v2). PERMANOVA -> Kruskal-Wallis
#
# Pregunta (contexto, según el asesor): ¿el comportamiento típico de los
# movimientos difiere entre unidades? Se usan TODOS los eventos y las
# métricas crudas: log10 del área, de H/L y de la elongación, estandarizadas
# y con distancia euclidiana.
#
# Factores: geología, geomorfología y cobertura, cada una con los grupos que
# entraron a la GWR y con las unidades sin agrupar, más el tipo de movimiento.
#   - Unidades sin agrupar: nombres limpiados (duplicados por escritura),
#     "Nube" como dato faltante y solo las unidades con >= MIN_EV eventos.
#   - PERMANOVA de un factor (9.999 permutaciones) + PERMDISP; Holm entre
#     los 7 factores.
#   - PERMANOVA secuencial "grupo + unidad": el término de la unidad mide si
#     las unidades difieren DENTRO de su grupo (lo que la GWR no vio).
#     Holm entre las 3 capas.
#   - Kruskal-Wallis por métrica con p por permutaciones (coin), épsilon²
#     y Holm entre las 3 métricas de cada factor. Se interpreta solo en los
#     factores con PERMANOVA significativa (columna PERMANOVA_sig).
#   - Medianas por clase en log10 sin estandarizar.
#
# Entrada: F5_eventos.csv (Parte 0 v2)
# Salidas en fase5/parteA: A_permanova.csv, A_anidada.csv, A_kruskal.csv,
#                          A_medianas.csv, A_unidades_excluidas.csv
# =====================================================================

paquetes <- c("vegan", "coin")
for (p in paquetes) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({ library(vegan); library(coin) })

# ------------------------- CONFIGURACIÓN -----------------------------
BASE    <- "D:/Users/tasal/Desktop/Tesis de grado"
ENTRADA <- file.path(BASE, "fase5", "Fase5_csv", "F5_eventos.csv")
SALIDA  <- file.path(BASE, "fase5", "parteA")

METRICAS <- c(Area = "LOG_AREA", HL = "LOG_HL", Elongacion = "LOG_ELONG")

# nombre, columna, columna del grupo (NA si el factor ya es un grupo)
FACTORES <- data.frame(
  nombre = c("Geologia (grupos)", "Geologia (unidades)",
             "Geomorfologia (grupos)", "Geomorfologia (unidades)",
             "Cobertura (grupos)", "Cobertura (unidades)", "Tipo de movimiento"),
  col    = c("GEOLOGIA", "GEO_UNI", "GEOMORFOLOGIA", "GEOM_UNI",
             "COBERTURA", "COB_UNI", "TIPO"),
  grupo  = c(NA, "GEOLOGIA", NA, "GEOMORFOLOGIA", NA, "COBERTURA", NA),
  stringsAsFactors = FALSE
)
SIN_AGRUPAR <- c("GEO_UNI", "GEOM_UNI", "COB_UNI")
EXCLUIR     <- c("Nube")      # clases que no son unidades
MIN_EV      <- 10             # eventos mínimos por clase (decisión de autor)
N_PERM      <- 9999
ALFA        <- 0.05
set.seed(2026)

# ------------------------- LECTURA Y LIMPIEZA ------------------------
dir.create(SALIDA, showWarnings = FALSE, recursive = TRUE)
d <- read.csv(ENTRADA, fileEncoding = "UTF-8-BOM", check.names = FALSE,
              stringsAsFactors = FALSE)
names(d) <- trimws(sub("^\ufeff", "", names(d)))

faltan <- setdiff(c(METRICAS, FACTORES$col), names(d))
if (length(faltan)) {
  stop("Faltan columnas: ", paste(faltan, collapse = ", "),
       "\nDisponibles: ", paste(names(d), collapse = ", "))
}

limpiar <- function(x) {
  x <- trimws(gsub("\\s+", " ", as.character(x)))
  x[x %in% c("", "NA", "NaN", "None")] <- NA
  x <- gsub("incisadas", "incisados", x, fixed = TRUE)
  x <- gsub("suelo residuales", "suelo residual", x, fixed = TRUE)
  # mismas unidades escritas con distinta mayúscula: se usa la forma más frecuente
  frec <- table(x)
  nom  <- names(frec)
  if (!length(nom)) return(x)
  etiq <- vapply(split(nom, tolower(nom)), function(v) v[which.max(frec[v])], "")
  x <- unname(etiq[tolower(x)])
  x[x %in% EXCLUIR] <- NA
  x
}
for (cc in SIN_AGRUPAR) d[[cc]] <- limpiar(d[[cc]])

Y_all <- scale(as.matrix(d[, METRICAS]))
colnames(Y_all) <- names(METRICAS)
cat("Eventos:", nrow(d), "\n")

# ------------------------- FUNCIONES ---------------------------------
preparar <- function(col) {
  g  <- d[[col]]
  ok <- !is.na(g) & complete.cases(Y_all)
  tb <- table(g[ok])
  chicas <- names(tb)[tb < MIN_EV]
  excl <- data.frame(Factor = rep(col, length(chicas)), Clase = chicas,
                     Eventos = as.integer(tb[chicas]))
  ok <- ok & !(g %in% chicas)
  list(ok = ok, g = droplevels(factor(g[ok])), excl = excl,
       sin_clase = sum(is.na(d[[col]])))
}

redondear <- function(x) {
  num <- vapply(x, is.numeric, logical(1))
  x[num] <- lapply(x[num], function(v) signif(v, 3))
  x
}

# ------------------------- CÁLCULO -----------------------------------
perm <- list(); anid <- list(); kw <- list(); med <- list(); excl <- list()

for (i in seq_len(nrow(FACTORES))) {
  f  <- FACTORES[i, ]
  p  <- preparar(f$col)
  y  <- Y_all[p$ok, , drop = FALSE]
  g  <- p$g
  D  <- dist(y)
  excl[[f$col]] <- p$excl
  
  # PERMANOVA de un factor y homogeneidad de dispersiones
  a  <- adonis2(D ~ g, data = data.frame(g = g), permutations = N_PERM)
  pb <- permutest(betadisper(D, g), permutations = N_PERM)
  perm[[f$col]] <- data.frame(
    Factor = f$nombre, n = length(g), Clases = nlevels(g), Clase_min = min(table(g)),
    Excluidos = sum(p$excl$Eventos) + p$sin_clase,
    gl = a$Df[1], R2 = a$R2[1], F = a$F[1], p = a$`Pr(>F)`[1],
    Disp_F = pb$tab$F[1], Disp_p = pb$tab$`Pr(>F)`[1]
  )
  
  # PERMANOVA secuencial grupo + unidad (solo unidades sin agrupar)
  if (!is.na(f$grupo)) {
    gr <- factor(d[[f$grupo]][p$ok])
    an <- adonis2(D ~ gr + g, data = data.frame(gr = gr, g = g),
                  permutations = N_PERM, by = "terms")
    anid[[f$col]] <- data.frame(
      Factor = f$nombre, n = length(g),
      R2_grupo = an["gr", "R2"], p_grupo = an["gr", "Pr(>F)"],
      gl_unidad = an["g", "Df"], R2_unidad_dentro = an["g", "R2"],
      F_unidad = an["g", "F"], p_unidad_dentro = an["g", "Pr(>F)"]
    )
  }
  
  # Kruskal-Wallis por métrica (p por permutaciones)
  for (m in colnames(y)) {
    k <- kruskal_test(v ~ g, data = data.frame(v = y[, m], g = g),
                      distribution = approximate(nresample = N_PERM))
    H <- as.numeric(statistic(k))
    kw[[paste(f$col, m)]] <- data.frame(
      Factor = f$nombre, Metrica = m, H = H, p = as.numeric(pvalue(k)),
      eps2 = H / (length(g) - 1)
    )
  }
  
  # Medianas por clase en log10 sin estandarizar
  crudo <- d[p$ok, METRICAS, drop = FALSE]
  names(crudo) <- paste0("med_log_", names(METRICAS))
  mm <- aggregate(crudo, by = list(Clase = g), FUN = median)
  med[[f$col]] <- cbind(Factor = f$nombre, n = as.vector(table(g)[as.character(mm$Clase)]), mm)
}

perm <- do.call(rbind, perm); rownames(perm) <- NULL
perm$p_Holm <- p.adjust(perm$p, "holm")

anid <- do.call(rbind, anid); rownames(anid) <- NULL
anid$p_unidad_Holm <- p.adjust(anid$p_unidad_dentro, "holm")

kw <- do.call(rbind, kw); rownames(kw) <- NULL
kw$p_Holm <- ave(kw$p, kw$Factor, FUN = function(x) p.adjust(x, "holm"))
kw$PERMANOVA_sig <- kw$Factor %in% perm$Factor[perm$p_Holm < ALFA]

med  <- do.call(rbind, med); rownames(med) <- NULL
excl <- do.call(rbind, excl); rownames(excl) <- NULL

write.csv(perm, file.path(SALIDA, "A_permanova.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(anid, file.path(SALIDA, "A_anidada.csv"),   row.names = FALSE, fileEncoding = "UTF-8")
write.csv(kw,   file.path(SALIDA, "A_kruskal.csv"),   row.names = FALSE, fileEncoding = "UTF-8")
write.csv(med,  file.path(SALIDA, "A_medianas.csv"),  row.names = FALSE, fileEncoding = "UTF-8")
write.csv(excl, file.path(SALIDA, "A_unidades_excluidas.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# ------------------------- IMPRESIÓN ---------------------------------
cat("\n================ PERMANOVA de un factor + dispersión ================\n")
print(redondear(perm), row.names = FALSE)
cat("\n================ PERMANOVA secuencial: unidades dentro de su grupo ================\n")
print(redondear(anid), row.names = FALSE)
cat("\n================ Kruskal-Wallis por métrica ================\n")
print(redondear(kw), row.names = FALSE)
cat("\n================ Clases excluidas por el corte de", MIN_EV, "eventos ================\n")
if (nrow(excl)) {
  print(data.frame(Clases = tapply(excl$Eventos, excl$Factor, length),
                   Eventos = tapply(excl$Eventos, excl$Factor, sum)))
}
cat("\nMedianas por clase en A_medianas.csv\nArchivos guardados en:", SALIDA, "\n")