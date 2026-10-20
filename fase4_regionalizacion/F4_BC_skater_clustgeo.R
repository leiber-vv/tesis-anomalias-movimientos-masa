# =============================================================================
# FASE 4 - PARTE B (R): SKATER, CALIDAD DE LAS REGIONES Y ESTABILIDAD POR BOOTSTRAP
# -----------------------------------------------------------------------------
# 1. Lee las coronas (eventos_f4) y el grafo de la Parte A (grafo_f4_tabla).
# 2. Solucion de referencia: SKATER (spdep) para k = 2 a 15, con minimo de 20
#    eventos por region y criterio de suma de cuadrados (Assuncao et al., 2006).
#    Para cada k: pseudo-F, silhouette, R2 por variable, Ward sin restriccion
#    y, si se pide, acuerdo (ARI) con la solucion de ArcGIS de la Parte B1.
# 3. Bootstrap: en cada replica remuestrea los eventos, reconstruye el grafo con
#    las reglas de la Parte A (Delaunay, aristas <= umbral, puentes minimos) y
#    corre SKATER desde k = 2 hasta el primer k en que la region mas pequena
#    llega al tamano minimo. Es REANUDABLE: guarda cada replica al terminarla.
#    Si los datos o los parametros cambian, las replicas viejas se apartan solas
#    (se renombran) y el bootstrap empieza de cero. Para detenerlo: Esc en RStudio,
#    o crear DETENER.txt en la carpeta fase4 (termina la replica en curso).
# 4. Resumen: Jaccard por region (Hennig, 2007), ARI y estabilidad por evento;
#    aplica la regla declarada para elegir k, guarda el k elegido para las
#    Partes C y D (F4_R_k_elegido.csv) e imprime las cifras para el texto.
# =============================================================================

## ---------------- CONFIGURACION ----------------
CARPETA      <- "D:/Users/tasal/Desktop/Tesis de grado/fase4"
GDB          <- file.path(CARPETA, "fase4.gdb")
CAPA_EV      <- "eventos_f4"
TABLA_VEC    <- "grafo_f4_tabla"
COMPARAR_ARCGIS <- FALSE  # TRUE solo si se volvio a correr la Parte B1 en ArcGIS con los datos actuales
ETIQ_ARCGIS  <- file.path(CARPETA, "F4_B1_etiquetas.csv")   # solucion de ArcGIS (solo para comparar)
PARAMETROS_A <- file.path(CARPETA, "F4_A_parametros.csv")   # lo escribe la Parte A
VARIABLES    <- c("LOG_AREA", "LOG_HL", "LOG_ELONG",
                  "PENDIENTE", "LOG_DIST_ESTRUCTURA", "LOG_DIST_DRENAJE", "DENS_LIN_1000")
FISICAS      <- c("PENDIENTE", "LOG_DIST_ESTRUCTURA", "LOG_DIST_DRENAJE", "DENS_LIN_1000")
K_REF        <- 2:15      # solucion de referencia
K_BOOT       <- NULL      # NULL = de 2 al primer k en que la region minima llega a TAM_MIN
TAM_MIN      <- 20
UMBRAL_M     <- NULL      # NULL = el mismo de la Parte A (se lee de F4_A_parametros.csv)
B            <- 100       # replicas (valor por defecto de Hennig, 2007)
ESQUEMA      <- "clasico" # "clasico" = con reemplazo | "submuestreo" = 80 % sin reemplazo (sensibilidad)
FRACCION_SUB <- 0.8
SEMILLA      <- 2026
UMBRAL_JACCARD <- 0.75
SOLO_RESUMEN <- FALSE     # TRUE = no correr replicas, solo resumir lo guardado
## ------------------------------------------------

for (p in c("sf", "spdep", "deldir", "cluster")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({ library(sf); library(spdep); library(cluster) })

CHECKPOINT <- file.path(CARPETA, sprintf("F4_R_boot_replicas_%s.csv", ESQUEMA))
DETENER    <- file.path(CARPETA, "DETENER.txt")

## ---------------- funciones ----------------
# Criterio de SKATER: suma de cuadrados de las desviaciones a la media del grupo
ssd <- function(data, id) sum(scale(data[id, , drop = FALSE], scale = FALSE)^2)

# nb (lista de vecinos de spdep) a partir de pares de indices; se simetriza
nb_desde_pares <- function(N, a, b) {
  v <- split(c(b, a), factor(c(a, b), levels = seq_len(N)))
  nb <- lapply(v, function(x) if (length(x)) sort(unique(as.integer(x))) else 0L)
  names(nb) <- NULL
  class(nb) <- "nb"
  attr(nb, "region.id") <- as.character(seq_len(N))
  attr(nb, "sym") <- TRUE
  nb
}

# SKATER secuencial: un corte por paso, de k = 2 hasta max(ks); soluciones anidadas
skater_secuencia <- function(nb, Z, ks) {
  n <- nrow(Z)
  costos <- lapply(nbcosts(nb, Z), function(x) pmax(x, 1e-12))
  lw <- nb2listw(nb, glist = costos, style = "B")
  arbol <- suppressWarnings(mstree(lw))
  out <- list(); res <- NULL
  for (k in 2:max(ks)) {
    res <- tryCatch(
      skater(if (is.null(res)) arbol[, 1:2] else res, Z, ncuts = 1,
             crit = TAM_MIN, vec.crit = rep(1, n), method = ssd),
      error = function(e) NULL)
    if (is.null(res)) break
    if (k %in% ks) out[[as.character(k)]] <- res$groups
  }
  out
}

# Numeracion estable de las regiones: 1 = la mas grande (empates: menor ID). Asi los numeros
# no cambian entre corridas y los archivos de distintas corridas se pueden cruzar.
canonica <- function(g, ids) {
  tam <- table(g); minid <- tapply(ids, g, min)
  orden <- order(-as.numeric(tam), as.numeric(minid[names(tam)]))
  mapa <- setNames(seq_along(orden), names(tam)[orden])
  as.integer(mapa[as.character(g)])
}

pseudo_f_r2 <- function(Z, g) {
  g <- as.integer(factor(g)); G <- max(g); n <- nrow(Z)
  medias <- rowsum(Z, g) / as.vector(table(g))
  w_j <- colSums((Z - medias[g, , drop = FALSE])^2)
  t_j <- colSums(sweep(Z, 2, colMeans(Z))^2)
  list(f = ((sum(t_j) - sum(w_j)) / (G - 1)) / (sum(w_j) / (n - G)), r2 = 1 - w_j / t_j)
}

ari <- function(a, b) {
  tab <- table(a, b); comb <- function(x) x * (x - 1) / 2
  s_ij <- sum(comb(tab)); s_a <- sum(comb(rowSums(tab))); s_b <- sum(comb(colSums(tab)))
  esp <- s_a * s_b / comb(length(a)); mx <- (s_a + s_b) / 2
  if (mx == esp) 1 else (s_ij - esp) / (mx - esp)
}

# Reglas de la Parte A sobre los eventos distintos U: Delaunay, aristas <= UMBRAL_M
# y puentes minimos (Kruskal) hasta conectar el grafo. Devuelve pares de eventos.
grafo_eventos <- function(U) {
  xyU <- xy[U, , drop = FALSE]
  tri <- suppressWarnings(deldir::deldir(xyU[, 1], xyU[, 2]))
  ar <- unique(cbind(pmin(tri$delsgs$ind1, tri$delsgs$ind2), pmax(tri$delsgs$ind1, tri$delsgs$ind2)))
  # Coronas que la triangulacion deja fuera (misma coordenada que otra): se unen a la mas cercana,
  # con la misma regla de la Parte A
  dentro <- sort(unique(c(ar[, 1], ar[, 2])))
  fuera <- setdiff(seq_along(U), dentro)
  for (a in fuera) {
    d2 <- colSums((t(xyU[dentro, , drop = FALSE]) - xyU[a, ])^2)
    b <- dentro[which.min(d2)]
    ar <- rbind(ar, c(min(a, b), max(a, b)))
  }
  forzada <- seq_len(nrow(ar)) > (nrow(ar) - length(fuera))
  lg <- sqrt(rowSums((xyU[ar[, 1], , drop = FALSE] - xyU[ar[, 2], , drop = FALSE])^2))
  corta <- lg <= UMBRAL_M | forzada
  padre <- seq_along(U)
  raiz <- function(x) { while (padre[x] != x) { padre[x] <<- padre[padre[x]]; x <- padre[x] }; x }
  for (i in which(corta)) { ra <- raiz(ar[i, 1]); rb <- raiz(ar[i, 2]); if (ra != rb) padre[ra] <- rb }
  puente <- rep(FALSE, nrow(ar))
  cand <- which(!corta); cand <- cand[order(lg[cand])]
  for (i in cand) {
    ra <- raiz(ar[i, 1]); rb <- raiz(ar[i, 2])
    if (ra != rb) { padre[ra] <- rb; puente[i] <- TRUE }
  }
  sel <- ar[corta | puente, , drop = FALSE]
  list(pares = cbind(U[sel[, 1]], U[sel[, 2]]), n_puentes = sum(puente))
}

correr_replica <- function(b) {
  set.seed(SEMILLA + b)
  filas_idx <- sort(if (ESQUEMA == "clasico") sample.int(n, n, replace = TRUE)
                    else sample.int(n, round(FRACCION_SUB * n)))
  U <- sort(unique(filas_idx))
  gr <- grafo_eventos(U)
  N <- length(filas_idx)
  filas_de <- split(seq_len(N), filas_idx)          # filas (copias) de cada evento
  A <- integer(0); Bv <- integer(0)
  for (p in seq_len(nrow(gr$pares))) {              # cada copia hereda los vecinos del evento
    e <- expand.grid(filas_de[[as.character(gr$pares[p, 1])]], filas_de[[as.character(gr$pares[p, 2])]])
    A <- c(A, e[[1]]); Bv <- c(Bv, e[[2]])
  }
  for (rs in filas_de) if (length(rs) > 1) {        # las copias de un evento quedan unidas
    e <- expand.grid(rs, rs); e <- e[e[[1]] != e[[2]], ]
    A <- c(A, e[[1]]); Bv <- c(Bv, e[[2]])
  }
  sol <- skater_secuencia(nb_desde_pares(N, A, Bv), scale(X[filas_idx, , drop = FALSE]), K_BOOT)
  primera <- vapply(filas_de, `[`, integer(1), 1)
  ev_idx <- as.integer(names(filas_de))
  out <- do.call(rbind, lapply(names(sol), function(k)
    data.frame(b = b, k = as.integer(k), ID = ids[ev_idx], ETIQUETA = sol[[k]][primera])))
  list(out = out, n_puentes = gr$n_puentes)
}

# Huella de los datos y parametros: si cambia, las replicas guardadas no sirven
huella <- function() {
  tf <- tempfile(fileext = ".csv")
  write.csv(data.frame(ID = ids, X = round(xy[, 1], 2), Y = round(xy[, 2], 2), round(X, 6)), tf, row.names = FALSE)
  cat(paste("UMBRAL_M", UMBRAL_M, "TAM_MIN", TAM_MIN, "ESQUEMA", ESQUEMA, "FRACCION_SUB", FRACCION_SUB,
            "SEMILLA", SEMILLA, "K_BOOT", paste(K_BOOT, collapse = "-"),
            "VARIABLES", paste(VARIABLES, collapse = "-")), file = tf, append = TRUE)
  unname(tools::md5sum(tf))
}

resumir <- function(boot, etiq, ks) {
  tk <- list(); treg <- list(); tev <- list()
  for (k in ks) {
    sub <- boot[boot$k == k, ]
    if (!nrow(sub)) next
    refk <- setNames(etiq[[paste0("k", k)]], as.character(etiq$ID))
    jac <- list(); aris <- numeric(0); co_s <- numeric(0); co_n <- numeric(0)
    for (bb in unique(sub$b)) {
      g <- sub[sub$b == bb, ]
      lr <- unname(refk[as.character(g$ID)]); lb <- g$ETIQUETA
      aris <- c(aris, ari(lr, lb))
      for (r in unique(lr)) {
        A <- lr == r
        jac[[as.character(r)]] <- c(jac[[as.character(r)]],
                                    max(vapply(unique(lb), function(cc) sum(A & lb == cc) / sum(A | lb == cc), numeric(1))))
        idx <- which(A)
        if (length(idx) > 1) {
          cnt <- table(lb[idx])
          v <- (as.numeric(cnt[as.character(lb[idx])]) - 1) / (length(idx) - 1)
          id_c <- as.character(g$ID[idx])
          co_s[id_c] <- ifelse(is.na(co_s[id_c]), 0, co_s[id_c]) + v
          co_n[id_c] <- ifelse(is.na(co_n[id_c]), 0, co_n[id_c]) + 1
        }
      }
    }
    medias <- vapply(jac, mean, numeric(1))
    tam_ref <- table(refk)
    treg[[length(treg) + 1]] <- data.frame(
      k = k, region = as.integer(names(jac)), n_eventos = as.integer(tam_ref[names(jac)]),
      jaccard_medio = medias, jaccard_de = vapply(jac, sd, numeric(1)),
      disuelta_pct = vapply(jac, function(v) 100 * mean(v <= 0.5), numeric(1)),
      replicas = vapply(jac, length, integer(1)))
    tk[[length(tk) + 1]] <- data.frame(
      k = k, replicas = length(unique(sub$b)), jaccard_min = min(medias), jaccard_medio = mean(medias),
      regiones_menos_0.6 = sum(medias < 0.6),
      regiones_0.6_0.75 = sum(medias >= 0.6 & medias < UMBRAL_JACCARD), ARI_medio = mean(aris))
    tev[[length(tev) + 1]] <- data.frame(ID = as.integer(names(co_s)), k = k,
                                         estabilidad = co_s / co_n[names(co_s)], replicas = co_n[names(co_s)])
  }
  list(k = do.call(rbind, tk), regiones = do.call(rbind, treg), eventos = do.call(rbind, tev))
}

## ---------------- datos ----------------
cat(strrep("=", 78), "\nFASE 4 - PARTE B (R): SKATER, CALIDAD Y ESTABILIDAD\n", strrep("=", 78), "\n", sep = "")
cat("Decisiones registradas:\n",
    "  1. SKATER (spdep; Assuncao et al., 2006) con criterio de suma de cuadrados y minimo de ",
    TAM_MIN, " eventos por region,\n     sobre el grafo de la Parte A (un solo grafo para el valle).\n",
    "  2. Las 7 variables entran estandarizadas (z) y con igual peso.\n",
    "  3. Referencia: k = ", min(K_REF), " a ", max(K_REF), "; calidad con pseudo-F, silhouette, R2 por variable",
    " y Ward sin restriccion.\n",
    "  4. Bootstrap ", ESQUEMA, ": ", B, " replicas, desde k = 2 hasta el primer k con la region minima en ",
    TAM_MIN, " eventos;\n     el grafo se reconstruye en cada replica con las reglas de la Parte A.\n",
    "  5. Regla: el mayor k con todas las regiones >= ", TAM_MIN, " eventos y Jaccard medio >= ",
    UMBRAL_JACCARD, " (Hennig, 2007).\n", sep = "")

if (is.null(UMBRAL_M)) {
  if (!file.exists(PARAMETROS_A)) stop("No encuentro ", PARAMETROS_A, ": corre primero la Parte A actualizada.")
  par_a <- read.csv(PARAMETROS_A)
  UMBRAL_M <- as.numeric(par_a$valor[par_a$parametro == "UMBRAL_M"])
  if (as.numeric(par_a$valor[par_a$parametro == "RESPETAR_RIO"]) == 1)
    warning("La Parte A se corrio con RESPETAR_RIO = True; las replicas de R no filtran por el rio.")
}
ev  <- st_read(GDB, layer = CAPA_EV, quiet = TRUE)
xy  <- st_coordinates(ev)[, 1:2]
ev  <- st_drop_geometry(ev)
ids <- ev$ID
X   <- as.matrix(ev[, VARIABLES])
n   <- nrow(X)
tab <- suppressWarnings(st_read(GDB, layer = TABLA_VEC, quiet = TRUE))
nb  <- nb_desde_pares(n, match(tab$ID, ids), match(tab$NID, ids))
cat("\nDATOS\n  Eventos:", n, "| por ladera:", paste(names(table(ev$LADERA)), table(ev$LADERA), sep = ": ", collapse = " | "),
    "\n  Grafo de la Parte A:", sum(card(nb)) / 2, "aristas | componentes:", n.comp.nb(nb)$nc, "\n")

# Control: el grafo que R reconstruye con todos los eventos debe coincidir con el de la Parte A
g_all <- grafo_eventos(seq_len(n))
clave <- function(a, b) paste(pmin(a, b), pmax(a, b))
k_r <- unique(clave(ids[g_all$pares[, 1]], ids[g_all$pares[, 2]]))
k_a <- unique(clave(tab$ID, tab$NID))
comunes <- length(intersect(k_r, k_a))
cat("  Umbral de vecindad:", UMBRAL_M, "m\n")
cat("  Control del grafo (R vs Parte A): aristas R", length(k_r), "| Parte A", length(k_a),
    "| en comun", comunes, "| puentes R", g_all$n_puentes,
    if (comunes == length(k_a) && comunes == length(k_r)) "-> COINCIDEN\n" else "-> NO COINCIDEN\n")
if (comunes < 0.98 * max(length(k_a), length(k_r)))
  stop("El grafo de R difiere del de la Parte A: revisa que UMBRAL_M y la capa eventos_f4 sean los de la ultima corrida de la Parte A.")

## ---------------- solucion de referencia ----------------
Z  <- scale(X)
D  <- dist(Z)
hc <- hclust(D, method = "ward.D2")
t0 <- Sys.time()
set.seed(SEMILLA)
ref <- skater_secuencia(nb, Z, K_REF)
cat("\nSKATER de referencia:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")
arc <- NULL
if (COMPARAR_ARCGIS && file.exists(ETIQ_ARCGIS)) {
  arc <- read.csv(ETIQ_ARCGIS, fileEncoding = "UTF-8-BOM")
  if (!setequal(arc$ID, ids)) {
    cat("AVISO: F4_B1_etiquetas.csv tiene otros eventos (es de una corrida anterior); no se compara con ArcGIS.\n")
    arc <- NULL
  }
}
filas <- list()
etiq <- data.frame(ID = ids, LADERA = ev$LADERA)
for (k in K_REF) {
  g <- ref[[as.character(k)]]
  if (is.null(g)) next
  if (length(unique(g)) < k) {
    cat("  k =", k, ": SKATER no pudo hacer mas cortes con el minimo de", TAM_MIN, "eventos; se omite\n")
    next
  }
  g <- canonica(g, ids)
  tam <- table(g)
  pf  <- pseudo_f_r2(Z, g)
  s   <- mean(silhouette(as.integer(factor(g)), D)[, "sil_width"])
  sw  <- mean(silhouette(cutree(hc, k), D)[, "sil_width"])
  ac  <- paste0("k", k)
  ari_arc <- if (!is.null(arc) && ac %in% names(arc)) ari(g, arc[[ac]][match(ids, arc$ID)]) else NA
  filas[[length(filas) + 1]] <- data.frame(
    k = k, regiones = length(tam), tam_min = min(tam), regiones_menores_20 = sum(tam < TAM_MIN),
    regiones_ambas_laderas = sum(tapply(ev$LADERA, g, function(x) length(unique(x)) > 1)),
    pseudo_F = pf$f, silhouette = s, silhouette_ward = sw, costo_restriccion = sw - s,
    R2_fisicas = mean(pf$r2[FISICAS]), R2_respuesta = mean(pf$r2[setdiff(VARIABLES, FISICAS)]),
    ARI_vs_ArcGIS = ari_arc, as.list(setNames(pf$r2, paste0("R2_", VARIABLES))))
  etiq[[ac]] <- g
}
res_ref <- do.call(rbind, filas)
write.csv(res_ref, file.path(CARPETA, "F4_R_referencia_k.csv"), row.names = FALSE)
write.csv(etiq, file.path(CARPETA, "F4_R_etiquetas.csv"), row.names = FALSE)
cat("\nRESULTADOS POR k (referencia, spdep)\n")
print(transform(res_ref[, c("k", "regiones", "tam_min", "regiones_menores_20", "regiones_ambas_laderas",
                            "pseudo_F", "silhouette", "silhouette_ward", "R2_fisicas", "R2_respuesta",
                            "ARI_vs_ArcGIS")],
                pseudo_F = round(pseudo_F, 1), silhouette = round(silhouette, 3),
                silhouette_ward = round(silhouette_ward, 3), R2_fisicas = round(R2_fisicas, 2),
                R2_respuesta = round(R2_respuesta, 2), ARI_vs_ArcGIS = round(ARI_vs_ArcGIS, 2)),
      row.names = FALSE)
cat("\nR2 POR VARIABLE\n")
print(round(res_ref[, c("k", paste0("R2_", VARIABLES))], 2), row.names = FALSE)
if (is.null(K_BOOT)) {
  k_tope <- res_ref$k[res_ref$tam_min <= TAM_MIN]
  K_BOOT <- 2:(if (length(k_tope)) min(k_tope) else max(res_ref$k))
}
cat("\nk que pasan al bootstrap:", paste(range(K_BOOT), collapse = " a "),
    "(desde ese k la region mas pequena tiene", TAM_MIN, "eventos)\n")

## ---------------- bootstrap ----------------
HUELLA_ARCH <- sub("\\.csv$", "_huella.txt", CHECKPOINT)
h_actual <- huella()
if (file.exists(CHECKPOINT)) {
  h_guardada <- if (file.exists(HUELLA_ARCH)) readLines(HUELLA_ARCH, warn = FALSE)[1] else ""
  if (!identical(h_guardada, h_actual)) {
    viejo <- sub("\\.csv$", format(Sys.time(), "_anterior_%Y%m%d_%H%M%S.csv"), CHECKPOINT)
    file.rename(CHECKPOINT, viejo)
    cat("\nAVISO: las replicas guardadas eran de otros datos o parametros. Se apartaron en\n ",
        basename(viejo), "\n  y el bootstrap empieza de cero.\n")
  }
}
writeLines(h_actual, HUELLA_ARCH)
if (SOLO_RESUMEN && !file.exists(CHECKPOINT))
  stop("No hay replicas de los datos actuales: corre con SOLO_RESUMEN <- FALSE.")
if (!SOLO_RESUMEN) {
  if (file.exists(DETENER)) file.remove(DETENER)
  hechas <- if (file.exists(CHECKPOINT)) unique(read.csv(CHECKPOINT)$b) else integer(0)
  pendientes <- setdiff(seq_len(B), hechas)
  cat("\nBOOTSTRAP | replicas guardadas:", length(hechas), "| pendientes:", length(pendientes), "\n")
  t0 <- Sys.time()
  for (j in seq_along(pendientes)) {
    if (file.exists(DETENER)) { cat("DETENIDO a pedido (DETENER.txt). Al volver a correr, continua.\n"); break }
    b <- pendientes[j]
    r <- correr_replica(b)
    if (!is.null(r$out) && nrow(r$out))
      write.table(r$out, CHECKPOINT, sep = ",", row.names = FALSE,
                  col.names = !file.exists(CHECKPOINT), append = file.exists(CHECKPOINT))
    seg <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("  replica %3d lista | puentes: %d | %.1f s por replica | faltan ~%.0f min\n",
                b, r$n_puentes, seg / j, seg / j * (length(pendientes) - j) / 60))
    flush.console()
  }
}

## ---------------- resumen ----------------
if (file.exists(CHECKPOINT)) {
  boot <- read.csv(CHECKPOINT)
  rs <- resumir(boot, etiq, K_BOOT)
  tk <- merge(rs$k, res_ref[, c("k", "regiones", "regiones_menores_20", "pseudo_F", "silhouette")], by = "k")
  tk$cumple_regla <- tk$regiones_menores_20 == 0 & tk$jaccard_min >= UMBRAL_JACCARD
  write.csv(tk, file.path(CARPETA, sprintf("F4_R_estabilidad_k_%s.csv", ESQUEMA)), row.names = FALSE)
  write.csv(rs$regiones, file.path(CARPETA, sprintf("F4_R_jaccard_regiones_%s.csv", ESQUEMA)), row.names = FALSE)
  write.csv(rs$eventos, file.path(CARPETA, sprintf("F4_R_estabilidad_eventos_%s.csv", ESQUEMA)), row.names = FALSE)
  cat("\nESTABILIDAD POR k\n")
  print(transform(tk, jaccard_min = round(jaccard_min, 2), jaccard_medio = round(jaccard_medio, 2),
                  ARI_medio = round(ARI_medio, 2), pseudo_F = round(pseudo_F, 1),
                  silhouette = round(silhouette, 3)), row.names = FALSE)
  if (any(tk$cumple_regla)) {
    k_fin <- max(tk$k[tk$cumple_regla])
    cat("\nk ELEGIDO POR LA REGLA:", k_fin, "\n")
    rr <- rs$regiones[rs$regiones$k == k_fin, ]; print(round(rr[order(rr$region), ], 2), row.names = FALSE)
  } else {
    k_fin <- tk$k[which.max(tk$jaccard_min)]
    cat("\nNingun k cumple la regla. Criterio de respaldo (definido despues de ver los resultados):\n",
        "el k con mayor Jaccard minimo -> k =", k_fin, "\n")
    rr <- rs$regiones[rs$regiones$k == k_fin, ]; print(round(rr[order(rr$region), ], 2), row.names = FALSE)
  }
  criterio <- if (any(tk$cumple_regla)) "regla" else "respaldo (mayor Jaccard minimo)"
  if (ESQUEMA == "clasico")
    write.csv(data.frame(k = k_fin, criterio = criterio, esquema = ESQUEMA),
              file.path(CARPETA, "F4_R_k_elegido.csv"), row.names = FALSE)
  
  # ---- cifras para el texto (las que usa la seccion de resultados)
  rk <- res_ref[res_ref$k == k_fin, ]; tkk <- tk[tk$k == k_fin, ]
  reg <- rs$regiones[rs$regiones$k == k_fin, ]; reg <- reg[order(reg$region), ]
  est_k <- rs$eventos[rs$eventos$k == k_fin, ]
  cat_h <- cut(reg$jaccard_medio, c(-Inf, 0.6, UMBRAL_JACCARD, Inf), right = FALSE,
               labels = c("no confiable", "limites dudosos", "estable"))
  pf_loc <- rk$pseudo_F >= max(res_ref$pseudo_F[abs(res_ref$k - k_fin) == 1])
  r2v <- unlist(rk[paste0("R2_", VARIABLES)]); names(r2v) <- VARIABLES
  cat("\nCIFRAS PARA EL TEXTO (k =", k_fin, "; criterio:", criterio, ")\n")
  cat(sprintf("  Silhouette SKATER: maxima %.3f (k = %d); negativa desde k = %s | Ward: maxima %.3f\n",
              max(res_ref$silhouette), res_ref$k[which.max(res_ref$silhouette)],
              if (any(res_ref$silhouette < 0)) min(res_ref$k[res_ref$silhouette < 0]) else "-",
              max(res_ref$silhouette_ward)))
  cat(sprintf("  Costo de la restriccion con k = %d: %.3f (silhouette Ward %.3f vs SKATER %.3f)\n",
              k_fin, rk$costo_restriccion, rk$silhouette_ward, rk$silhouette))
  cat(sprintf("  R2 medio: fisicas %.2f | comportamiento %.2f\n", rk$R2_fisicas, rk$R2_respuesta))
  cat("  R2 por variable (de mayor a menor):",
      paste(sprintf("%s %.2f", names(sort(r2v, decreasing = TRUE)), sort(r2v, decreasing = TRUE)), collapse = " | "), "\n")
  cat(sprintf("  Pseudo-F: %.1f (%s maximo local) | Jaccard minimo %.2f | medio %.2f | ARI %.2f\n",
              rk$pseudo_F, if (pf_loc) "es" else "no es", tkk$jaccard_min, tkk$jaccard_medio, tkk$ARI_medio))
  cat(sprintf("  ARI con k = 2: %.2f\n", tk$ARI_medio[tk$k == 2]))
  cat("  Regiones:", paste(sprintf("%d (%d ev.; Jaccard %.2f; %s; disuelta en %.0f %%)", reg$region, reg$n_eventos,
                                   reg$jaccard_medio, cat_h, reg$disuelta_pct), collapse = "\n           "), "\n")
  cat(sprintf("  Eventos con estabilidad >= 0,75: %.0f %%\n", 100 * mean(est_k$estabilidad >= 0.75)))
  cat("  Regiones con eventos de ambas laderas:", rk$regiones_ambas_laderas, "de", rk$regiones, "\n")
  
  cat("\nGuardado en", CARPETA, ": F4_R_referencia_k, F4_R_etiquetas, F4_R_estabilidad_k,",
      "F4_R_jaccard_regiones, F4_R_estabilidad_eventos (.csv) y F4_R_k_elegido.csv\n")
}
cat(strrep("=", 78), "\n")

# =============================================================================
# FASE 4 - PARTE C (R): CLUSTGEO COMO CONTRASTE DE SKATER
# -----------------------------------------------------------------------------
# ClustGeo (Chavent et al., 2018) agrupa con una restriccion espacial BLANDA:
# mezcla una distancia en atributos (D0) y una distancia espacial (D1) con un
# peso alfa. Aqui D1 sale del MISMO grafo de la Parte A (1 - adyacencia), asi
# que la unica diferencia con SKATER es restriccion dura frente a blanda.
# 1. Calcula Q0 y Q1 (choicealpha) para alfa = 0 a 1 y elige alfa con la regla
#    declarada: el alfa que maximiza Q0norm + Q1norm (mejor compromiso entre
#    homogeneidad en atributos y cohesion de vecindad).
# 2. Para cada alfa: ARI con SKATER, % de eventos en fragmentos no contiguos,
#    silhouette y regiones con eventos de ambas laderas.
# 3. Detalle del alfa elegido: tabla cruzada SKATER x ClustGeo.
# =============================================================================

ETIQ_SKATER <- file.path(CARPETA, "F4_R_etiquetas.csv")   # solucion de referencia de la Parte B

K           <- 5
RANGO_ALFA  <- seq(0, 1, 0.1)


for (p in c("sf", "spdep", "cluster", "ClustGeo")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({ library(sf); library(spdep); library(cluster); library(ClustGeo) })

if (is.null(K)) {
  arch_k <- file.path(CARPETA, "F4_R_k_elegido.csv")
  if (!file.exists(arch_k)) stop("No encuentro F4_R_k_elegido.csv: corre primero la Parte B (esquema clasico).")
  K <- read.csv(arch_k)$k[1]
}

ari <- function(a, b) {
  tab <- table(a, b); comb <- function(x) x * (x - 1) / 2
  s_ij <- sum(comb(tab)); s_a <- sum(comb(rowSums(tab))); s_b <- sum(comb(colSums(tab)))
  esp <- s_a * s_b / comb(length(a)); mx <- (s_a + s_b) / 2
  if (mx == esp) 1 else (s_ij - esp) / (mx - esp)
}

# Eventos que quedan fuera de la parte conectada mas grande de su region (segun el grafo)
fragmentados <- function(g, nb) {
  frag <- rep(FALSE, length(g))
  for (r in unique(g)) {
    idx <- which(g == r); en <- g == r
    comp <- integer(length(g)); cc <- 0L
    for (s in idx) if (comp[s] == 0L) {
      cc <- cc + 1L; comp[s] <- cc; pila <- s
      while (length(pila)) {
        v <- pila[length(pila)]; pila <- pila[-length(pila)]
        for (w in nb[[v]]) if (w > 0 && en[w] && comp[w] == 0L) { comp[w] <- cc; pila <- c(pila, w) }
      }
    }
    mayor <- which.max(tabulate(comp[idx]))
    frag[idx[comp[idx] != mayor]] <- TRUE
  }
  frag
}

cat(strrep("=", 78), "\nFASE 4 - PARTE C (R): CLUSTGEO COMO CONTRASTE DE SKATER\n", strrep("=", 78), "\n", sep = "")
cat("Decisiones registradas:\n",
    "  1. ClustGeo (Chavent et al., 2018) con k = ", K, ", el mismo de SKATER.\n",
    "  2. D0: distancia euclidiana con las 7 variables estandarizadas (las mismas de SKATER).\n",
    "  3. D1: 1 - adyacencia del grafo de la Parte A (restriccion blanda sobre la misma vecindad).\n",
    "  4. Regla para alfa (declarada antes de ver resultados): el alfa que maximiza Q0norm + Q1norm.\n",
    "  5. Comparacion con SKATER: ARI, % de eventos en fragmentos no contiguos y silhouette.\n", sep = "")

## ---------------- datos ----------------
ev  <- st_drop_geometry(st_read(GDB, layer = CAPA_EV, quiet = TRUE))
ids <- ev$ID
n   <- nrow(ev)
Z   <- scale(as.matrix(ev[, VARIABLES]))
tab <- suppressWarnings(st_read(GDB, layer = TABLA_VEC, quiet = TRUE))
a <- match(tab$ID, ids); b <- match(tab$NID, ids)
nb <- lapply(split(c(b, a), factor(c(a, b), levels = seq_len(n))), function(x) sort(unique(as.integer(x))))
names(nb) <- NULL; class(nb) <- "nb"; attr(nb, "region.id") <- as.character(ids); attr(nb, "sym") <- TRUE
etq <- read.csv(ETIQ_SKATER)
if (!setequal(etq$ID, ids))
  stop("F4_R_etiquetas.csv no corresponde a los eventos actuales de eventos_f4: corre primero la Parte B.")
skt <- etq[[paste0("k", K)]][match(ids, etq$ID)]

D0 <- dist(Z)
A  <- nb2mat(nb, style = "B"); diag(A) <- 1
D1 <- as.dist(1 - A)
cat("\nDATOS\n  Eventos:", n, "| aristas del grafo:", sum(card(nb)) / 2,
    "| SKATER k =", K, "con", length(unique(skt)), "regiones\n")
cat("  Control: eventos de SKATER en fragmentos no contiguos =", sum(fragmentados(skt, nb)), "(debe ser 0)\n")

## ---------------- eleccion de alfa ----------------
cr <- choicealpha(D0, D1, range.alpha = RANGO_ALFA, K = K, graph = FALSE)
Q  <- as.matrix(cr$Q); Qn <- as.matrix(cr$Qnorm)
filas <- list(); grupos <- list()
for (i in seq_along(RANGO_ALFA)) {
  g <- cutree(hclustgeo(D0, D1, alpha = RANGO_ALFA[i]), K)
  grupos[[i]] <- g
  fr <- fragmentados(g, nb)
  filas[[i]] <- data.frame(
    alfa = RANGO_ALFA[i], Q0 = Q[i, 1], Q1 = Q[i, 2], Q0norm = Qn[i, 1], Q1norm = Qn[i, 2],
    suma_norm = Qn[i, 1] + Qn[i, 2], ARI_vs_SKATER = ari(g, skt),
    pct_fragmentados = 100 * mean(fr), tam_min = min(table(g)),
    silhouette = mean(silhouette(g, D0)[, "sil_width"]),
    regiones_ambas_laderas = sum(tapply(ev$LADERA, g, function(x) length(unique(x)) > 1)))
}
tabla <- do.call(rbind, filas)
i_sel <- which.max(tabla$suma_norm)
alfa  <- tabla$alfa[i_sel]
cg    <- grupos[[i_sel]]
s_skt <- mean(silhouette(skt, D0)[, "sil_width"])

cat("\nCLUSTGEO PARA CADA ALFA (k =", K, ")\n")
print(transform(tabla, Q0 = round(Q0, 3), Q1 = round(Q1, 3), Q0norm = round(Q0norm, 2),
                Q1norm = round(Q1norm, 2), suma_norm = round(suma_norm, 2),
                ARI_vs_SKATER = round(ARI_vs_SKATER, 2), pct_fragmentados = round(pct_fragmentados, 1),
                silhouette = round(silhouette, 3)), row.names = FALSE)
cat("\nALFA ELEGIDO POR LA REGLA:", alfa, "\n")
cat("  ARI ClustGeo vs SKATER:", round(tabla$ARI_vs_SKATER[i_sel], 2),
    "| eventos en fragmentos no contiguos:", sum(fragmentados(cg, nb)),
    sprintf("(%.1f %%)", tabla$pct_fragmentados[i_sel]),
    "\n  Silhouette: ClustGeo", round(tabla$silhouette[i_sel], 3), "| SKATER", round(s_skt, 3), "\n")
cat("\nTABLA CRUZADA (filas: SKATER, columnas: ClustGeo)\n")
print(table(SKATER = skt, ClustGeo = cg))

write.csv(tabla, file.path(CARPETA, sprintf("F4_R_clustgeo_alfa_k%d.csv", K)), row.names = FALSE)
write.csv(data.frame(ID = ids, SKATER = skt, CLUSTGEO = cg, FRAGMENTO = as.integer(fragmentados(cg, nb))),
          file.path(CARPETA, sprintf("F4_R_clustgeo_k%d.csv", K)), row.names = FALSE)
cat("\nGuardado en ", CARPETA, ": F4_R_clustgeo_alfa_k", K, ".csv y F4_R_clustgeo_k", K, ".csv\n", sep = "")
cat(strrep("=", 78), "\n")