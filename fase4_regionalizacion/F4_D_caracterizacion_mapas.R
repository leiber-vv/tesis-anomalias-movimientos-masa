# =============================================================================
# FASE 4 - PARTE D (R): CARACTERIZACION DE LAS REGIONES Y MAPAS
# -----------------------------------------------------------------------------
# 1. Caracteriza las regiones de SKATER (k elegido en la Parte B): tamano y laderas, estabilidad
#    del bootstrap clasico (Jaccard y categoria de Hennig), silhouette, perfil
#    de las 7 variables, unidades geologicas, geomorfologicas y de cobertura
#    dominantes, y composicion por tipo, fuente, anomalia y completitud.
# 2. V de Cramer corregida (Bergsma, 2013) entre regiones y cada variable
#    categorica, para SKATER y para ClustGeo. Descriptiva: sin valor p, porque
#    regiones y unidades son contiguas por construccion.
# 3. Mapas (PNG) de SKATER y ClustGeo: puntos y poligonos de Thiessen recortados
#    a 1 km alrededor de los eventos, mas una comparacion lado a lado. Los colores
#    de ClustGeo se alinean con la region de SKATER con la que mas coincide.
# 4. Guarda todo en CSV y las capas en F4_regiones.gpkg (para armar figuras en ArcGIS).
# =============================================================================

## ---------------- CONFIGURACION ----------------
CARPETA     <- "D:/Users/tasal/Desktop/Tesis de grado/fase4"
GDB         <- file.path(CARPETA, "fase4.gdb")
CAPA_EV     <- "eventos_f4"
RIO_DSN     <- file.path(CARPETA, "fase4.gdb")   # archivo donde esta el eje del rio (el mismo de la Parte A)
RIO_CAPA    <- "Rio_medellin"
K           <- NULL      # NULL = el k elegido en la Parte B (F4_R_k_elegido.csv); o un numero para forzarlo
ETIQ_SKATER <- file.path(CARPETA, "F4_R_etiquetas.csv")
JACCARD     <- file.path(CARPETA, "F4_R_jaccard_regiones_clasico.csv")
ESTAB_EV    <- file.path(CARPETA, "F4_R_estabilidad_eventos_clasico.csv")
VALLE_DSN    <- "D:/Users/tasal/Desktop/Tesis de grado/fase2/area de estudio.shp"      # archivo (.shp, .gpkg o .gdb) con las margenes del valle: la misma capa de los mapas anteriores
VALLE_CAPA   <- NULL      # nombre de la capa dentro de ese archivo (para .shp puede quedar en NULL)
VENTANA_DSN  <- NULL      # opcional: archivo (.gdb o .gpkg) con la ventana de laderas >= 10 grados
VENTANA_CAPA <- NULL      # opcional: nombre de esa capa de poligonos
RADIO_THIESSEN <- 1000    # recorte de los poligonos de Thiessen (escala de vecindad del grafo)
VARIABLES   <- c("LOG_AREA", "LOG_HL", "LOG_ELONG",
                 "PENDIENTE", "LOG_DIST_ESTRUCTURA", "LOG_DIST_DRENAJE", "DENS_LIN_1000")
CRUCE       <- c("GEOLOGIA", "GEOMORFOLOGIA", "COBERTURA")
COMPOSICION <- c("TIPO", "GRUPO", "ANOMALIA", "BAJO_UMBRAL")
ETIQ_VAR    <- c(LOG_AREA = "log \u00e1rea", LOG_HL = "log H/L", LOG_ELONG = "log elongaci\u00f3n",
                 PENDIENTE = "Pendiente", LOG_DIST_ESTRUCTURA = "log dist. a estructuras",
                 LOG_DIST_DRENAJE = "log dist. a drenajes", DENS_LIN_1000 = "Densidad de lineamientos")
COLORES     <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2", "#D55E00", "#F0E442")
## ------------------------------------------------

for (p in c("sf", "cluster", "ggplot2")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
suppressPackageStartupMessages({ library(sf); library(cluster); library(ggplot2) })
sf_use_s2(FALSE)

## ---------------- funciones ----------------
canonica <- function(g, ids) {          # 1 = region mas grande (empates: menor ID)
  tam <- table(g); minid <- tapply(ids, g, min)
  orden <- order(-as.numeric(tam), as.numeric(minid[names(tam)]))
  mapa <- setNames(seq_along(orden), names(tam)[orden])
  as.integer(mapa[as.character(g)])
}

cramer_v <- function(x, y) {            # V de Cramer con correccion de sesgo (Bergsma, 2013)
  ok <- !is.na(x) & !is.na(y) & x != "" & y != ""
  tab <- table(x[ok], y[ok]); tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  n <- sum(tab); r <- nrow(tab); k <- ncol(tab)
  if (r < 2 || k < 2) return(NA)
  chi2 <- suppressWarnings(chisq.test(tab, correct = FALSE)$statistic)
  phi2c <- max(0, chi2 / n - (k - 1) * (r - 1) / (n - 1))
  rc <- r - (r - 1)^2 / (n - 1); kc <- k - (k - 1)^2 / (n - 1)
  as.numeric(sqrt(phi2c / min(kc - 1, rc - 1)))
}

hennig <- function(j) ifelse(j >= 0.75, "estable", ifelse(j >= 0.6, "patron, limites dudosos", "no confiable"))

permutaciones <- function(v) {
  if (length(v) <= 1) return(list(v))
  do.call(c, lapply(seq_along(v), function(i) lapply(permutaciones(v[-i]), function(p) c(v[i], p))))
}

alinear <- function(cg, skt) {          # renumera ClustGeo segun la region de SKATER con mayor coincidencia
  tab <- table(cg, skt)
  if (nrow(tab) != ncol(tab) || nrow(tab) > 8) return(cg)
  mejor <- NULL; mx <- -1
  for (p in permutaciones(seq_len(ncol(tab)))) {
    s <- sum(tab[cbind(seq_len(nrow(tab)), p)])
    if (s > mx) { mx <- s; mejor <- p }
  }
  mapa <- setNames(as.integer(colnames(tab)[mejor]), rownames(tab))
  as.integer(mapa[as.character(cg)])
}

thiessen_por_region <- function(pts, campo, area) {
  env <- st_as_sfc(st_bbox(st_buffer(pts, 2 * RADIO_THIESSEN)))
  v <- st_sf(geometry = st_collection_extract(st_voronoi(st_union(pts), envelope = env), "POLYGON"))
  v$fila <- vapply(st_intersects(v, pts), function(x) if (length(x)) x[1] else NA_integer_, integer(1))
  v <- v[!is.na(v$fila), ]
  v$REGION <- pts[[campo]][v$fila]
  v <- suppressWarnings(st_intersection(v, area))
  aggregate(v["REGION"], by = list(REGION_ID = v$REGION), FUN = min)
}

tema <- theme_minimal(base_size = 10) +
  theme(axis.text = element_text(size = 7), panel.grid = element_line(color = "grey92"))
tema_mapa <- theme_void(base_size = 10) +
  theme(plot.title = element_blank(), plot.subtitle = element_blank(),
        panel.border = element_rect(fill = NA, colour = "#B9B5AB", linewidth = 0.6),
        legend.text = element_text(size = 11), legend.title = element_text(size = 12),
        legend.key.size = unit(0.6, "cm"),
        plot.margin = margin(6, 6, 6, 6), strip.text = element_text(face = "bold", size = 10))

# Muestra la figura en el panel Plots de RStudio y la guarda en PNG
mostrar_guardar <- function(g, archivo, ancho, alto) {
  print(g)
  ggsave(file.path(CARPETA, archivo), g, width = ancho, height = alto, dpi = 300)
}

## ---------------- datos ----------------
cat(strrep("=", 78), "\nFASE 4 - PARTE D (R): CARACTERIZACION DE LAS REGIONES Y MAPAS\n", strrep("=", 78), "\n", sep = "")
ev  <- st_read(GDB, layer = CAPA_EV, quiet = TRUE)
ids <- ev$ID
n   <- nrow(ev)
rio <- tryCatch(st_transform(st_zm(st_read(RIO_DSN, layer = RIO_CAPA, quiet = TRUE)), st_crs(ev)),
                error = function(e) NULL)
if (is.null(rio)) cat("ATENCION: no se encontro la capa del rio (", RIO_CAPA, "en", RIO_DSN,
                      "); los mapas saldran sin el rio. Revisa RIO_DSN y RIO_CAPA.\n")
valle <- NULL
if (!is.null(VALLE_DSN)) {
  valle <- if (is.null(VALLE_CAPA)) st_read(VALLE_DSN, quiet = TRUE) else st_read(VALLE_DSN, layer = VALLE_CAPA, quiet = TRUE)
  valle <- st_transform(st_zm(valle), st_crs(ev))
}
if (is.null(K)) {
  arch_k <- file.path(CARPETA, "F4_R_k_elegido.csv")
  if (!file.exists(arch_k)) stop("No encuentro F4_R_k_elegido.csv: corre primero la Parte B (esquema clasico).")
  K <- read.csv(arch_k)$k[1]
}
CLUSTGEO <- file.path(CARPETA, sprintf("F4_R_clustgeo_k%d.csv", K))
etq <- read.csv(ETIQ_SKATER)
cgf <- read.csv(CLUSTGEO)
if (!setequal(etq$ID, ids) || !setequal(cgf$ID, ids))
  stop("Los resultados de las Partes B o C no corresponden a los eventos actuales: corre primero B y C.")
ev$SKATER    <- canonica(etq[[paste0("k", K)]][match(ids, etq$ID)], ids)
ev$CLUSTGEO  <- alinear(cgf$CLUSTGEO[match(ids, cgf$ID)], ev$SKATER)
ev$FRAGMENTO <- cgf$FRAGMENTO[match(ids, cgf$ID)]
est <- read.csv(ESTAB_EV); est <- est[est$k == K, ]
ev$ESTABILIDAD <- est$estabilidad[match(ids, est$ID)]
Z <- scale(as.matrix(st_drop_geometry(ev)[, VARIABLES]))
D0 <- dist(Z)
cat("Eventos:", n, "| regiones SKATER:", length(unique(ev$SKATER)),
    "| regiones ClustGeo:", length(unique(ev$CLUSTGEO)), "\n")

## ---------------- caracterizacion ----------------
evd <- st_drop_geometry(ev)
reg <- sort(unique(evd$SKATER))
tam <- as.integer(table(factor(evd$SKATER, levels = reg)))
car <- data.frame(REGION = reg, eventos = tam,
                  ladera_1 = as.integer(table(factor(evd$SKATER[evd$LADERA == 1], levels = reg))),
                  ladera_2 = as.integer(table(factor(evd$SKATER[evd$LADERA == 2], levels = reg))))
jac <- read.csv(JACCARD); jac <- jac[jac$k == K, ]
car$jaccard <- jac$jaccard_medio[match(car$REGION, jac$region)]   # misma numeracion: 1 = la mas grande
if (!all(jac$n_eventos[match(car$REGION, jac$region)] == car$eventos))
  stop("El Jaccard guardado no corresponde a estas regiones: vuelve a correr la Parte B.")
car$categoria_hennig <- hennig(car$jaccard)
car$estabilidad_eventos <- as.numeric(tapply(evd$ESTABILIDAD, evd$SKATER, mean, na.rm = TRUE)[as.character(reg)])
sil <- silhouette(evd$SKATER, D0)
car$silhouette <- as.numeric(tapply(sil[, "sil_width"], evd$SKATER, mean)[as.character(reg)])
for (v in VARIABLES) car[[paste0("mediana_", v)]] <- as.numeric(tapply(evd[[v]], evd$SKATER, median)[as.character(reg)])
for (v in CRUCE) if (v %in% names(evd)) {
  dom <- t(sapply(reg, function(r) {
    t <- sort(table(evd[[v]][evd$SKATER == r]), decreasing = TRUE)
    c(names(t)[1], round(100 * t[1] / sum(t)))
  }))
  car[[paste0(v, "_dominante")]] <- dom[, 1]
  car[[paste0(v, "_pct")]] <- as.numeric(dom[, 2])
}
comp <- list()
for (v in COMPOSICION) if (v %in% names(evd)) {
  x <- ifelse(is.na(evd[[v]]) | evd[[v]] == "", "(vacio)", as.character(evd[[v]]))
  tb <- table(factor(evd$SKATER, levels = reg), x)
  pct <- round(100 * prop.table(tb, 1), 1)
  for (cat_ in colnames(pct)) car[[paste0(v, "_", cat_, "_pct")]] <- as.numeric(pct[, cat_])
  comp[[v]] <- rbind(pct, Valle = round(100 * prop.table(table(x)), 1)[colnames(pct)])
}
write.csv(car, file.path(CARPETA, sprintf("F4_D_caracterizacion_k%d.csv", K)), row.names = FALSE)

cat("\nCARACTERIZACION DE LAS REGIONES SKATER (k =", K, "; 1 = la mas grande)\n")
print(transform(car[, c("REGION", "eventos", "ladera_1", "ladera_2", "jaccard", "categoria_hennig",
                        "estabilidad_eventos", "silhouette")],
                jaccard = round(jaccard, 2), estabilidad_eventos = round(estabilidad_eventos, 2),
                silhouette = round(silhouette, 3)), row.names = FALSE)
cat("\nMEDIANAS DE LAS VARIABLES (en las unidades guardadas en eventos_f4)\n")
med <- car[, c("REGION", paste0("mediana_", VARIABLES))]; names(med) <- c("REGION", VARIABLES)
med <- rbind(med, data.frame(REGION = "Valle", as.list(sapply(evd[, VARIABLES], median))))
print(cbind(REGION = med$REGION, round(med[, VARIABLES], 2)), row.names = FALSE)
cat("\nUNIDADES DOMINANTES (% de los eventos de la region)\n")
for (v in CRUCE) if (paste0(v, "_dominante") %in% names(car))
  cat(sprintf("  %-14s %s\n", v, paste(sprintf("R%d: %s (%d %%)", car$REGION, car[[paste0(v, "_dominante")]],
                                                car[[paste0(v, "_pct")]]), collapse = " | ")))
for (v in names(comp)) { cat("\nCOMPOSICION POR", v, "(% por region)\n"); print(comp[[v]]) }

## ---------------- V de Cramer ----------------
vc <- data.frame(variable = c(CRUCE, COMPOSICION)[c(CRUCE, COMPOSICION) %in% names(evd)])
vc$V_SKATER   <- sapply(vc$variable, function(v) cramer_v(evd$SKATER, evd[[v]]))
vc$V_CLUSTGEO <- sapply(vc$variable, function(v) cramer_v(evd$CLUSTGEO, evd[[v]]))
write.csv(vc, file.path(CARPETA, sprintf("F4_D_cramer_v_k%d.csv", K)), row.names = FALSE)
cat("\nV DE CRAMER CORREGIDA (0 = sin relacion, 1 = correspondencia total; descriptiva)\n")
print(transform(vc, V_SKATER = round(V_SKATER, 2), V_CLUSTGEO = round(V_CLUSTGEO, 2)), row.names = FALSE)

## ---------------- perfil de las regiones ----------------
perfil <- do.call(rbind, lapply(reg, function(r)
  data.frame(REGION = paste("Regi\u00f3n", r), variable = VARIABLES,
             z = colMeans(Z[evd$SKATER == r, , drop = FALSE]))))
perfil$variable <- factor(ETIQ_VAR[perfil$variable], levels = rev(ETIQ_VAR[VARIABLES]))
g_perfil <- ggplot(perfil, aes(REGION, variable, fill = z)) + geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", z)), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       name = "Media z") +
  labs(x = NULL, y = NULL, title = sprintf("Perfil de las regiones de SKATER (k = %d)", K),
       subtitle = "Media z de cada variable; positivo = sobre el promedio del valle") + tema
mostrar_guardar(g_perfil, "F4_D_perfil_regiones.png", 8, 4.5)

## ---------------- mapas ----------------
## Despues de la primera corrida, esta seccion se puede volver a correr sola (seleccionarla y
## Ctrl + Enter) para ajustar colores, tamanos o textos y ver el cambio en el panel Plots.
area <- st_union(st_buffer(ev, RADIO_THIESSEN))
if (!is.null(VENTANA_DSN) && !is.null(VENTANA_CAPA)) {
  ven <- st_transform(st_read(VENTANA_DSN, layer = VENTANA_CAPA, quiet = TRUE), st_crs(ev))
  area <- st_intersection(area, st_union(st_make_valid(ven)))
}
th_s <- thiessen_por_region(ev, "SKATER", area)
th_c <- thiessen_por_region(ev, "CLUSTGEO", area)
col_k <- setNames(COLORES[seq_len(max(reg))], reg)
ev$CLASE_EST <- cut(ev$ESTABILIDAD, c(-Inf, 0.5, 0.75, Inf), right = FALSE,
                    labels = c("menor que 0,5", "0,5 a 0,75", "0,75 o m\u00e1s"))
marco <- st_bbox(if (!is.null(valle)) valle else st_buffer(ev, RADIO_THIESSEN))
encuadre <- coord_sf(xlim = marco[c("xmin", "xmax")], ylim = marco[c("ymin", "ymax")], datum = NA)
asp <- as.numeric((marco["ymax"] - marco["ymin"]) / (marco["xmax"] - marco["xmin"]))
ALTO_1 <- 5.2 * asp + 1.0      # alto de un mapa de 7 pulgadas de ancho, sin espacio en blanco
ALTO_2 <- 4.4 * asp + 1.2      # alto del panel comparativo de 11 pulgadas de ancho
capa_valle <- if (!is.null(valle)) geom_sf(data = valle, fill = NA, color = "grey30", linewidth = 0.5) else NULL
capa_rio <- if (!is.null(rio)) geom_sf(data = rio, color = "#3B6FB6", linewidth = 0.6) else NULL

g1 <- ggplot() + capa_valle + capa_rio +
  geom_sf(data = ev, aes(color = factor(SKATER), shape = CLASE_EST), size = 1.7, stroke = 0.7) +
  scale_color_manual(values = col_k, name = "Regi\u00f3n") +
  scale_shape_manual(values = setNames(c(4, 1, 16), levels(ev$CLASE_EST)),
                     name = "Estabilidad\n(bootstrap)", drop = FALSE) + 
  encuadre + labs(title = sprintf("SKATER, k = %d", K),
                  subtitle = "Coronas por regi\u00f3n; forma = estabilidad en el bootstrap") + tema_mapa +
  theme(legend.position = "inside", legend.position.inside = c(0.98, 0.02),
        legend.justification = c(1, 0), legend.box = "horizontal",
        legend.background = element_rect(fill = "white", colour = NA),
        legend.margin = margin(4, 6, 4, 6)) +
  guides(color = guide_legend(override.aes = list(size = 3.5)),
         shape = guide_legend(override.aes = list(size = 3.5)))

g2 <- ggplot() + geom_sf(data = th_s, aes(fill = factor(REGION_ID)), color = "white", linewidth = 0.1) +
  capa_valle + capa_rio + geom_sf(data = ev, size = 0.3) +
  scale_fill_manual(values = col_k, name = "Regi\u00f3n") +
  encuadre + labs(title = sprintf("SKATER, k = %d", K),
                  subtitle = sprintf("Pol\u00edgonos de Thiessen recortados a %d m de los eventos", RADIO_THIESSEN)) + tema_mapa

g3 <- ggplot() + capa_valle + capa_rio +
  geom_sf(data = ev, aes(color = factor(CLUSTGEO),
                         shape = factor(FRAGMENTO, levels = 0:1, labels = c("contiguo", "fragmento"))),
          size = 1.7, stroke = 0.7) +
  scale_color_manual(values = col_k, name = "Regi\u00f3n") + scale_shape_manual(values = c(16, 4), name = NULL) +
  encuadre + labs(title = sprintf("ClustGeo, k = %d", K),
                  subtitle = "Colores alineados con SKATER; x = fuera de la parte contigua") + tema_mapa
g4 <- ggplot() + geom_sf(data = th_c, aes(fill = factor(REGION_ID)), color = "white", linewidth = 0.1) +
  capa_valle + capa_rio + geom_sf(data = ev, size = 0.3) +
  scale_fill_manual(values = col_k, name = "Regi\u00f3n") +
  encuadre + labs(title = sprintf("ClustGeo, k = %d", K),
                  subtitle = sprintf("Pol\u00edgonos de Thiessen recortados a %d m de los eventos", RADIO_THIESSEN)) + tema_mapa
mostrar_guardar(g1, "F4_D_mapa_skater_puntos.png", 7, ALTO_1)
mostrar_guardar(g2, "F4_D_mapa_skater_thiessen.png", 7, ALTO_1)
mostrar_guardar(g3, "F4_D_mapa_clustgeo_puntos.png", 7, ALTO_1)
mostrar_guardar(g4, "F4_D_mapa_clustgeo_thiessen.png", 7, ALTO_1)

comp_th <- rbind(transform(th_s, metodo = "SKATER"), transform(th_c, metodo = "ClustGeo"))
comp_th$metodo <- factor(comp_th$metodo, levels = c("SKATER", "ClustGeo"))
g5 <- ggplot() + geom_sf(data = comp_th, aes(fill = factor(REGION_ID)), color = "white", linewidth = 0.1) +
  capa_valle + capa_rio + scale_fill_manual(values = col_k, name = "Regi\u00f3n") + facet_wrap(~metodo) +
  encuadre + labs(title = sprintf("Regiones de SKATER y ClustGeo (k = %d)", K),
                  subtitle = "Colores de ClustGeo alineados con la regi\u00f3n de SKATER con la que m\u00e1s coincide") + tema_mapa
mostrar_guardar(g5, "F4_D_mapa_comparacion.png", 11, ALTO_2)

## ---------------- capas para ArcGIS ----------------
gpkg <- file.path(CARPETA, "F4_regiones.gpkg")
st_write(ev[, c("ID", "LADERA", "SKATER", "CLUSTGEO", "FRAGMENTO", "ESTABILIDAD")], gpkg,
         layer = "eventos_regiones", delete_layer = TRUE, quiet = TRUE)
st_write(th_s, gpkg, layer = "thiessen_skater", delete_layer = TRUE, quiet = TRUE)
st_write(th_c, gpkg, layer = "thiessen_clustgeo", delete_layer = TRUE, quiet = TRUE)

cat("\nGuardado en ", CARPETA, ":\n",
    "  Tablas: F4_D_caracterizacion_k", K, ".csv, F4_D_cramer_v_k", K, ".csv\n",
    "  Figuras: F4_D_perfil_regiones.png, F4_D_mapa_skater_puntos.png, F4_D_mapa_skater_thiessen.png,\n",
    "           F4_D_mapa_clustgeo_puntos.png, F4_D_mapa_clustgeo_thiessen.png, F4_D_mapa_comparacion.png\n",
    "  Capas: F4_regiones.gpkg (eventos_regiones, thiessen_skater, thiessen_clustgeo)\n", sep = "")
cat(strrep("=", 78), "\n")
