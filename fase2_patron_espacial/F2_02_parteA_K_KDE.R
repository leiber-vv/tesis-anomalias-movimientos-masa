# =====================================================================
#  FASE 2 - PARTE A: ¿DONDE SE AGRUPAN LOS EVENTOS Y A QUE ESCALA?
#
#  1. Lee las coronas y las dos ventanas de laderas (10 y 4,7 grados).
#  2. Define el conjunto de la fase: las coronas dentro de la ventana de 10.
#  3. KDE del inventario integrado con el radio elegido por los datos
#     (validacion cruzada por verosimilitud) y correccion de borde.
#     Se exporta como raster (.asc) para hacer el mapa en ArcGIS.
#  4. Funcion L de Ripley en cinco corridas:
#       principal       integrado, ventana 10, con reactivaciones
#       sensibilidad 1  integrado, ventana 10, sin reactivaciones
#       sensibilidad 2  integrado, ventana 4,7, con reactivaciones
#       sesgo           mapeo y registros institucionales, ventana 10
#     Banda puntual (para ver distancias) y banda global (para decidir),
#     hasta 3.000 m, interpretando desde 100 m.
#
#  Como correrlo: un bloque a la vez (cursor dentro + Ctrl+Alt+T).
# =====================================================================


# ---- BLOQUE 0: paquetes, carpeta y parametros ----------------------------

library(sf)
library(spatstat.geom)
library(spatstat.random)
library(spatstat.explore)

carpeta <- "D:/Users/tasal/Desktop/Tesis de grado/fase2"
setwd(carpeta)

N_SIM <- 999        # simulaciones
NRANK <- 25         # banda puntual del 95 % con 999 simulaciones
R_MIN <- 100        # desde aqui se interpreta (por debajo: reactivaciones)
R_MAX <- 3000       # distancia maxima
CELDA_KDE <- 50     # tamano de pixel del mapa de densidad, en metros
set.seed(2026)


# ---- BLOQUE 1: leer y revisar los datos -----------------------------------

coronas <- st_zm(st_read("coronas_fase2.shp", quiet = TRUE))
v10     <- st_read("ventana_10grados.shp", quiet = TRUE)     # 10g = 10 grados
v47     <- st_read("ventana_4_7grados.shp", quiet = TRUE)

# Si tus columnas se llaman distinto, renombralas aqui. Ejemplo:
# names(coronas)[names(coronas) == "Reactivac"] <- "reactiv"
faltan <- setdiff(c("ID_POL", "grupo", "Reactivaci"), names(coronas))
if (length(faltan) > 0) stop("No encuentro estas columnas: ", paste(faltan, collapse = ", "),
                             "\n  Columnas disponibles: ", paste(names(coronas), collapse = ", "))

# Se normalizan los textos por su primera letra, asi no importan tildes,
# mayusculas ni espacios. La reactivacion se acepta como texto ("Si", "Sí",
# "SI") o como numero (1 = si, 0 = no).
primera <- function(x) tolower(substr(trimws(as.character(x)), 1, 1))
es_si   <- function(x) primera(x) %in% c("s", "y", "1", "t")
coronas$grupo_n <- factor(ifelse(primera(coronas$grupo) == "m", "Mapeo", "Institucional"))
coronas$es_reac <- es_si(coronas$Reactivaci)

cat("Coronas:", nrow(coronas), "\n")
print(table(grupo = coronas$grupo_n))
cat("\nValores originales del campo de reactivacion:\n")
print(table(coronas$reactiv, useNA = "ifany"))
cat("Interpretados como reactivacion (TRUE = si):\n")
print(table(reactivacion = coronas$es_reac))


# ---- BLOQUE 2: ventanas y conjunto de la fase -----------------------------

a_ventana <- function(v) {
  g <- st_union(st_make_valid(st_zm(st_geometry(st_transform(v, st_crs(coronas))))))
  g <- st_make_valid(st_simplify(g, preserveTopology = TRUE, dTolerance = 5))
  as.owin(g)
}
W10 <- a_ventana(v10)
W47 <- a_ventana(v47)
cat(sprintf("Area ventana 10 grados:  %.1f km2\n", area(W10) / 1e6))
cat(sprintf("Area ventana 4,7 grados: %.1f km2\n", area(W47) / 1e6))

xy <- st_coordinates(coronas)
dentro10 <- inside.owin(xy[, 1], xy[, 2], W10)
dentro47 <- inside.owin(xy[, 1], xy[, 2], W47)
cat("\nCoronas fuera de la ventana de 10 grados: ", sum(!dentro10), "\n  ID:",
    paste(coronas$ID_POL[!dentro10], collapse = ", "), "\n")
cat("Coronas fuera de la ventana de 4,7 grados:", sum(!dentro47), "\n  ID:",
    paste(coronas$ID_POL[!dentro47], collapse = ", "), "\n")

patron <- function(sel, W) {
  X <- ppp(xy[sel, 1], xy[sel, 2], window = W,
           marks = data.frame(ID_POL = coronas$ID_POL[sel],
                              grupo  = coronas$grupo_n[sel],
                              reac   = coronas$es_reac[sel]))
  unitname(X) <- c("metro", "metros")
  X
}
P10 <- patron(dentro10, W10)      # conjunto de la fase
P47 <- patron(dentro47, W47)

corridas <- list(
  "Principal (10, con reactivaciones)"   = unmark(P10),
  "Sin reactivaciones (10)"              = unmark(subset(P10, !reac)),
  "Ventana de 4,7 grados"                = unmark(P47),
  "Mapeo (10)"                           = unmark(subset(P10, grupo == "Mapeo")),
  "Registros institucionales (10)"       = unmark(subset(P10, grupo == "Institucional"))
)
cat("\n")
for (k in names(corridas)) cat(sprintf("  %-36s %4d eventos\n", k, npoints(corridas[[k]])))

# Registro del conjunto de la fase (lo usa tambien la parte B)
write.csv(data.frame(ID_POL = coronas$ID_POL, dentro_10g = dentro10, dentro_4_7g = dentro47),
          "conjunto_fase2.csv", row.names = FALSE)
saveRDS(list(W10 = W10, W47 = W47), "ventanas_laderas.rds")


# ---- BLOQUE 3: KDE del inventario integrado -------------------------------
# El radio (sigma, desviacion del nucleo gaussiano) se elige por validacion
# cruzada: se prueba cada radio prediciendo la posicion de cada evento con
# los demas, y se queda el que mejor predice.

X_kde <- unmark(P10)
sigma_ppl    <- bw.ppl(X_kde)
sigma_diggle <- bw.diggle(X_kde)
cat(sprintf("Radio elegido (validacion cruzada por verosimilitud): %.0f m\n", sigma_ppl))
cat(sprintf("Referencia (criterio de Diggle):                      %.0f m\n", sigma_diggle))
cat(sprintf("Radio equivalente si lo corres en ArcGIS (nucleo cuartico): %.0f m\n",
            sigma_ppl * sqrt(8)))

marco <- Frame(W10)
xs <- seq(marco$xrange[1] + CELDA_KDE / 2, marco$xrange[2] + CELDA_KDE, by = CELDA_KDE)
ys <- seq(marco$yrange[1] + CELDA_KDE / 2, marco$yrange[2] + CELDA_KDE, by = CELDA_KDE)
kde <- density(X_kde, sigma = sigma_ppl, edge = TRUE, xy = list(x = xs, y = ys))
kde <- eval.im(pmax(kde, 0) * 1e6)          # eventos por km2

fig_kde <- function() {
  par(mar = c(1, 1, 3, 4))
  plot(kde, main = sprintf("Densidad de eventos (radio %.0f m)", sigma_ppl),
       ribargs = list(las = 1))
  plot(W10, add = TRUE, border = "grey40", lwd = 0.5)
  points(X_kde, pch = 16, cex = 0.25)
}
fig_kde()

# Exportar como raster ASCII para ArcGIS, con la misma proyeccion de las coronas
exportar_asc <- function(im, archivo) {
  v <- as.matrix(im)[nrow(im$v):1, ]        # primera fila = la de mas al norte
  v[is.na(v)] <- -9999
  cabecera <- c(paste("ncols", ncol(v)), paste("nrows", nrow(v)),
                paste("xllcorner", im$xrange[1]), paste("yllcorner", im$yrange[1]),
                paste("cellsize", im$xstep), "NODATA_value -9999")
  writeLines(cabecera, archivo)
  v <- matrix(formatC(v, format = "f", digits = 5), nrow = nrow(v))
  write.table(v, archivo, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
  invisible(file.copy("coronas.prj", sub("\\.asc$", ".prj", archivo), overwrite = TRUE))
}
exportar_asc(kde, "kde_integrado_10g.asc")
png("kde_integrado_10g.png", width = 1600, height = 1800, res = 220); fig_kde(); invisible(dev.off())


# ---- BLOQUE 4: funcion L en las cinco corridas (tarda varios minutos) -----

r_vec <- seq(0, R_MAX, length.out = 301)        # cada 10 m
env_puntual <- list(); env_global <- list()
for (k in names(corridas)) {
  cat("Calculando", k, "...\n")
  env_puntual[[k]] <- envelope(corridas[[k]], Lest, correction = "translate", r = r_vec,
                               nsim = N_SIM, nrank = NRANK, savefuns = TRUE, verbose = FALSE)
  # La banda global reutiliza las mismas 999 simulaciones
  env_global[[k]]  <- envelope(env_puntual[[k]], global = TRUE,
                               ginterval = c(R_MIN, R_MAX), verbose = FALSE)
}
saveRDS(list(puntual = env_puntual, global = env_global), "envolventes_L_parteA.rds")
cat("Listo.\n")


# ---- BLOQUE 5: resumen --------------------------------------------------------

tramos <- function(r, cond) {
  cond[is.na(cond)] <- FALSE
  if (!any(cond)) return("ninguno")
  rl <- rle(cond); fin <- cumsum(rl$lengths); ini <- fin - rl$lengths + 1
  paste(sprintf("%.0f a %.0f m", r[ini[rl$values]], r[fin[rl$values]]), collapse = "; ")
}

resumen <- data.frame()
sink("resultados_parteA.txt")
cat("FASE 2 - PARTE A\n")
cat(sprintf("KDE: radio por validacion cruzada = %.0f m (Diggle = %.0f m)\n", sigma_ppl, sigma_diggle))
cat(sprintf("Funcion L: %d simulaciones | interpretacion de %d a %d m\n\n", N_SIM, R_MIN, R_MAX))
for (k in names(env_puntual)) {
  d <- as.data.frame(env_puntual[[k]]); g <- as.data.frame(env_global[[k]])
  d <- d[d$r >= R_MIN, ]
  p <- dclf.test(env_puntual[[k]], rinterval = c(R_MIN, R_MAX))$p.value
  dif <- d$obs - d$r; i <- which.max(dif)
  cat("==", k, "==\n")
  cat(sprintf("  Eventos: %d | intensidad: %.3f eventos/km2\n",
              npoints(corridas[[k]]), intensity(corridas[[k]]) * 1e6))
  cat(sprintf("  Prueba global DCLF (%d-%d m): p = %.3f\n", R_MIN, R_MAX, p))
  cat("  Fuera de la banda GLOBAL por encima:", tramos(g$r, g$obs > g$hi), "\n")
  cat("  Fuera de la banda GLOBAL por debajo:", tramos(g$r, g$obs < g$lo), "\n")
  cat("  Fuera de la banda PUNTUAL por encima:", tramos(d$r, d$obs > d$hi), "\n")
  cat("  Fuera de la banda PUNTUAL por debajo:", tramos(d$r, d$obs < d$lo), "\n")
  cat(sprintf("  Maximo de L(r) - r: %.0f m, a %.0f m\n\n", dif[i], d$r[i]))
  resumen <- rbind(resumen, data.frame(corrida = k, n = npoints(corridas[[k]]), p_DCLF = p,
                                       max_L_menos_r = dif[i], r_del_maximo = d$r[i]))
}
sink()
cat(readLines("resultados_parteA.txt"), sep = "\n")
write.csv(resumen, "resumen_parteA.csv", row.names = FALSE)


# ---- BLOQUE 6: figuras (en pantalla y en PNG) ------------------------------

# 6a. Resultado principal, con las dos bandas
fig_principal <- function() {
  k <- names(env_puntual)[1]
  d <- as.data.frame(env_puntual[[k]]); g <- as.data.frame(env_global[[k]])
  ylim_ok <- range(c(d$obs - d$r, d$lo - d$r, d$hi - d$r, g$lo - g$r, g$hi - g$r), na.rm = TRUE)
  par(mar = c(4.5, 4.5, 2, 1))
  plot(NA, xlim = c(0, R_MAX), ylim = ylim_ok,
       xlab = "Distancia r (m)", ylab = "L(r) - r (m)", main = "Inventario integrado")
  polygon(c(d$r, rev(d$r)), c(d$lo - d$r, rev(d$hi - d$r)), col = "grey85", border = NA)
  lines(g$r, g$hi - g$r, lty = 2, col = "grey30"); lines(g$r, g$lo - g$r, lty = 2, col = "grey30")
  abline(h = 0, lty = 3, col = "red"); abline(v = R_MIN, lty = 3, col = "grey50")
  lines(d$r, d$obs - d$r, lwd = 2)
  legend("topright", bty = "o", bg = "white", box.col = "white", cex = 0.8,
         legend = c("Observada", "Banda puntual 95 %", "Banda global 95 %", "Inicio de interpretacion"),
         lwd = c(2, 8, 1, 1), lty = c(1, 1, 2, 3), col = c("black", "grey85", "grey30", "grey50"))
}

# 6b y 6c. Curvas observadas superpuestas
fig_superpuestas <- function(claves, titulo, colores) {
  cur <- lapply(env_puntual[claves], as.data.frame)
  yl  <- range(sapply(cur, function(d) range(d$obs - d$r, na.rm = TRUE)))
  par(mar = c(4.5, 4.5, 2, 1))
  plot(NA, xlim = c(0, R_MAX), ylim = yl, xlab = "Distancia r (m)", ylab = "L(r) - r (m)", main = titulo)
  abline(h = 0, lty = 3, col = "red"); abline(v = R_MIN, lty = 3, col = "grey50")
  for (i in seq_along(cur)) lines(cur[[i]]$r, cur[[i]]$obs - cur[[i]]$r, col = colores[i], lwd = 2)
  legend("topleft", legend = claves, col = colores, lwd = 2, bty = "n", cex = 0.8)
}
fig_sensibilidad <- function() fig_superpuestas(names(env_puntual)[1:3], "Sensibilidad",
                                                 c("black", "#E69F00", "#009E73"))
fig_sesgo <- function() fig_superpuestas(names(env_puntual)[c(1, 4, 5)], "Sesgo de registro",
                                          c("black", "#1f77b4", "#d62728"))

# En pantalla
fig_principal(); fig_sensibilidad(); fig_sesgo()

# En archivo
png("L_principal.png", width = 1600, height = 1200, res = 220); fig_principal(); invisible(dev.off())
png("L_sensibilidad.png", width = 1600, height = 1200, res = 220); fig_sensibilidad(); invisible(dev.off())
png("L_sesgo.png", width = 1600, height = 1200, res = 220); fig_sesgo(); invisible(dev.off())

cat("Archivos guardados en:", getwd(), "\n")
