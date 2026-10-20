#' =====================================================================
#'  # FASE 2 - PARTE A: ¿DONDE SE AGRUPAN LOS EVENTOS Y A QUE ESCALA?
#'  
#' Estenescript se encarga de responde la pregunta ónde se concentran los eventos en el valle
#' usa las coronas generadas por el script previo en python.
#' 
#' El análisis se limita a las laderas,
#'  con una ventana formada por las zonas de pendiente igual o mayor a 10°.
#'  En el fondo plano del valle casi no hay movimientos en masa
#'  razon por la cual se evitan estas zonas para no generar un agrupamiento aparente spor por estar en las laderas. 
#'  El conjunto de la fase son las coronas que caen dentro de esa ventana. Una segunda ventana, de 4,7°,
#'   sirve como prueba de sensibilidad, así como tambíen dejarlo sin ventana.
#'   
#' El script hace dos cosas:
#'
#' 1. Prepara la estimación de densidad por kernel (KDE): elige el radio con los datos, 
#' lo convierte al núcleo que usa ArcGIS y calcula un ráster de corrección de borde. 
#' El mapa de densidad se hace en ArcGIS.
#' 
#' 2. Calcula la función L de Ripley en cinco corridas y compara cada una con 999 patrones de puntos al azar
#' 
# =====================================================================

#' =====================================================================
#' # BLOQUE 0: paquetes, carpeta y parametros
#' 
#' ### Qué hace
#' Carga los paquetes, fija la carpeta de trabajo y los parámetros de las simulaciones,
#'imprime las versiones y fija la semilla de los números aleatorios.
#' 
#'  
#' set.seed(2026) fija el punto de partida de los números aleatorios.
#'  Las líneas de versiones imprimen la de R y la de cada paquete instalado 
#' 
#' ====================================================================



library(sf)
library(spatstat.geom)
library(spatstat.random)
library(spatstat.explore)

carpeta <- "D:/Users/tasal/Desktop/Tesis de grado/fase2"
setwd(carpeta)

N_SIM <- 999        # simulaciones
NRANK <- 25         # banda puntual del 95 % con 999 simulaciones
NRANK_GLOBAL <- 50  # banda global del 95 % con 999 simulaciones
R_MIN <- 100        # desde aqui se interpreta (por debajo: posibles reactivaciones)
R_MAX <- 3000       # distancia maxima
set.seed(2026)



# ======================================================
#' # BLOQUE 1: lectura y revisión de datos 
#' 
#' Lee coronas.shp, encuentra las columnas de grupo 
#' y de reactivación aunque tengan otro nombre,
#'y crea dos columnas normalizadas: el grupo (grupo_n) 
#'y si el evento es una reactivación (es_reac). Imprime los conteos para revisarlos.
#' 
#' 
#' Im prime los conteos por grupo, 
#' los valores originales del campo de reactivación
#'  y cómo quedaron interpretados.
#' 
#' ====================================================

coronas <- st_zm(st_read("coronas.shp", quiet = TRUE))

# Busca cada columna por su nombre; si no existe con ese nombre exacto,
# busca una que empiece igual (por ejemplo "Reactivac" o "REACTIV").
buscar <- function(nombre, prefijo) {
  if (nombre %in% names(coronas)) return(nombre)
  cand <- grep(paste0("^", prefijo), names(coronas), ignore.case = TRUE, value = TRUE)
  if (length(cand) == 0) stop("No encuentro la columna '", nombre, "'.\n  Columnas disponibles: ",
                              paste(names(coronas), collapse = ", "))
  cat("Se usa la columna", cand[1], "como", nombre, "\n")
  cand[1]
}
if (!"ID_POL" %in% names(coronas)) stop("No encuentro la columna ID_POL.\n  Columnas disponibles: ",
                                        paste(names(coronas), collapse = ", "))
col_grupo <- buscar("grupo", "grup")
col_reac  <- buscar("reactivaci", "reac")

# Se normalizan los textos por su primera letra, asi no importan tildes,
# mayusculas ni espacios. La reactivacion se acepta como texto ("Si", "Sí",
# "SI") o como numero (1 = si, 0 = no).
primera <- function(x) tolower(substr(trimws(as.character(x)), 1, 1))
es_si   <- function(x) primera(x) %in% c("s", "y", "1", "t")
coronas$grupo_n <- factor(ifelse(primera(coronas[[col_grupo]]) == "m", "Mapeo", "Institucional"))
coronas$es_reac <- es_si(coronas[[col_reac]])

cat("Coronas:", nrow(coronas), "\n")
print(table(grupo = coronas$grupo_n))
cat("\nValores originales del campo de reactivacion:\n")
print(table(coronas[[col_reac]], useNA = "ifany"))
cat("Interpretados como reactivacion (TRUE = si):\n")
print(table(reactivacion = coronas$es_reac))




# =============================================================================
#' # BLOQUE 2: ventanas y conjunto de la fase 
#'
#' Lee las dos ventanas de laderas como máscaras de píxeles de 20 m,
#' define qué coronas caen dentro de cada una y arma los cinco patrones
#' de puntos de las corridas. Guarda el conjunto de la fase y las ventanas
#' para la parte B. 
#'
#' ** Las ventanas se crean desde ArcGIS Pro, con un filtro de pendiente en la misma capa**
#'
#' Las ventanas se leen como mascaras de pixeles (raster ASCII de ArcGIS).
#' Un poligono de laderas tiene miles de partes y huecos, y convertirlo al
#' formato de spatstat desborda la memoria de R; la mascara no tiene ese problema.
#===============================================================================

leer_mascara <- function(archivo, resolucion = 20) {
  if (!file.exists(archivo)) stop("No encuentro '", archivo, "' en ", getwd())
  ext_ <- tolower(tools::file_ext(archivo))
  if (ext_ %in% c("asc", "txt")) {                 # raster ASCII de ArcGIS
    cab <- trimws(readLines(archivo, n = 6))
    v   <- setNames(as.numeric(sub("^\\S+\\s+", "", cab)), tolower(sub("\\s.*$", "", cab)))
    nc <- v[["ncols"]]; nr <- v[["nrows"]]; cs <- v[["cellsize"]]
    x0 <- if (!is.na(v["xllcorner"])) v[["xllcorner"]] else v[["xllcenter"]] - cs / 2
    y0 <- if (!is.na(v["yllcorner"])) v[["yllcorner"]] else v[["yllcenter"]] - cs / 2
    nd <- if (!is.na(v["nodata_value"])) v[["nodata_value"]] else -9999
    m  <- matrix(scan(archivo, skip = 6, quiet = TRUE), nrow = nr, ncol = nc, byrow = TRUE)
    dentro <- !is.na(m) & m != nd
    return(owin(xrange = c(x0, x0 + nc * cs), yrange = c(y0, y0 + nr * cs),
                mask = dentro[nr:1, ]))
  }
  if (!requireNamespace("terra", quietly = TRUE)) install.packages("terra")
  if (ext_ %in% c("shp", "gpkg")) {                # ventana guardada como poligono
    pol <- terra::vect(archivo)
    r   <- terra::rasterize(pol, terra::rast(pol, resolution = resolucion), field = 1)
  } else {                                         # .tif, .img o carpeta de raster de ArcGIS
    r <- terra::rast(archivo)[[1]]
  }
  m <- terra::as.matrix(r, wide = TRUE)            # fila 1 = norte
  dentro <- !is.na(m) & m != 0                     # 0 o sin dato = fuera de la ventana
  e <- as.vector(terra::ext(r))                    # xmin, xmax, ymin, ymax
  owin(xrange = e[1:2], yrange = e[3:4], mask = dentro[nrow(m):1, ])
}
W10 <- leer_mascara("ventana_10grados.shp")
W47 <- leer_mascara("ventana_4_7grados.shp")
cat(sprintf("Area ventana 10 grados:  %.1f km2\n", area(W10) / 1e6))
cat(sprintf("Area ventana 4,7 grados: %.1f km2\n", area(W47) / 1e6))
W10 <- leer_mascara("ventana_10grados.shp")      # 10g = 10 grados
W47 <- leer_mascara("ventana_4_7grados.shp")
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

# ==============================================================================
#' # BLOQUE 3: radio y apoyos para el KDE en ArcGIS 
#' 
#' **¿Que hace?**
#' Calcula en R lo que ArcGIS no puede calcular solo para el mapa de densidad:
#' el radio del núcleo elegido por los datos, su equivalente para el núcleo de ArcGIS,
#' un ráster de corrección de borde, los cortes de clase del mapa y la capa de puntos de entrada.
#' Deja los parámetros en kde_para_arcgis.txt.
#' 
#' a) bw.ppl elige el ancho de banda $\sigma$ de un núcleo gaussiano por validación cruzada de verosimilitud. 
#' Para cada $\sigma$ candidato estima la densidad en cada evento sin usar ese evento $\hat{\lambda}_{-i}(x_i)$ y calcula 
#' $$\mathrm{CV}(\sigma) = \sum_{i} \log \hat{\lambda}_{-i}(x_i) - \int_{W} \hat{\lambda}(u) \, du,$$
#' y se queda con el $\sigma$ que da el mayor valor, es decir, 
#' el que hace más probable el patrón observado.
#' 
#' b) **Radio para ArcGIS** La herramienta Densidad kernel de ArcGIS usa un núcleo cuártico,
#'  y no uno gaussiano.
#'  Para que los dos suavicen lo mismo se igualan sus varianzas.
#'  
#' c) **Factor de borde** Cerca del borde de la ventana,
#'  que es el límite del área de estudio y el del fondo del valle,
#'  parte del núcleo de cada evento cae fuera de las laderas. 
#'  ArcGIS reparte ahí una parte de la densidad, así que la subestima en el borde.
#'  El factor de borde de cada píxel es la fracción del núcleo, centrado en ese píxel, 
#'  que cae dentro de la ventana; al dividir el mapa de ArcGIS por este factor se corrige esa pérdida.
#' 
# ==============================================================================================


CELDA_KDE <- 50                       # tamano de pixel del mapa, en metros
X_kde <- unmark(P10)

# a) Radio. Se calcula en pixeles de 100 m para aligerar la carga

sigma_ppl    <- bw.ppl(X_kde, eps = 100)
sigma_diggle <- NA                    ++


# b) ArcGIS usa un nucleo cuartico, no gaussiano. El radio que produce el
#    mismo suavizado es sigma x raiz de 8.
radio_arcgis <- sigma_ppl * sqrt(8)

# c) Factor de borde: para cada pixel, que fraccion del nucleo cae dentro de
#    la ventana. Cerca de los bordes (del AMVA y del fondo del valle) ArcGIS
#    reparte parte de la densidad fuera de las laderas y la subestima;
#    dividir el mapa de ArcGIS por este factor lo corrige.

marco <- Frame(W10)
xs <- seq(marco$xrange[1] + CELDA_KDE / 2, marco$xrange[2] + CELDA_KDE, by = CELDA_KDE)
ys <- seq(marco$yrange[1] + CELDA_KDE / 2, marco$yrange[2] + CELDA_KDE, by = CELDA_KDE)
Wm <- as.mask(W10, xy = list(x = xs, y = ys))
nk <- ceiling(radio_arcgis / CELDA_KDE); u <- seq(-nk, nk) * CELDA_KDE
d2 <- outer(u^2, u^2, "+")
nucleo <- ifelse(d2 <= radio_arcgis^2, (1 - d2 / radio_arcgis^2)^2, 0)
nucleo <- im(nucleo / (sum(nucleo) * CELDA_KDE^2), xcol = u, yrow = u)
borde <- as.im(convolve.im(as.im(Wm, value = 1, na.replace = 0), nucleo), W = Wm)
borde$v[!Wm$m] <- NA                  # solo dentro de la ventana


# d) Cortes de clase: multiplos de la densidad media en las laderas.
#    Son los mismos para cualquier mapa (integrado, mapeo, institucional).
dens_media <- npoints(X_kde) / area(W10) * 1e6        # eventos por km2 de ladera
cortes <- dens_media * c(0.5, 1, 2, 4, 8)


# e) Puntos del conjunto de la fase, para usarlos como entrada en ArcGIS
pts <- coronas[dentro10, "ID_POL"]
pts$grupo   <- as.character(coronas$grupo_n[dentro10])
pts$reactiv <- as.integer(coronas$es_reac[dentro10])
st_write(pts, "coronas_fase2_10g.shp", delete_dsn = file.exists("coronas_fase2_10g.shp"), quiet = TRUE)

# Exportar el factor de borde como raster ASCII con la proyeccion de las coronas
exportar_asc <- function(im, archivo) {
  v <- as.matrix(im)[nrow(im$v):1, ]        # primera fila = la de mas al norte
  v[is.na(v)] <- -9999
  writeLines(c(paste("ncols", ncol(v)), paste("nrows", nrow(v)),
               paste("xllcorner", im$xrange[1]), paste("yllcorner", im$yrange[1]),
               paste("cellsize", im$xstep), "NODATA_value -9999"), archivo)
  v <- matrix(formatC(v, format = "f", digits = 5), nrow = nrow(v))
  write.table(v, archivo, append = TRUE, row.names = FALSE, col.names = FALSE, quote = FALSE)
  invisible(file.copy("coronas.prj", sub("\\.asc$", ".prj", archivo), overwrite = TRUE))
}
exportar_asc(borde, "borde_10g.asc")


# Todo lo que se necesita para ArcGIS, en un archivo de texto
sink("kde_para_arcgis.txt")
cat("PARAMETROS PARA EL KDE EN ARCGIS (ventana de 10 grados)\n\n")
cat(sprintf("Eventos: %d (archivo coronas_fase2_10g.shp)\n", npoints(X_kde)))
cat(sprintf("Radio elegido por los datos (sigma gaussiano): %.0f m\n", sigma_ppl))
cat(sprintf("Search radius para ArcGIS (nucleo cuartico):   %.0f m\n", radio_arcgis))
cat(sprintf("Output cell size:                              %d m\n", CELDA_KDE))
cat(sprintf("Densidad media en las laderas: %.3f eventos/km2\n\n", dens_media))
cat("Cortes de clase (eventos/km2), multiplos de la densidad media:\n")
cat(sprintf("  0,5x = %.3f | 1x = %.3f | 2x = %.3f | 4x = %.3f | 8x = %.3f\n",
            cortes[1], cortes[2], cortes[3], cortes[4], cortes[5]))
sink()
cat(readLines("kde_para_arcgis.txt"), sep = "\n")




# ===============================================================================
#' # BLOQUE 4: funcion L en las cinco corridas (tarda varios minutos) 
#' 
#' **¿Qué hace?**
#' 
#' Calcula la función L de Ripley de cada corrida y
#'  la compara con 999 patrones de puntos repartidos al azar en la misma ventana. 
#' 
#' La función K de Ripley mide cuántos eventos hay, en promedio, 
#' a menos de una distancia r de cada evento, dividido por la densidad media.
#' Si los eventos estuvieran repartidos al azar $K(r) = \pi r^2$.
#' La función L la transforma para que el azar sea una recta, $L(r)= \sqrt{K(r)/\pi}$
#' que bajo azar es $L(r)= r$ y se grafica $L(r) - r$
#' 
#' L(r) con correccion de borde por traslacion. La covarianza de la ventana se
#' calcula una sola vez por ventana; asi cada simulacion toma fracciones de
#' segundo en vez de minutos.
#' 
# ==============================================================================

r_vec <- seq(0, R_MAX, length.out = 301)
en_int <- r_vec >= R_MIN & r_vec <= R_MAX       # tramo de interpretacion
cov_ventana <- list(W10 = setcov(as.mask(W10, eps = 50)),
                    W47 = setcov(as.mask(W47, eps = 50)))

L_trans <- function(X, gW) {
  n <- npoints(X); A <- area(Window(X))
  cp <- closepairs(X, R_MAX, twice = TRUE, what = "all")
  w  <- edge.Trans(dx = cp$dx, dy = cp$dy, W = Window(X), paired = TRUE, gW = gW)
  o  <- order(cp$d); acum <- c(0, cumsum(w[o]))
  sqrt(A / (n * (n - 1)) * acum[findInterval(r_vec, cp$d[o]) + 1] / pi)
}

envolvente_L <- function(X, gW) {
  obs  <- L_trans(X, gW)
  # una simulacion a la vez, para no acumular 999 patrones en memoria
  sims <- sapply(seq_len(N_SIM), function(i) L_trans(runifpoint(npoints(X), Window(X)), gW))
  ord  <- apply(sims, 1, sort)                              # cada columna: un r
  lo   <- ord[NRANK, ]; hi <- ord[N_SIM + 1 - NRANK, ]      # banda puntual
  # banda global: maxima desviacion respecto al azar (L = r) en el tramo
  desv  <- apply(abs(sims[en_int, ] - r_vec[en_int]), 2, max)
  u     <- sort(desv, decreasing = TRUE)[NRANK_GLOBAL]
  # prueba DCLF: suma de desvios cuadrados en el tramo
  t_obs <- sum((obs[en_int] - r_vec[en_int])^2)
  t_sim <- colSums((sims[en_int, ] - r_vec[en_int])^2)
  list(tabla = data.frame(r = r_vec, obs = obs, lo = lo, hi = hi,
                          glo = ifelse(en_int, r_vec - u, NA), ghi = ifelse(en_int, r_vec + u, NA)),
       p_DCLF = (1 + sum(t_sim >= t_obs)) / (N_SIM + 1))
}

ventana_de <- c("W10", "W10", "W47", "W10", "W10")      # ventana de cada corrida
resultados_L <- list()
for (i in seq_along(corridas)) {
  k <- names(corridas)[i]
  cat("Calculando", k, "...\n")
  resultados_L[[k]] <- envolvente_L(corridas[[k]], cov_ventana[[ventana_de[i]]])
}
saveRDS(resultados_L, "envolventes_L_parteA.rds")
cat("Listo.\n")



# ===========================================================================
#' # BLOQUE 5: resumen  
#' 
#' 
#' Resume cada corrida en resultados_parteA.txt y resumen_parteA.csv:
#'  número de eventos, densidad, valor p de la prueba DCLF,
#'   tramos de distancia donde la curva sale de cada banda y máximo de $L(r) - r$  
#' 
#' 
# ============================================================================


if (exists("resultados_L")) {                       # version con ventanas como mascara
  tablas <- lapply(resultados_L, function(x) x$tabla)
  p_dclf <- sapply(resultados_L, function(x) x$p_DCLF)
} else {                                            # version con ventanas como poligono
  tablas <- lapply(names(env_puntual), function(k) {
    d <- as.data.frame(env_puntual[[k]])
    g <- as.data.frame(envelope(env_puntual[[k]], global = TRUE, ginterval = c(R_MIN, R_MAX),
                                nrank = 50, verbose = FALSE))   # banda global del 95 %
    d$glo <- g$lo[match(d$r, g$r)]; d$ghi <- g$hi[match(d$r, g$r)]
    d[, c("r", "obs", "lo", "hi", "glo", "ghi")]
  })
  names(tablas) <- names(env_puntual)
  p_dclf <- sapply(env_puntual, function(e) dclf.test(e, rinterval = c(R_MIN, R_MAX))$p.value)
}

tramos <- function(r, cond) {
  cond[is.na(cond)] <- FALSE
  if (!any(cond)) return("ninguno")
  rl <- rle(cond); fin <- cumsum(rl$lengths); ini <- fin - rl$lengths + 1
  paste(sprintf("%.0f a %.0f m", r[ini[rl$values]], r[fin[rl$values]]), collapse = "; ")
}

resumen <- data.frame()
sink("resultados_parteA.txt")
cat("FASE 2 - PARTE A\n")
cat(sprintf("KDE: radio por validacion cruzada = %.0f m (ArcGIS: %.0f m)\n", sigma_ppl, radio_arcgis))
cat(sprintf("Funcion L: %d simulaciones | bandas del 95 %% | interpretacion de %d a %d m\n\n",
            N_SIM, R_MIN, R_MAX))
for (k in names(tablas)) {
  di <- tablas[[k]]; di <- di[di$r >= R_MIN & di$r <= R_MAX, ]
  dif <- di$obs - di$r; i <- which.max(dif)
  cat("==", k, "==\n")
  cat(sprintf("  Eventos: %d | intensidad: %.3f eventos/km2\n",
              npoints(corridas[[k]]), intensity(corridas[[k]]) * 1e6))
  cat(sprintf("  Prueba global DCLF (%d-%d m): p = %.3f\n", R_MIN, R_MAX, p_dclf[[k]]))
  cat("  Fuera de la banda GLOBAL por encima:", tramos(di$r, di$obs > di$ghi), "\n")
  cat("  Fuera de la banda GLOBAL por debajo:", tramos(di$r, di$obs < di$glo), "\n")
  cat("  Fuera de la banda PUNTUAL por encima:", tramos(di$r, di$obs > di$hi), "\n")
  cat("  Fuera de la banda PUNTUAL por debajo:", tramos(di$r, di$obs < di$lo), "\n")
  cat(sprintf("  Maximo de L(r) - r: %.0f m, a %.0f m\n\n", dif[i], di$r[i]))
  resumen <- rbind(resumen, data.frame(corrida = k, n = npoints(corridas[[k]]), p_DCLF = p_dclf[[k]],
                                       max_L_menos_r = dif[i], r_del_maximo = di$r[i]))
}
sink()
cat(readLines("resultados_parteA.txt"), sep = "\n")
write.csv(resumen, "resumen_parteA.csv", row.names = FALSE)



# ==================================================================================

#' # BLOQUE 6: figuras 
#' 
#' ## Qué hace
#' 
#' Dibuja tres figuras de $L(r) - r$ , las muestra en el panel Plots de RStudio, 
#' las guarda en PNG y escribe todas las curvas en curvas_L_parteA.csv.
#' 
# =================================================================================

if (exists("resultados_L")) {
  tablas <- lapply(resultados_L, function(x) x$tabla)
} else {
  tablas <- lapply(names(env_puntual), function(k) {
    d <- as.data.frame(env_puntual[[k]])
    # banda global del 95 % (rango 50 de 999), con las mismas simulaciones guardadas
    g <- as.data.frame(envelope(env_puntual[[k]], global = TRUE, ginterval = c(R_MIN, R_MAX),
                                nrank = 50, verbose = FALSE))
    d$glo <- g$lo[match(d$r, g$r)]; d$ghi <- g$hi[match(d$r, g$r)]
    d[, c("r", "obs", "lo", "hi", "glo", "ghi")]
  })
  names(tablas) <- names(env_puntual)
}
k_ppal <- names(tablas)[1]; k_sens <- names(tablas)[1:3]; k_sesgo <- names(tablas)[c(1, 4, 5)]

rango_y <- function(claves) range(unlist(lapply(tablas[claves], function(d)
  c(d$obs - d$r, d$lo - d$r, d$hi - d$r, d$glo - d$r, d$ghi - d$r))), na.rm = TRUE)

# Un panel: curva observada con su banda puntual (gris) y su banda global (punteada)
panel_L <- function(k, ylim, titulo) {
  d <- tablas[[k]]
  plot(NA, xlim = c(0, R_MAX), ylim = ylim, xlab = "Distancia r (m)", ylab = "L(r) - r (m)",
       main = titulo, cex.main = 1)
  polygon(c(d$r, rev(d$r)), c(d$lo - d$r, rev(d$hi - d$r)), col = "grey85", border = NA)
  lines(d$r, d$ghi - d$r, lty = 2, col = "grey30"); lines(d$r, d$glo - d$r, lty = 2, col = "grey30")
  abline(h = 0, lty = 3, col = "red"); abline(v = R_MIN, lty = 3, col = "grey50")
  lines(d$r, d$obs - d$r, lwd = 2)
}
leyenda_bandas <- function() legend("right", inset= c(0.08,0), bty = "o", bg = "white", box.col = "white", cex = 0.65,
                                    legend = c("Patrón observado", "Banda puntual 95 %", "Banda global 95 %", "100 m"),
                                    lwd = c(2, 8, 1, 1), lty = c(1, 1, 2, 3), col = c("black", "grey85", "grey30", "grey50"))

# Un panel con varias curvas observadas superpuestas
panel_superpuestas <- function(claves, colores, titulo, ylim, etiquetas= claves) {
  plot(NA, xlim = c(0, R_MAX), ylim = ylim, xlab = "Distancia r (m)", ylab = "L(r) - r (m)",
       main = titulo, cex.main = 1)
  abline(h = 0, lty = 3, col = "red"); abline(v = R_MIN, lty = 3, col = "grey50")
  for (i in seq_along(claves)) lines(tablas[[claves[i]]]$r,
                                     tablas[[claves[i]]]$obs - tablas[[claves[i]]]$r,
                                     col = colores[i], lwd = 2)
  legend("bottomright",inset= c(0.02,0.20), legend = etiquetas, col = colores, lwd = 2, bty = "n", cex = 0.75)
}

# 6a. Resultado principal
fig_principal <- function() {
  par(mfrow = c(1, 1), mar = c(4.5, 4.5, 2, 1))
  panel_L(k_ppal, rango_y(k_ppal), "Ventana de 10°")
  leyenda_bandas()
}

# 6b. Sensibilidad: cada corrida con sus propias bandas, en la misma escala
fig_sensibilidad <- function() {
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2.5, 1))
  yl <- rango_y(k_sens)
  yl_curvas <- range(unlist(lapply(tablas[k_sens], function(d) d$obs - d$r)), na.rm = TRUE)
  panel_L(k_sens[1], yl, "(a) Ventana de 10°"); leyenda_bandas()
  panel_L(k_sens[2], yl, "(b) Ventana de 10°, sin reactivaciones")
  panel_L(k_sens[3], yl, "(c) Ventana de 4,7°")
  panel_superpuestas(k_sens, c("black", "#E69F00", "#009E73"), "(d) Comparación", yl_curvas,
                     etiquetas = c("Principal (10°)", "Sin reactivaciones (10°)", "Ventana de 4,7°"))
  par(mfrow = c(1, 1))
}

# 6c. Sesgo de registro: cada fuente con sus bandas, y las tres juntas
fig_sesgo <- function() {
  par(mfrow = c(1, 3), mar = c(4.5, 4.5, 2.5, 1))
  yl <- rango_y(k_sesgo)
  yl_curvas <- range(unlist(lapply(tablas[k_sesgo], function(d) d$obs - d$r)), na.rm = TRUE)
  panel_L(k_sesgo[2], yl, "(a) Mapeo sobre imágenes"); leyenda_bandas()
  panel_L(k_sesgo[3], yl, "(b) Registros institucionales")
  panel_superpuestas(k_sesgo, c("black", "#1f77b4", "#d62728"), "(c) Comparación", yl_curvas,
                     etiquetas = c("Inventario integrado", "Mapeo sobre imágenes",
                                   "Registros institucionales"))
  par(mfrow = c(1, 1))
}


# En pantalla
fig_principal(); fig_sensibilidad(); fig_sesgo()
  
# En archivo
png("L_principal.png", width = 1600, height = 1200, res = 220); fig_principal(); invisible(dev.off())
png("L_sensibilidad.png", width = 2600, height = 2000, res = 220); fig_sensibilidad(); invisible(dev.off())
png("L_sesgo.png", width = 3300, height = 1100, res = 220); fig_sesgo(); invisible(dev.off())

# Tabla con todas las curvas y bandas
write.csv(do.call(rbind, lapply(names(tablas), function(k) cbind(corrida = k, tablas[[k]]))),
          "curvas_L_parteA.csv", row.names = FALSE)
cat("Archivos guardados en:", getwd(), "\n")