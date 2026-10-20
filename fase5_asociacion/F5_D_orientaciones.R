# =====================================================================
#' #  FASE 5 - Parte D:. Alineación de los movimientos con las fallas
#' 
#' Esta parte final tiene el objetivo de responder: 
#' ¿los eventos anómalos están más alineados con las fallas que el resto del inventario?
#' 
#' Inicialmente se midió el ángulo entre el eje mayor y la orientación local de la falla más cercana de cada evento
#' 
#' * `ANG_REL_S`  entre −90 y 90°: la diferencia con signo.
#' * `ANG_REL` entre 0 y 90°: su valor absoluto 
#' 
#' **Pruebas**
#' 
#' |Prueba|Que trata de resolver|
#'  |---|---|
#'  |$U^{2}$ de Watson de dos muestras ( Watson, 1962)| Si al distribución  del ángulo relativo de los anómalos difiere de la del resto.|
#'  |Mann-Whitney unilateral sobre ANG_REL|Los eventos anomalos tienen ángulos menores, es decir, son más paralelos a la falla  |
#'  |Kuiper de una muestra|¿El ángulo relativo de todos los eventos se aparta de orientaciones al azar?|
#'  |U^{2}$ eventos frente a tramos de falla |Existe una diferencia ente la orientación general de los eventos y la fallas|
#' 
#' ** $U^{2}$ es una prueba hecha para el círculo: no depende del punto desde el que se empiece a medir el ángulo, a diferencia de Kolmogorov-Smirnov.**
#' 
#' 
#' Se realizan tres corridas por tres contextos distintos:
#' 
#' 1. eventos con elongación de al menos 1,2, cuya orientación está definida.
#' 
#' 2. mismos eventos pero solo a menos de 1 km de la falla.
#' 
#' 3. todos los eventos, sin ningún filtro.
#' 
#' 
#' **Datos de entrada**
#' 
#' F5_eventos.csv y F5_estructuras_seg.csv de la parte 0. 
#' 
#' 
#' 

# =====================================================================

#' # Configuración 
#' 
#' 
#' Reúne las rutas, los nombres de capas y campos, y los parámetros de la fase. Es lo único que hay que editar.



#-------------------------------------------------------------------------------

#Rutas
# Carpeta del proyecto, archivos csv y carpeta de resultados

BASE      <- "D:/Users/tasal/Desktop/Tesis de grado"
EVENTOS   <- file.path(BASE, "fase5", "Fase5_csv", "F5_eventos.csv")
SEGMENTOS <- file.path(BASE, "fase5", "Fase5_csv", "F5_estructuras_seg.csv")
SALIDA    <- file.path(BASE, "fase5", "parteD")



E_MIN       <- 1.2       # Elongación mínima de la corrida principal.
DIST_CERCA  <- 1000      # Distancia de la sensibilidad.
N_PERM      <- 9999      # Permutaciones y simulaciones de cada prueba.

# Clases del histograma y pétalos de las rosetas
ANCHO_BIN   <- 15        
ANCHO_PETALO <- 10      


COL_FALLA <- "grey35"; COL_NORM <- "grey65"; COL_ANOM <- "#C0392B"
set.seed(2026)


#===============================================================================

#' # Lectura
#' 
#' ### Que hace
#' 
#' lee los eventos, toma ANOM (o la arma con ANOMALIA), 
#' comprueba las columnas y deja los eventos con ángulo relativo.
#' 
#' 
#' 

#-------------------------------------------------------------------------------
leer <- function(ruta) {
  x <- read.csv(ruta, fileEncoding = "UTF-8-BOM", check.names = FALSE, stringsAsFactors = FALSE)
  names(x) <- trimws(sub("^\ufeff", "", names(x)))
  x
}
dir.create(SALIDA, showWarnings = FALSE, recursive = TRUE)
d <- leer(EVENTOS)
if (!"ANOM" %in% names(d)) d$ANOM <- as.integer(tolower(d$ANOMALIA) != "normal")
faltan <- setdiff(c("ANOM", "ORIENT_AX", "ELONG_RECT", "ANG_REL", "ANG_REL_S", "D_FALLA"), names(d))
if (length(faltan)) stop("Faltan columnas: ", paste(faltan, collapse = ", "))
d <- d[!is.na(d$ANG_REL_S), ]
an <- d$ANOM == 1
cat("Eventos con ángulo relativo:", nrow(d), "| anomalías:", sum(an), "\n")



#================================================================================


#' # Funciones 
#' 
#' |Función|Que hace|
#'  |---|---|
#'  |a_rad2(th)| duplica un ángulo axial y lo pasa a radianes.|
#'  |eje_medio(th)|promedia los cosenos y senos de los ángulos duplicados. El ángulo de ese promedio, dividido por 2, es el eje medio |
#'  | u2_stat(a, b) |Calcula el $U^{2}$ de Watson para dos muestras de ángulos|
#'  |u2_prueba(th1, th2) |duplica los ángulos, saca $U^{2}$  y la compara con 9999 repartos al azar de los ángulos entre las dos muestras. |
#' |kuiper_V(r)| Calcula el estadistico V de Kuiper frente a la distribución uniforme en el círculo|
#' |kuiper_prueba(th) | Compara la v de los ángulos duplicados con la de 9999 muestras de ángulos uniformes del mismo tamaño. |
#' 
#' 
#' 
#' ** La $U^{2}$ de Watson**
#' 
#' La función recorre la muestra combinada ordenada y calcula en cada punto $d_k$,
#'  la diferencia entre las distribuciones acumuladas de las dos muestra:
#'  
#'  $$ U^2 = \frac{n m}{N^2} \left[ \sum_{k} d_k^2 - \frac{1}{N} \left( \sum_{k} d_k \right)^2 \right],$$
#'  
#' Con n y m datos en cada muestra y N su suma. Restar la media de las $d_k$ genera un resultado
#' independiente del punto del círculo desde el que se empiece a medir.
#' 
#' 

#------------------------------------------------------------------------------------------


a_rad2 <- function(th) ((2 * th) %% 360) * pi / 180      # axial -> ángulo duplicado

eje_medio <- function(th) {                              # eje medio (0-180) y R
  r <- a_rad2(th); C <- mean(cos(r)); S <- mean(sin(r))
  c(eje = ((atan2(S, C) * 180 / pi) / 2) %% 180, R = sqrt(C^2 + S^2))
}

u2_stat <- function(a, b) {                              # U² de Watson, dos muestras
  n <- length(a); m <- length(b); N <- n + m
  t  <- sort(c(a, b))
  dk <- findInterval(t, sort(a)) / n - findInterval(t, sort(b)) / m
  (n * m / N^2) * (sum(dk^2) - sum(dk)^2 / N)
}
u2_prueba <- function(th1, th2) {
  a <- a_rad2(th1); b <- a_rad2(th2); obs <- u2_stat(a, b)
  z <- c(a, b); n <- length(a)
  sim <- replicate(N_PERM, { i <- sample.int(length(z), n); u2_stat(z[i], z[-i]) })
  c(U2 = obs, p = (1 + sum(sim >= obs)) / (N_PERM + 1))
}

kuiper_V <- function(r) {                                # Kuiper, una muestra vs uniforme
  u <- sort((r %% (2 * pi)) / (2 * pi)); n <- length(u); i <- seq_len(n)
  max(i / n - u) + max(u - (i - 1) / n)
}
kuiper_prueba <- function(th) {
  r <- a_rad2(th); V <- kuiper_V(r)
  sim <- replicate(N_PERM, kuiper_V(runif(length(r), 0, 2 * pi)))
  c(V = V, p = (1 + sum(sim >= V)) / (N_PERM + 1))
}



#======================================================================================

#' # Alineación: anómalos frente al resto
#' 
#' 
#' ### Qué hace
#' 
#' 
#' Corre la comparación principal en las tres corridas y las guarda
#' 
#' 
#' ** Como lo hace**
#' 
#' 1. Separa los ángulos de los anómalos y del resto: con signo para las pruebas circulares
#'  y en valor absoluto para Mann-Whitney.
#'  
#' 2. Calcula $U^{2}$ entre los dos grupos, el Mann-Whitney unilateral 
#' y la prueba de Kuiper con todos los eventos de la corrida.
#' 
#' 3. Reporta la mediana del ángulo y el porcentaje de ángulos de hasta 30° 
#' en cada grupo, y su eje medio relativo.
#' 
#' 

#---------------------------------------------------------------------------------------


filtros <- list(
  "Principal (E >= 1,2)"            = d$ELONG_RECT >= E_MIN,
  "E >= 1,2 y <= 1 km de la falla"  = d$ELONG_RECT >= E_MIN & d$D_FALLA <= DIST_CERCA,
  "Sin filtro"                       = rep(TRUE, nrow(d))
)
res <- list()
for (nf in names(filtros)) {
  s  <- filtros[[nf]]
  a1 <- d$ANG_REL_S[s & an];  a0 <- d$ANG_REL_S[s & !an]
  b1 <- d$ANG_REL[s & an];    b0 <- d$ANG_REL[s & !an]
  u  <- u2_prueba(a1, a0)
  mw <- suppressWarnings(wilcox.test(b1, b0, alternative = "less"))
  k  <- kuiper_prueba(c(a1, a0))
  e1 <- eje_medio(a1); e0 <- eje_medio(a0)
  res[[nf]] <- data.frame(
    Corrida = nf, n_anom = length(a1), n_resto = length(a0),
    Mediana_ang_anom = median(b1), Mediana_ang_resto = median(b0),
    Pct30_anom = 100 * mean(b1 <= 30), Pct30_resto = 100 * mean(b0 <= 30),
    Rel_medio_anom = ((e1[["eje"]] + 90) %% 180) - 90, R_anom = e1[["R"]],
    Rel_medio_resto = ((e0[["eje"]] + 90) %% 180) - 90, R_resto = e0[["R"]],
    U2 = u[["U2"]], p_U2 = u[["p"]], p_MW_mas_paralelos = mw$p.value,
    Kuiper_todos_V = k[["V"]], Kuiper_todos_p = k[["p"]]
  )
}
res <- do.call(rbind, res); rownames(res) <- NULL
write.csv(res, file.path(SALIDA, "D_alineacion.csv"), row.names = FALSE, fileEncoding = "UTF-8")




#======================================================================================

#' # Contexto: orientación general
#' 
#' 
#' Si existe el archivo de tramos, compara la orientación de los ejes de los eventos (con elongación de al menos 1,2) 
#' con la de los tramos de falla: el eje medio y el $R$ de cada uno, el ángulo entre los dos ejes medios 
#' y la $U^{2}$  entre las dos muestras.
#' 



ctx <- NULL
fallas <- NULL
if (file.exists(SEGMENTOS)) {
  sg <- leer(SEGMENTOS)
  fallas <- sg$ORIENT_AX[!is.na(sg$ORIENT_AX)]
  ev <- d$ORIENT_AX[d$ELONG_RECT >= E_MIN]
  ef <- eje_medio(fallas); ee <- eje_medio(ev)
  u  <- u2_prueba(ev, fallas)
  dif <- abs(ee[["eje"]] - ef[["eje"]]) %% 180
  ctx <- data.frame(n_eventos = length(ev), eje_eventos = ee[["eje"]], R_eventos = ee[["R"]],
                    n_tramos = length(fallas), eje_fallas = ef[["eje"]], R_fallas = ef[["R"]],
                    dif_ejes = min(dif, 180 - dif), U2 = u[["U2"]], p_U2 = u[["p"]])
  write.csv(ctx, file.path(SALIDA, "D_contexto.csv"), row.names = FALSE, fileEncoding = "UTF-8")
}


# ================================================================================

#' # Figura del ángulo relativo
#' 
#' Prepara la proporción de eventos por clase de 15° en la corrida principal, para anómalos y resto.
#' ibuja el panel (a) con barras, y la línea punteada de lo esperado al azar (la misma proporción en cada clase), 
#' y el panel (b) con las dos distribuciones acumuladas, la diagonal del azar y el resultado de la $U^{2}$
#' 
#' 
#' 
#' 

sp  <- filtros[[1]]
r1  <- d$ANG_REL[sp & an]; r0 <- d$ANG_REL[sp & !an]
cortes <- seq(0, 90, by = ANCHO_BIN)
h1 <- hist(pmin(r1, 89.999), breaks = cortes, plot = FALSE)$counts / length(r1)
h0 <- hist(pmin(r0, 89.999), breaks = cortes, plot = FALSE)$counts / length(r0)

fig_alineacion <- function() {
  op  <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1), mgp = c(2.6, 0.7, 0), las = 1,
             cex.main = 1, cex.lab = 0.95, cex.axis = 0.85)
  dec <- options(OutDec = ",")
  on.exit({ par(op); options(dec) })
  etiquetas <- paste0(head(cortes, -1), "-", cortes[-1])
  barplot(rbind(h0, h1), beside = TRUE, names.arg = etiquetas, col = c(COL_NORM, COL_ANOM),
          border = NA, ylim = c(0, max(c(h0, h1)) * 1.3), cex.names = 0.7, xlab = "Ángulo con la falla más cercana (°)",
          ylab = "Proporción de eventos", main = "(a) Distribución del ángulo")
  abline(h = ANCHO_BIN / 90, lty = 3)
  legend("topright", bty = "n", cex = 0.75, fill = c(COL_NORM, COL_ANOM, NA), border = NA,
         lty = c(NA, NA, 3), legend = c(sprintf("Resto (n = %d)", length(r0)),
                                        sprintf("Anómalos (n = %d)", length(r1)), "Orientación al azar"))
  plot(ecdf(r0), do.points = FALSE, verticals = TRUE, col = COL_NORM, lwd = 1.5,
       xlim = c(0, 90), xlab = "Ángulo con la falla más cercana (°)", ylab = "Proporción acumulada",
       main = "(b) Distribución acumulada")
  plot(ecdf(r1), do.points = FALSE, verticals = TRUE, col = COL_ANOM, lwd = 1.5, add = TRUE)
  abline(0, 1 / 90, lty = 3)
  legend("bottomright", bty = "n", cex = 0.8,
         legend = paste0("U² = ", format(round(res$U2[1], 3), nsmall = 3),
                         "; p = ", format(round(res$p_U2[1], 2), nsmall = 2)))
}


#================================================================================

#' Figura de las rosetas
#' 
#' Dos funciones:
#' 
#' 1. `rosa(th, tit, col, rmax, cortes_p)` dibuja una roseta. Cuenta los ejes por pétalo de 10° entre 0 y 180°, 
#' dibuja cada pétalo dos veces, en su dirección y en la opuesta, porque un eje apunta a los dos lados,
#'  y traza el eje medio con línea discontinua.
#' 
#' 2. `fig_rosetas()` pone lado a lado las rosetas de los tramos de falla, del resto de eventos y de los anómalos. Las tres usan la misma escala
#' 

#-----------------------------------------------------------------------------------------

rosa <- function(th, tit, col, rmax, cortes_p) {
  h <- hist(th %% 180, breaks = cortes_p, plot = FALSE)$counts / length(th)
  plot.new(); plot.window(c(-1.2, 1.2) * rmax, c(-1.2, 1.2) * rmax, asp = 1)
  aa <- seq(0, 2 * pi, length.out = 200)
  for (k in c(0.5, 1)) lines(k * rmax * sin(aa), k * rmax * cos(aa), col = "grey85")
  r <- sqrt(h)
  for (i in seq_along(h)) for (off in c(0, 180)) {
    a <- (seq(cortes_p[i], cortes_p[i + 1], length.out = 12) + off) * pi / 180
    polygon(c(0, r[i] * sin(a)), c(0, r[i] * cos(a)), col = col, border = "white", lwd = 0.5)
  }
  em <- eje_medio(th); ae <- em[["eje"]] * pi / 180
  segments(-rmax * sin(ae), -rmax * cos(ae), rmax * sin(ae), rmax * cos(ae), lty = 2, lwd = 1.2)
  text(0, 1.13 * rmax, "N", font = 2, cex = 0.9)
  title(main = tit, line = 0.8)
  mtext(paste0("n = ", length(th), "; eje medio = ", round(em[["eje"]]), "°; R = ",
               format(round(em[["R"]], 2), nsmall = 2)), side = 1, line = 0.3, cex = 0.7)
}
fig_rosetas <- function() {
  if (is.null(fallas)) return(invisible(NULL))
  cortes_p <- seq(0, 180, by = ANCHO_PETALO)
  grupos <- list(list(fallas, "(a) Tramos de falla", COL_FALLA),
                 list(d$ORIENT_AX[sp & !an], "(b) Resto de eventos", COL_NORM),
                 list(d$ORIENT_AX[sp & an], "(c) Anómalos", COL_ANOM))
  rmax <- sqrt(max(sapply(grupos, function(g) max(hist(g[[1]] %% 180, breaks = cortes_p,
                                                       plot = FALSE)$counts / length(g[[1]])))))
  op  <- par(mfrow = c(1, 3), mar = c(1.8, 0.5, 2.2, 0.5), cex.main = 1)
  dec <- options(OutDec = ",")
  on.exit({ par(op); options(dec) })
  for (g in grupos) rosa(g[[1]], g[[2]], g[[3]], rmax, cortes_p)
}

fig_alineacion()
png(file.path(SALIDA, "D_alineacion.png"), width = 2400, height = 1100, res = 300); fig_alineacion(); invisible(dev.off())
if (!is.null(fallas)) {
  fig_rosetas()
  png(file.path(SALIDA, "D_rosetas.png"), width = 2700, height = 1050, res = 300); fig_rosetas(); invisible(dev.off())
}

# ------------------------- IMPRESIÓN ---------------------------------
redondear <- function(v) {
  num <- vapply(v, is.numeric, logical(1))
  v[num] <- lapply(v[num], function(z) signif(z, 3))
  v
}
cat("\n================ Alineación con la falla más cercana: anómalos vs resto ================\n")
print(redondear(res), row.names = FALSE)
if (!is.null(ctx)) {
  cat("\n================ Contexto: orientación global de eventos y fallas ================\n")
  print(redondear(ctx), row.names = FALSE)
}
cat("\nArchivos guardados en:", SALIDA, "\n")

