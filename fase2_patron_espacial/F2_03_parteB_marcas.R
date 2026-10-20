# =====================================================================
#' # FASE 2 - ETAPA 3: ¿LOS EVENTOS CERCANOS SE PARECEN EN COMPORTAMIENTO?
#'
#' **Que hace**
#' 
#' Con este script se busca responder a la segunda incognita de la FASE 2. Estudiar 
#' si los eventos cercanos entre sí se parecen en área, movilidad y forma,
#'  y si los de un mismo tipo de movimiento tienden a estar juntos. 
#'  Usa el mismo conjunto de eventos y la misma ventana de laderas de 10° de la parte A.
#' 
#' Hace tres análisis: 
#' 
#' 1. Marcas numéricas ($log_{10}$ del área, $H/L$ y elongación): 
#' función de correlación de marcas y variograma de marcas, 
#' para el inventario integrado, que es el resultado principal, 
#' y para el mapeo sobre imágenes, como contraste frente al sesgo de registro.
#' 
#' 2. Marca categórica (tipo de movimiento: deslizamiento, flujo y complejo): 
#' función de conexión de marcas normalizada de cada tipo consigo mismo.
#' 
#' Variograma experimental de H/L, en todas las direcciones y en cuatro direcciones,
#'  con el paquete `gstat`.
#'  
#' En cada uno se plantea el etiquetado aleatorio: en donde 
#' los puntos se quedan donde están y se barajan los valores entre ellos. 
#' Así se prueba si los valores dependen de la posición, 
#' sin mezclarlo con el agrupamiento de los puntos
#' 
#' 
#' **Salidas**
#' 
#' Este script entrega en la misma carpeta de la fase 2. 
#' Las, curvas y bandas de las marcas numéricas, 
#' las curvas y bandas de los tipos de movimiento, y el valor p conjunto.
#' Los variogramas de H/L con sus bandas.
#' El Resumen de resultados
#' y las figuras en pantalla y en png.
#' 
#' **Como reproducirlo**
#' 
#' Primero se debe correr la parte A, que crea conjunto_fase2.csv y ventanas_laderas.rds en la misma carpeta.
#' Se debe verificar la ruta de la carpeta en la variable carpeta del bloque 0
#' Ejecutar los bloques en orden 
#' 
#' `set.seed(2026)` fija los números aleatorios; para obtener exactamente los mismos resultados

# =====================================================================

#' # BLOQUE 0: paquetes, carpeta y parametros 
#' 
#' Carga los paquetes, fija la carpeta de trabajo y los parámetros, 
#' imprime las versiones y fija la semilla de los números aleatorios.
#' 
#' sf lee el shapefile; spatstat.geom, spatstat.random y spatstat.explore
#'  manejan los patrones de puntos y calculan las funciones de marcas,
#'   y gstat calcula el variograma.
#' 
# ===========================================================================


library(sf)
library(spatstat.geom)
library(spatstat.random)
library(spatstat.explore)
library(gstat)

carpeta <- "D:/Users/tasal/Desktop/Tesis de grado/fase2" 
setwd(carpeta)

N_SIM        <- 999          # permutaciones
NRANK        <- 25           # banda puntual del 95 %
NRANK_GLOBAL <- 50           # banda global del 95 %
R_MIN        <- 100          # desde aqui se interpreta
R_MAX        <- 3000         # distancia maxima
CORRECCION   <- "none"       # correccion de borde
set.seed(2026)

r_vec  <- seq(0, R_MAX, length.out = 301)     # cada 10 m
en_int <- r_vec >= R_MIN & r_vec <= R_MAX     # tramo de interpretacion



# ===========================================================================
#' # BLOQUE 1: Lectura de los datos y búsqueda de columnas
#' 
#' **Qué hace**
#' 
#' Lee las coronas y el conjunto de la fase,
#'  encuentra las columnas que necesita aunque tengan otro nombre,
#'   arma una tabla con los valores de cada evento y revisa que no falten datos.
#' 
#' Imprime cuántos eventos hay de cada tipo y avisa si alguna marca tiene valores vacíos o menores o iguales a cero.
#'  La función de correlación de marcas multiplica los valores de los pares, 
#'  y con valores negativos no se puede interpretar.
#' 
# ===========================================================================

coronas <- st_zm(st_read("coronas.shp", quiet = TRUE))
conjunto <- read.csv("conjunto_fase2.csv")


buscar <- function(prefijo, nombre) {
  cand <- grep(paste0("^", prefijo), names(coronas), ignore.case = TRUE, value = TRUE)
  if (length(cand) == 0) return(NA)
  cat(sprintf("  %-22s -> columna %s\n", nombre, cand[1]))
  cand[1]
}
cat("Columnas usadas:\n")
col_grupo <- buscar("grup", "grupo")
col_hl    <- buscar("h_?l", "H/L")
col_elong <- buscar("elong", "elongacion")
col_tipo  <- buscar("tipo", "tipo de movimiento")
col_log   <- buscar("log", "log del area")
col_area  <- if (is.na(col_log)) buscar("area", "area (se le saca log10)") else NA
faltan <- c(grupo = col_grupo, HL = col_hl, elongacion = col_elong, tipo = col_tipo)
if (anyNA(faltan) || (is.na(col_log) && is.na(col_area)))
  stop("Falta alguna columna. Columnas disponibles: ", paste(names(coronas), collapse = ", "))

primera <- function(x) tolower(substr(trimws(as.character(x)), 1, 1))
tipo_normal <- function(x) {
  t <- tolower(trimws(as.character(x)))
  ifelse(startsWith(t, "flu"), "Flujo",
  ifelse(startsWith(t, "des"), "Deslizamiento",
  ifelse(startsWith(t, "com"), "Complejo",
  ifelse(startsWith(t, "ca"),  "Caida",
  ifelse(startsWith(t, "rep"), "Reptacion", NA)))))
}

datos <- data.frame(
  ID_POL = coronas$ID_POL,
  grupo  = ifelse(primera(coronas[[col_grupo]]) == "m", "Mapeo", "Institucional"),
  logA   = if (!is.na(col_log)) as.numeric(coronas[[col_log]]) else log10(as.numeric(coronas[[col_area]])),
  HL     = as.numeric(coronas[[col_hl]]),
  elong  = as.numeric(coronas[[col_elong]]),
  tipo   = tipo_normal(coronas[[col_tipo]])
)
xy <- st_coordinates(coronas)

# Mismo conjunto que la parte A
dentro <- datos$ID_POL %in% conjunto$ID_POL[conjunto$dentro_10g]
cat(sprintf("\nEventos en el conjunto de la fase: %d de %d\n", sum(dentro), nrow(datos)))
cat("Tipo de movimiento (conjunto de la fase):\n")
print(table(datos$tipo[dentro], useNA = "ifany"))
for (m in c("logA", "HL", "elong")) {
  v <- datos[[m]][dentro]
  if (anyNA(v)) cat("ATENCION:", sum(is.na(v)), "eventos sin valor de", m, "\n")
  if (any(v <= 0, na.rm = TRUE)) cat("ATENCION: hay valores <= 0 en", m,
                                     "(la correlacion de marcas necesita valores positivos)\n")
}

# ===========================================================================
#' # Bloque 2. Patrones de puntos
#' 
#' **¿Que hace?**
#' 
#' Arma los patrones de puntos de spatstat con las marcas de cada evento:
#' el integrado, el del mapeo sobre imágenes y el de los tres tipos de movimiento principales.
#' Imprime cuántos eventos tiene cada patrón
#' 
#'  Las 6 caídas y los eventos sin tipo quedan fuera, porque son muy pocos.
#' 
# =========================================================================== 

W10 <- readRDS("ventanas_laderas.rds")$W10          # la ventana de la parte A
if (W10$type == "mask" && CORRECCION != "none") {
  CORRECCION <- "none"
  cat("La ventana es de pixeles: se usa correccion 'none' (con 'translate' tardaria horas).\n")
}

X_todos <- ppp(xy[dentro, 1], xy[dentro, 2], window = W10,
               marks = datos[dentro, c("ID_POL", "grupo", "logA", "HL", "elong", "tipo")])
unitname(X_todos) <- c("metro", "metros")
conjuntos <- list("Integrado" = X_todos,
                  "Mapeo"     = subset(X_todos, grupo == "Mapeo"))

TIPOS  <- c("Deslizamiento", "Flujo", "Complejo")
X_tipo <- subset(X_todos, tipo %in% TIPOS)
marks(X_tipo) <- factor(marks(X_tipo)$tipo, levels = TIPOS)

MARCAS <- c(logA = "Log del área", HL = "H/L", elong = "Elongación")
for (cj in names(conjuntos)) cat(sprintf("  %-10s %4d eventos\n", cj, npoints(conjuntos[[cj]])))
cat(sprintf("  %-10s %4d eventos (%s)\n", "Tipos", npoints(X_tipo), paste(TIPOS, collapse = ", ")))

# =========================================================================
#' # BLOQUE 3: marcas numericas: correlación y variograma de marcas 
#' 
#' **¿Qué hace?**
#' 
#' Para cada conjunto y cada marca , calcula la función de correlación de marcas 
#' y el variograma de marcas del patrón observado y de 999 copias con los valores barajados.
#' Con eso arma las bandas y la prueba DCLF.
#' 
#' 
#' 
#' **Nota sobre la correccion de borde**: en el etiquetado aleatorio los puntos no
#' se mueven, asi que la curva observada y las 999 permutadas se calculan con
#' exactamente los mismos pares. La correccion de borde pondera esos pares igual
#' en todas las curvas: no cambia la comparacion, solo el estimador. Ademas, estas
#' funciones son cocientes y el efecto de borde se cancela en buena parte. En
#' pruebas, la curva observada con y sin correccion difirio menos de 2 % a 100 m
#' y menos de 0,1 % desde 500 m. Por eso se usa "none", que es mucho mas rapido.
#'
# =========================================================================


valor <- function(f) as.data.frame(f)[[fvnames(f, ".y")]]

bandas <- function(obs, sims) {
  ref  <- rowMeans(sims, na.rm = TRUE)
  kesimo <- function(x, k, dec) { x <- sort(x, decreasing = dec); if (length(x) < k) NA else x[k] }
  lo   <- apply(sims, 1, kesimo, k = NRANK, dec = FALSE)
  hi   <- apply(sims, 1, kesimo, k = NRANK, dec = TRUE)
  desv <- apply(abs(sims[en_int, ] - ref[en_int]), 2, max, na.rm = TRUE)
  u    <- sort(desv, decreasing = TRUE)[NRANK_GLOBAL]
  t_obs <- sum((obs[en_int] - ref[en_int])^2, na.rm = TRUE)
  t_sim <- colSums((sims[en_int, ] - ref[en_int])^2, na.rm = TRUE)
  list(tabla = data.frame(r = r_vec, obs = obs, ref = ref, lo = lo, hi = hi,
                          glo = ifelse(en_int, ref - u, NA), ghi = ifelse(en_int, ref + u, NA)),
       p_DCLF = (1 + sum(t_sim >= t_obs)) / (N_SIM + 1))
}

res_marcas <- list()
for (cj in names(conjuntos)) for (m in names(MARCAS)) {
  X <- unmark(conjuntos[[cj]]); v <- marks(conjuntos[[cj]])[[m]]
  cat(sprintf("%s - %s ", cj, MARCAS[m])); flush.console()
  Xm <- setmarks(X, v)
  obs_c <- valor(markcorr(Xm, r = r_vec, correction = CORRECCION))
  obs_v <- valor(markvario(Xm, r = r_vec, correction = CORRECCION))
  sims_c <- sims_v <- matrix(NA_real_, length(r_vec), N_SIM)
  for (s in seq_len(N_SIM)) {
    Y <- setmarks(X, sample(v))                 # los mismos puntos, valores barajados
    sims_c[, s] <- valor(markcorr(Y, r = r_vec, correction = CORRECCION))
    sims_v[, s] <- valor(markvario(Y, r = r_vec, correction = CORRECCION))
    if (s %% 100 == 0) { cat("."); flush.console() }
  }
  res_marcas[[paste(cj, m, "correlacion", sep = "|")]] <- bandas(obs_c, sims_c)
  res_marcas[[paste(cj, m, "variograma",  sep = "|")]] <- bandas(obs_v, sims_v)
  saveRDS(res_marcas, "marcas_parteB.rds")      # guarda lo avanzado
  cat(" listo\n")
}

# ====================================================================================
#' #  Marca categórica: tipo de movimiento
#' 
#' ¿Que hace?
#' 
#' Prueba si los eventos de un mismo tipo de movimiento tienden a estar cerca entre sí, 
#' con la función de conexión de marcas normalizada,
#'  para cada tipo y para los tres en conjunto.
#' 
#' La función de conexión normalizada de un tipo t 
#' es la probabilidad de que dos eventos separados por una distancia r sean ambos de tipo t,
#' dividida por la que tendrían si los tipos estuvieran repartidos al azar
#' (la proporción de ese tipo al cuadrado). Vale 1 sin estructura; por encima de 1,
#' los eventos de ese tipo se juntan entre sí.
#' 
#' 
# ====================================================================================
conexion <- function(Y, tipo) valor(markconnect(Y, tipo, tipo, r = r_vec,
                                                correction = CORRECCION, normalise = TRUE))
obs_t  <- sapply(TIPOS, function(ti) conexion(X_tipo, ti))
sims_t <- array(NA_real_, c(length(r_vec), length(TIPOS), N_SIM), dimnames = list(NULL, TIPOS, NULL))
cat("Tipos de movimiento "); flush.console()
for (s in seq_len(N_SIM)) {
  Y <- setmarks(X_tipo, sample(marks(X_tipo)))
  for (ti in TIPOS) sims_t[, ti, s] <- conexion(Y, ti)
  if (s %% 100 == 0) { cat("."); flush.console() }
}
cat(" listo\n")
res_tipos <- lapply(TIPOS, function(ti) bandas(obs_t[, ti], sims_t[, ti, ]))
names(res_tipos) <- TIPOS

# Prueba global conjunta: ¿los tipos, en conjunto, se agrupan con los de su tipo?
ref_t <- apply(sims_t, c(1, 2), mean, na.rm = TRUE)
t_obs <- sum((obs_t[en_int, ] - ref_t[en_int, ])^2, na.rm = TRUE)
t_sim <- apply(sims_t, 3, function(S) sum((S[en_int, ] - ref_t[en_int, ])^2, na.rm = TRUE))
p_tipos_conjunta <- (1 + sum(t_sim >= t_obs)) / (N_SIM + 1)
saveRDS(list(tipos = res_tipos, p_conjunta = p_tipos_conjunta), "tipos_parteB.rds")

# =================================================================================

#' # BLOQUE 5: variograma experimental y direccional de H/L 
#' 
#' Calcula el variograma experimental de $H/L$,
#'  primero en todas las direcciones y luego en cuatro 
#'  (norte–sur, noreste–suroeste, este–oeste y noroeste–sureste), 
#' y lo compara con 999 permutaciones de los valores.
#' 
#' El variograma experimental agrupa los pares de eventos en clases de distancia de 200 m 
#' y calcula, en cada clase h. 
#' 
# ===================================================================================


datos_v <- data.frame(x = coords(X_todos)$x, y = coords(X_todos)$y, HL = marks(X_todos)$HL)
ANCHO <- 200                                     # intervalos de 200 m
DIRECCIONES <- c(0, 45, 90, 135)

variograma_con_banda <- function(alpha = 0) {    # alpha = 0 solo: omnidireccional
  calc <- function(d) variogram(HL ~ 1, locations = ~x + y, data = d,
                                cutoff = R_MAX, width = ANCHO, alpha = alpha)
  obs  <- calc(datos_v)
  sims <- replicate(N_SIM, { d2 <- datos_v; d2$HL <- sample(d2$HL); calc(d2)$gamma })
  obs$lo  <- apply(sims, 1, quantile, 0.025)
  obs$hi  <- apply(sims, 1, quantile, 0.975)
  obs$ref <- rowMeans(sims)
  obs$p_global <- NA
  for (dd in unique(obs$dir.hor)) {
    i <- obs$dir.hor == dd & obs$dist >= R_MIN
    t_obs <- sum((obs$gamma[i] - obs$ref[i])^2)
    t_sim <- colSums((sims[i, , drop = FALSE] - obs$ref[i])^2)
    obs$p_global[obs$dir.hor == dd] <- (1 + sum(t_sim >= t_obs)) / (N_SIM + 1)
  }
  obs
}
cat("Variograma omnidireccional...\n");  vario_omni <- variograma_con_banda()
cat("Variogramas direccionales...\n");   vario_dir  <- variograma_con_banda(DIRECCIONES)
saveRDS(list(omni = vario_omni, direccional = vario_dir), "variogramas_parteB.rds")


# =============================================================================

#' # BLOQUE 6: resumen
#' 
#' ¿Que hace?
#' 
#' Escribe los resultados en un .txt , y muestra en la consola,
#'  cómo leer cada función y, para cada análisis, el valor p y 
#'  los tramos de distancia donde la curva sale de las bandas.
#' 
#' 
# ================================================================================


tramos <- function(r, cond) {
  cond[is.na(cond)] <- FALSE
  if (!any(cond)) return("ninguno")
  rl <- rle(cond); fin <- cumsum(rl$lengths); ini <- fin - rl$lengths + 1
  paste(sprintf("%.0f a %.0f m", r[ini[rl$values]], r[fin[rl$values]]), collapse = "; ")
}
resumir <- function(res, titulo) {
  d <- res$tabla; d <- d[en_int, ]
  cat("==", titulo, "==\n")
  cat(sprintf("  Prueba global DCLF (%d-%d m): p = %.3f\n", R_MIN, R_MAX, res$p_DCLF))
  cat("  Fuera de la banda GLOBAL por encima:", tramos(d$r, d$obs > d$ghi), "\n")
  cat("  Fuera de la banda GLOBAL por debajo:", tramos(d$r, d$obs < d$glo), "\n")
  cat("  Fuera de la banda PUNTUAL por encima:", tramos(d$r, d$obs > d$hi), "\n")
  cat("  Fuera de la banda PUNTUAL por debajo:", tramos(d$r, d$obs < d$lo), "\n\n")
}

sink("resultados_parteB.txt")
cat("FASE 2 - PARTE B\n")
cat(sprintf("%d permutaciones | bandas del 95 %% | interpretacion de %d a %d m | correccion: %s\n\n",
            N_SIM, R_MIN, R_MAX, CORRECCION))
cat("COMO LEER:\n")
cat("  Correlacion de marcas (referencia 1): por ENCIMA = vecinos con valores altos a la vez.\n")
cat("  Variograma de marcas (referencia = varianza): por DEBAJO = vecinos parecidos.\n")
cat("  Conexion de tipos (referencia 1): por ENCIMA = los de ese tipo se juntan entre si.\n\n")
cat("---------------- MARCAS NUMERICAS ----------------\n\n")
for (k in names(res_marcas)) {
  p <- strsplit(k, "\\|")[[1]]
  resumir(res_marcas[[k]], sprintf("%s | %s | %s", p[1], MARCAS[p[2]], p[3]))
}
cat("---------------- TIPO DE MOVIMIENTO ----------------\n\n")
cat(sprintf("Prueba global conjunta de los tres tipos: p = %.3f\n\n", p_tipos_conjunta))
for (ti in TIPOS) resumir(res_tipos[[ti]], paste(ti, "con", tolower(ti)))
cat("---------------- VARIOGRAMA DE H/L (gstat) ----------------\n\n")
cat(sprintf("Omnidireccional: p global = %.3f\n", vario_omni$p_global[1]))
cat("  Por debajo de la banda (vecinos parecidos):",
    tramos(vario_omni$dist, vario_omni$gamma < vario_omni$lo & vario_omni$dist >= R_MIN), "\n")
for (dd in DIRECCIONES) {
  v <- vario_dir[vario_dir$dir.hor == dd, ]
  cat(sprintf("Direccion %3d grados: p global = %.3f\n", dd, v$p_global[1]))
  cat("  Por debajo de la banda:", tramos(v$dist, v$gamma < v$lo & v$dist >= R_MIN), "\n")
}
sink()
cat(readLines("resultados_parteB.txt"), sep = "\n")


# =============================================================================
#' # BLOQUE 7: figuras 
#' 
#' Dibuja cinco figuras, las muestra en el panel Plots de RStudio y las guarda en PNG.
#' 
#' 1. Correlación de marcas: una fila por marca y una columna por conjunto (seis paneles).
#' 
#' 2. Variograma de marcas, con la misma organización.
#' 
#' 3. Conexión de cada tipo consigo mismo (tres paneles)
#' 
#' 4. Variograma omnidireccional de H/L.
#' 
#' 5. Variograma de H/L en las cuatro direcciones
#' 
# ==============================================================================

leyenda_inferior <- function(textos, lwd, lty, col) {
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  legend("bottom", legend = textos, lwd = lwd, lty = lty, col = col,
         horiz = TRUE, bty = "n", cex = 0.8, xpd = NA)
}
leyenda_bandas <- function() leyenda_inferior(
  c("Curva de los eventos", "Referencia sin estructura", "Azar: banda puntual 95 %",
    "Azar: banda global 95 %", "Inicio de interpretación (100 m)"),
  lwd = c(2, 1, 8, 1, 1), lty = c(1, 3, 1, 2, 3), col = c("black", "red", "grey85", "grey30", "grey50"))

panel_banda <- function(d, titulo, ylab) {
  yl <- range(unlist(d[en_int, c("obs", "lo", "hi", "glo", "ghi")]), na.rm = TRUE)
  plot(NA, xlim = c(0, R_MAX), ylim = yl, xlab = "Distancia r (m)", ylab = ylab,
       main = titulo, cex.main = 1)
  polygon(c(d$r, rev(d$r)), c(d$lo, rev(d$hi)), col = "grey85", border = NA)
  lines(d$r, d$ghi, lty = 2, col = "grey30"); lines(d$r, d$glo, lty = 2, col = "grey30")
  lines(d$r, d$ref, lty = 3, col = "red"); abline(v = R_MIN, lty = 3, col = "grey50")
  lines(d$r, d$obs, lwd = 2)
}

fig_marcas <- function(funcion, ylab) {
  par(mfrow = c(3, 2), mar = c(4.5, 4.5, 2.5, 1), oma = c(3, 0, 0, 0))
  letras <- letters[1:6]; i <- 0
  for (m in names(MARCAS)) for (cj in names(conjuntos)) {
    i <- i + 1
    panel_banda(res_marcas[[paste(cj, m, funcion, sep = "|")]]$tabla,
                sprintf("(%s) %s · %s", letras[i], MARCAS[m], cj), ylab)
  }
  leyenda_bandas()
}
fig_correlacion  <- function() fig_marcas("correlacion", "Correlación de marcas")
fig_vario_marcas <- function() fig_marcas("variograma", "Variograma de marcas")

fig_tipos <- function() {
  par(mfrow = c(1, 3), mar = c(4.5, 4.5, 2.5, 1), oma = c(3, 0, 0, 0))
  for (i in seq_along(TIPOS))
    panel_banda(res_tipos[[TIPOS[i]]]$tabla,
                sprintf("(%s) %s con %s", letters[i], TIPOS[i], tolower(TIPOS[i])),
                "Conexión de marcas normalizada")
  leyenda_bandas()
}

panel_vario <- function(v, titulo) {
  plot(v$dist, v$gamma, type = "n", ylim = range(c(v$gamma, v$lo, v$hi)),
       xlab = "Distancia (m)", ylab = "Semivarianza de H/L", main = titulo, cex.main = 1)
  polygon(c(v$dist, rev(v$dist)), c(v$lo, rev(v$hi)), col = "grey85", border = NA)
  lines(v$dist, v$ref, lty = 3, col = "red"); abline(v = R_MIN, lty = 3, col = "grey50")
  lines(v$dist, v$gamma, type = "b", pch = 16)
}
leyenda_vario <- function() leyenda_inferior(
  c("Semivarianza de los eventos", "Referencia sin estructura", "Azar: banda 95 %",
    "Inicio de interpretación (100 m)"),
  lwd = c(2, 1, 8, 1), lty = c(1, 3, 1, 3), col = c("black", "red", "grey85", "grey50"))
fig_vario_omni <- function() {
  par(mfrow = c(1, 1), mar = c(4.5, 4.5, 2.5, 1), oma = c(3, 0, 0, 0))
  panel_vario(vario_omni, "Variograma de H/L (omnidireccional)"); leyenda_vario()
}
fig_vario_dir <- function() {
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2.5, 1), oma = c(3, 0, 0, 0))
  nombres <- c("0" = "norte-sur", "45" = "noreste-suroeste", "90" = "este-oeste", "135" = "noroeste-sureste")
  for (i in seq_along(DIRECCIONES)) {
    dd <- DIRECCIONES[i]
    panel_vario(vario_dir[vario_dir$dir.hor == dd, ],
                sprintf("(%s) %d grados (%s)", letters[i], dd, nombres[as.character(dd)]))
  }
  leyenda_vario()
}

# En pantalla (flechas del panel Plots para pasar de una a otra)
fig_correlacion(); fig_vario_marcas(); fig_tipos(); fig_vario_omni(); fig_vario_dir()

# En archivo
png("correlacion_marcas_parteB.png", width = 2200, height = 2600, res = 220); fig_correlacion(); invisible(dev.off())
png("variograma_marcas_parteB.png", width = 2200, height = 2600, res = 220); fig_vario_marcas(); invisible(dev.off())
png("conexion_tipos_parteB.png", width = 3300, height = 1200, res = 220); fig_tipos(); invisible(dev.off())
png("variograma_HL_parteB.png", width = 1600, height = 1300, res = 220); fig_vario_omni(); invisible(dev.off())
png("variograma_HL_direccional_parteB.png", width = 2200, height = 2000, res = 220); fig_vario_dir(); invisible(dev.off())

cat("Archivos guardados en:", getwd(), "\n")
