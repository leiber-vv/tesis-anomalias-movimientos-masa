# =============================================================================
#' # Fase 6 - Parte B: Validación e incertidumbre
#'
#' Este script tiene el objetivo de responder: ¿las 19 anomalías de la Fase 3 corresponden a procesos reales o a artefactos del inventario, y qué tan estable es su identificación?
#'
#' Tiene cinco partes:
#'
#' | Parte | Qué hace |
#' |---|---|
#' | A | Prueba si las anomalías se concentran por fuente, por completitud o en los eventos anteriores al DEM de 2007 (Fisher). |
#' | B | Resume la revisión en imágenes: la lectura de cada anomalía y el desfase entre el año del evento y la imagen usada. |
#' | C | Sensibilidad: cambia una decisión a la vez (h del MCD, umbral, corte adaptativo y covarianza por tercil de pendiente). |
#' | D | Remuestreo: 500 submuestras sin reemplazo del 80 %, reajustando la GWR mixta y el MCD en cada una. |
#' | E | Arma la Tabla 4 y marca las anomalías confirmadas: proceso real y estables. |
#'
#' **Datos de entrada**
#'
#' | Archivo | Contenido |
#' |---|---|
#' | `F3_fase3_integrado.csv` | Fase 3: residuales estudentizados, distancia robusta, anomalía, fuente y completitud de cada evento. |
#' | `F3_S2_tabla_modelo.csv` | Fase 3: variables de la GWR mixta. |
#' | `F6_lecturas.csv` | Lectura de cada anomalía y fecha de la imagen usada en Google Earth Pro. Sale de `F6_planilla_revision.xlsx`. |
#' | `tabla_inventario.xls` | Tabla de atributos de los 453 polígonos, de donde sale el año de cada evento. |
#'
#' **Salidas** (carpeta `salidas`)
#'
#' | Archivo | Contenido |
#' |---|---|
#' | `F6_resumen.txt` | Las cifras que se usan en la redacción. |
#' | `F6_C_sensibilidad.csv` | Si cada evento es anómalo en cada escenario de sensibilidad. |
#' | `F6_C_tercil.csv` | Varianzas y distancia riemanniana de cada tercil de pendiente. |
#' | `F6_D_reproduccion.csv` | Residuales reproducidos frente a los de la Fase 3. |
#' | `F6_D_remuestreo.csv` | Veces que cada evento entró a una submuestra, veces que resultó anómalo y su P. |
#' | `F6_D_remuestreo_nanom.csv` | Número de anomalías en cada submuestra. |
#' | `F6_tabla4.csv` | La Tabla 4 de la tesis. |
#'
#' **Requisitos:** GWmodel y robustbase, los mismos de la Fase 3, y readxl solo si el inventario se lee en `.xls`. Tarda de 3 a 4 minutos.


#==============================================================================

#' # Paquetes
#'
#' GWmodel da los pesos de la GWR y robustbase el MCD.

#------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GWmodel)
  library(robustbase)
})


#==============================================================================

#' # Configuración
#'
#' Reúne las rutas y los parámetros de la fase. Es lo único que hay que editar.

#------------------------------------------------------------------------------

CARPETA_F3 <- "D:/Users/tasal/Desktop/Tesis de grado/Fase 3"   # CSV de la Fase 3
CARPETA_F6 <- "D:/Users/tasal/Desktop/Tesis de grado/Fase 6"   # F6_lecturas.csv
SALIDA     <- file.path(CARPETA_F6, "salidas")

ARCH_INTEGRADO <- file.path(CARPETA_F3, "F3_fase3_integrado.csv")
ARCH_MODELO    <- file.path(CARPETA_F3, "F3_S2_tabla_modelo.csv")
ARCH_LECTURAS  <- file.path(CARPETA_F6, "F6_lecturas.csv")
ARCH_INVENT    <- file.path(CARPETA_F6, "tabla_inventario.xls")   # o "tabla_inventario.csv"
ANIO_DEM       <- 2007                # microzonificacion sismica (AMVA, 2007)

H_BASE   <- 0.75                      # fraccion del MCD de la Fase 3
H_ALT    <- c(0.50, 0.90)             # sensibilidad de h
BW       <- c(AREA = 129, HL = 228, ELONG = 290)  # anchos de banda (vecinos) de la Fase 3, por AICc
B_PERM   <- 999                       # permutaciones de la prueba por tercil
B_SUB    <- 500                       # submuestras del remuestreo
FRAC_SUB <- 0.80                      # fraccion de eventos en cada submuestra
CORTE_P  <- 0.75                      # estabilidad (mismo corte de la Fase 4)
SEMILLA  <- 2026


#==============================================================================

#' # Reporte y funciones comunes
#'
#' | Función | Qué hace |
#' |---|---|
#' | `msg(...)` | Imprime un texto y lo agrega a `F6_resumen.txt`. |
#' | `leer(f)` | Lee un CSV quitando la marca invisible del comienzo del archivo (BOM). |
#' | `z(v)` | Estandariza un vector. |
#' | `mcd(X, h, semilla)` | Corre el MCD con una semilla fija, para que cada llamada dé lo mismo, y devuelve el generador aleatorio a su estado anterior para no alterar las permutaciones ni las submuestras. |
#' | `dist_rob(X, m)` | Distancia robusta de cada evento al centro del MCD. |
#'
#' `UMBRAL_975` y `UMBRAL_99` son los cortes de anomalía de la Fase 3: la raíz de los cuantiles 0,975 y 0,99 de una $\chi^2$ con 3 grados de libertad.

#------------------------------------------------------------------------------

dir.create(SALIDA, showWarnings = FALSE, recursive = TRUE)
RESUMEN <- file.path(SALIDA, "F6_resumen.txt")
cat("FASE 6 - VALIDACION E INCERTIDUMBRE\n", format(Sys.time()), "\n\n", file = RESUMEN)
msg <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n")
  cat(txt, "\n", file = RESUMEN, append = TRUE)
}

leer <- function(f) read.csv(f, fileEncoding = "UTF-8-BOM", check.names = FALSE)
z    <- function(v) (v - mean(v)) / sd(v)
UMBRAL_975 <- sqrt(qchisq(0.975, 3))
UMBRAL_99  <- sqrt(qchisq(0.990, 3))

# MCD con semilla fija para que cada llamada sea reproducible. Guarda y restaura el
# estado del generador aleatorio para no alterar las permutaciones ni las submuestras.
mcd <- function(X, h = H_BASE, semilla = 1) {
  if (exists(".Random.seed", envir = globalenv())) {
    estado <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", estado, envir = globalenv()))
  }
  set.seed(semilla)
  covMcd(X, alpha = h)
}
dist_rob <- function(X, m) sqrt(mahalanobis(X, m$center, m$cov))


#==============================================================================

#' # Datos
#'
#' ### Qué hace
#'
#' Une en una tabla los resultados de la Fase 3, las variables del modelo y el año de cada evento.
#'
#' **Cómo lo hace**
#'
#' 1. Lee los cuatro archivos y deja solo los eventos con `EN_MODELO = 1`.
#' 2. Toma el año del inventario. La columna se busca como `A.o` para que la ñ de «Año» no cause problemas de codificación.
#' 3. Une todo por el identificador y arma `R`, la matriz con los tres residuales estudentizados de la Fase 3, y `anom`, que marca las 19 anomalías.

#------------------------------------------------------------------------------

f3  <- leer(ARCH_INTEGRADO)
mod <- leer(ARCH_MODELO)
mod <- mod[mod$EN_MODELO == 1, ]
lec <- leer(ARCH_LECTURAS)
inv <- if (grepl("\\.xlsx?$", ARCH_INVENT)) as.data.frame(readxl::read_excel(ARCH_INVENT)) else leer(ARCH_INVENT)
anio <- inv[, grep("^A.o$", names(inv))[1]]                 # columna Año del inventario
f3$ANIO <- anio[match(f3$ID, inv$ID_evento)]

d <- merge(f3, mod[, c("ID", "LOG_AREA", "LOG_HL", "LOG_ELONG", "PENDIENTE",
                       "LOG_DIST_ESTRUCTURA", "LOG_DIST_DRENAJE", "DENS_LIN_1000",
                       "GEOLOGIA", "GEOMORFOLOGIA", "COBERTURA")], by = "ID")
d <- d[order(d$ID), ]
n <- nrow(d)
R <- as.matrix(d[, c("RE_AREA", "RE_HL", "RE_ELONG")])   # residuales estudentizados de la Fase 3
anom <- d$ANOMALIA != "normal"
msg("Eventos: ", n, " | anomalias de la Fase 3: ", sum(anom),
    " (", sum(d$ANOMALIA == "fuerte"), " fuertes)")


#==============================================================================

#' # Parte A. Concentración de las anomalías
#'
#' ### Qué hace
#'
#' Con la prueba exacta de Fisher, prueba si las anomalías son más frecuentes en una fuente (mapeo o institucional), entre los eventos bajo el umbral de completitud o entre los eventos anteriores al DEM de 2007.
#'
#' El DEM es el de la microzonificación sísmica (AMVA, 2007). En un evento anterior, la topografía del DEM ya incluye la cicatriz o el depósito, así que H y L se miden sobre terreno modificado. Al final lista las anomalías anteriores al DEM.

#------------------------------------------------------------------------------

msg("\n--- PARTE A. Concentracion de las anomalias ---")
for (v in c("GRUPO", "BAJO_UMBRAL")) {
  t <- table(d[[v]], anom)
  p <- fisher.test(t)$p.value
  msg(v, ": ", paste(sprintf("%s %d/%d", rownames(t), t[, "TRUE"], rowSums(t)), collapse = "; "),
      " | Fisher p = ", sprintf("%.3f", p))
}
# Eventos anteriores al DEM: su topografia ya incluye la cicatriz o el deposito, de modo que
# H y L se miden sobre terreno modificado.
con_anio <- !is.na(d$ANIO)
previo <- d$ANIO < ANIO_DEM
t_dem <- table(previo[con_anio], anom[con_anio])
msg(sprintf("Anteriores al DEM (%d): %d anomalias de %d; posteriores: %d de %d; sin ano: %d | Fisher p = %.3f",
            ANIO_DEM, t_dem["TRUE", "TRUE"], sum(t_dem["TRUE", ]), t_dem["FALSE", "TRUE"], sum(t_dem["FALSE", ]),
            sum(!con_anio), fisher.test(t_dem)$p.value))
msg("Anomalias anteriores al DEM: ", paste(d$ID[anom & con_anio & previo], collapse = ", "))


#==============================================================================

#' # Parte B. Lectura de las 19 anomalías
#'
#' Cada anomalía se revisó en las imágenes históricas de Google Earth Pro, respondiendo seis preguntas en la planilla de revisión:
#'
#' | Columna | Pregunta | Respuestas posibles |
#' |---|---|---|
#' | P1 | ¿Es un solo evento? | Sí, No, No se puede ver |
#' | P2 | ¿Incluye arranque y depósito? | Sí, No, No se puede ver |
#' | P3 | ¿Incluye algo ajeno (cauce, vía, otro depósito)? | Sí, No, No se puede ver |
#' | P4 | ¿La corona está en la parte alta de la cicatriz? | Sí, No, No se puede ver |
#' | P5 | Si hay un error, ¿empuja la métrica en el sentido de la anomalía? | Sí, No, No aplica (no hay error) |
#' | P6 | ¿Se ve en la imagen el rasgo que la hace anómala? | Sí, No, No se puede ver |
#'
#' P1 a P4 revisan el polígono y sirven para responder P5. La planilla asigna la lectura con una fórmula, en este orden:
#'
#' 1. P5 = Sí: **error de inventario**.
#' 2. Evento bajo el umbral de completitud: **sesgo de registro**.
#' 3. P6 = Sí: **proceso real**.
#' 4. Polígono digitalizado desde un punto de la fuente: **sesgo de registro**.
#' 5. Ninguna de las anteriores: **indeterminada**.
#'
#' ### Qué hace
#'
#' 1. Comprueba que `F6_lecturas.csv` tenga exactamente las 19 anomalías.
#' 2. Cruza el perfil de cada anomalía con su lectura.
#' 3. Compara el año de la imagen usada con el año del evento: cuántas revisiones usaron una imagen del mismo año o del siguiente, cuáles no tuvieron imagen útil y cuáles tuvieron un desfase mayor.

#------------------------------------------------------------------------------

msg("\n--- PARTE B. Lectura de las 19 anomalias ---")
stopifnot(setequal(lec$ID, d$ID[anom]))
lb <- merge(lec, d[, c("ID", "PERFIL")], by = "ID")
tb <- table(lb$PERFIL, lb$lectura)
msg(paste(capture.output(print(addmargins(tb))), collapse = "\n"))
# Fecha de la imagen usada frente al ano del evento
li <- merge(lec, d[, c("ID", "ANIO")], by = "ID")
li$desfase <- as.integer(substr(li$fecha_imagen, 1, 4)) - li$ANIO
msg(sprintf("Imagen del mismo ano o del siguiente: %d | sin imagen util: %d (%s) | desfase mayor: %s",
            sum(li$desfase %in% 0:1), sum(is.na(li$desfase)), paste(li$ID[is.na(li$desfase)], collapse = ", "),
            paste(sprintf("%d (%d anos)", li$ID[!is.na(li$desfase) & li$desfase > 1],
                          li$desfase[!is.na(li$desfase) & li$desfase > 1]), collapse = ", ")))


#==============================================================================

#' # Parte C. Sensibilidad
#'
#' ### Qué hace
#'
#' Cambia una decisión de la Fase 3 a la vez, sobre los mismos residuales, y mira qué anomalías se mantienen.
#'
#' | Escenario | Cambio |
#' |---|---|
#' | `h0.50` y `h0.90` | Fracción del MCD de 0,50 y de 0,90 en lugar de 0,75. |
#' | `u99` | Umbral del cuantil 0,99 en lugar del 0,975 (solo quedan las fuertes). |
#' | `adapt` | Corte adaptativo de Filzmoser et al. (2005), de la Fase 3. |
#' | `tercil` | Una covarianza para cada tercil de pendiente en lugar de una global. |
#'
#' `escenarios` cuenta en cuántos de los cinco escenarios sigue siendo anómalo cada evento. Con h = 0,75 el script comprueba que salen las mismas anomalías de la Fase 3.
#'
#' **Covarianza por tercil de pendiente**
#'
#' La diferencia entre la covarianza de un tercil ($B$) y la global ($A$) se mide con la distancia riemanniana afín-invariante (Pennec et al., 2006):
#'
#' $$d(A, B) = \sqrt{\sum_i \log^2 \lambda_i} ,$$
#'
#' donde $\lambda_i$ son los valores propios de $A^{-1}B$. Vale 0 si las dos matrices son iguales. El estadístico suma las distancias de los tres terciles, y la p sale de 999 permutaciones de la etiqueta de tercil. Después recalcula la distancia robusta de cada evento con la covarianza de su tercil.

#------------------------------------------------------------------------------

msg("\n--- PARTE C. Sensibilidad ---")
sens <- data.frame(ID = d$ID, base = anom)

# C.1 Fraccion h del MCD
for (h in c(H_ALT[1], H_BASE, H_ALT[2])) {
  a <- dist_rob(R, mcd(R, h)) > UMBRAL_975
  sens[[sprintf("h%.2f", h)]] <- a
  msg(sprintf("h = %.2f n: %d anomalias; se mantienen %d de %d; nuevas: %s",
              h, sum(a), sum(a & anom), sum(anom), paste(d$ID[a & !anom], collapse = ", ")))
}
stopifnot(all(sens$h0.75 == anom))   # h = 0,75 reproduce la Fase 3

# C.2 Umbral chi2(3; 0,99) y corte adaptativo de Filzmoser et al. (2005), de la Fase 3
sens$u99   <- d$ANOMALIA == "fuerte"
sens$adapt <- d$ADAPTATIVO == 1
msg("Umbral 0,99: ", sum(sens$u99), " | corte adaptativo: ", sum(sens$adapt))

# C.3 Covarianza por tercil de pendiente
#     Diferencia con la covarianza global: suma de distancias riemannianas
#     afin-invariantes (Pennec et al., 2006) d(A,B) = sqrt(sum(log(eig(A^-1 B))^2)),
#     con 999 permutaciones de la etiqueta de tercil.
q   <- quantile(d$PENDIENTE, c(1/3, 2/3))
ter <- cut(d$PENDIENTE, c(-Inf, q, Inf), labels = c("bajo", "medio", "alto"))
airm <- function(A, B) sqrt(sum(log(Re(eigen(solve(A, B), only.values = TRUE)$values))^2))
cov_glob <- mcd(R)$cov
estad <- function(et) sum(sapply(levels(et), function(k) airm(cov_glob, mcd(R[et == k, ])$cov)))
S_obs <- estad(ter)
set.seed(SEMILLA)
S_perm <- replicate(B_PERM, estad(factor(sample(as.character(ter)), levels = levels(ter))))
p_ter <- (sum(S_perm >= S_obs) + 1) / (B_PERM + 1)

rd_ter <- numeric(n)
tab_ter <- NULL
for (k in levels(ter)) {
  i <- ter == k
  mk <- mcd(R[i, ])
  rd_ter[i] <- dist_rob(R[i, ], mk)
  tab_ter <- rbind(tab_ter, data.frame(tercil = k,
                                       pendiente = sprintf("%.1f-%.1f", min(d$PENDIENTE[i]), max(d$PENDIENTE[i])),
                                       eventos = sum(i), anomalias_F3 = sum(anom & i),
                                       var_area = diag(mk$cov)[1], var_HL = diag(mk$cov)[2], var_elong = diag(mk$cov)[3],
                                       d_riemann = airm(cov_glob, mk$cov)))
}
tab_ter <- rbind(tab_ter, data.frame(tercil = "global", pendiente = "", eventos = n, anomalias_F3 = sum(anom),
                                     var_area = diag(cov_glob)[1], var_HL = diag(cov_glob)[2],
                                     var_elong = diag(cov_glob)[3], d_riemann = NA))
sens$tercil <- rd_ter > UMBRAL_975
write.csv(tab_ter, file.path(SALIDA, "F6_C_tercil.csv"), row.names = FALSE)
msg("Terciles de pendiente (cortes ", paste(round(q, 1), collapse = " y "), " grados):")
msg(paste(capture.output(print(tab_ter, digits = 3, row.names = FALSE)), collapse = "\n"))
msg(sprintf("Estadistico = %.3f; p (permutaciones) = %.3f", S_obs, p_ter))
msg(sprintf("Con la covarianza de cada tercil: %d anomalias; se mantienen %d de %d (tercil bajo: %d de %d); salen: %s; nuevas: %s",
            sum(sens$tercil), sum(sens$tercil & anom), sum(anom),
            sum(sens$tercil & anom & ter == "bajo"), sum(anom & ter == "bajo"),
            paste(d$ID[anom & !sens$tercil], collapse = ", "),
            paste(d$ID[!anom & sens$tercil], collapse = ", ")))

sens$escenarios <- rowSums(sens[, c("h0.50", "h0.90", "u99", "adapt", "tercil")])
sens$tercil_pendiente <- ter
write.csv(sens, file.path(SALIDA, "F6_C_sensibilidad.csv"), row.names = FALSE)


#==============================================================================

#' # Parte D. Remuestreo
#'
#' ### Qué hace
#'
#' Mide qué tan estable es cada anomalía frente a la muestra. Repite la GWR mixta y el MCD en 500 submuestras sin reemplazo del 80 % de los eventos (Politis et al., 1999); sin reemplazo, un evento no se repite dentro de su propio vecindario. Los anchos de banda se multiplican por 0,8, porque con menos eventos el mismo número de vecinos abarcaría más área.
#'
#' **Nota sobre la GWR mixta**
#'
#' `gwr.mixed` no entrega la diagonal de la matriz sombrero, que hace falta para estudentizar los residuales. Por eso el script calcula la matriz sombrero exacta del modelo mixto (Fotheringham et al., 2002), con los pesos bisquare adaptativos de `gw.weight`:
#'
#' $$S = S_b + (I - S_b) X_2 \left[X_2^\top (I - S_b)^\top (I - S_b) X_2\right]^{-1} X_2^\top (I - S_b)^\top (I - S_b) ,$$
#'
#' donde $S_b$ es la matriz sombrero de la GWR local con el intercepto y las cuatro continuas, y $X_2$ son las indicadoras de las categóricas con coeficiente global. Los residuales se estudentizan como en `gwr.basic`: $e_i / \sqrt{\sigma^2 (1 - S_{ii})}$, con $\sigma^2 = RSS / (n - 2\,\mathrm{tr}(S) + \mathrm{tr}(S^\top S))$.
#'
#' | Función | Qué hace |
#' |---|---|
#' | `mixed_hat(X1, X2, D, bw)` | Calcula la matriz sombrero de la GWR mixta. |
#' | `estudentizar(y, S)` | Calcula los residuales estudentizados. |
#' | `ref_mas_frecuente(v)` | Toma como referencia la clase más frecuente de cada categórica. |
#' | `ajustar(idx, escala_bw)` | Ajusta las tres dimensiones con los eventos `idx` y devuelve sus residuales estudentizados. |
#'
#' **Cómo lo hace**
#'
#' 1. **Verificación (D.1).** Ajusta con todos los eventos y lo compara con la Fase 3: correlación y diferencia máxima de los residuales, y cuántas de las 19 anomalías se reproducen.
#' 2. **Submuestras (D.2).** En cada submuestra ajusta la GWR mixta y el MCD y marca las anomalías.
#' 3. Para cada evento calcula $P$, la fracción de las submuestras que lo incluyeron en que resultó anómalo. Es estable si $P \geq 0{,}75$, el corte de la Fase 4.

#------------------------------------------------------------------------------

msg("\n--- PARTE D. Remuestreo ---")

# Matriz sombrero exacta de la GWR mixta (Fotheringham et al., 2002):
#   S = Sb + (I - Sb) X2 [X2'(I-Sb)'(I-Sb)X2]^-1 X2'(I-Sb)'(I-Sb)
#   Sb: matriz sombrero de la GWR local con X1 (intercepto + 4 continuas)
#   X2: indicadoras de las categoricas con coeficiente global
mixed_hat <- function(X1, X2, D, bw) {
  m <- nrow(X1)
  W <- gw.weight(D, bw, kernel = "bisquare", adaptive = TRUE)
  Sb <- t(sapply(seq_len(m), function(i) {
    XtW <- t(X1 * W[, i])
    X1[i, ] %*% solve(XtW %*% X1, XtW)
  }))
  IS <- diag(m) - Sb
  X3 <- IS %*% X2
  Sb + IS %*% X2 %*% solve(crossprod(X3), t(X3) %*% IS)
}
estudentizar <- function(y, S) {
  e <- as.vector(y - S %*% y)
  m <- length(y)
  s2 <- sum(e^2) / (m - 2 * sum(diag(S)) + sum(S * S))
  e / sqrt(s2 * (1 - diag(S)))
}
ref_mas_frecuente <- function(v) { v <- factor(v); relevel(v, ref = names(which.max(table(v)))) }

CONT <- c("PENDIENTE", "LOG_DIST_ESTRUCTURA", "LOG_DIST_DRENAJE", "DENS_LIN_1000")
X2_todo <- model.matrix(~ ref_mas_frecuente(d$GEOLOGIA) + ref_mas_frecuente(d$GEOMORFOLOGIA) +
                          ref_mas_frecuente(d$COBERTURA) + ref_mas_frecuente(d$TIPO) +
                          ref_mas_frecuente(d$GRUPO))[, -1]
RESP <- c(AREA = "LOG_AREA", HL = "LOG_HL", ELONG = "LOG_ELONG")
D_todo <- gw.dist(dp.locat = cbind(d$X, d$Y))

ajustar <- function(idx, escala_bw = 1) {
  X1 <- cbind(1, apply(as.matrix(d[idx, CONT]), 2, z))
  X2 <- X2_todo[idx, , drop = FALSE]
  X2 <- X2[, colSums(X2) > 0, drop = FALSE]
  D  <- D_todo[idx, idx]
  sapply(names(RESP), function(v) {
    S <- mixed_hat(X1, X2, D, round(escala_bw * BW[[v]]))
    estudentizar(z(d[idx, RESP[[v]]]), S)
  })
}

# D.1 Verificacion: reproducir la Fase 3 con todos los eventos
RE_rep <- ajustar(seq_len(n))
rd_rep <- dist_rob(RE_rep, mcd(RE_rep))
rep_df <- data.frame(ID = d$ID, RE_rep, RE_F3 = R, RD_rep = rd_rep, RD_F3 = d$RD)
write.csv(rep_df, file.path(SALIDA, "F6_D_reproduccion.csv"), row.names = FALSE)
for (j in 1:3) msg(sprintf("Reproduccion %s: r = %.5f; max |dif| = %.3f",
                           names(RESP)[j], cor(RE_rep[, j], R[, j]), max(abs(RE_rep[, j] - R[, j]))))
msg(sprintf("Anomalias reproducidas: %d; coinciden con la Fase 3: %d de %d",
            sum(rd_rep > UMBRAL_975), sum(rd_rep > UMBRAL_975 & anom), sum(anom)))

# D.2 Submuestras sin reemplazo (Politis et al., 1999); anchos escalados a FRAC_SUB
k <- round(FRAC_SUB * n)
incluido <- marcado <- integer(n)
n_anom <- integer(B_SUB)
set.seed(SEMILLA)
t0 <- Sys.time()
for (b in seq_len(B_SUB)) {
  s <- sort(sample.int(n, k))
  REs <- ajustar(s, escala_bw = FRAC_SUB)
  mc <- covMcd(REs, alpha = H_BASE)
  a <- dist_rob(REs, mc) > UMBRAL_975
  incluido[s] <- incluido[s] + 1L
  marcado[s]  <- marcado[s] + a
  n_anom[b]   <- sum(a)
  if (b %% 50 == 0) cat(sprintf("  submuestra %d de %d (%s)\n", b, B_SUB, format(Sys.time() - t0, digits = 3)))
}
rem <- data.frame(ID = d$ID, RD = d$RD, ANOMALIA = d$ANOMALIA,
                  incluido = incluido, marcado = marcado, P = marcado / pmax(incluido, 1))
write.csv(rem, file.path(SALIDA, "F6_D_remuestreo.csv"), row.names = FALSE)
write.csv(data.frame(submuestra = seq_len(B_SUB), n_anomalias = n_anom),
          file.path(SALIDA, "F6_D_remuestreo_nanom.csv"), row.names = FALSE)

Pa <- rem$P[anom]
msg(sprintf("Submuestras de %d eventos: mediana de %d anomalias (%.1f %%); intervalo del 90 %%: %d a %d; esperadas por azar: %.1f",
            k, median(n_anom), 100 * median(n_anom) / k,
            quantile(n_anom, 0.05), quantile(n_anom, 0.95), 0.025 * k))
msg(sprintf("Anomalias con P = 1: %d | con P >= %.2f: %d | inestables: %s (P de %.2f a %.2f)",
            sum(Pa == 1), CORTE_P, sum(Pa >= CORTE_P),
            paste(rem$ID[anom & rem$P < CORTE_P], collapse = ", "),
            min(Pa[Pa < CORTE_P]), max(Pa[Pa < CORTE_P])))
msg(sprintf("Eventos normales con P >= 0,5: %d (P maxima: %.2f)",
            sum(!anom & rem$P >= 0.5), max(rem$P[!anom])))


#==============================================================================

#' # Parte E. Tabla 4 y confirmación
#'
#' ### Qué hace
#'
#' Reúne para cada anomalía su fuente, tipo, año, perfil, distancia robusta, escenarios de sensibilidad, $P$ y lectura. Una anomalía se confirma si su lectura es proceso real y además es estable ($P \geq 0{,}75$). Guarda `F6_tabla4.csv` e imprime cuántas anomalías hay de cada lectura, los errores de inventario por tercil de pendiente y las confirmadas.

#------------------------------------------------------------------------------

msg("\n--- PARTE E. Tabla 4 ---")
t4 <- merge(d[anom, c("ID", "GRUPO", "TIPO", "ANIO", "PERFIL", "RD", "ANOMALIA")],
            sens[, c("ID", "escenarios", "tercil_pendiente")], by = "ID")
t4 <- merge(t4, rem[, c("ID", "P")], by = "ID")
t4 <- merge(t4, lec, by = "ID")
t4$estable    <- t4$P >= CORTE_P
t4$confirmada <- t4$lectura == "Proceso real" & t4$estable
t4 <- t4[order(-t4$RD), ]
t4$RD <- round(t4$RD, 2); t4$P <- round(t4$P, 2)
write.csv(t4, file.path(SALIDA, "F6_tabla4.csv"), row.names = FALSE)
msg(paste(capture.output(print(t4[, c("ID", "GRUPO", "TIPO", "ANIO", "RD", "escenarios", "P", "lectura", "confirmada")],
                               row.names = FALSE)), collapse = "\n"))
msg("Lecturas: ", paste(names(table(t4$lectura)), table(t4$lectura), sep = " = ", collapse = "; "))
msg("Errores de inventario por tercil de pendiente: ",
    paste(names(table(t4$tercil_pendiente[t4$lectura == "Error de inventario"])),
          table(t4$tercil_pendiente[t4$lectura == "Error de inventario"]), sep = " = ", collapse = "; "))
msg("Confirmadas: ", paste(t4$ID[t4$confirmada], collapse = ", "))
msg("\nListo. Resultados en: ", SALIDA)
