# Comportamiento anómalo de movimientos en masa en el Valle de Aburrá
El trabajo busca los movimientos en masa cuyo comportamiento (área, movilidad H/L y elongación) se aparta de lo esperado según sus condicionantes. Está dividido en siete fases, cada una en su carpeta. Cada código empieza con una explicación de su objetivo, sus datos de entrada, sus salidas y cómo correrlo.
 
## Requisitos
 
- **ArcGIS Pro 3.7**, con su entorno de Python (`arcpy`, `numpy`, `pandas`, `scipy`, `matplotlib`), para los notebooks (`.ipynb`). Algunos pasos usan la extensión Spatial Analyst.
- **R 4.6.1** y RStudio para los scripts (`.R`). Los paquetes principales son `spatstat` y `gstat` (Fase 2), `GWmodel`, `robustbase` y `spdep` (Fases 3 y 6), `ClustGeo` (Fase 4) y `vegan` y `coin` (Fase 5). La mayoría de los scripts instala los que falten.

## Orden del trabajo
 
Cada fase usa las salidas de las anteriores. En general, la preparación espacial se hace en ArcGIS Pro (notebooks) y la estadística en R.
 
**Fase 0. Completitud del inventario.** Ajusta la relación entre el área y la frecuencia de los eventos y calcula el umbral de completitud de cada fuente (mapeo propio e institucional). Los eventos con área menor que ese umbral quedan marcados en el inventario; esa marca acompaña a cada evento hasta la Fase 6.
 
**Fase 1. Métricas de movilidad.** Mide en el DEM la altura (H) y la longitud (L) de cada evento, entre la corona y el pie, y calcula el H/L.
 
**Fase 2. Patrón espacial.** Define la corona de cada evento, el punto que lo representa en el resto del trabajo, y estudia cómo se distribuyen las coronas y si los eventos vecinos se parecen. De aquí sale la escala de 1 km que usa la Fase 4. La carpeta `exploratorio` guarda versiones anteriores que no hacen falta para reproducir los resultados.
 
**Fase 3. Comportamiento esperado y anomalías.** Revisa la muestra, estima con GWR mixta lo que se espera de cada evento según sus condicionantes y marca como anomalías contextuales los eventos cuyos residuales se apartan más (distancia robusta con MCD). Al final integra todo en `F3_fase3_integrado` (CSV y GeoPackage), la entrada de las fases 4 a 6.
 
**Fase 4. Regionalización.** Arma en ArcGIS Pro el grafo de vecindad de las coronas y, en R, divide el valle en regiones con SKATER, mide su estabilidad con bootstrap, las compara con ClustGeo y las caracteriza.
 
**Fase 5. Asociación con los condicionantes.** Un notebook une en `F5_eventos.csv` los resultados anteriores, las unidades sin agrupar y la orientación de cada evento frente a la falla más cercana. Los scripts de R prueban si el comportamiento cambia entre unidades y si las anomalías se asocian con algún condicionante, con la distancia a las fallas o con su orientación.
 
**Fase 6. Validación.** Un notebook arma una planilla y un archivo KML para revisar las 19 anomalías en imágenes. Después de la revisión, el script de R mide la concentración de las anomalías, su sensibilidad a las decisiones del método y su estabilidad con remuestreo, y arma la Tabla 4.
 
## Procesos que no están en los códigos publicados
 
Algunos pasos se hicieron con herramientas de ArcGIS Pro, a mano o con códigos que no están en este repositorio:
 
- **Fase 1.** El recorte del DEM a los polígonos (Extraer por máscara), la elongación con el rectángulo mínimo (Geometría mínima de delimitación) y los logaritmos de las métricas (Calcular campo). El notebook describe cada paso.
- **Fase 2.** El ráster de pendiente del valle, las ventanas de laderas con pendiente de al menos 10° y 4,7°, y el mapa final de densidad kernel, hechos con herramientas de ArcGIS Pro.
- **Fase 3.** La tabla del modelo (`F3_S2_tabla_modelo.csv`) se armó con un script aparte. Esa tabla:
  - une las métricas de la Fase 1 con los condicionantes de cada corona, a saber, la clase mayoritaria de geología, geomorfología y cobertura en un círculo de 18 m (con las clases ya agrupadas), la pendiente, la distancia a fallas y a drenajes y la densidad de lineamientos en 1 km;
  - marca con `EN_MODELO = 1` los 427 eventos del modelo. Quedan fuera, entre otros, las caídas y los eventos sin tipo, porque son muy pocos para estimar su efecto.
  Los mapas de los coeficientes locales también se hicieron con un script aparte, a partir de `F3_S3_resultados.csv`.
- **Fase 4.** La capa del río Medellín se copió a `fase4.gdb` para dibujarla en los mapas.
- **Fase 6.** La tabla de atributos del inventario se exportó a `tabla_inventario.xls`. Cada anomalía se revisó a mano en Google Earth Pro con la planilla y el KML, respondiendo seis preguntas sobre el polígono y la imagen. Las lecturas se guardaron en `F6_lecturas.csv`, con las columnas `ID`, `lectura` y `fecha_imagen`.
## Datos
 
Las capas de entrada (inventario, DEM y condicionantes) no están en este repositorio. Cada código indica en su encabezado las capas y tablas que usa.
