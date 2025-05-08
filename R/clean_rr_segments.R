#' Limpieza y segmentación de intervalos RR de registros ECG
#'
#' La función `clean_rr_segments()` permite procesar la señal de intervalos RR que se ha obtenido a partir
#' de los registros ECG (electrocardiogramas) para su posterior análisis, eliminando artefactos, latidos ectópicos,
#' valores atípicos y tramos de señal de baja calidad
#'
#' @section Objetivo principal:
#' El objetivo principal de la función es reconstruir una señal continua y fisiológicamente plausible por medio de los siguientes pasos:
#' \itemize{
#'   \item \strong{Exclusión de anotaciones no latido y ectópicas.}
#'   \item \strong{Identificación de latidos anómalos} (outliers) mediante estadística local.
#'   \item \strong{Detección de saltos bruscos} en la señal RR que puedan alterar el análisis.
#'   \item \strong{Evaluación de la calidad de la señal por épocas}, marcando segmentos inestables.
#'   \item \strong{Segmentación de tramos válidos} y eliminación de huecos prolongados (≥ 5 s).
#'   \item \strong{Interpolación de huecos} mediante vecinos válidos, usando interpolación cúbica tipo PCHIP.
#' }
#'
#'
#' @param rr Vector numérico. Son los intervalos RR en muestras (tiempo entre latidos consecutivos).
#' @param rt Vector numérico. Tiempos acumulados en muestras desde el inicio del registro.
#' @param ann Vector numérico. Códigos ASCII de tipo de latido (ej. 4='N')
#' @param p Valor numérico. Parámetro de sensibilidad para la detección de outliers.
#' @param w Entero. Tamaño de la ventana local para la detección de outliers (número de latidos).
#' @param fs Entero. Frecuencia de muestreo de los electocardiogramas. 1000 Hz.
#' @param figsavefullpath Ruta opcional para guardar las figuras del diagnóstico.
#'
#' @return La función devuelve una lista con 3 elementos:
#'  \describe{
#'   \item{rr_new_all}{Vector numérico. Señal RR procesada, limpia y con interpolaciones aplicadas.}
#'   \item{rt_new_all}{Vector numérico. Tiempos correspondientes a cada intervalo RR en la señal limpia.}
#'   \item{seg_new}{Matriz numérica. Cada fila representa el índice de inicio y fin de un segmento válido.}
#' }
#' @section Reglas de limpieza aplicadas:
#' \itemize{
#'   \item Se excluyen anotaciones no latido (ej. '+', '|', 's').
#'   \item Se eliminan latidos ectópicos (códigos distintos de 78) y valores fuera de rango fisiológico.
#'   \item Se detectan outliers por ventana deslizante con umbral estadístico.
#'   \item Se invalidan épocas completas si hay menos de 15 latidos válidos, más de 10\% de ruido,
#'         o menos del 90\% del tiempo útil.
#'   \item Se segmenta la señal en tramos continuos válidos, excluyendo huecos ≥ 5 segundos.
#'   \item Se interpolan huecos breves usando interpolación cúbica PCHIP (`signal::interp1()`).
#' }
#'
#' @section Requisitos:
#' Necesita el paquete \pkg{signal} para la función \code{interp1()}.
#'
#' @seealso \code{\link{read_ann_file}} para la lectura inicial de archivos .ANN,
#' y \code{\link{signal::interp1}} para interpolación de tipo pchip.
#'
#' @examples
#' \dontrun{
#' # Simulación de una señal con latidos válidos y artefactos
#' rr <- c(800, 805, 810, NA, NA, 780, 1600, 800, 790, 810)
#' rt <- cumsum(rr, na.rm = TRUE)
#' rt <- cumsum(ifelse(is.na(rr), 0, rr))
#' ann <- rep(78, length(rr))  # todos normales
#' ann[7] <- 20  # un latido ectópico
#'
#' resultado <- clean_rr_segments(rr, rt, ann, p = 2, w = 5, fs = 1000)
#' plot(resultado$rt_new_all, resultado$rr_new_all, type = "l")
#' }
#'
#' @importFrom signal interp1
#' @export
clean_rr_segments <- function(rr, rt, ann, p, w, fs, figsavefullpath = NULL) {
  numd <- 10        # Número mínimo de datos válidos alrededor de un hueco para poder interpolar
  least_min <- 4    # Duración mínima calculada en minutos de latidos normales requeridos

  # Función auxiliar: reemplaza valores de un vector x entre índices i[1] y i[2] con nuevos valores x1
  replace <- function(x, i, x1) {
    if (i[2] >= i[1] && i[1] <= length(x) && i[2] <= length(x)) {
      before <- if (i[1] > 1) x[1:(i[1] - 1)] else c()
      after  <- if (i[2] < length(x)) x[(i[2] + 1):length(x)] else c()
      return(c(before, x1, after))
    } else {
      return(x)  # si los índices son inválidos, no se reemplaza nada
    }
  }

  # Función auxiliar: detecta intervalos contiguos donde rr es NA
  get_nan_intervals <- function(rr) {
    isnan_vec <- as.integer(is.na(rr))
    ss <- which(diff(c(0, isnan_vec)) == 1)  # inicio de huecos
    ee <- which(diff(c(isnan_vec, 0)) == -1) # fin de huecos
    return(list(start = ss, end = ee))
  }

  # Función para detectar outliers estadísticos en la señal RR, ignorando ectópicos
  deoutlier <- function(a, idx_ecto, w, p) {
    wl <- floor(w / 2)  # mitad del tamaño de la ventana deslizante para que sea simétrica
    ser <- 1:length(a) # índices originales

    # Si hay latidos ectópicos (idx_ecto) los vamos a eliminar temporalmente del análisis
    if (length(idx_ecto) > 0) {
      a <- a[-idx_ecto]          # quitamos los valores de la serie
      ser <- ser[-idx_ecto]      # quitamos los índices para mantender la correspondencia
    }
    # muy poca señal para evaluar con menos de 10 puntos
    if (length(a) < 10) return(numeric(0))

    # Se añaden valores al principio y al final para evitar errores en los bordes al aplicar la ventana
    a <- c(a[1:wl], a, a[(length(a) - wl + 1):length(a)])

    # Inicio de un vector para marcar qué posiciones son outliers (0=normal, 1=outlier)
    bad <- rep(0, length(a))

    # Para cada punto, calcula si se desvía mucho de la mediana local
    for (i in (wl + 1):(length(a) - wl)) {
      M <- median(a[(i - wl):(i + wl)], na.rm = TRUE) # mediana local
      ST <- sd(a[(i - wl):(i + wl)], na.rm = TRUE) # desviación estándar local

      # Si el valor se aleja más de p*ST de la mediana se marcará como outlier
      if (abs(a[i] - M) > ST * p) {
        bad[i] <- 1
      }
    }

    # Elimina los bordes añadidos artificialmente y recuperamos los índices originales
    bad <- bad[(wl + 1):(length(bad) - wl)]
    badidx <- ser[bad == 1] # devuelve los índices de los outliers detectados
    return(badidx)
  }

  # Función usada para interpolar tramos RR que tienen valores faltantes (NA)
  # Usamos la interpolación cúbica tipo "pchip" basada en valores vecinos válidos
  intp_rr <- function(rr, rt, numd) {
    # Eliminar NAs al inicio/final de la señal
    while (is.na(rr[1])) { rr <- rr[-1]; rt <- rt[-1] }
    while (is.na(rr[length(rr)])) { rr <- rr[-length(rr)]; rt <- rt[-length(rt)] }

    # Detectar huecos (intervalos de NA)
    nan_intervals <- get_nan_intervals(rr) # Usamos la función anterior para devolver una lista con los índices de inicio y fin de cada hueco
    ss <- nan_intervals$start # posición de inicio de hueco
    ee <- nan_intervals$end # posición de fin de hueco

    # Une una matriz de huecos (inicio y fin de cada tramo NA)
    irr <- if (length(ss) == 0) matrix(numeric(0), ncol = 2) else cbind(ss, ee)

    idx <- 1:length(rr) # índice global
    rslt <- rr  # copia de rr que será modificada con la interpolación
    rt_new <- rt # copia de rt ajustada
    r <- ceiling(numd / 2) # tamaño inicial de margen

    # Si hay huecos que necesitan ser rellenados
    if (nrow(irr) > 0) {
      # Se recorre cada intervalo con NAs para interpolar
      for (ii in nrow(irr):1) {
        # Selecciona una ventana alrededor del hueco (r vecinos a izquierda y derecha)
        iintv <- intersect(idx, (irr[ii, 1] - r):(irr[ii, 2] + r))
        intv <- rr[iintv] # extraemos los valores dentro de esta ventana

        # Asegura que haya suficientes valores válidos
        while (sum(!is.na(intv)) < numd) {
          r <- r + 1
          iintv <- intersect(idx, (irr[ii, 1] - r):(irr[ii, 2] + r))
          intv <- rr[iintv]
        }

        # Estimar duración del hueco e interpolar
        mu <- mean(intv, na.rm = TRUE) # media de RR local
        intvdur <- rt[irr[ii, 2] + 1] - rt[irr[ii, 1] - 1] # duración del hueco en tiempo
        numi <- round(intvdur / mu) - 1 # número estimado de RR que caben en ese hueco

        if (numi > 0) {
          # Realizamos la interpolación
          x <- iintv[!is.na(intv)] # índices de puntos válidos

          # Corrige desplazamiento de índices si la interpolación cambia longitud
          x[x > irr[ii, 2]] <- x[x > irr[ii, 2]] + numi - (irr[ii, 2] - irr[ii, 1] + 1)

          y <- intv[!is.na(intv)] # valores de RR válidos

          # Nuevos índices donde se insertarán los valores interpolados
          xi <- irr[ii, 1]:(irr[ii, 1] + numi - 1)

          # Interpolación PCHIP con signal::interp1 (cúbica, forma suave)
          yi <- signal::interp1(x, y, xi, method = "pchip")

          # Reemplazo de los NA en rr por los valores interpolados
          rslt <- replace(rslt, c(irr[ii, 1], irr[ii, 2]), yi)

          # Crea nuevos tiempos correspondientes en rt, espaciados uniformemente
          rt_new <- replace(rt_new, c(irr[ii, 1], irr[ii, 2]), rt[irr[ii, 1] - 1] + (1:numi) * (intvdur / (numi + 1)))
        } else {
          rslt <- replace(rslt, c(irr[ii, 1], irr[ii, 2]), numeric(0))
          rt_new <- replace(rt_new, c(irr[ii, 1], irr[ii, 2]), numeric(0))
        }
      }
    }
    # Devuelve la señal RR interpolada y el nuevo vector de tiempos
    return(list(rslt = rslt, rt_new = rt_new))
  }

  # Primer paso: Eliminar anotaciones no deseadas (no son latidos: '+', 's', '|')
  # ASCII: '+' = 43, 's' = 115, '|' = 124
  nonbeats <- which(ann %in% c(43, 115, 124))
  if (length(nonbeats) > 0) {
    ann <- ann[-nonbeats]
    rt  <- rt[-nonbeats]
    rr  <- rr[-nonbeats]
  }

  # Segundo paso: Recalcular RR por seguridad (diferencias entre tiempos)
  rr <- c(rr[1], diff(rt))

  # Tercer paso: Comprobar si hay suficiente señal de calidad (>= 4 min de latidos normales 'N')
  if (sum(rr[ann == 78], na.rm = TRUE) < 0.8 * least_min * 60 * fs) {
    return(list(rr_new_all = numeric(0), rt_new_all = numeric(0), seg_new = numeric(0)))
  }

  # Cuarto paso: Detectar ectópicos por tipo de latido o RR muy fuera de rango fisiológico
  idx1 <- which(ann != 78)  # 78 = 'N'
  if (length(idx1) > 0 && idx1[length(idx1)] == length(ann)) idx1 <- idx1[-length(idx1)]
  idx2 <- idx1 + 1
  idx_ecto <- union(union(idx1, idx2),
                    union(which(rr < fs * 60 / 220), which(rr > fs * 60 / 30)))

  # Quinto paso: Detectar latidos estadísticamente aberrantes (outliers)
  badidx <- deoutlier(rr, idx_ecto, w, p)

  # Sexto paso: Detectar cambios abruptos entre latidos consecutivos (saltos fuertes)
  para_diffRR <- 500
  badidx2 <- which(abs(diff(rr)) > para_diffRR) + 1

  # Ajustar índice del latido conflictivo anterior
  for (iib in seq_along(badidx2)) {
    thisbeat <- badidx2[iib]
    if (thisbeat <= 2) next
    if ((rr[thisbeat - 1] > rr[thisbeat] && rr[thisbeat - 1] > rr[thisbeat - 2]) ||
        (rr[thisbeat - 1] < rr[thisbeat] && rr[thisbeat - 1] < rr[thisbeat - 2])) {
      badidx2[iib] <- badidx2[iib] - 1
    }
  }

  badidx_total <- union(badidx, badidx2)

  # Séptimo paso: Descarta épocas enteras si son demasiado ruidosas o inestables
  prmt_epoch_length <- 30 # duraciónd e cada época en segundos
  prmt_num_N_beat <- 15 # número mínimo de latidos válidos esperados por épocas
  prmt_prop_bad_beat <- 0.1 # máximo 10 latidos con valores faltantes o malos
  prmt_prop_normal_time <- 0.9 # tiempo útil

  # Copiamos rr y marcamos como NA los latidos ectópicos y los detectados como erróneos
  rr_tmp <- rr
  rr_tmp[idx_ecto] <- NA
  rr_tmp[badidx_total] <- NA

  # Vector de igual longitud que rr: marcará con 1 los índices de épocas malas
  idx_badepoch <- rep(0, length(rr))

  # Número total de épocas en el registro, basado en tiempo total y duración de época
  Nepoch <- ceiling(rt[length(rt)] / fs / prmt_epoch_length)

  # Recorremos cada época
  for (iepoch in 1:Nepoch) {
    # Selección de los índices de latidos cuya marca de tiempo caen sobre esta época
    idx <- which(rt >= (iepoch - 1) * fs * prmt_epoch_length & rt < iepoch * fs * prmt_epoch_length)
    # cortamos los rr de esa época y evaluamos su calidad
    rr_cut <- rr_tmp[idx]

    # Condiciones para marcar una época mala
    if (length(rr_cut) < prmt_num_N_beat ||
        sum(is.na(rr_cut)) / length(rr_cut) > prmt_prop_bad_beat ||
        sum(rr_cut, na.rm = TRUE) / fs / prmt_epoch_length < prmt_prop_normal_time) {
      # Marcamos todos esos latidos como pertenecientes a una época no confiable
      idx_badepoch[idx] <- 1
    }
  }

  # Paso 8: Reasignación de NA en rr para limpiar la señal según lo detectado

  rr[badidx_total] <- NA  # quita valores etiquetados como outliers
  rr[idx_ecto] <- NA  # quita ectópicos
  rr[idx_badepoch == 1] <- NA # quita épocas enteras marcadas como malas

  # Creamos vector lógico para encontrar inicios y finales de bloques con NA
  isnan <- is.na(rr)

  # Detectamos los bordes de cada bloque NA
  ss <- which(diff(c(FALSE, isnan)) == 1)  # inicios de NA
  ee <- which(diff(c(isnan, FALSE)) == -1) # finales de NA

  # Si no hay huecos: señal continua. Si los hay: matriz con inicios y finales de huecos
  segNaN <- if (length(ss) == 0 || length(ee) == 0) {
    matrix(numeric(0), ncol = 2)
  } else {
    cbind(ss, ee)
  }

  # Ahora identificamos las partes entre huecos que tienen duración suficiente
  prmt_seg_def_sec <- 5  # Se requiere que el hueco dure al menos 5 s para generar nuevos segmentos

  if (nrow(segNaN) > 0) {
    intervs <- rt[ee] - rt[ss]  # duración de cada hueco

    # Identificamos huecos suficientemente largos para separar segmentos
    gap_idx <- which(intervs >= prmt_seg_def_sec * fs)

    if (length(gap_idx) > 0) {
      # Construimos inicio y fin de cada segmento limpio entre huecos largos
      seg_start <- c(1, ee[gap_idx] + 1)
      seg_end <- c(ss[gap_idx] - 1, length(rr))
      seg <- cbind(seg_start, seg_end)

      # Eliminamos segmentos vacíos (mal definidos)
      seg <- seg[seg[, 2] > seg[, 1], , drop = FALSE]
    } else {
      seg <- matrix(c(1, length(rr)), nrow = 1)  # si no hay huecos largos, solo 1 segmento
    }
  } else {
    seg <- matrix(c(1, length(rr)), nrow = 1)    # sin huecos: todo el RR es 1 solo segmento
  }


  # Paso 9: Interpolación por segmentos válidos detectados

  rr_new_all <- numeric(0)  # Vector final concatenado de RR interpolados
  rt_new_all <- numeric(0)  # Vector de tiempos correspondiente
  seg_new <- matrix(1, nrow = nrow(seg), ncol = 2)  # Nueva tabla de segmentos
  curr_seg_st <- 1 # posición inicial del siguiente segmento interpolado

  # Recorrer cada uno de los segmentos definidos como válidos
  for (iseg in 1:nrow(seg)) {
    # Extraemos el tramo rr y rt del segmento
    rr_seg <- rr[seg[iseg, 1]:seg[iseg, 2]]
    rt_seg <- rt[seg[iseg, 1]:seg[iseg, 2]]

    # Interpolamos huecos internos dentro del segmento
    intp_result <- intp_rr(rr_seg, rt_seg, numd)
    rr_new <- intp_result$rslt
    rt_new <- intp_result$rt_new

    # Concatenamos los resultados al vector final
    rr_new_all <- c(rr_new_all, rr_new)
    rt_new_all <- c(rt_new_all, rt_new)

    # Guardamos el índice de inicio y fin de este segmento en la señal reconstruida
    seg_new[iseg, ] <- c(curr_seg_st, curr_seg_st + length(rr_new) - 1)

    # Actualizamos el contador de posición
    curr_seg_st <- curr_seg_st + length(rr_new)
  }


  # Resultado final
  return(list(rr_new_all = rr_new_all, rt_new_all = rt_new_all, seg_new = seg_new))
}
