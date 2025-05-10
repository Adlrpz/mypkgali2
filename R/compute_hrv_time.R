#' @title Cálculo de métricas de variabilidad de la frecuencia cardíaca (HRV) en el dominio del tiempo
#' @description
#' La función `compute_hrv_time()` permite calcular las métricas estándar de la variabilidad de la frecuencia cardíaca a partir de una serie de intervalos RR ya limpios y validados.
#'
#'
#' @section Objetivo principal:
#' \itemize{
#'   \item Calcular los indicadores cásicos de HRV en registros de electrocardigramas preprocesados
#'    }
#'
#' @section Detalles:
#' Entre las métricas que se calculan en el análisis encontramos:
#' #' \enumerate{
#'   \item Cálculo de las diferencias sucesivas entre intervalos RR (`diff(rr_ms)`).
#'   \item Estimación de RMSSD (raíz del promedio de diferencias cuadradas).
#'   \item Proporciones pNN50 y pNN20 (porcentaje de diferencias mayores a 50 o 20 ms).
#'   \item Desviación estándar global de RR (`SDRR`, equivalente a SDNN).
#'   \item Componentes SD1 y SD2 del gráfico de dispersión de Poincaré.
#'   \item Cociente SDRatio = SD1 / SD2, como medida del balance autonómico.
#' }
#'
#' @param rr_ms Vector numérico con los intervalos RR en milisegundos o ms, ya filtrados y validados
#'
#' @return Una lista con las siguientes métricas redondeadas a 4 decimales:
#' \describe{
#'   \item{Mean}{Media de los intervalos RR en ms.}
#'   \item{SDRR}{Desviación estándar de los RR (ms). Representa variabilidad global.}
#'   \item{RMSSD}{Raíz cuadrada del promedio de las diferencias sucesivas (ms). Indicador de variabilidad rápida.}
#'   \item{pNN50}{Porcentaje de pares de RR con diferencia > 50 ms.}
#'   \item{pNN20}{Porcentaje de pares de RR con diferencia > 20 ms.}
#'   \item{SD1}{Componente transversal del gráfico de Poincaré. Variabilidad a corto plazo.}
#'   \item{SD2}{Componente longitudinal del gráfico de Poincaré. Variabilidad a largo plazo.}
#'   \item{SDRatio}{Relación SD1 / SD2. Representa el balance entre componentes vagales y simpáticas.}
#' }
#' @seealso \code{\link{clean_rr_segments}} para preparar la señal RR antes de este análisis.
#'
#' @examples
#' \dontrun{
#' # Simulación de una serie RR artificial en milisegundos
#' rr <- c(800, 810, 795, 805, 820, 790, 780, 810, 800)
#' result <- compute_hrv_time(rr)
#' print(result)
#' }
#'
#' @importFrom stats sd
#'
#' @export
compute_hrv_time <- function(rr_ms) {

  # Primer paso: Cálculo de las diferencias sucesivas entre intervalos RR
  rr_diff <- diff # cómo varía la duración entre latidos adyacentes

  # Segundo paso: Cálculo de SD1
  # Se calcula la desviación estándar de las diferencias
  SDSD <- sd(rr_diff, na.rm = TRUE)
  # Es una medida de las dispersión perpendicular a la línea de identidad en el gráfico de Pointcaré
  SD1 <- sqrt(0.5 * SDSD^2)

  # Tercer paso: Cálculo de SDRR
  # Es la desviación estándar de todos los intervalos RR
  SDRR <- sd(rr_ms, na.rm = TRUE) # variabilidad global de la señal

  # Cuarto paso: Cálculo de SD2
  SD2 <- sqrt(2 * SDRR^2 - 0.5 * SDSD^2) # componente longitudinal en el gráfico de Pointcaré

  # Quinto paso: Cociente entre SD1 y SD2
  SDRatio <- SD1 / SD2

  # Sexto paso: Cálculo Root Mean Square of Successive Differences
  # Raíz cuadrada de la media de las diferencias sucesivas elevadas al cuadrado
  RMSSD <- sqrt(mean(rr_diff^2, na.rm = TRUE))

  # Séptimo paso: Cálculo pNN50
  # Porcentaje de diferencias sucesivas mayores a 50 ms.
  pNN50 <- mean(abs(rr_diff) > 50, na.rm = TRUE) * 100

  # Octavo paso: pNN20
  # Porcentaje de las diferencias sucesivas mayores a 20 ms
  pNN20 <- mean(abs(rr_diff) > 20, na.rm = TRUE) * 100

  # Noveno paso: Media de los intervalos RR
  Mean <- mean(rr_ms, na.rm = TRUE)


  # Se devuelve una lista con los resultados redondedos a 4 decimales
  return(list(
    Mean = round(Mean, 4),
    SDRR = round(SDRR, 4),
    RMSSD = round(RMSSD, 4),
    pNN50 = round(pNN50, 4),
    pNN20 = round(pNN20, 4),
    SD1 = round(SD1, 4),
    SD2 = round(SD2, 4),
    SDRatio = round(SDRatio, 4)
  ))
}
