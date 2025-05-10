#' Función para el análisis espectral de la variabilidad de la frecuencia cardíaca en el dominio de la frecuencia
#' @description
#' La función `compute_hrv_frecuence()` realiza un análisis espectral de la variabilidad de la frecuencia cardíaca (HRV)
#' usando la transformada rápida de Fourier (FFT) para calcular la potencia en diferentes bandas de frecuencia estandarizadas
#' asociadas a componentes del sistema nervioso autónomo.
#'
#' @section Objetivo principal:
#' \itemize{
#'   \item Obtener potencias espectrales en bandas ULF, VLF, LF y HF.
#'   \item Calcular el cociente LF/HF como estimador del balance simpático-parasimpático.
#'}
#'
#'
#' @section Detalles:
#' \enumerate{
#'   \item Interpolación de la señal RR sobre una malla equiespaciada mediante \code{interp1()} tipo PCHIP.
#'   \item Eliminación de la media para centrar la señal (variabilidad pura).
#'   \item Cálculo de autocorrelación para suavizar y estabilizar el espectro.
#'   \item Aplicación de ventana de Hann para mitigar efectos de borde.
#'   \item Cálculo de la FFT y conversión a espectro de potencia.
#'   \item Integración espectral por bandas estándar: ULF, VLF, LF, HF.
#' }
#' @param t1 Vector numérico. Tiempos acumulados (en milisegundos) de los latidos válidos.
#' @param b1 Vector numérico. Intervalos RR (en milisegundos) correspondientes a los tiempos `t1`.
#' @param resample_index resample_index Intervalo de remuestreo deseado (en milisegundos). Por defecto, 125 ms (≈ 8 Hz).
#' @param fs Frecuencia de muestreo original del ECG (en Hz). Por defecto son 1000 Hz.
#'
#' @return Lista con:
#' \describe{
#'   \item{P_ULF}{Potencia en la banda ultra baja frecuencia (<0.003 Hz).}
#'   \item{P_VLF}{Potencia en la banda muy baja frecuencia (0.003–0.04 Hz).}
#'   \item{P_LF}{Potencia en la banda baja frecuencia (0.04–0.15 Hz).}
#'   \item{P_HF}{Potencia en la banda alta frecuencia (0.15–0.4 Hz).}
#'   \item{P_TP}{Potencia total (0–0.4 Hz).}
#'   \item{P_TPP}{Potencia total de la señal autocorrelada.}
#'   \item{LF_HF_ratio}{Cociente entre bandas LF y HF.}
#'   \item{fdatx}{Magnitudes del espectro de frecuencia (FFT).}
#'   \item{freq}{Vector de frecuencias en Hz correspondiente a `fdatx`.}
#' }
#'
#' @section Requisitos:
#' Esta función requiere el paquete \pkg{signal} para las funciones \code{interp1()} y \code{conv()}.
#'
#' @seealso \code{\link[signal]{interp1}} para interpolación cúbica, y \code{\link{compute_hrv_time}} para análisis HRV en el dominio del tiempo.
#'
#' @examples
#' \dontrun{
#' # Simulación de señal RR:
#' t <- seq(0, 60000, by = 1000)  # Tiempos en ms (1 latido por segundo)
#' rr <- rep(1000, length(t))     # Intervalos RR constantes de 1000 ms
#'
#' result <- HRVf2_R(t, rr)
#' plot(result$freq, result$fdatx, type = "l", xlab = "Frecuencia (Hz)", ylab = "Amplitud")
#' }
#'
#' @importFrom signal interp1 conv
#' @export
compute_hrv_frecuence <- function(t1, b1, resample_index = 125, fs = 1000) {

  # Primer paso: Generar una malla de tiempo equiespaciada new_t_RR con paso de resample_index en ms
  new_t_RR <- seq(t1[1], t1[length(t1)], by = resample_index) # FFT necesita una señal muestreada a un tiempo uniforme

  # Interpolación cúbica tipo PCHIP de los intervalos RR sobre la nueva malla
  b2 <- interp1(t1, b1, new_t_RR, method = "pchip")

  # Longitud de la señal una vez remuestreada
  f_n <- length(b2)

  # Segundo paso: Se elimina la media para centrar la señal
  da_t <- b2 - mean(b2)

  # Tercer paso: Autocorrelación
  # Se calula la autocorrelación de la señal centrada, usando convolución con su reverso
  cdat <- signal::conv(da_t, rev(da_t))

  # Se normaliza
  cdat <- cdat / c(1:f_n, (f_n - 1):1)

  # Cuarto paso: Aplicar la ventana de Hann para reducir los artefactos por los bordes
  L <- round(f_n / 4) # longitud de la ventana

  # Eje simétrico de la ventana centrado en 0
  tn <- seq(-L + 1, L - 1) / (L * 2)

  # Ventana de Hann con forma coseno^2
  han_w <- cos(pi * tn)^2

  # Se aplica la ventana en la posición central del vector autocorrelado
  han_w2 <- rep(0, length(cdat))
  han_start <- (length(cdat) - length(han_w)) %/% 2 + 1
  han_w2[han_start:(han_start + length(han_w) - 1)] <- han_w
  cdat <- cdat * han_w2

  # Quinto paso: Transformada rápida de Fourier o FFT y normalización
  Resample_cdat <- c(cdat, 0)
  fdatx <- abs(fft(Resample_cdat)) / f_n

  # Sexto paso: Definición de las bandas espectrales
  max_frequency <- fs / resample_index # Frecuencia de muestreo de la señal RR remuestreada en Hz

  ULFband <- 1:(ceiling(0.003 * (f_n * 2) / max_frequency) - 1)
  VLFband <- ceiling(0.003 * (f_n * 2) / max_frequency):(ceiling(0.04 * (f_n * 2) / max_frequency) - 1)
  LFband  <- ceiling(0.04 * (f_n * 2) / max_frequency):(ceiling(0.15 * (f_n * 2) / max_frequency) - 1)
  HFband  <- ceiling(0.15 * (f_n * 2) / max_frequency):(ceiling(0.4 * (f_n * 2) / max_frequency) - 1)
  TPband  <- 1:(ceiling(0.4 * (f_n * 2) / max_frequency) - 1)

  # Séptimo paso: Cálculo de las potencias de cada banda
  # Se calcula la "potencia" espectral como la suma de las amplitudes en cada banda
  P_TPP <- sum(fdatx)
  P_TP  <- sum(fdatx[TPband])
  P_ULF <- sum(fdatx[ULFband])
  P_VLF <- sum(fdatx[VLFband])
  P_LF  <- sum(fdatx[LFband])
  P_HF  <- sum(fdatx[HFband])
  LHratio <- P_LF / P_HF

  # Octavo paso: Preparación de espectro y frecuencia para visualización
  # Se selecciona el rango de frecuencias hasta 0.5 Hz
  displayband <- 1:(ceiling(0.5 * (f_n * 2) / max_frequency) - 1)
  displaybandf <- displayband * max_frequency / (f_n * 2)

  # Noveno paso: Devolver los resultados como una lista
  return(list(
    P_ULF = P_ULF,
    P_VLF = P_VLF,
    P_LF = P_LF,
    P_HF = P_HF,
    P_TP = P_TP,
    P_TPP = P_TPP,
    LF_HF_ratio = LHratio,
    fdatx = fdatx[displayband],
    freq = displaybandf
  ))
}
