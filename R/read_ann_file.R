#' Leer y decodificar archivos .ANN de ECG
#'
#' La función `read_ann_file()` nos permite leer archivos binarios `.ANN` que contienen las anotaciones limpias de latidos cardíacos procedentes de señales de electocardiogramas o ECG.
#' Los archivos presentan la siguiente estructura: comienzan con 12 bytes de cabecera, terminan con 14 bytes de cierre, y contienen bloques intermedios de 8 bytes por cada latido.
#'
#' Cada bloque de 8 bytes codifica:
#' \itemize{
#'   \item El intervalo RR, es decir, tiempo entre latidos consecutivos, en muestras (típicamente a 1000 Hz).
#'   \item El tipo de latido, codificado en formato ASCII.
#' }
#'
#' La función realiza las siguientes operaciones:
#' \enumerate{
#'   \item Lee el archivo binario completo.
#'   \item Elimina los bytes de cabecera y pie.
#'   \item Extrae los bloques de 8 bytes por latido.
#'   \item Calcula los intervalos RR acumulados (RR_time).
#'   \item Interpreta los códigos de tipo de latido (ANNtype) según su valor ASCII.
#' }
#'
#' Algunos ejemplos de tipo de latido son:
#' \tabular{llll}{
#'   \strong{Código ASCII} \tab \strong{Símbolo} \tab \strong{Significado clínico} \tab \strong{Tratamiento} \cr
#'   4 \tab N \tab Latido normal del seno \tab Usar para HRV \cr
#'   20 \tab V \tab PVC (latido ventricular prematuro) \tab Outlier, excluir \cr
#'   36 \tab S \tab Latido supraventricular ectópico \tab Excluir \cr
#'   40 \tab E \tab Escape ventricular \tab Incluir si relevante \cr
#'   88 \tab * \tab Artefacto/anotación especial \tab Ignorar \cr
#'   (y otros)
#' }
#'
#' @param filepath Ruta al directorio donde se encuentra el archivo `.ANN`.
#' @param filename Nombre del archivo `.ANN` (por ejemplo, `"Thoth_20230208220000s.ANN"`).
#'
#' @return Una lista con:
#' \describe{
#'   \item{RR}{Vector numérico con los intervalos RR en milisegundos.}
#'   \item{RR_time}{Vector acumulado del tiempo total transcurrido desde el inicio del registro.}
#'   \item{ANNtype}{Vector de caracteres con el tipo de latido interpretado (por ejemplo, 'N', 'V', '*').}
#' }
#'
#' @examples
#' \dontrun{
#' # Suponiendo que tienes un archivo ANN en "datos/Thoth_20230208220000s.ANN"
#' resultado <- read_ann_file(filepath = "datos", filename = "Thoth_20230208220000s.ANN")
#'
#' # Ver los primeros 10 intervalos RR
#' head(resultado$RR, 10)
#'
#' # Ver los tipos de los primeros latidos
#' head(resultado$ANNtype, 10)
#' }
#'
#' @export

read_ann_file <- function(filepath, filename) {
  full_path <- file.path(filepath, filename)

  # Primer paso: leer archivo como binario
  con <- file(full_path, "rb")
  bytes <- readBin(con, what = "raw", n = file.info(full_path)$size)
  close(con)

  # Segundo paso: saltar cabecera y pie
  bytes <- bytes[13:(length(bytes) - 14)]

  n_beats <- floor(length(bytes) / 8)
  byte_matrix <- matrix(as.integer(bytes[1:(n_beats * 8)]), nrow = 8)

  # Tercer paso: Verificar campos inesperados
  if (any(byte_matrix[6, ] > 0) || any(byte_matrix[8, ] > 0)) {
    stop("Formato de archivo ANN inesperado.")
  }

  # Cuarto paso: RR en muestras
  RR <- byte_matrix[4, ] +
    byte_matrix[5, ] * 256 +
    byte_matrix[2, ] * 256^2 +
    byte_matrix[3, ] * 256^3

  # Quinto paso: suma acumulada
  RR_time <- cumsum(RR)  # Tiempo acumulado en muestras

  # Sexto paso: Mapeamos del tipo de latido
  ann_code <- byte_matrix[7, ]
  ANNtype <- rep("X", n_beats)

  ANNtype[ann_code == 112] <- "+"
  ANNtype[ann_code == 120] <- "?"
  ANNtype[ann_code ==  40] <- "E"
  ANNtype[ann_code ==   4] <- "N"
  ANNtype[ann_code ==  36] <- "S"
  ANNtype[ann_code ==  52] <- "Q"
  ANNtype[ann_code ==  20] <- "V"
  ANNtype[ann_code == 164] <- "r"
  ANNtype[ann_code ==  72] <- "s"
  ANNtype[ann_code ==  32] <- "A"
  ANNtype[ann_code ==  88] <- "*"

  return(list(RR = RR, RR_time = RR_time, ANNtype = ANNtype))
}
