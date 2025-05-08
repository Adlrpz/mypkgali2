\name{datos_ANN}
\alias{datos_ANN}
\docType{data}
\title{Archivos .ANN incluidos en el paquete}
\description{
Este paquete contiene un conjunto de archivos reales de anotaciones de ECG (\code{.ANN}) procedentes del estudio ONTIME, en el que participaron 229 pacientes con registros ambulatorios de hasta 3,5 días, capturados con dispositivos Holter a 1000 Hz.

Los archivos han sido corregidos manualmente tras la detección automática de ondas R, eliminando artefactos y latidos ectópicos. Cada archivo corresponde a una hora de registro limpio y está organizado por paciente, por ejemplo:

\code{inst/extdata/ECG_A0191_PM/Thoth_20230208220000s.ANN}

Cada archivo se puede acceder desde el paquete utilizando la función \code{system.file()}.

La función \code{\link{read_ann_file}} permite decodificar estos archivos y extraer:
\itemize{
  \item Intervalos RR (en milisegundos)
  \item Tiempo acumulado (\code{RR_time})
  \item Tipo de latido (\code{ANNtype})
}

El formato de cada archivo .ANN es el siguiente:
\itemize{
  \item 12 bytes de cabecera (reservados)
  \item Bloques de 8 bytes por cada latido:
    \itemize{
      \item Bytes 2–5 codifican el intervalo RR (en muestras a 1000 Hz)
      \item Byte 7 codifica el tipo de latido en formato ASCII
      \item Bytes 6 y 8 deben ser cero
    }
  \item 14 bytes de cierre donde acabaría el archivo
}
}

\usage{
# Solo informativo, no se accede con data()
}

\format{
Archivos binarios externos, accedidos mediante rutas como:

\code{
  file <- system.file("extdata/ECG_A0191_PM", "Thoth_20230208190000.ANN", package = "mypkgali2")
}
}

\seealso{
\code{\link{read_ann_file}} para leer los archivos ANN.
}
