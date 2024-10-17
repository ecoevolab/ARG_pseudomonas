# script para filtrar archivo de snvs con 33 genomas 
setwd("/Users/itzyyp/Documents/lab_sur/ARG/metadata/snvs/33genomas_bacter/")
x = read.table("snps_combined.catalog.tsv", head = TRUE)

library(dplyr)
#realizar columna de alternativo
x <- x %>% 
  mutate(alternativo = substr(snp_id, nchar(snp_id), nchar(snp_id)))
## borrar del nombre de fila 
x <- x %>% 
  mutate(snp_id = substr(snp_id, 1, nchar(snp_id) - 3))
## realizar columna de referencia
x <- x %>% 
  mutate(referencia = substr(snp_id, nchar(snp_id), nchar(snp_id)))
## borrar
x <- x %>% 
  mutate(snp_id = substr(snp_id, 1, nchar(snp_id) - 3))
## posicion
library(stringr)
#borra todo lo que esta depsues de ||, el ?>= es para decir despues de esos sin tomarlos en cuenta 
x <- x %>% 
  mutate(posicion = str_extract(snp_id, "(?<=\\|\\|).*"))
## borrar lo delimitadores 
x <- x %>% 
  mutate(snp_id = str_remove(snp_id, "(\\|\\|).*"))
# columna de contig
x <- x %>%
  mutate(contig = str_extract(snp_id, "(?<=\\|)[^|]*(?=\\|)")) # Extrae lo que está entre el primer y segundo '|'
## borrar lo delimitadores 
x <- x %>%
  mutate(contig = str_extract(snp_id, "(?<=\\|)[^|]*(?=\\|)")) %>% # Extrae lo que está entre el primer y segundo '|'
  select(-snp_id)
#ordenar columnas, tiene 86650 filas
x <- x[, c(37,1:36)]

#filtrar para eliminar filas donde todos tiene alternativos
i <- 1
while (i <= nrow(x)) {
  suma <- sum(x[i, 2:34])
  if (suma == 33) {
    x <- x[-i, ]
    # No incrementamos i aquí para que en la próxima iteración i apunte a la siguiente fila
  } else {
    # Solo incrementamos i si no eliminamos una fila
    i <- i + 1
  }
}

column_names <- colnames(x)
num_cols <- ncol(x)
selected_columns <- x[, 2:(num_cols-3)]
selected_names <- colnames(selected_columns)
rm(selected_columns)
names = substr(selected_names, 6,12)

#anotar para cada genoma los alternativos con los cuales cuenta
num = 2 
for (nombre in names) {  
  assign(nombre, x[,c(1,num,35,37)])
  num = num+1
}

#eliminar los 0 y quedarme con los que tengan 1
for (nombre in names) {
  assign(nombre, get(nombre)[get(nombre)[, 2] != 0, ])
}

#cambiar nombre de columna para cuanod lo juntemos con demás, poder identificar
#a qn pertenece esa columna
for (nombre in names) {
  name = get(nombre)
  colnames(name)[colnames(name) == "alternativo"] <- paste("alternativo", nombre, sep = "_")
  assign(nombre, name)
}

for (nombre in names) {
  # Obtener el dataframe actual
  data = get(nombre)
  # Ordenar el dataframe por la columna 'contig'
  data_sorted = data[order(as.numeric(gsub("\\D", "", data$contig))), ]
  # Reasignar el dataframe ordenado a la variable con el mismo nombre
  assign(nombre, data_sorted)
}

###
nombres= c()
for (nombre in names) {
  b = 1
  unique_nom = unique(get(nombre)$contig)
  u = 1
  while (u <= nrow(get(nombre))){
    assign(paste(nombre, unique_nom[b], sep = "_"), filter(get(nombre), contig == unique_nom[b]))
    nombres = c(nombres, paste(nombre, unique_nom[b], sep = "_"))
    u = u + 1
    if (u <= nrow(get(nombre)) && get(nombre)$contig[u] != unique_nom[b]){
      b = b + 1
    }
  }
}

nombres = unique(nombres)
long = length(unique_nom)

# Inicialización de vectores dinámicamente
for (name in unique_nom) {
  assign(name, c())  # Crear un vector vacío con el nombre de cada contig
}

# hacer lista de nombres
for (contig in nombres) {
  for (w in 1:long) {
    if (substr(contig, start = 8, stop = 19) == unique_nom[w]) {
      name_contig <- unique_nom[w]
      current_vector <- get(name_contig)
      updated_vector <- c(current_vector, contig)
      assign(name_contig, updated_vector)
    }
  }
}

for (contig in unique_nom) {
  cont = get(contig)
  if (length(cont) > 0) {
    # Inicializar el primer vector para el merge
    merged_data <- get(cont[1])
    for (q in 2:length(cont)) {
      val = get(cont[q])
      merged_data <- merge(merged_data, val, all = TRUE)
    }
    # Reemplazar NA con N, pues bacter no acepta "-"
    merged_data[is.na(merged_data)] <- "-"
    # Asignar el resultado final del merge al nombre original
    assign(contig, merged_data)
  } else {
    # Si no hay vectores en 'cont', asignar un dataframe vacío
    assign(contig, data.frame(matrix(ncol = 0, nrow = 0)))
  }
}

# eliminar posiciones que no quiero
for (contig in unique_nom) {
  # Obtener el dataframe resultante del nombre de 'contig'
  merged_data <- get(contig)
  # Verificar si el dataframe no es NULL y tiene columnas
  if (!is.null(merged_data) && ncol(merged_data) > 0) {
    # Usar dplyr::select para eliminar columnas que empiezan con "snvs_" y la columna "contig"
    merged_data <- merged_data %>%
      select(-starts_with("snvs_"), -contig)
    
    # Asignar el dataframe limpio al nombre original
    assign(contig, merged_data)
  }
}

unique_nom <- setdiff(unique_nom, "UUNS01000009")
carpeta <- "/Users/itzyyp/Documents/lab_sur/ARG/metadata/snvs/secuencias/"
for (contig in unique_nom) {
  # Obtener el dataframe resultante del nombre de 'contig'
  merged_data <- get(contig)
  # Verificar si el dataframe no es NULL y tiene columnas
  if (!is.null(merged_data) && ncol(merged_data) > 0) {
    # Eliminar la columna de posiciones
    merged_data <- merged_data %>%
      select(-starts_with("posicion"))
    # Crear el contenido del archivo FASTA
    fasta_lines <- c()
    # Iterar sobre las columnas restantes
    for (col_name in colnames(merged_data)) {
      # Extraer los datos de la columna
      column_data <- merged_data[[col_name]]
      # Crear el contenido del header y los datos de la columna
      fasta_header <- paste0(">", col_name)
      fasta_content <- paste(column_data, collapse = "")
      # Añadir al contenido total
      fasta_lines <- c(fasta_lines, fasta_header, fasta_content)
    }
    # Nombre del archivo para este contig
    archivo <- paste0(carpeta, contig, ".fasta")
    writeLines(fasta_lines, archivo)
  }
}

