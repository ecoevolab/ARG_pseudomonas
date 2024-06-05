# script unir tabla de metadatos con archivo correspondiente
setwd("/Users/itzyyp/Documents/lab_sur/ARG/metadata/")
library(tidyverse)
metadata = read.csv("metadata_1524.csv")
plate = metadata
# eliminar todas las columnas, menos donde viene el id de la placa 
plate <- plate[ ,-c(2:7)]
plate <- gsub("plate", "", plate)
plate <- trimws(plate)
# eliminar todo lo que viene despues del punto
plate = stringr::str_extract(plate, pattern = "^(.+?)")
cepas = read.table("nombrescepas")
#quitar nombre de la carpeta para solamente tener el nombre de la cepa 
cepas$V1 = gsub("PRJEB24450/", "", cepas$V1)
#quitar segunda columna
cepas <- cepas[,-c(2)]
#convertir en dataframe
cepas <- data.frame(cepas = cepas)
placa = read_lines("pcepas")
#convertir en dataframe
placa <- data.frame(strain_name=placa)
# juntar en una sola tabla, se uso cbind ya que solamente lo queremos juntar y 
# no buscamos que coincidan en algo 
strainskarasov <- cbind(cepas, placa)

#intercambiar nombre dentro de columna strain_name para poder asociar con otra tabla
metadata$strain_name <- gsub("plate", "p", metadata$strain_name)
#quitar informacion adicional y quedarnos solamnete con nombre de placa
strainskarasov$strain_name <- gsub("Pseudomonas sp. 286 isolate", "", strainskarasov$strain_name)
strainskarasov$strain_name <- gsub("genome assembly", "", strainskarasov$strain_name)
strainskarasov$cepas <- gsub(".fasta", "", strainskarasov$cepas)

#eliminar espacios en blanco para poder emparejar 
metadata$strain_name <- trimws(metadata$strain_name)
strainskarasov$strain_name <- trimws(strainskarasov$strain_name)
strainskarasov$cepas <- trimws(strainskarasov$cepas)

# crear nueva tabla donde este toda la informacion
metadatacomp = metadata %>% 
  left_join(strainskarasov, by = "strain_name")
metadatacomp = cbind(metadatacomp, plate)

d <- sapply(metadatacomp[, 3:6], unique)

library(dplyr)

# Crear una lista para almacenar las tablas filtradas
tablas_filtradas <- list()

# Obtener los valores únicos de las columnas relevantes
unique_plants <- unique(metadatacomp$plant_num)
unique_sites <- unique(metadatacomp$site_collection)
unique_dates <- unique(metadatacomp$date_collection)
unique_leaves <- unique(metadatacomp$leaf_num)
unique_OTU <- unique(metadatacomp$OTU_assignment)

# Recorrer los valores únicos de cada columna
for (planta in unique_plants) {
  for (site in unique_sites) {
    for (date in unique_dates) {
      for (leaf in unique_leaves) {
        for (OTU in unique_OTU) {
          # Filtrar el dataframe 'metadatacomp' con los valores actuales de cada columna
          filtrada <- metadatacomp %>%
            filter(
              plant_num == planta,
              site_collection == site,
              date_collection == date,
              leaf_num == leaf, 
              OTU_assignment == OTU
              )
          # Crear un nombre de tabla combinando los valores de los criterios
          nombre_tabla <- paste(planta, site, date, leaf, OTU, sep = "_")
          if (nrow(filtrada) > 0 && ncol(filtrada) > 0) {
            tablas_filtradas[[nombre_tabla]] <- filtrada
          }
        }
      }
    }
  }
}



# Imprimir los nombres de las tablas filtradas para verificar
print(names(tablas_filtradas))

#ftable(site_collection ~ plant_num + leaf_num, metadatacomp)

library(tidyverse)
x = metadatacomp %>%
  group_by(site_collection, OTU_assignment, plant_num, date_collection, leaf_num) %>%
  summarise(n = n(),
            .groups = 'drop') %>%
  arrange(desc(n))

## Tablas que elegimos 
tablas_elegidas <- list()

# Iterar sobre las tablas y aplicar el filtrado
for (nombre_tabla in names(tablas_filtradas)) {
  tabla <- tablas_filtradas[[nombre_tabla]]
  
  # Verificar condiciones
  if (all(tabla$site_collection == "Eyach",
          tabla$date_collection == "11_12_2015",
          tabla$OTU_assignment == "OTU1", 
          tabla$plant_num %in% c("137", "155", "107", "99", "134", "182"))) {
    tablas_elegidas[[nombre_tabla]] <- tabla
  }
}
# Imprimir nombres de las tablas filtradas
names(tablas_elegidas)

num = c("3","1","6","4","2","5")

# Agregra columna de a que grupo pertenece
for (i in seq_along(tablas_elegidas)) {
  tablas_elegidas[[i]]$group <- num[i]
}

#hacer una lista con los nombres de los genomas que queremos
tabla_genomas = list()
nombre = c()

for (i in seq_along(tablas_elegidas)) {
  tabla_genomas[[i]] <- tablas_elegidas[[i]]$cepas
  nombre[i] = tablas_elegidas[[i]]$plant_num
  
}
names(tabla_genomas) <- nombre

## tabla sin site_collection "Eyach"
metadatasineyach = metadatacomp[metadatacomp$site_collection != "Eyach",]

library(tidyverse)
sineyach = metadatasineyach %>%
  group_by(site_collection, OTU_assignment, plant_num, date_collection, leaf_num) %>%
  summarise(n = n(),
            .groups = 'drop') %>%
  arrange(desc(n))

## buscar genomas con n5o grande 

nombre_buscado <- "UUNS01"
resultados <- list()

for (i in seq_along(tabla_genomas)) {
  if (nombre_buscado %in% tabla_genomas[[i]]) {
    resultados[[names(tabla_genomas)[i]]] <- tabla_genomas[[i]][tabla_genomas[[i]] == nombre_buscado]
  }
}

# Imprimir resultados
print(resultados)
