#### MATCH DE FPKM CON CLINICAL Y BIOSPECIMEN METADATA PARA CREAR SUBCONJUNTOS DE AD Y HEALTHIES ###

# Libraries
pacman::p_load("dplyr", 
               'vroom', 
               'stringr')



# Unión de las metadata (clinical/biospecimen) ---------------------------------
#Nombrar a mi data 
clinical_metadata <- vroom(file = 'ROSMAP_clinical.csv')
biospecimen_metadata <- vroom(file = 'ROSMAP_biospecimen_metadata.csv')


#combinar metadata en una sola
cli_bio_metadata <- left_join(x = clinical_metadata, 
                              y = biospecimen_metadata, 
                              by = "individualID")

#guardar la tabla en csv
vroom_write(cli_bio_metadata,
           file = "cli_bio_metadata.csv", 
           delim = ",")





# Unión de los FPKM normalizados -----------------------------------------------
#Nombrar mis datos de RNAseq
RNAseq_normalized_1_6 <- vroom(file = 'ROSMAP_RNAseq_FPKM_gene_plates_1_to_6_normalized.tsv')
RNAseq_normalized_7_8 <- vroom(file = 'ROSMAP_RNAseq_FPKM_gene_plates_7_to_8_normalized.tsv')

#Juntar los datos de expresion en un solo formato
RNAseq_normalized_1_8 <- left_join(x = RNAseq_normalized_1_6, 
                              y = RNAseq_normalized_7_8, 
                              by = c("gene_id" = "gene_id"))

#Remover data inoportuna
RNAseq_normalized_1_8 <- RNAseq_normalized_1_8 %>% 
  dplyr::select(-c("tracking_id.x", "tracking_id.y",
                   "492_120515_6", "492_120515_7")) 

#Quitar las versiones de los ensembl 
RNAseq_normalized_1_8$gene_id <- strtrim(RNAseq_normalized_1_8$gene_id, 15)

#Guardar la tabla 
vroom_write(RNAseq_normalized_1_8,
            file = "RNAseq_normalized_1_8.csv", 
            delim = ",")





# Creación de los subconjuntos de metadata de AD y healthies en mujeres --------
#Leer la data
metadata <- vroom ( file = "cli_bio_metadata.csv")
FPKM <- vroom( file = "RNAseq_normalized_1_8.csv")


#Simplificar la metadata para quedarme unicamente con los datos RNAseq en todas las mujeres
metadata_women_RNAseq<- metadata %>% 
  filter(assay == 'rnaSeq', 
         organ == 'brain',
         msex == "0")

#Crear una matriz unicamente con los datos que me interesan de la metadata 
metadata_women <- metadata_women_RNAseq[,c(18,19,3,15,16)] #18-individualID, 19-specimenID, 3-Sex, 15-CERADsc 16-cogdx


#Crear subconjunto de AD y healthies en mujeres segun la metadata 

# cogdx = Diagnostico clinico más probable en el momento de la muerte

# 1 = NCI no cognitive impairment
# 2 = MCI mild cognitive impairment and no other cause of CI
# 3 = MCI mild cognitive impairment and other cause of CI
# 4 = Alzheimer's dementia and no other cause of CI (Prob AD)
# 5 = Alzheimer's dementia and other cause of CI (Poss AD)
# 6 = Other dementia


# Ceradsc = Medida semicuantitativa de placas neuriticas para evualuación y diagnostico post mortem 

# 1 = Definite AD / Yes
# 2 = Probable AD / Yes
# 3 = Possible AD / No
# 4 = No AD / No


#Si queremos aquellas con alzheimer segun los criterios anteriores
women_AD <- metadata_women %>% 
  filter(ceradsc == '1' | ceradsc == '2' , 
         cogdx == '4')

#Si queremos sanas segun los criterios anteriores
women_healthy <- metadata_women%>% 
  filter(ceradsc == '4', 
         cogdx == '1')

#Guardar tabla
vroom_write(metadata_women,
            file = "metadata_women.csv",
            delim = ",")
vroom_write(women_AD,
            file = "women_AD.csv",
            delim = ",")
vroom_write(women_healthy,
            file = "women_healthy.csv",
            delim = ",")




# Creación de las FPKM matrix  --------------------------------------------------
# Se crea una matriz de los specimenID presentes en la matriz de FPKM y de la metadata de mujeres
colnames_woID <- substr(colnames(FPKM),1, nchar(colnames(FPKM))-2)
colnames_woID[1] <- "gene_id"


FPKM_matrix_women <- FPKM[, (colnames_woID %in% metadata_women$specimenID)] %>% 
  mutate(gene_id = FPKM$gene_id, .before = 1) %>%
  rename_at(-1, ~str_sub(., end = -3))  

FPKM_matrix_AD <-FPKM[, (colnames_woID %in% women_AD$specimenID)] %>% 
  mutate(gene_id = FPKM$gene_id, .before = 1) %>%
  rename_at(-1, ~str_sub(., end = -3))

FPKM_matrix_healthy <-FPKM[, (colnames_woID %in% women_healthy$specimenID)] %>% 
  mutate(gene_id = FPKM$gene_id, .before = 1) %>%
  rename_at(-1, ~str_sub(., end = -3))

#Guardar tablas
vroom_write(FPKM_matrix_AD,
            file = "FPKM_matrix_AD.csv",
            delim = ",")
vroom_write(FPKM_matrix_healthy,
            file = "FPKM_matrix_healthy.csv",
            delim = ",")


####### FIN ########