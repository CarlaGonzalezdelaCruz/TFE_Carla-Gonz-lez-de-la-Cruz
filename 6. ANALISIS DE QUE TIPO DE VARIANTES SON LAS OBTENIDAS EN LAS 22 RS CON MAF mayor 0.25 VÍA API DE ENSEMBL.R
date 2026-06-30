## 6. ANALISIS DE QUE TIPO DE VARIANTES SON LAS OBTENIDAS EN LAS 22 RS CON MAF > 0.25 VÍA API DE ENSEMBL
```r

library(httr2)
library(jsonlite)
library(openxlsx)

#22 rsIDs originales
rsids <- c("rs759174", "rs28946889", "rs4148326", "rs4663971", "rs4148328", 
           "rs4148329", "rs6717546", "rs6719561", "rs2003363", "rs56133230", 
           "rs10203853", "rs10209214", "rs6728940", "rs6746002", "rs1500475", 
           "rs6728520", "rs6735414", "rs10202032", "rs10199512", "rs4663336", 
           "rs6723936", "rs62648950")

url_vep <- "https://grch37.rest.ensembl.org/vep/human/id?hgvs=1"
body_data <- list(ids = rsids)

#Consulta a la API
message("1. Consultando Ensembl VEP (Extrayendo Dataset Máximo)...")
req <- request(url_vep) %>%
  req_headers(`Content-Type` = "application/json", `Accept` = "application/json") %>%
  req_body_raw(jsonlite::toJSON(body_data, auto_unbox = TRUE), "application/json") %>%
  req_timeout(20) # Evita bloqueos en el servidor

res <- req_perform(req)

#Procesamiento de los campos
if (resp_status(res) == 200) {
  raw_data <- resp_body_string(res) %>% fromJSON(flatten = TRUE)
  
  message("2. Desempaquetando variables globales y consecuencias de transcritos...")
  lista_total <- list()
  contador <- 1
  
  for (i in 1:nrow(raw_data)) {
    variante <- raw_data[i, ]
    
    # --- Datos Nivel Variante ---
    rsid       <- variante$id
    cromosoma  <- if("seq_region_name" %in% names(variante)) variante$seq_region_name else NA
    pos_inicio <- if("start" %in% names(variante)) variante$start else NA
    pos_fin    <- if("end" %in% names(variante)) variante$end else NA
    hebra      <- if("strand" %in% names(variante)) variante$strand else NA
    localizacion_global <- paste0("chr", cromosoma, ":", pos_inicio, "-", pos_fin, " (", ifelse(hebra == 1, "+", "-"), ")")
    
    # Bloque de transcritos
    consecuencias <- variante$transcript_consequences[[1]]
    
    if (!is.null(consecuencias) && nrow(consecuencias) > 0) {
      for (j in 1:nrow(consecuencias)) {
        tx <- consecuencias[j, ]
        
        # Extracción y mapeo de los campos que has solicitado 
        tx_id        <- if("transcript_id" %in% names(tx)) tx$transcript_id else NA
        coding_impact<- if("consequence_terms" %in% names(tx)) paste(unlist(tx$consequence_terms), collapse = ", ") else NA
        gen          <- if("gene_symbol" %in% names(tx)) tx$gene_symbol else NA
        
        # Nomenclaturas HGVS puras de Ensembl
        hgvs_coding  <- if("hgvsc" %in% names(tx)) tx$hgvsc else "No disponible/No codificante"
        hgvs_protein <- if("hgvsp" %in% names(tx)) tx$hgvsp else "No modifica proteína"
        
        # Datos moleculares del transcrito
        splice_dist  <- if("splice_distance" %in% names(tx)) tx$splice_distance else "N/A"
        tsl          <- if("tsl" %in% names(tx)) tx$tsl else NA
        mane_select  <- if("mane_select" %in% names(tx)) tx$mane_select else "N/A"
        appris       <- if("appris" %in% names(tx)) tx$appris else "N/A"
        
        # Accesiones a bases de datos de proteínas externas
        uniprot_acc  <- if("swissprot" %in% names(tx)) {
          paste(unlist(tx$swissprot), collapse = ", ")
        } else if ("uniprot_isoform" %in% names(tx)) {
          paste(unlist(tx$uniprot_isoform), collapse = ", ")
        } else {
          "N/A"
        }
        
        biotype      <- if("biotype" %in% names(tx)) tx$biotype else NA
        
        # Guardar todo estructurado en la fila
        lista_total[[contador]] <- data.frame(
          RSID = rsid,
          Location_hg19 = localizacion_global,
          Chrom = cromosoma,
          Pos_Start = pos_inicio,
          Transcript = tx_id,
          Gene = gen,
          Coding_Impact = coding_impact,
          HGVS_Coding = hgvs_coding,
          HGVS_Protein = hgvs_protein,
          Splice_Distance = splice_dist,
          TSL = tsl,
          MANE_Select = mane_select,
          APPRIS = appris,
          UniProt_Accession = uniprot_acc,
          Biotype = biotype,
          stringsAsFactors = FALSE
        )
        contador <- contador + 1
      }
    }
  }
  
  tabla_maxima <- do.call(rbind, lista_total)
  
  #Exporotar a Excel
  message("3. Escribiendo reporte masivo en Excel...")
  wb <- createWorkbook()
  addWorksheet(wb, "Anotación Completa VEP")
  
  #Estilo
  estilo_top <- createStyle(fontName = "Calibri", fontColour = "#FFFFFF", 
                            textDecoration = "bold", bgFill = "#1F4E79", halign = "center")
  
  writeData(wb, "Anotación Completa VEP", tabla_maxima, headerStyle = estilo_top)
  setColWidths(wb, "Anotación Completa VEP", cols = 1:ncol(tabla_maxima), widths = "auto")
  
  saveWorkbook(wb, "Anotacion_Completa_Ensembl_UGT1A.xlsx", overwrite = TRUE)
  
  message("¡PROCESO COMPLETADO!")
  message("Revisa tu directorio de trabajo. Se ha generado: 'Anotacion_Completa_Ensembl_UGT1A.xlsx'")
  
} else {
  stop("El servidor de Ensembl no ha respondido de forma correcta. Código: ", resp_status(res))
}
```
