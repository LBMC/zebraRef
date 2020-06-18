#' Dre_genes
#' 
#' A gene ID table for *D. rerio*.
#' 
#' This table was built directly from the ensembl biomaRt.
#' 
#' 
#' @docType data
#' 
#' @section Format:
#' A dataframe with the folowing columns :
#' 
#'  - `ens_id`: Ensembl ID (*e.g* `ENSDARG00000069301`)
#'  - `transcript_id`: Transcript ID (*e.g* `ENSDART00000100744`)
#'  - `public_name`: Gene name (*e.g* `tmem177`)
#'  - `transcript_length`: Transcript length from 5'UTR start to 3'UTR end (*e.g* `2099`)
#' 
#' 
#' @source \href{https://www.ensembl.org/index.html}{ensembl biomaRt} 
#' 
"Dre_genes"