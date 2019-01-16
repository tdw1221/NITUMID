
#this file is for utility functions including:
#signature genes finding and data cleaning
#results wrapping up and output
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

#A function that finds signature genes positions in raw data frame/data matrix,
#It also outputs those signature genes that are not present
#this is for fomulating the right matrix for NMF


#' Find signature genes' positions in gene expression data
#'
#' `Signature_Match` takes the input gene expression data and a vector of signature gene symbols, return signature
#' genes' positions in data matrix, as well as those not-matched signature genes.
#'
#' @param gene_expression this should be a numeric matrix whose columns are samples and rows are genes,
#'                        its row names should be genes symbol
#' @param signature_genes this should be a vector of string/character speciying signature gene symbols
#'                        default value is our generic signature genes list
#' @param signature_genes_alias this optional argument allows you to give signature genes' alias, default is NULL.
#'                              this also should be a vector of string/character, same length as signature_genes, those who don't have an alias should
#'                              be set as NA.
#' @return 'Signature_Match' returns a list consists of two elements: 'matched_index' is a vector same length as signature genes
#'          its each entry specifying the corresponding signature genes' positions in 'gene_expression', for those genes don't present in
#'          data matrix, its corresponding entries are left with 1221; 'missing_row_index' is a vector indicating which genes are not found
#'          in 'gene_expression'
#' @examples data is needed otherwise it won't run, the second example won't work because the alias vector has one non-character entry
#' (having NAs in signature_genes_alias is fine)
#'  Signature_Match(gene_expression_matrix,signature_genes=c("MIA","EGFR") )
#'  \dontrun{
#'  Signature_Match(gene_expression_matrix,signature_genes=c("MIA","EGFR"),signature_genes_alias=c(NA,2) )
#'  }
#'  @export
Signature_Match <- function(gene_expression, signature_genes=as.character(signature_marker_melanoma$gene_symbol), signature_genes_alias=NULL){

  #check if the gene expression matrix and signature_genes are correctly formatted
  if (!is.matrix(gene_expression) ){
    stop("input gene_expression should be a numeric matrix!")
  }
  if ( !is.character(signature_genes)){
    stop("input signature_genes should be a character vector")
  }
  if (!is.null(signature_genes_alias)){
    if ( !is.character(signature_genes_alias)){
      stop("input signature_genes_alias should be a character vector")
    }
  }

  #match signature genes' position
  primary_match <- match(as.character(signature_marker_melanoma$gene_symbol),rownames(gene_expression))
  if (!is.null(signature_genes_alias)){
  second_match <- match(final_marker_v5_geneSymbol,rownames(gene_expression))
  primary_match[which(is.na(primary_match))] <- second_match[which(is.na(primary_match))]
  }
  missing_row_index <- which(is.na(primary_match) )
  #need to put a random number for missing signatures otherwise the NMF function will fail
  primary_match[which(is.na(primary_match) )] <- 1221
  #return
  return(list(matched_index = primary_match, missing_row_index=missing_row_index))
}

