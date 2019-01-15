
#this file is for package required data including:
#information for all cell types' marker genes
#trichotomous matrix A
#mean expression level for all signature genes

# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'



#' Siganture genes' information for immune cells and melanoma
#'
#' @description A dataset contains all 53 signature genes' information for 10 immune cell types and
#' melanoma, including their gene symbols, corresponding primary and secondary cell types
#' and supporting evidence for each gene being signature. Notice that the sigantures are
#' subjective to changes between different tumor types. For a detailed introduction to primary and
#' secondary cell types, please refer to our paper.
#'
#' @name signature_marker_melanoma
#' @docType data
#' @format this is a dataframe with 53 rows and 4 variables:
#' \describe {
#'   \item{Cell}{to which cell type is this gene a primary signature}
#'   \item{Secondary_cells}{to which cell type(s) is this gene a secondary signature}
#'   \item{gene_symbol}{the gene symbol}
#'   \item{Supporting_Evidence}{links or publications that justify each gene as a siganture gene}
#' }
NULL



#' Trichotomous immune cells and melanoma
#'
#' @description This is a trichotomous guide matrix (53 by 11) corresponding to the 10 immune cell types and melanoma
#' each row represents a gene. For the i-j entry, it is encoded as 1 when this gene i is NOT a signature gene to
#' celltype j; and it is encoded as 0.5 it is a secodnary siganture gene for the j cell type; and 0 if it is a primary
#' signature gene for j-th cell type. This encoding is counter-intuitive (if you think it should be the other way around)
#' is because it is easier to take this form for penalization. PLEASE PAY ATTENTION THAT HERE THE A IS NOT GUIDE MATRIX A, 
#' but 1-A 
#'
#' @name A_melanoma_v5
#' @docType data
NULL


#' Siganture genes' mean expression in training data
#'
#' @description This vector is the mean expression across all training data of 10 immune cell types, it is
#' to provide information on general relative gene expression level of all signature genes, which is used
#' in some normalizationg circumstances, see our paper for more details
#'
#' @name row_mean_v5
#' @docType data
#' @format this is a numeric vector length 53
NULL

#' Example of cell types matching table
#'
#' @description This is a dataframe that serves as reference when testing NITUMID's performance on datasets that with
#' different cell types specification. For example, if you have a melanoma dataset with known cell types proportions
#' for T cells, DC, Macrophages, Mature B cells (MB), NK cells and other cells together, then it does not match NITUMID's 11-celltype
#' setup. In this case, you could use this dataframe to specify the matching relationships between
#' your known cell types and NITUMID's cell types. This dataframe has two column, `origina_index` marks our 11 cell types
#' from 1 to 11, you are not supposed to modify that; you can modify `destin_index` accordingly as DC-1, CD4-2, CD8-2, Macrophages-3,
#' MB-4, NK-5, Monocytes-6, Plasma-6, MAST-6, Eosinophil-6, Melanoma-6. Basically, we merged CD4+ T cell and CD8+ T cell into one cell
#' type, and all our cell types after monocytes into the 6-th cell type.
#'\describe {
#'   \item{origina_index}{Indexed NITUMID's 11 cell types, from 1 to 11}
#'   \item{destin_index}{Indexed customized cell types }
#' }
#' @name cell_structure_example
#' @docType data
NULL

#' Example CD8+ T cell gene expression data
#'
#' @description This is an example dataset of 20 CD8+ T cells
#' @name ematrix_GSE6740_CD8
#' @docType data
NULL

#' Example bulk melanoma gene expressiond ata
#'
#' @description This is an example dataset of 119 bulk melanoma samples
#' @name Nivolumab_rna_sub
#' @docType data
NULL

#' Example bulk melanoma gene expression with known fractions
#'
#' @description This is an example dataset of 4 bulk melanoma samples with known fractions
#' @name real_mix_gene_expr
#' @docType data
NULL

#' Bulk melanoma samples' fractions
#'
#' @description This is the true fractions table for `real_mix_gene_expr`
#' @name real_H_mix
#' @docType data
NULL

#' Bulk melanoma samples' cell structure
#'
#' @description This is the cell structure table for `real_mix_gene_expr`
#' @name real_H_mix
#' @docType data
NULL



