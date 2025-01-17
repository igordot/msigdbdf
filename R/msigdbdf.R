#' Retrieve the MSigDB gene sets
#'
#' Retrieve a data frame of MSigDB gene sets and their member genes.
#' Starting with release 2022.1, MSigDB was split into human and mouse resources, each one provided in the approved gene symbols of its respective species.
#' The MSigDB versioning convention is in the format `Year.Release.Species`.
#' The species referenced in this function is the one specified in the release version.
#'
#' @param target_species Species abbreviation for human or mouse databases (`"HS"` or `"MM"`).
#'
#' @return A tibble (a data frame with class [`tbl_df`]) of gene sets with one gene per row.
#'
#' @importFrom dplyr arrange inner_join
#'
#' @export
msigdbdf <- function(target_species = c("HS", "MM")) {
  target_species <- toupper(target_species)
  target_species <- match.arg(target_species)
  if (target_species == "HS") {
    mdb <- dplyr::inner_join(gene_set_members_hs, gene_set_details_hs, by = "gs_id")
  }
  if (target_species == "MM") {
    mdb <- dplyr::inner_join(gene_set_members_mm, gene_set_details_mm, by = "gs_id")
  }

  mdb <- dplyr::arrange(
    mdb,
    .data$gs_id,
    .data$db_gene_symbol,
    .data$db_ensembl_gene,
    .data$source_gene
  )

  return(mdb)
}
