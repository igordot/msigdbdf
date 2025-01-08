#' Retrieve the gene sets data frame
#'
#' Retrieve a data frame of gene sets and their member genes.
#' Starting with release 2022.1, MSigDB is split into human and mouse resources.
#' Each one is provided in the approved gene symbols of its respective species.
#' The versioning convention of MSigDB is in the format `Year.Release.Species`.
#'
#' @param species Species abbreviation for human or mouse databases (`"Hs"` or `"Mm"`).
#'
#' @return A tibble (a data frame with class [`tbl_df`]) of gene sets with one gene per row.
#'
#' @importFrom dplyr inner_join
#'
#' @export
msigdbrdata <- function(species = c("Hs", "Mm")) {
  species <- match.arg(species)
  if (species == "Hs") {
    gs <- dplyr::inner_join(gene_set_details_hs, gene_set_members_hs, by = "gs_id")
  }
  if (species == "Mm") {
    gs <- dplyr::inner_join(gene_set_details_mm, gene_set_members_mm, by = "gs_id")
  }

  return(gs)
}
