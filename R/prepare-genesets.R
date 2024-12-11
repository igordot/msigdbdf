#' Generate a table of gene set details
#'
#' Convert the MSigDB SQLite database tables to a single table of gene set information.
#'
#' @param x A list of data frames returned by `read_msigdb_sqlite()`.
#'
#' @returns A data frame with gene set details.
#'
#' @import dplyr stringr tibble tidyr
prepare_genesets <- function(x) {
  if (!is.list(x)) {
    stop("Input must be a list of data frames")
  }

  # Combine core information about each gene set with details
  mgs <- inner_join(x[["gene_set"]], x[["gene_set_details"]], by = c("id" = "gene_set_id"))

  # Add publication information
  mgs <- left_join(mgs, x[["publication"]], by = c("publication_id" = "id"))

  # Extract collection and subcollection information
  mgs <- separate_wider_delim(
    mgs,
    cols = .data$collection_name,
    delim = ":",
    names = c("gs_collection", "gs_subcollection"),
    too_few = "align_start",
    too_many = "merge",
    cols_remove = TRUE
  )

  # Select and rename the relevant columns
  mgs <- select(
    mgs,
    gs_id = .data$systematic_name,
    gs_name = .data$standard_name,
    gs_description = .data$description_brief,
    gs_pmid = .data$PMID,
    gs_geoid = .data$GEO_id,
    gs_url = .data$external_details_URL,
    .data$gs_collection,
    .data$gs_subcollection
  )

  # Removed archived gene sets (not relevant in 2024)
  # mgs <- filter(mgs, gs_cat != "ARCHIVED")

  # Replace NA values with empty strings
  mgs <- replace_na(
    mgs,
    list(
      gs_subcollection = "",
      gs_description = "",
      gs_pmid = "",
      gs_geoid = "",
      gs_url = ""
    )
  )

  # Clean up the final table
  mgs <- distinct(mgs)
  mgs <- arrange(mgs, .data$gs_name, .data$gs_id)

  # Check that the final table contains all the initial gene sets
  if (nrow(mgs) != nrow(x[["gene_set"]])) {
    stop("Some gene sets were lost during merging")
  }
  if (!identical(sort(mgs$gs_id), sort(x[["gene_set_details"]]$systematic_name))) {
    stop("Some gene sets were changed during merging")
  }
  if (any(is.na(mgs$gs_id))) {
    stop("NAs in gs_id")
  }
  if (any(is.na(mgs$gs_name))) {
    stop("NAs in gs_name")
  }
  if (any(is.na(mgs$gs_collection))) {
    stop("NAs in gs_collection")
  }

  return(mgs)
}
