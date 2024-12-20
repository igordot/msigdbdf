#' Generate a table of gene set details
#'
#' Convert the MSigDB SQLite database tables to a single table of gene set information.
#'
#' @param x A list of data frames returned by `msigdb_sqlite()`.
#'
#' @returns A data frame with gene set details.
#'
#' @export
gene_set_details <- function(x) {
  if (!is.list(x)) {
    stop("Input must be a list of data frames")
  }
  if (nrow(x$MSigDB) != 1) {
    stop("MSigDB data frame must have one row")
  }

  # Combine core information about each gene set with details
  mgs <- inner_join(x$gene_set, x$gene_set_details, by = c("id" = "gene_set_id"))

  # Add publication information
  mgs <- left_join(mgs, x$publication, by = c("publication_id" = "id"))

  # Extract collection and subcollection information
  mgs <- separate_wider_delim(
    mgs,
    cols = "collection_name",
    delim = ":",
    names = c("gs_collection", "gs_subcollection"),
    too_few = "align_start",
    too_many = "merge",
    cols_remove = TRUE
  )

  # Select and rename the relevant columns
  mgs <- select(
    mgs,
    "gs_collection",
    "gs_subcollection",
    gs_id = "systematic_name",
    gs_name = "standard_name",
    gs_description = "description_brief",
    gs_source_species = "source_species_code",
    gs_pmid = "PMID",
    gs_geoid = "GEO_id",
    gs_url = "external_details_URL"
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

  # Check that standard_name is always present
  if (length(x$MSigDB$version_name) != 1) {
    stop("msigdb does not have exactly one entry for the specified version")
  }

  # Add MSigDB database information
  mgs$db_version <- x$MSigDB$version_name
  mgs$db_target_species <- x$MSigDB$target_species_code

  # Clean up the final table
  mgs <- distinct(mgs)
  mgs <- arrange(mgs, .data$gs_name, .data$gs_id)

  # Check that the final table seems reasonable
  if (ncol(mgs) != 11) {
    stop("Missing columns")
  }
  if (nrow(mgs) != nrow(x$gene_set)) {
    stop("Some gene sets were lost during merging")
  }
  if (!identical(sort(mgs$gs_id), sort(x$gene_set_details$systematic_name))) {
    stop("Some gene sets were altered during merging")
  }
  if (any(is.na(mgs$gs_id))) {
    stop("NAs in column gs_id")
  }
  if (any(is.na(mgs$gs_name))) {
    stop("NAs in column gs_name")
  }
  if (any(is.na(mgs$gs_collection))) {
    stop("NAs in column gs_collection")
  }

  return(mgs)
}
