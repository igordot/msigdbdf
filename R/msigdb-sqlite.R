#' Read MSigDB SQLite database
#'
#' Download the MSigDB SQLite database and extract the relevant tables as data frames.
#' Each database file holds one MSigDB release for one resource (human or mouse).
#'
#' @param x MSigDB version, such as `2023.1.Hs`.
#'
#' @returns A list of data frames.
#'
#' @references MSigDB SQLite database documentation: <https://docs.gsea-msigdb.org/#MSigDB/MSigDB_SQLite_Database/>
#'
#' @importFrom dplyr filter select tbl
#' @importFrom tibble as_tibble
msigdb_sqlite <- function(x) {
  if (!requireNamespace("DBI", quietly = TRUE)) {
    stop("Package 'DBI' must be installed to use this function")
  }
  if (!requireNamespace("RSQLite", quietly = TRUE)) {
    stop("Package 'RSQLite' must be installed to use this function")
  }

  # Define MSigDB download variables
  # https://data.broadinstitute.org/gsea-msigdb/msigdb/release/
  mdb_version <- x
  mdb_db <- str_glue("msigdb_v{mdb_version}.db")
  mdb_zip <- str_glue("{mdb_db}.zip")
  url_base <- "https://data.broadinstitute.org/gsea-msigdb/msigdb"
  mdb_zip_url <- str_glue("{url_base}/release/{mdb_version}/{mdb_zip}")

  # Check that the database version and the resulting URL are valid
  if (requireNamespace("RCurl", quietly = TRUE)) {
    if (!RCurl::url.exists(mdb_zip_url)) {
      stop("The MSigDB SQLite file URL does not exist: ", mdb_zip_url)
    }
  }

  # Download the MSigDB SQLite file
  options(timeout = 100)
  download.file(url = mdb_zip_url, destfile = mdb_zip)
  unzip(mdb_zip)

  # Check MSigDB SQLite file size in bytes
  # utils:::format.object_size(file.size(mdb_db), units = "auto")

  # Open database connection to SQLite file and extract tables as tibbles
  # https://docs.gsea-msigdb.org/#MSigDB/MSigDB_SQLite_Database/
  db <- DBI::dbConnect(RSQLite::SQLite(), dbname = mdb_db, flags = RSQLite::SQLITE_RO)

  db_list <- list()

  db_list$MSigDB <- tibble::as_tibble(dplyr::tbl(db, "MSigDB"))
  db_list$MSigDB <- dplyr::filter(db_list$MSigDB, .data$version_name == mdb_version)

  # Check that the MSigDB table is one row (corresponding to the selected version)
  if (nrow(db_list$MSigDB) != 1) {
    stop("The MSigDB table does not have exactly one entry for the specified version")
  }

  # The gene_set table holds the core information about each gene set
  # Columns: id, standard_name, collection_name, tags, license_code
  db_list$gene_set <- tibble::as_tibble(dplyr::tbl(db, "gene_set"))
  db_list$gene_set <- dplyr::select(
    db_list$gene_set,
    .data$id,
    .data$standard_name,
    .data$collection_name
  )

  # Check that the gene_set table has a reasonable number of rows
  if (nrow(db_list$gene_set) < 10000) {
    stop("The gene_set table has too few entries")
  }

  # Check that standard_name is always present
  if (any(is.na(db_list$gene_set$standard_name))) {
    stop("Missing standard_name in gene_set")
  }

  # The gene_set_details table gives a variety of additional details for each gene set
  # Columns: gene_set_id, added_in_MSigDB_id, description_brief, description_full, systematic_name, exact_source, external_details_URL, source_species_code, primary_namespace_id, second_namespace_id, num_namespaces, publication_id, GEO_id, contributor, contrib_organization, changed_in_MSigDB_id, changed_reason
  db_list$gene_set_details <- tibble::as_tibble(dplyr::tbl(db, "gene_set_details"))
  db_list$gene_set_details <- dplyr::select(
    db_list$gene_set_details,
    .data$gene_set_id,
    .data$description_brief,
    .data$description_full,
    .data$systematic_name,
    .data$external_details_URL,
    .data$source_species_code,
    .data$publication_id,
    .data$GEO_id
  )

  # Check that systematic_name is always present
  if (any(is.na(db_list$gene_set_details$systematic_name))) {
    stop("Missing systematic_name in gene_set_details")
  }

  # The gene_symbol table holds the canonical information for the genes, including the official symbol and the NCBI (formerly Entrez) Gene ID
  # The namespace_id is constant as all symbols are mapped into the same namespace
  # Columns: id, symbol, NCBI_id, namespace_id
  db_list$gene_symbol <- tibble::as_tibble(dplyr::tbl(db, "gene_symbol"))
  db_list$gene_symbol <- dplyr::select(
    db_list$gene_symbol,
    .data$id,
    .data$symbol,
    .data$NCBI_id
  )

  # The namespace table allow us to identify the mapping info associated with each gene_symbol
  # Columns: id, label, species_code
  db_list$namespace <- tibble::as_tibble(dplyr::tbl(db, "namespace"))

  # The source_member table contains original gene set member identifiers
  # The gene_symbol_id column gives the mapping to our uniformly mapped gene symbols
  # Columns: id, source_id, gene_symbol_id, namespace_id
  db_list$source_member <- tibble::as_tibble(dplyr::tbl(db, "source_member"))

  # The gene_set_source_member is for joining source_member identifiers
  # Columns: gene_set_id, source_member_id
  db_list$gene_set_source_member <- tibble::as_tibble(dplyr::tbl(db, "gene_set_source_member"))

  # The publication and author tables associate publication info to gene sets
  # Columns: id, title, PMID, DOI, URL
  db_list$publication <- tibble::as_tibble(dplyr::tbl(db, "publication"))
  db_list$publication <- dplyr::select(
    db_list$publication,
    .data$id,
    .data$PMID
  )

  # The collection table holds the information for each MSigDB Collection
  # Columns: id, collection_name, full_name, description, parent_collection_id, gparent_collection_id, ggparent_collection_id
  db_list$collection <- tibble::as_tibble(dplyr::tbl(db, "collection"))
  db_list$collection <- dplyr::select(
    db_list$collection,
    .data$id,
    .data$collection_name,
    .data$full_name,
    .data$description
  )

  # Close database connection
  DBI::dbDisconnect(db)

  # Delete the downloaded file
  file.remove(mdb_zip)
  file.remove(mdb_db)

  return(db_list)
}
