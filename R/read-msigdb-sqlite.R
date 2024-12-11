#' Read MSigDB SQLite database
#'
#' Download the MSigDB SQLite database and extract the relevant tables as data frames.
#'
#' @param x MSigDB version, such as `2023.1.Hs`.
#'
#' @returns A list of data frames.
#'
#' @import dplyr stringr tibble utils
read_msigdb_sqlite <- function(x) {
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
  mdb_url_base <- "https://data.broadinstitute.org/gsea-msigdb/msigdb"
  mdb_zip_url <- str_glue("{mdb_url_base}/release/{mdb_version}/{mdb_zip}")

  # Download the MSigDB SQLite file
  options(timeout = 100)
  download.file(url = mdb_zip_url, destfile = mdb_zip)
  unzip(mdb_zip, exdir = ".")

  # Check MSigDB SQLite file size in bytes
  # utils:::format.object_size(file.size(mdb_db), units = "auto")

  # Open database connection to SQLite file and extract tables as tibbles
  # https://docs.gsea-msigdb.org/#MSigDB/MSigDB_SQLite_Database/
  db <- DBI::dbConnect(RSQLite::SQLite(), dbname = mdb_db, flags = RSQLite::SQLITE_RO)

  db_list <- list()

  # The gene_set table holds the core information about each gene set
  # Columns: id, standard_name, collection_name, tags, license_code
  gene_set <- as_tibble(tbl(db, "gene_set"))
  gene_set <- select(
    gene_set,
    .data$id,
    .data$standard_name,
    .data$collection_name
  )
  db_list[["gene_set"]] <- gene_set

  # Check that standard_name is always present
  if (any(is.na(gene_set_details$standard_name))) {
    stop("Missing standard_name in gene_set")
  }

  # The gene_set_details table gives a variety of additional details for each gene set
  # Columns: gene_set_id, added_in_MSigDB_id, description_brief, description_full, systematic_name, exact_source, external_details_URL, source_species_code, primary_namespace_id, second_namespace_id, num_namespaces, publication_id, GEO_id, contributor, contrib_organization, changed_in_MSigDB_id, changed_reason
  gene_set_details <- as_tibble(tbl(db, "gene_set_details"))
  gene_set_details <- select(
    gene_set_details,
    .data$gene_set_id,
    .data$description_brief,
    .data$description_full,
    .data$systematic_name,
    .data$external_details_URL,
    .data$source_species_code,
    .data$publication_id,
    .data$GEO_id
  )
  db_list[["gene_set_details"]] <- gene_set_details

  # Check that the systematic_name is always present
  if (any(is.na(gene_set_details$systematic_name))) {
    stop("Missing systematic_name in gene_set_details")
  }

  # The gene_symbol table holds the canonical information for the genes, including the official symbol and the NCBI (formerly Entrez) Gene ID
  # The namespace_id is constant as all symbols are mapped into the same namespace
  # Columns: id, symbol, NCBI_id, namespace_id
  gene_symbol <- as_tibble(tbl(db, "gene_symbol"))
  gene_symbol <- select(
    gene_symbol,
    .data$id,
    .data$symbol,
    .data$NCBI_id
  )
  gene_symbol$NCBI_id <- as.integer(gene_symbol$NCBI_id)
  db_list[["gene_symbol"]] <- gene_symbol

  # The namespace table allow us to identify the mapping info associated with each gene_symbol
  # Columns: id, label, species_code
  namespace <- as_tibble(tbl(db, "namespace"))
  db_list[["namespace"]] <- namespace

  # The source_member table contains original gene set member identifiers
  # The gene_symbol_id column gives the mapping to our uniformly mapped gene symbols
  # Columns: id, source_id, gene_symbol_id, namespace_id
  source_member <- as_tibble(tbl(db, "source_member"))
  db_list[["source_member"]] <- source_member

  # The gene_set_source_member is for joining source_member identifiers
  # Columns: gene_set_id, source_member_id
  gene_set_source_member <- as_tibble(tbl(db, "gene_set_source_member"))
  db_list[["gene_set_source_member"]] <- gene_set_source_member

  # The publication and author tables associate publication info to gene sets
  # Columns: id, title, PMID, DOI, URL
  publication <- as_tibble(tbl(db, "publication"))
  publication <- select(
    publication,
    .data$id,
    .data$PMID
  )
  db_list[["publication"]] <- publication

  # The collection table holds the information for each MSigDB Collection
  # Columns: id, collection_name, full_name, description, parent_collection_id, gparent_collection_id, ggparent_collection_id
  collection <- as_tibble(tbl(db, "collection"))
  collection <- select(
    collection,
    .data$id,
    .data$collection_name,
    .data$full_name,
    .data$description
  )
  db_list[["collection"]] <- collection

  # Close database connection to MSigDB SQLite
  DBI::dbDisconnect(db)

  # Delete the downloaded file
  file.remove(mdb_zip)
  file.remove(mdb_db)

  return(db_list)
}
