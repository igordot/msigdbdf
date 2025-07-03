# Exporting data for msigdbr

# Load packages
library(msigdbdf)
library(glue)

# Save human data
mdb_hs <- as.data.frame(msigdbdf::msigdbdf(target_species = "HS"))
version_hs <- unique(mdb_hs$db_version)
stopifnot(length(version_hs) == 1)
saveRDS(
  mdb_hs,
  file = glue("msigdb.{version_hs}.rds"),
  compress = "xz"
)

# Save mouse data
mdb_mm <- as.data.frame(msigdbdf::msigdbdf(target_species = "MM"))
version_mm <- unique(mdb_mm$db_version)
stopifnot(length(version_mm) == 1)
saveRDS(
  mdb_mm,
  file = glue("msigdb.{version_mm}.rds"),
  compress = "xz"
)
