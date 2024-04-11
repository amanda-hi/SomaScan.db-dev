# ---------- Load packages and set up -------------#
suppressPackageStartupMessages({
    library(AnnotationDbi)
    library(AnnotationForge)
    library(DBI)
    library(dplyr)
    library(human.db0)
    library(org.Hs.eg.db)
    library(readxl)
    library(tidyr)
    library(withr)
})

regenerateDB <- function(file, dest, ver, overwrite = TRUE) {
    temp_db <- tempdir()
    withr::with_dir(temp_db,
                    AnnotationForge::makeDBPackage("HUMANCHIP_DB",
                                                   affy = FALSE,
                                                   prefix = "SomaScan",
                                                   fileName = file,
                                                   baseMapType = "eg", # Using Entrez Gene ("eg") IDs
                                                   version = ver,
                                                   manufacturer = "SomaLogic",
                                                   chipName = "SomaScan",
                                                   manufacturerUrl = "https://somalogic.com/somascan-platform/",
                                                   author = "Amanda Hiser"))
    file.copy(from = file.path(temp_db, "SomaScan.sqlite"), to = dest, 
              overwrite = overwrite)
    unlink(temp_db, recursive = TRUE)
}

# The 11k version of the SomaScan assay menu encompasses all analytes
# from the previous menus, and as such represents a "superset"
# of all currently available assay analytes and their unique identifiers
superset <- "v5.0"


# --------- Read in SomaScan menus ---------------#
# Starting w/ copies of the most up-to-date versions of the extended annots 
# available to SomaScan customers.
# Beginning on row 9 will skip the file headers
ann_content <- list(
    v4.0 = read_xlsx("SomaScan_V4.0_5K_Annotated_Content_20210616.xlsx", 
                     skip = 8),
    v4.1 = read_xlsx("SomaScan_V4.1_7K_Annotated_Content_20240214.xlsx", 
                     skip = 8),
    v5.0 = read_xlsx("SomaScan_V5.0_11K_Annotated_Content_20240214.xlsx", 
                     skip = 8)
)

# Input data for annotation package creation must have 2 fields:
# 1. Central key, i.e. the proprietary identifier/manufacturer's probe ID 
#    (the SomaScan SeqId)
# 2. A gene mapping identifier from an external database. Valid external IDs 
#    include GenBank, Unigene, Entrez, or RefSeq ID
# Per annotation files above, currently available IDs = UniProt & Entrez, but 
# Entrez is the only ChipDb-compatible ID


#--------- Address missing identifiers ---------------#
# Final table can't contain missing values for either of the 2 required fields 
# (see https://support.bioconductor.org/p/122523/#122973 or 
# https://github.com/Bioconductor/AnnotationForge/issues/39). 
# Missing IDs will pose a problem if not dealt with or removed. Some of these 
# missing values will be manually curated and resolved in later revisions of 
# the menu
missing_entrezIds <- lapply(ann_content, function(.x) {
    .x %>%
        dplyr::filter(is.na(`Entrez Gene ID`) | is.null(`Entrez Gene ID`)) %>%
        dplyr::select(SeqId:`UniProt ID`, contains("Entrez"))
})

lapply(missing_entrezIds, nrow)

# Take a look at the analytes that are missing Entrez Gene IDs,
# majority belong to "Fc_MOUSE" or "No Protein" category
lapply(missing_entrezIds, function(x) {
    count(x, `Target Full Name`) %>%
        arrange(desc(n)) %>%
        head(n = 10)
})

# Any analytes missing an Entrez ID *must* be removed before making the final 
# ChipDb
cleaned_tables <- lapply(ann_content, function(x) {
    dplyr::select(x, SeqId, `Entrez Gene ID`, `UniProt ID`, 
                  `Target Full Name`) %>%
        as.data.frame() %>%
        filter(!is.na(`Entrez Gene ID`))
})

# Making sure all NA entries were removed - should be FALSE for all data sets
lapply(cleaned_tables, function(x) any(is.na(x$`Entrez Gene ID`)))

# Checking dimensions of final, cleaned data frames
lapply(cleaned_tables, dim)


#------- Address protein complexes ---------------#
# Many SeqIds are associated with a protein complex, rather than a single 
# protein. This is indicated by the use of a "|" character to delineate that 
# multiple external IDs associated with the SeqId in question.
# In theory, the first ID is the most unique, and subsequent IDs are potential 
# members of the protein complex, but all are 
# theorized to be associated (at some level) with the SeqId
prot_complexes <- lapply(cleaned_tables, function(x) {
    filter(x, grepl("\\|", `Entrez Gene ID`))
})

lapply(prot_complexes, nrow) # Count of multimapping instances

cleaned_tables <- lapply(cleaned_tables, function(x) {
    x %>%
        group_by(SeqId) %>%
        tidyr::separate_rows(`Entrez Gene ID`, sep = "\\|") %>%
        # Clean up extra characters from the split
        mutate(across(`Entrez Gene ID`, ~gsub(",", "", .))) %>% 
        mutate(across(`Entrez Gene ID`, ~trimws(., which = "both"))) %>%
        ungroup()
})

# Make sure everything was cleaned up; should be 0 for each dataset
lapply(cleaned_tables, function(x) {
    mergedRowCount <- x$`Entrez Gene ID`[which(grepl("\\|", x$`Entrez Gene ID`))]
    stopifnot(length(mergedRowCount) == 0)
})


# ----------Create final data frames --------------#
# These data frames will be used to create package data objects & 
# database tables
pkg_menus <- lapply(cleaned_tables, function(x) x$SeqId)
pkg_seed <- cleaned_tables$v5.0[,1:2]
colnames(pkg_seed) <- colnames <- c("probe_id", "entrez_id")

# Must be saved to file, which will be used to generate the annotation 
# database. Per AnnotationForge docs, this file must NOT have column headers
f <- "SomaScan_BiocAnnoPkg_seqIDs_EntrezIDs.txt"
write.table(pkg_seed, f, sep = "\t", col.names = FALSE, row.names = FALSE, 
            quote = FALSE)
write.table(missing_entrezIds[[superset]], file = "SomaScan_BiocAnnoPkg_missingEntrezIDs.txt",
            sep = "\t", quote = FALSE)


#----------- Generate SQLite db --------------#
# This database will require a "HUMANCHIP_DB" schema.
# Note: The human.db0 package MUST be installed, or the AnnotationForge 
#.      functions below won't work
#grep("HUMAN", available.dbschemas(), value = T)

# Metadata information required to populate a human ChipDb
pkgmeta <- c(DBSCHEMA = "HUMANCHIP_DB",
             ORGANISM = "Human", 
             SPECIES = "Homo sapiens", 
             MANUFACTURER = "SomaLogic", 
             CHIPNAME = "SomaScan", 
             MANUFACTURERURL = "https://somalogic.com/somascan-platform/")

# This package was originally created using the function call below.
# Note: this only needed to be run ONCE, during initial package creation
# AnnotationForge::makeDBPackage("HUMANCHIP_DB",
#                                affy = FALSE,
#                                prefix = "SomaScan",
#                                fileName = f,
#                                # Using Entrez Gene IDs
#                                baseMapType = "eg",
#                                # Pre-submission pkg version required by Bioc
#                                version = "0.99.0",
#                                manufacturer = "SomaLogic",
#                                chipName = "SomaScan",
#                                author = "Amanda Hiser")

# When new SomaScan menus are released, the SQLite db must be regenerated to 
# incorporate new annotations. AnnotationForge::populateDB() does not 
# correctly link the PROBEID for ChipDBs. 
# This requires a hacky workaround to re-create the database, ignoring other 
# newly generated package components
regenerateDB(file.path("/home/rstudio/SomaScan.db-dev", f), 
             dest = ".", ver = "1.0.0")


# ---------- Create database connection--------------#
con <- dbConnect(dbDriver("SQLite"), "SomaScan.sqlite")


#----------- Modify metadata --------------#
# This step suggested by the AnnotationForge vignette "Making New Annotation 
# Packages", section 7
dbGetQuery(con, "SELECT * FROM metadata") # Check current metadata
dbSendQuery(con, paste("UPDATE metadata SET value='SomaScan.db' WHERE",
                       "value='AnnotationDbi'"))

superset <- cleaned_tables[[grep(superset, names(cleaned_tables))]][[1]]

#----------- Add tables to database ------------------#
# Create a dataframe to delineate which analytes are available in each 
# version of the menu
menu_table <- data.frame(probe_id = superset,
                         v4_0 = ifelse(superset %in% cleaned_tables$v4.0[[1]], 1, 0),
                         v4_1 = ifelse(superset %in% cleaned_tables$v4.1[[1]], 1, 0),
                         v5_0 = ifelse(superset %in% cleaned_tables$v5.0[[1]], 1, 0))

# These won't be identical due to probe ID duplicates, 
# but should be in the same ballpark
apply(menu_table[,-1], 2, sum)

# Adding above data as a table to the SQLite database
dbGetQuery(con, 
           paste("CREATE Table map_menu", 
                 "(probe_id TEXT, v4_0 INTEGER, v4_1 INTEGER, v5_0 INTEGER)")
)

dbAppendTable(con, "map_menu", menu_table)

# Adding target full name table to SQLite database
seqId2targName <- dplyr::select(cleaned_tables$v5.0, 
                                probe_id = "SeqId", 
                                target_full_name = "Target Full Name")

# These should be identical
identical(pkg_seed$probe_id, seqId2targName$probe_id) 

dbGetQuery(con, paste("CREATE Table target_names", 
                      "(probe_id TEXT, target_full_name TEXT)"))
dbAppendTable(con, "target_names", seqId2targName)


# ------- Close database connection ------#
dbDisconnect(con)
