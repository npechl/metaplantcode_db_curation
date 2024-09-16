


# load libraries -------------------

library(data.table)
library(stringr)

library(seqinr)

# tax files --------------------

fls = "dbs/BOLD/20240510/" |> list.files(full.names = TRUE, pattern = "tsv") |> sort()
fas = "dbs/BOLD/20240510/" |> list.files(full.names = TRUE, pattern = "fasta") |> sort()

clean_bold <- function(path) {
    
    x = path |> fread(sep = "\t", quote = "", fill = TRUE, showProgress = FALSE)
    
    colnames(x)[1] = "seqid"
    
    return(x)
    
}

# workflow -----------------------

dir.create("dbs/BOLD/20240510/r-curation", showWarnings = FALSE)

for(i in seq_along(fls)) {
    
    marker = fls[i] |> basename() |> str_sub(1, -5) |> str_split("_") |> lapply(function(x) x[length(x)]) |> unlist()
    
    q = fls[i] |> clean_bold()
    q = q[which(kingdom == "Plantae")]
    
    q = q[, c(
        "seqid",
        "taxid", "taxon_name", "taxon_rank",
        "kingdom", "phylum", "class", "order", "family", "genus", "species",
        "sampleid",                "fieldid",                 "museumid",                "record_id",              
        "specimenid",              "processid_minted_date",   "bin_uri",                 "bin_created_date",        "collection_code",        
        "inst",                    "subfamily",              
        "tribe",                   "subspecies",              "species_reference",      
        "identification",          "identification_method",   "identification_rank",     "identified_by",           "identifier_email",       
        "taxonomy_notes",          "sex",                     "reproduction",            "life_stage",              "short_note",
        "notes",                   "voucher_type",            "tissue_type",             "specimen_linkout",        "associated_specimens",
        "associated_taxa",         "collectors",              "collection_date_start",   "collection_date_end",     "collection_event_id",    
        "collection_time",         "collection_notes",        "geoid",                   "country/ocean",           "country_iso",            
        "province/state",          "region",                  "sector",                  "site",                    "site_code",              
        "coord",                   "coord_accuracy",          "coord_source",            "elev",                    "elev_accuracy",          
        "depth",                   "depth_accuracy",          "habitat",                 "sampling_protocol",       "nuc",                    
        "nuc_basecount",           "insdc_acs",               "funding_src",             "marker_code",             "primers_forward",        
        "primers_reverse",         "sequence_run_site",       "sequence_upload_date",    "bold_recordset_code_arr"
    ), with = FALSE]
    
    f = fas[i] |> read.fasta()
    
    names(f) = f |> names() |> str_split_i("\\|", 1)
    
    f = f[q$seqid]
    
    names(f) = paste0(names(f), "_", marker)
    q$seqid = paste0(q$seqid, "_", marker)
    
    write.fasta(f, names(f) |> as.list(), file.out = paste0("dbs/BOLD/20240510/r-curation/", marker, ".fasta"), nbchar = 120)
    fwrite(q, paste0("dbs/BOLD/20240510/r-curation/", marker, ".tax"), row.names = FALSE, quote = FALSE, sep = "\t")
}


# -----------------------------------------

# library(bold)
# library(data.table)
# 
# its  = bold_seqspec(marker = "ITS") |> setDT()
# its1 = bold_seqspec(marker = "ITS1") |> setDT()
# its2 = bold_seqspec(marker = "ITS2") |> setDT()
# trnl = bold_seqspec(marker = "trnL") |> setDT()
# 
# x = bold_seq()









