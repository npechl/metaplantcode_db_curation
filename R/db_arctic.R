


# load libraries -------------------

library(data.table)
library(stringr)

library(insect)

# helper functions --------------------

extract_arctic <- function(path) {
    
    
    reads <- path |> readLines()
    reads <- reads[which(str_sub(reads, 1, 1) == ">")] |> str_sub(2, -1)
    
    seq_ids <- reads |> str_split_i("\\ ", 1)
    
    reads = reads |> str_split(pattern = "\\ ", n = 2) |> lapply(function(q) q[2]) |> unlist()
    
    reads = reads |>
        str_split("\\;") |>
        lapply(str_squish) |>
        lapply(function(q) {
            q = ifelse(str_detect(q, "\\="), q, paste0("info=", q))
            
            arc_names = q |> str_split_i("\\=", 1) 
            arc_names = ifelse(duplicated(arc_names), paste0(arc_names, arc_names |> rowid(prefix = "_")), arc_names)
            
            
            q = q |> str_split_i("\\=", 2) |> t() |> as.data.table()
            
            colnames(q) = arc_names
            
            return(q)
            
        })
    
    names(reads) = seq_ids
    
    reads = reads |> rbindlist(idcol = "seqID", use.names = TRUE, fill = TRUE)
    
    return(reads)
    
}

# workflow -----------------------

p   <- "dbs/ARCTIC/trnL.fasta"
tax <- extract_arctic(p)


taxonomy_ncbi = taxonomy()

l = get_lineage(tax$taxid |> as.numeric(), db = taxonomy_ncbi) |>
    lapply(t) |>
    lapply(as.data.table) |>
    rbindlist(idcol = "index", use.names = TRUE, fill = TRUE)

l = l[, c("index", "kingdom", "phylum", "class", "order", "family", "genus", "species"), with = FALSE]

tax$family  = NULL
tax$genus   = NULL
tax$species = NULL

tax = cbind(tax, l[, 2:ncol(l)])

colnames(tax)[c(1, 5, 6)] = c("seqid", "taxon_name", "taxon_rank")

tax$taxon_name = ifelse(is.na(tax$taxon_name), tax$species, tax$taxon_name)
tax$taxon_rank = ifelse(is.na(tax$taxon_rank), "species", tax$taxon_rank)


fwrite(tax, paste0(dirname(p), "/trnL.tax"), row.names = FALSE, quote = FALSE, sep = "\t")
