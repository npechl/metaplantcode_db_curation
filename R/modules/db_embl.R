




# load libraries -------------------

library(data.table)
library(stringr)

library(insect)

dir.create("dbs/EMBL/r-curation")


extract_embl <- function(path) {
    
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

p   <- "dbs/EMBL/trnL.fasta"
tax <- extract_embl(p)



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

colnames(tax)[c(1, 4, 9)] = c("seqid", "taxon_name", "taxon_rank")

tax$taxon_name = ifelse(is.na(tax$taxon_name), tax$species, tax$taxon_name)
tax$taxon_rank = ifelse(is.na(tax$taxon_rank), "species", tax$taxon_rank)

tax = tax[which(kingdom == "Viridiplantae"), c(
    "seqid", "kingdom", "phylum", "class", "order", "family", "genus", "species",
    "taxon_name",  "taxid", "taxon_rank",
    "count",           "merged_taxid",     "family_name",                                              
    "reverse_match",   "scientific_name",                                          
    "forward_error",   "forward_tm",       "genus_name",                                              
    "seq_length_ori",  "forward_match",    "reverse_tm",                                              
    "reverse_error",   "strand",           "info",                                                    
    "info_2",          "info_3",           "info_4",                                                  
    "info_5",          "info_6",           "info_7",                                                  
    "Hypnum hamulosum partial tRNA-Leu gene, specimen_voucher",
    "Corynebacterium doosanense CAU 212 ",
    "Vibrio vulnificus NBRC 15645 "                           
                                                       
                                                       
    
), with = FALSE]

f = read.fasta("dbs/EMBL/trnL.fasta")
f = f[tax$seqid]

write.fasta(f, names(f) |> as.list(), "dbs/EMBL/r-curation/trnL.fasta", nbchar = 120)
fwrite(tax, "dbs/EMBL/r-curation/trnL.tax", row.names = FALSE, quote = FALSE, sep = "\t")

