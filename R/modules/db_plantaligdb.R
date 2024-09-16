





# load libraries -------------------

library(data.table)
library(stringr)

library(seqinr)

library(insect)

# tax files --------------------

fls = "dbs/PLANTALIGDB/" |> list.files(full.names = TRUE, pattern = "tax") |> sort()
fas = "dbs/PLANTALIGDB/" |> list.files(full.names = TRUE, pattern = "fasta") |> sort()

clean_tax <- function(path) {
    
    x = path |> readLines()
    x = x[2:length(x)]
    
    x = data.table(
        "species" = x |> str_split_i("\\t", 1) |> str_squish(),
        "seqid"   = x |> str_split_i("\\t", 2) |> str_squish()
    )
    
    return(x)
    
}


clean_fasta <- function(path) {
    
    x = path |> read.fasta()
    
    x = x |> lapply(function(q) q[which(q != "-")])
    
    names(x) = names(x) |> 
        str_remove_all("\\(|\\)|modified|reversed") |>
        str_split("\\-|\\_") |> 
        lapply(function(q) q[which(q != "")]) |>
        lapply(function(q) q[length(q)]) |>
        unlist()
    
    return(x)
}

# workflow -----------------------

dir.create("dbs/PLANTALIGDB/r-curation", showWarnings = FALSE)

taxonomy_ncbi = taxonomy()

for(i in seq_along(fls)) {
    
    marker = fls[i] |> basename() |> str_split_i("\\.", 1)
    
    t = fls[i] |> clean_tax()
    # f = fas[i] |> clean_fasta()
    # 
    # t = t[which(seqid %in% names(f))]
    # f = f[t$seqid]
    
    taxids = get_taxID(t$species, taxonomy_ncbi)
    
    
    l = get_lineage(taxids, db = taxonomy_ncbi) |>
        lapply(t) |>
        lapply(as.data.table) |>
        rbindlist(idcol = "index", use.names = TRUE, fill = TRUE)
    
    l = l[, c("index", "kingdom", "phylum", "class", "order", "family", "genus", "species"), with = FALSE]
    
    
    
    tax = cbind(t[l$index, 2], l[, 2:ncol(l)])
    
    tax = tax[which(!is.na(species) & !is.na(seqid))]
    

    f = fas[i] |> clean_fasta()
    
    tax = tax[which(tax$seqid %in% names(f))]
    
    f = f[tax$seqid] 

    
    write.fasta(f, names(f) |> as.list(), file.out = paste0("dbs/PLANTALIGDB/r-curation/", marker, ".fasta"), nbchar = 120)
    fwrite(tax, paste0("dbs/PLANTALIGDB/r-curation/", marker, ".tax"), row.names = FALSE, quote = FALSE, sep = "\t")
}










