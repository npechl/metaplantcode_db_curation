




# load libraries -------------------

library(data.table)
library(stringr)

library(seqinr)

# tax files --------------------

dir.create("dbs/PLANITS/20200329/r-curation", showWarnings = FALSE)

fls = "dbs/PLANITS/20200329/" |> list.files(full.names = TRUE, pattern = "tsv") |> sort()
fas = "dbs/PLANITS/20200329/" |> list.files(full.names = TRUE, pattern = "fasta") |> sort()

clean_planits <- function(path) {
    
    x = path |> fread(sep = "\t", quote = "", fill = TRUE, header = FALSE)
    
    colnames(x) = c("seqid", "taxonomy")
    
    y = x$taxonomy |> str_split("\\;", simplify = TRUE) |> as.data.table()
    
    y$kingdom = ifelse(y$V1 == "Chlorophyta", "Viridiplantae", "Plantae")
    
    colnames(y) = c("phylum", "class", "order", "family", "genus", "species", "kingdom")
    
    x$taxonomy = NULL
    
    x = cbind(x, y)
    
    x = x[, c("seqid", "kingdom", "phylum", "class", "order", "family", "genus", "species"), with = FALSE]
    
    x = x |> lapply(function(q) { ifelse(q == "NA", NA, q) }) |> setDT()
    
    x$taxon_name = ifelse(
        !is.na(x$species), x$species,
        ifelse(
            !is.na(x$genus), x$genus,
            ifelse(
                !is.na(x$family), x$family,
                ifelse(
                    !is.na(x$order), x$order,
                    ifelse(
                        !is.na(x$class), x$class,
                        ifelse(!is.na(x$phylum), x$phylum, x$kingdom)
                    )
                )
            )
        )
    )
    
    x$taxon_rank = ifelse(
        !is.na(x$species), "species",
        ifelse(
            !is.na(x$genus), "genus",
            ifelse(
                !is.na(x$family), "family",
                ifelse(
                    !is.na(x$order), "order",
                    ifelse(
                        !is.na(x$class), "class",
                        ifelse(!is.na(x$phylum), "phylum", "kingdom")
                    )
                )
            )
        )
    )
    
    return(x)
    
}

# workflow -----------------------

for(i in seq_along(fls)) {
    
    q = fls[i] |> clean_planits()
    f = fas[i] |> read.fasta()
    
    marker = fls[i] |> basename() |> str_split_i("\\.", 1)
    
    names(f) = paste0(names(f), "_", marker)
    q$seqid  = paste0(q$seqid, "_", marker)
    
    write.fasta(f, names(f) |> as.list(), file.out = paste0("dbs/PLANITS/20200329/r-curation/", marker, ".fasta"), nbchar = 120)
    fwrite(q, paste0("dbs/PLANITS/20200329/r-curation/", marker, ".tax"), row.names = FALSE, quote = FALSE, sep = "\t")
}
