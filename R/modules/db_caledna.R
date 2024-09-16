




# load libraries -------------------

library(data.table)
library(stringr)

# tax files --------------------

fls = "dbs/CALEDNA/" |> list.files(full.names = TRUE, pattern = "txt") |> sort()
fas = "dbs/CALEDNA/" |> list.files(full.names = TRUE, pattern = "fasta") |> sort()

dir.create("dbs/CALEDNA/r-curation", showWarnings = FALSE)

clean_caledna <- function(path) {
    
    x = path |> fread(sep = "\t", quote = "", fill = TRUE, header = FALSE)
    
    colnames(x) = c("seqid", "taxonomy")
    
    y = x$taxonomy |> str_split("\\;", simplify = TRUE) |> as.data.table()
    
    colnames(y) = c("kingdom", "phylum", "class", "order", "family", "genus", "species")
    
    x$taxonomy = NULL
    
    x = cbind(x, y)
    
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

library(seqinr)

for(i in seq_along(fls)) {
    
    marker = fls[i] |> basename() |> str_sub(1, -5)
    
    q = fls[i] |> clean_caledna()
    
    q = q[which(phylum %in% c("Streptophyta", "Chlorophyta", "Bacillariophyta", "Rhodophyta"))]
    
    f = fas[i] |> read.fasta()
    f = f[q$seqid]
    
    write.fasta(f, names(f) |> as.list(), file.out = paste0("dbs/CALEDNA/r-curation/", marker, ".fasta"), nbchar = 120)
    
    fwrite(q, paste0("dbs/CALEDNA/r-curation/", marker, ".tax"), row.names = FALSE, quote = FALSE, sep = "\t")
}
