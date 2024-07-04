




# load libraries -------------------

library(data.table)
library(stringr)

# tax files --------------------

fls = "dbs/PLANITS/20200329/" |> list.files(full.names = TRUE, pattern = "tsv")


clean_planits <- function(path) {
    
    x = path |> fread(sep = "\t", quote = "", fill = TRUE, header = FALSE)
    
    colnames(x) = c("seqid", "taxonomy")
    
    y = x$taxonomy |> str_split("\\;", simplify = TRUE) |> as.data.table()
    
    colnames(y) = c("phylum", "class", "order", "family", "genus", "species")
    
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

for(i in fls) {
    
    q = clean_planits(i)
    
    fwrite(q, paste0(str_sub(i, 1, -4), "tax"), row.names = FALSE, quote = FALSE, sep = "\t")
}
