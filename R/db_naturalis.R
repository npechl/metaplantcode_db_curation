
library(data.table)
library(stringr)

clean_naturalis <- function(path) {
    
    x = path |> fread(quote = "", fill = TRUE, sep = "\t", header = FALSE)
    
    colnames(x) = c("seqid", "taxonomy")
    
    y = x$taxonomy |> 
        str_split("\\;") |> 
        lapply(function(q) {
            
            q = q |> str_split("\\__", n = 2, simplify = TRUE)
            
            q = q[which(q[, 1] != "" & q[, 2] != ""), ]
            
            rank_names = q[, 1]
            
            q = q[, 2] |> t() |> as.data.table()
            
            colnames(q) = rank_names 
            
            return(q)
        }) |> 
        
        rbindlist(use.names = TRUE, fill = TRUE)
    
    x = cbind(x, y)
    
    x$taxonomy = NULL
    
    colnames(x) = c("seqid", "kingdom", "phylum", "class", "order", "family", "genus", "species")
    
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



p = "dbs/inhouse/NATURALIS/rbcL1-arch.tax"

tax = clean_naturalis(p)

fwrite(tax, paste0(dirname(p), "/rbcL1.tax"), row.names = FALSE, quote = FALSE, sep = "\t")







