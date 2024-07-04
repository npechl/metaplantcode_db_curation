
library(data.table)
library(stringr)

extract_davidlab <- function(path) {
    
    reads <- path |> readLines()
    reads <- reads[which(str_sub(reads, 1, 1) == ">")] |> str_sub(2, -1)
    
    reads = data.table(
        "seqid" = reads |> str_split_i("\\ ", 1),
        "organism" = reads |> str_split("\\ ", n = 2) |> lapply(function(q) q[2]) |> unlist()
    )
    
    return(reads)
}





# workflow -------------------------------------------



## davidlab ----------------------------
p   <- "dbs/inhouse/DAVIDLAB/trnL.fasta"
tax <- extract_davidlab(p)

taxonomy_ncbi = taxonomy()


taxids = get_taxID(tax$organism, taxonomy_ncbi)


l = get_lineage(taxids, db = taxonomy_ncbi) |>
    lapply(t) |>
    lapply(as.data.table) |>
    rbindlist(idcol = "index", use.names = TRUE, fill = TRUE)

l = l[, c("index", "kingdom", "phylum", "class", "order", "family", "genus", "species"), with = FALSE]

tax = cbind(tax, l[, 2:ncol(l)])

tax$taxon_name = ifelse(
    !is.na(tax$species), tax$species,
    ifelse(
        !is.na(tax$genus), tax$genus,
        ifelse(
            !is.na(tax$family), tax$family,
            ifelse(
                !is.na(tax$order), tax$order,
                ifelse(
                    !is.na(tax$class), tax$class,
                    ifelse(!is.na(tax$phylum), tax$phylum, tax$kingdom)
                )
            )
        )
    )
)

tax$taxon_rank = ifelse(
    !is.na(tax$species), "species",
    ifelse(
        !is.na(tax$genus), "genus",
        ifelse(
            !is.na(tax$family), "family",
            ifelse(
                !is.na(tax$order), "order",
                ifelse(
                    !is.na(tax$class), "class",
                    ifelse(!is.na(tax$phylum), "phylum", "kingdom")
                )
            )
        )
    )
)

fwrite(tax, paste0(dirname(p), "/trnL.tax"), row.names = FALSE, quote = FALSE, sep = "\t")






