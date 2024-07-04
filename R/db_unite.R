

library(data.table)
library(stringr)



tax = "dbs/UNITE/20240424/UNITE_public_all_21.04.2024_taxonomy.tsv" |> fread(header = TRUE)

tax[[2]] = NULL


colnames(tax) = c("seqid", "kingdom", "phylum", "class", "order", "family", "genus", "species")

for(i in 2:ncol(tax)) {
    
    tax[[i]] = tax[[i]] |> str_split_i("\\__", 2)
    
}


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



fwrite(tax, "dbs/UNITE/20240424/UNITE_public_all_21.04.2024_taxonomy.tax", row.names = FALSE, quote = FALSE, sep = "\t")


