

# pip install pyfaidx
# faidx UNITE_public_all_21.04.2024.fasta -g "k__Viridiplantae" > ITS.fasta
# grep -e ">" ITS.fasta > ITS-headers.txt

library(data.table)
library(stringr)

dir.create("dbs/UNITE/20240424/r-curation")


tax = "dbs/UNITE/20240424/ITS-headers.txt" |> readLines()

tax = tax |> str_split("\\|", simplify = TRUE) |> as.data.frame() |> setDT()

tax$V3 = NULL

t = tax$V2 |> str_split("\\;", simplify = TRUE) |> as.data.frame() |> setDT()

tax$V2 = NULL

tax = cbind(tax, t)

colnames(tax) = c("seqid", "kingdom", "phylum", "class", "order", "family", "genus", "species")

tax$seqid = tax$seqid |> str_sub(2, -1)

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



fwrite(tax, "dbs/UNITE/20240424/r-curation/ITS.tax", row.names = FALSE, quote = FALSE, sep = "\t")


