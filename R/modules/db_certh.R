


library(data.table)
library(stringr)
library(readxl)
library(insect)

dir.create("dbs/inhouse/CERTH/r-curation", showWarnings = FALSE)

extract_certh <- function(path) {
    
    reads = path |> read_xlsx() |> setDT()
    
    # reads$Most_probable_seq = NULL
    
    reads$ID_sequence  = reads$ID_sequence |> str_remove_all("\\n")
    reads$Marker       = reads$Marker |> str_remove_all("\\n")
    reads$Organism     = reads$Organism |> str_remove_all("\\n")
    reads$Taxonomy     = reads$Taxonomy |> str_remove_all("\\n")
    reads$Area         = reads$Area |> str_remove_all("\\n")
    
    return(reads)
    
}


## certh ----------------------------
p   <- "dbs/inhouse/CERTH/Master_File_Total.xlsx"
tax <- extract_certh(p)

taxonomy_ncbi = taxonomy()


taxids = get_taxID(tax$Organism, taxonomy_ncbi)


l = get_lineage(taxids, db = taxonomy_ncbi) |>
    lapply(t) |>
    lapply(as.data.table) |>
    rbindlist(idcol = "index", use.names = TRUE, fill = TRUE)

l = l[, c("index", "kingdom", "phylum", "class", "order", "family", "genus", "species"), with = FALSE]

tax = cbind(tax, l[, 2:ncol(l)])

colnames(tax)[1] = "seqid"


tax = tax[which(kingdom == "Viridiplantae")]


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

library(seqinr)


tax = tax |> split(tax$Marker)


for(i in names(tax)) {
    
    t = tax[[i]]
    
    write.fasta(
        sequences = t$Most_probable_seq |> str_to_lower() |> as.list(),
        names = t$seqid |> as.list(),
        file.out = paste0("dbs/inhouse/CERTH/r-curation/", i, ".fasta"),
        nbchar = 120
    )
    
    t$Most_probable_seq = NULL
    
    t = t[, c(
        "seqid",
        "kingdom", "phylum", "class", "order", "family", "genus", "species",
        "taxon_name", "taxon_rank",
        "Marker", "Organism", "Taxonomy", "Area"
    ), with = FALSE]
    
    fwrite(
        t, file = paste0("dbs/inhouse/CERTH/r-curation/", i, ".tax"),
        row.names = FALSE, quote = FALSE, sep = "\t"
    )
    
}
























