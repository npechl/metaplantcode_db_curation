
rm(list = ls())

# load libraries -------------------------

library(seqinr)

library(data.table)
library(stringr)

# define functions -------------------


read_taxonomy <- function(file, nReps = 10) {
    
    df = file |> fread()
    
    df = df[which(kingdom == "Plantae")]
    
    df$taxonomy = paste0(
        "k__", df$kingdom, ";",
        "p__", df$phylum, ";",
        "c__", df$class, ";",
        "o__", df$order, ";",
        "f__", df$family, ";",
        "g__", df$genus, ";",
        "s__", df$species, ";"
    )
    
    df = df[, c(
        "processid", 
        "taxonomy"
    ), with = FALSE]
    
    key = df[, by = taxonomy, head(.SD, 1)]
    key = key$processid |> sample(nReps)
    
    
    return(
        list(
            "taxonomy" = df,
            "rep_ids" = key
        )
    )
}

read_sequences <- function(file, taxa) {
    
    s = file |> read.fasta()
    
    s = s |> lapply(str_to_upper)
    
    names(s) = s |> names() |> str_split_i("\\|", 1)
    
    s = s[taxa$taxonomy$processid]
    
    s_reps = s[taxa$rep_ids]
    
    return(
        list(
            "sequences" = s,
            "reps"      = s_reps
        )
    )
    
}

write_taxonomy <- function(x, file_name) {
    
    fwrite(x, file_name, row.names = FALSE, col.names = FALSE, sep = "\t")
    
}

write_sequences <- function(sequences, reps, file_name) {
    
    write.fasta(sequences, names(sequences), file.out = paste0(file_name, "_seqs.fa"), nbchar = 3000)
    write.fasta(reps, names(reps), file.out = paste0(file_name, "_reps.fa"), nbchar = 3000)
    
}

# input files ----------------------

taxonomy_file <- "BOLD_20240510/BOLD_Public.10-May-2024_ITS1.tsv"
seqs_file     <- "BOLD_20240510/BOLD_Public.10-May-2024_ITS1.fasta"
out_name      <- "BOLD_20240510/BOLD_Public.10-May-2024_ITS1"

# run functions -----------------------

o1 = taxonomy_file |> read_taxonomy()
o2 = seqs_file |> read_sequences(o1)

write_taxonomy(o1$taxonomy, paste0(out_name, ".tax"))
write_sequences(o2$sequences, o2$reps, out_name)
