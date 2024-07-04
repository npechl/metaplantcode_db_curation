



library(data.table)
library(stringr)



h = readLines("UNITE_20240424_ITS/headers.txt")

h = h |> str_sub(2, -1)

m = h |> str_split("\\|", simplify = TRUE) |> as.data.frame() |> setDT()

t = m$V2 |> str_split("\\;", simplify = TRUE) |> as.data.table()

m = cbind(m[, c(1, 3), with = FALSE], t)


fwrite(
    m, "UNITE_20240424_ITS/UNITE_public_all_21.04.2024_taxonomy.tsv", 
    row.names = FALSE, quote = FALSE, sep = "\t"
)
