


# load libraries -------------------

library(data.table)
library(stringr)

# tax files --------------------

fls = "dbs/BOLD/20240510/" |> list.files(full.names = TRUE, pattern = "tsv")


clean_bold <- function(path) {
    
    x = path |> fread(sep = "\t", quote = "", fill = TRUE)
    
    colnames(x)[2] = "seqid"
    
    return(x)
    
}

# workflow -----------------------

for(i in fls) {
    
    q = clean_bold(i)
    
    fwrite(q, paste0(str_sub(i, 1, -4), "tax"), row.names = FALSE, quote = FALSE, sep = "\t")
}
