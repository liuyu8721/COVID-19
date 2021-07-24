Args <- commandArgs(TRUE)

require(stringr)

fasta <- read.table(paste0(Args[1], ".fasta"), col.names = "seq", quote = "", sep = "\n")
header <- word(fasta$seq[grep("^>", fasta$seq)], 2, 2, fixed("|"))
header.row <- grep("^>", fasta$seq)
header.end <- c(header.row[-1] - 1, nrow(fasta))
fasta.info <- data.frame(id = header,
                         start = header.row,
                         end = header.end,
                         stringsAsFactors = F)
write.table(fasta.info, paste0(Args[1], ".header"), col.names = F, row.names = F, quote = F, sep=",")
