#taxonomy:"Plasmodium [5820]" annotation:(type:signal evidence:manual) AND reviewed:yes
#111 proteins
#taxonomy:"Plasmodium [5820]" NOT annotation:(type:signal evidence:manual) AND reviewed:yes
#361 proteins

read_signalp41 <- function(connection) {
  all_lines <- readLines(connection)[-c(1L:2)]
  do.call(rbind, lapply(all_lines, function(i) {
    line <- strsplit(i, " ")[[1]]
    line <- line[!line == ""]
    res <- data.frame(sp.probability = line[10] == "Y",
                      sp.start = ifelse(line[10] == "Y", 1, NA),
                      sp.end = ifelse(line[10] == "Y", as.numeric(line[5]) - 1, NA))
    rownames(res) <- line[1]
    res
  }))
}

if(Sys.info()["nodename"] == "phobos" )
  pathway <- "/home/michal/Dropbox/signal-peptide2_data/"


library(signalHsmm)
library(seqinr)
library(hmeasure)

signal.hsmm1987 <- train_hsmm(read_uniprot(paste0(pathway, "pos_ultrahard_data.txt"), euk = TRUE),
                              aaaggregation)
plas_pos <- read.fasta("plas_pos.fasta")
plas_neg <- read.fasta("plas_neg.fasta")

system("/home/michal/signalp-4.1/signalp -t euk -f short plas_neg.fasta > plas_neg.short_out")
system("/home/michal/signalp-4.1/signalp -t euk -f short plas_pos.fasta > plas_pos.short_out")
predSignal <- rbind(read_signalp41("plas_pos.short_out"), 
                    read_signalp41("plas_neg.short_out"))

signal.hsmm1987_preds <- pred2df(predict(signal.hsmm1987, c(plas_pos, plas_neg)))


metrics <- do.call(rbind, lapply(list(signalP = predSignal, signalHsmm = signal.hsmm1987_preds), function(predictor)
  HMeasure(c(rep(1, length(plas_pos)), rep(0, length(plas_neg))), predictor[["sp.probability"]])[["metrics"]]))
#wow!