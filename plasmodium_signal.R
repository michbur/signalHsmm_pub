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

if(Sys.info()["nodename"] == "MICHALKOMP" )
  pathway <- "C:/Users/Michal/Dropbox/signal-peptide2_data/"


library(signalHsmm)
library(seqinr)
library(hmeasure)

# signal.plas <- train_hsmm(read_uniprot(paste0(pathway, "plasmodium_big.txt"), euk = TRUE),
#                           aaaggregation)

#train algotrithm
signal.hsmm1987 <- train_hsmm(read_uniprot(paste0(pathway, "pos_ultrahard_data.txt"), euk = TRUE),
                              aaaggregation)

plas_pos <- read_uniprot("plasmodium_pos.txt", euk = TRUE, what = "signal")
plas_neg <- read.fasta("plas_neg.fasta", seqtype = "AA")

short_plas_neg <- lapply(plas_neg, function(i)
  i[1L:ifelse(length(i) > 500, 500, length(i))])

write.fasta(c(plas_pos, short_plas_neg),
            c(names(plas_pos), names(plas_neg)),
            file = "plas_total.fasta")

write.fasta(plas_pos, names(plas_pos), file = "plasmodium_pos_filtered.fasta")

system("/home/michal/signalp-4.1/signalp -t euk -f short plas_neg.fasta > plas_neg.short_out")
system("/home/michal/signalp-4.1/signalp -t euk -f short plasmodium_pos_filtered.fasta > plas_pos.short_out")
predSignal <- rbind(read_signalp41("plas_pos.short_out"), 
                    read_signalp41("plas_neg.short_out"))

signal.hsmm1987_preds <- pred2df(predict(signal.hsmm1987, c(plas_pos, plas_neg)))
signal.hsmm_preds <- pred2df(run_signalHsmm(c(plas_pos, plas_neg)))

real_labels <- c(rep(1, length(plas_pos)), rep(0, length(plas_neg)))

metrics <- do.call(rbind, lapply(list(signalP = predSignal, signalHsmm = signal.hsmm1987_preds,
                                      signalHsmm1987 = signal.hsmm1987_preds), function(predictor)
                                        HMeasure(real_labels, predictor[["sp.probability"]])[["metrics"]]))

plas_bench <- c(plas_pos, plas_neg)
save(plas_bench, real_labels, file = "plasmodium_benchmark.RData")


#wow!

#library(pROC)
#auc(real_labels, signal.hsmm_preds[["sp.probability"]])
HMeasure(real_labels, signal.hsmm_preds[["sp.probability"]])[["metrics"]]

signalHSMM_proteins <- plas_pos[(signal.hsmm_preds[["sp.probability"]] > 0.5 & !predSignal[["sp.probability"]])[as.logical(real_labels)]]
regions <- sapply(signalHSMM_proteins,
                  function(single_protein) find_nhc(single_protein, attr(single_protein, "sig")))
save(signalHSMM_proteins, regions, file = "dane_dla_Piotra.RData")
