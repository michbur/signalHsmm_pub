#taxonomy:"Eukaryota [2759]" keyword:"Transit peptide [KW-0809]" keyword:"Chloroplast [KW-0150]" created:[20080000 TO 20140000] AND reviewed:yes
#taxonomy:"Eukaryota [2759]" keyword:"Transit peptide [KW-0809]" keyword:"Chloroplast [KW-0150]" created:[19800000 TO 20070000] AND reviewed:yes

library(signalHsmm)
test_dat <- read_uniprot("C:/Users/Michal/Dropbox/signal-peptide2_data/chloroplast_from2008.txt",
                         euk = TRUE, what = "transit")
train_dat <- read_uniprot("C:/Users/Michal/Dropbox/signal-peptide2_data/chloroplast_to2007.txt",
                         euk = TRUE, what = "transit")

#remove unproperly read seqs
train_dat <- train_dat[sapply(train_dat, function(i) length(attr(i, "sig"))) == 2]
test_dat <- test_dat[sapply(test_dat, function(i) length(attr(i, "sig"))) == 2]



signal_dat <- read_uniprot("C:/Users/Michal/Dropbox/signal-peptide2_data/pub_pos_train.txt",
                          euk = TRUE, what = "signal")

signal_model <- train_hsmm(signal_dat, aaaggregation, max_length = 32)
#fix bug in measure_region
transit_model <- train_hsmm(train_dat, aaaggregation, max_length = 81)

pred <- predict(transit_model, test_dat)
sum(round(pred2df(pred)[, "sp.probability"], 2) > 0.5)/length(test_dat)