require(ggplot2)

k_accuracies <- read.delim("~/Documents/accuracies/k_accuracies.txt", header=FALSE)
colnames(k_accuracies) = c("accuracy", "n_folds")
#ggplot(k_accuracies, aes(x=factor(n_folds), y=accuracy)) + geom_point() + geom_jitter(width = 0.1, height = 0.1)
#ggplot(k_accuracies, aes(x=factor(n_folds), y=accuracy)) + geom_point() + geom_jitter(width = 0.1, height = 0.1)

ggplot(kmer_accuracies, aes(x=factor(n_folds), y=accuracy)) + 
  geom_boxplot() + geom_point() + 
  geom_jitter(width = 0.1, height = 0) +
  theme_bw(base_size=18) +
  xlab("N-folds cross-validation") +
  ylim(0,1)


kmer_accuracies <- read.delim("~/Documents/accuracies/kmer_accuracies.txt", header=FALSE)
colnames(kmer_accuracies) = c("accuracy", "n_folds")

ggplot(kmer_accuracies, aes(x=factor(n_folds), y=accuracy)) + 
  geom_boxplot() + geom_point() + 
  geom_jitter(width = 0.1, height = 0) +
  theme_bw(base_size=18) +
  xlab("K-mer Length") +
  ylim(0,1)
