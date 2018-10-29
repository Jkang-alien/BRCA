library('Cairo')
library("trackViewer")
library(stringr)

data <- read.delim('BRCA_data.txt')
variant <- read.delim('table1.txt')
variant <- variant[-27,]
variant$Predictied.effect
position <- str_match(variant$Mutation.in.corresponding.cDNA, '[0-9]{1,4}')
position <- unique(floor(as.numeric(position)/3 + 1))
position
sample_BRCA1 <- GRanges("chr17", IRanges(position[seq(1,13)],
                                         width=rep(1, 13) ,
                                         names=variant$Predictied.effect[variant$Gene == 'BRCA1' & (duplicated(variant$Predictied.effect)==FALSE)]))

sample_BRCA1$score <- c(2,3,1,1,1,1,2,2,1,1,1,1,2)
sample_BRCA1$color <- c('white', 'white', 'white', 'white', 'black', 'white', 'white',
                        'white','white', 'white', 'white','white','white')

sample_BRCA2 <- GRanges("chr17", IRanges(position[seq(14,21)],
                        width=rep(1, 8),
                        names=variant$Predictied.effect[variant$Gene == 'BRCA2' & (duplicated(variant$Predictied.effect)==FALSE)]))
sample_BRCA2$score <- c(1,1,1,1,1,1,1,2)
sample_BRCA2$color <- c('black', 'white', 'white', 'white', 'black',
                        'white','white', 'white')
features_BRCA1 <- GRanges("chr17", IRanges(start = c(1, 24, 344, 648, 1662, 1757, 1863),
                                           width=c(0, 41, 164, 331, 62, 86, 0),
                                           names=c('','RING-finger',
                                                   'Serine-rich domain associated with BRCT',
                                                   'Ethylene insensitive 3',
                                                   'BRCA1 C Terminus (BRCT)',
                                                   'BRCA1 C Terminus (BRCT)', '')),
                          height = rep(0.05, 7),
                          fill=c( '#FFFFFF', "#FF8833", "#51C6E6", "#DFA32D", "#E495A5", "#E495A5", '#FFFFFF'))

features_BRCA2 <- GRanges("chr17", IRanges(start = c(1, 1002, 1212, 1421, 1517,
                                                     1696, 1837, 1972, 2051,
                                                     2479, 2670, 2831, 3052, 3418),
                                           width=c(0, 35, 35, 34, 34, 33, 33, 34, 34,
                                                   189, 130, 41, 139, 0),
                                           names=c('','BRCA2 repeat',
                                                   'BRCA2 repeat',
                                                   'BRCA2 repeat',
                                                   'BRCA2 repeat',
                                                   'BRCA2 repeat',
                                                   'BRCA2 repeat',
                                                   'BRCA2 repeat',
                                                   'BRCA2 repeat',
                                                   'BRCA2, helical',
                                                   "BRCA2, oligonucleotide/oligosaccharide-binding, domain 1",
                                                   'Tower',
                                                   'BRCA2, oligonucleotide/oligosaccharide-binding, domain 3',
                                                   '')),
                          height = rep(0.05, 14),
                          fill=c( '#FFFFFF', "#FF8833", "#FF8833","#FF8833","#FF8833","#FF8833","#FF8833","#FF8833","#FF8833",
                                  "#51C6E6", "#DFA32D", "#E495A5", "#BEBADA", '#FFFFFF'))


CairoPDF(file = 'lolliplot_BRCA1.pdf',
         width =7, height = 7, pointsize = 12)
lolliplot(sample_BRCA1, features_BRCA1,
          xaxis = c(1,400, 800, 1200, 1600, 1863),
          ylab = FALSE,
          yaxis=FALSE,
          legend = legend)
dev.off()

CairoPDF(file = 'lolliplot_BRCA2.pdf',
         width =7, height = 7, pointsize = 12)
lolliplot(sample_BRCA2, features_BRCA2,
          xaxis = c(1,500, 1000, 1500, 2000, 2500, 3000, 3418),
          ylab = FALSE,
          yaxis=FALSE,
          legend = legend)
dev.off()

CairoPNG(file = 'lolliplot_BRCA1.png',
         width =1400, height = 1400, pointsize = 36)
lolliplot(sample_BRCA1, features_BRCA1,
          xaxis = c(1,400, 800, 1200, 1600, 1863),
          ylab = FALSE,
          yaxis=FALSE,
          legend = legend)
dev.off()

CairoPNG(file = 'lolliplot_BRCA2.png',
         width =1400, height = 1400, pointsize = 36)
lolliplot(sample_BRCA2, features_BRCA2,
          xaxis = c(1,500, 1000, 1500, 2000, 2500, 3000, 3418),
          ylab = FALSE,
          yaxis=FALSE,
          legend = legend)
dev.off()

data <- data[c(-104, -110), ]
colnames(data)
summary(data$FIGO.stage)
data$Diganosis

#############################################
######## Baysian prevalence #################

theta <- seq(0,0.15, 0.001)


sample_size <- c(47, 103, 316, 235, 46)
somatic <- c(3, 5, 19, 11, 4)
non_somatic <- sample_size - somatic
sample_n <- sum(sample_size)
somatic_n <- sum(somatic)
non_somatic_n <- sum(non_somatic)
sample_size_cmc <- 88
somatic_n_cmc <- 2

CI <- quantile(rbeta(1000, somatic_n + somatic_n_cmc, 
                     non_somatic_n + 86), probs = c(0.025, 0.975))
CI
mean_posterior <- (somatic_n + somatic_n_cmc) / (somatic_n + somatic_n_cmc + non_somatic_n + 86 +1)
mean_posterior
plot(theta, dbeta(theta, somatic_n, non_somatic_n), 
     type = 'l', lty = 2, ylim = c(0,60),
     ylab = 'Probability density',
     xlab = 'Prevalence of somatic pathogenic mutation')
#polygon(CI[1], CI[2], col = 'red')
lines(theta, dbeta(theta, somatic_n + somatic_n_cmc, 
                   non_somatic_n + 86), lty=1)
lines(theta, dbeta(theta, somatic_n_cmc, 
                    86), lty=3)
legend("topright", legend = c("Posterior", "Prior", "Current"),
       lty = 1:3, xjust = 1, yjust = 1)



CI <- quantile(rbeta(1000, somatic_n_cmc, 
                     86), probs = c(0.025,0.5, 0.975))
CI
abline(v = CI[1], lyt = 1)
abline(v = CI[2], lty = 1)

(somatic / sample_size)*100
