library('Cairo')
library("trackViewer")
variant <- read.delim('table1.txt')
sample_BRCA1 <- GRanges("chr17", IRanges(c(1348, 1777, 1833, 265, 1210, 812, 661, 432, 1677, 1148, 1691, 308, 494, 265, 1210),
                                         width=rep(1, 15) ,
                                         names=variant$Predictied.effect[variant$Predictied.effect!=''&variant$Gene == 'BRCA1'][c(-3, -12, -15, -16)]))
sample_BRCA1$score <- c(1,1,2,1,1,1,1,1,1,2,3,1,1,1,1)
sample_BRCA1$color <- c('white', 'white', 'white', 'white', 'white', 'white', 'black',
                        'white','white', 'white', 'white','white','white','white','white')
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
sample_BRCA2 <- GRanges("chr17", IRanges(c(1248, 1200, 2494, 2722, 42, 1609, 2241, 467),
                                         width=rep(1, 8) ,
                                         names=variant$Predictied.effect[variant$Predictied.effect!=''&variant$Gene == 'BRCA2'][-8]))
sample_BRCA2$score <- c(1,1,1,1,1,1,2,1)
sample_BRCA2$color <- c('white', 'white', 'white', 'white', 'black',
                        'white','white', 'white')


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


coln