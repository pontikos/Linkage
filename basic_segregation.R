cases <- read.table('cases.txt',header=FALSE)[,1]
controls <- read.table('controls.txt',header=FALSE)[,1]

d <- read('t-EVA18942.traw')

ca <- d[,cases]
co <- d[,controls]

ca <- ca>0
co <- co>0

n.cases <- length(cases)
n.controls <- length(controls)

i <- which( rowSums(ca)==n.cases & rowSums(co)==0 )
d[i,c('CHR','POS', 'SNP',cases)]

i <- which( rowSums(ca)==0 & rowSums(co)==n.controls )
d[i,c('CHR','POS', 'SNP',controls)]
