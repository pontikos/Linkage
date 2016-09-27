#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(parallel, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
#suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(optparse, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))

#option_list <- list( 
#make_option(c("--chr"), help = "chromosome")
#)
#OptionParser(option_list=option_list) -> option.parser
#parse_args(option.parser) -> opt

# Find the right marker for Merlin is key
# keep SNPs seen in all individuals
# and that recombine 


m <- as.data.frame(fread('EVA18942.map',header=FALSE))
colnames(m) <- c('chr','snp','cm', 'pos')



# b) only SNPs with rs numbers


# d) snps at an interval of 0.1cm


# a) only keep polymorphic SNPs with no missingness
d <- read('t-EVA18942.traw')
rownames(d) <- d$SNP
print(length(i <- which(d$CHR %in% seq(1,22))))
d <- d[i,] 
X <- d[,grep('jelinkovi',colnames(d))]
x <- t(apply(X,1,function(x) { return( c('miss'=length(which(is.na(x))), 'wt'=length(which(x==0)), 'het'=length(which(x==1)), 'hom'=length(which(x==2))) )
}))
rownames(x) <- d$SNP 
x <- x[which(x[,'miss'] == 0),]
i <- which( x[,'het']>=2 & x[,'wt']>x[,'het'] & x[,'het']>x[,'hom'] )
d <- d[rownames(x[i,]),]
# exclude SNPs which do not recombine
d <- d[!duplicated( paste(d[,'CHR'],d[,'(C)M'],sep='-') ),]

# b) Mendelian inconsistent markers were discarded
# exclude SNPs with ME problems
ped <- read.csv('genotype-ped.csv')
ped$ID <- paste(ped$Family,ped$genotyped,sep='_')
ped$Father <- paste(ped$Family,ped$Father,sep='_')
ped$Mother <- paste(ped$Family,ped$Mother,sep='_')
rownames(ped) <- ped$ID
X <- d[,grep('jelinkovi',colnames(d))]
# parallel version
ok <- mclapply(1:nrow(X), function(i) {
    x <- X[i,]
    for (n in names(x)) {
            Mother <- ped[n,'Mother']
            Father <- ped[n,'Father']
            if (Mother %in% names(x) & Father %in% names(x)) {
                if (x[Mother] == 0 & x[Father] == 0  & x[n]!=0) { return(FALSE) }
                else if (x[Mother] == 0 & x[Father] == 0  & x[n]>0) { return(FALSE) }
                else if (x[Mother] == 2 & x[Father] == 2  & x[n]!=2){ return(FALSE) }
            }
    }
    return(TRUE)
}, mc.cores=20)
ok <- as.logical(unlist(ok))
d <- d[ok,]

selected.variants <- rownames(d) 
d <- as.data.frame(fread('clean3.map',header=FALSE))
colnames(d) <- c('chr','snp','cm','pos')
d <- d[grep('^rs', d$snp),]
rownames(d) <- d$snp 
length(selected.variants <- intersect(selected.variants,rownames(d))) 
d <- d[selected.variants,]

print('map files')
snps <- list()
for (chr in seq(1,22)) {
    print(chr)
    print(nrow(x <- d[which(d$chr==chr),]))
    x$cm <- x$cm-min(x$cm)
    print(nrow(x <- x[!duplicated(round(x$cm,0)),]))
    snps[[chr]] <- x
    write.table(snps[[chr]][,c('chr','snp','cm')],file=sprintf('clean_%s.map',chr),row.name=FALSE,quote=FALSE,sep='\t',col.name=FALSE)
}

print('dat files')
for (chr in seq(1,22)) {
    print(chr)
    x <- data.frame(A='M',Analysis=snps[[chr]]$snp)
    print(nrow(x))
    write.table(x,file=sprintf('clean_%s.dat',chr),row.name=FALSE,quote=FALSE,sep='\t')
    write.table(x$Analysis,file=sprintf('clean_%s.txt',chr),row.name=FALSE,quote=FALSE,sep='\t',col.name=FALSE)
}

for (chr in seq(1,22)) {
    print(system(sprintf('plink  --file linkage/final --extract clean_%s.txt --recode --out merlin_%s', chr, chr, chr)))
}

# add WGS data
for (chr in seq(1,22)) {
    print(system(sprintf('Rscript wgs_to_ped.R --sample %s --chr %s','JW15',chr)))
}

# run merlin
for (chr in seq(1,22)) {
cmd <- "merlin -p merlin_%s.ped -d  clean_%s.dat -m clean_%s.map --npl --exp --bits 28 --megabytes 9999 --tabulate --markerNames > chr%s.results"
system(sprintf(cmd,chr,chr,chr,chr))

}




