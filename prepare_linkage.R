#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(kinship2))

option_list <- list(
    make_option(c('--chr'), default='22', help=''),
    make_option(c('--linkage.markers'),  help='Markers to use for linkage.  Finding the right marker for Merlin is key. Ideally keep SNPs genotyped in all individuals and that recombine.'),
    make_options(c('--base.dir'),  help='Base dir'),
    make_option(c('--trim'), default='', help='Comma separated list of indviduals to remove from pedigree. Descendants and spouses are also removed.'),
    make_option(c('--skip'), default='', help='Comma separated list of indviduals to skip from pedigree. Their spouse is removed and their descendants are directly linked to their parents.')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

chr <- opt$chr
trim <- unlist(strsplit(opt$trim,','))
skip <- unlist(strsplit(opt$skip,','))
linkage.markers.file <- opt$linkage.markers
base.dir <- opt$base.dir


#output files for merlin
#merlin_%s.ped
#merlin_%s.dat
#merlin_%s.map

# Actually here just use the exome markers, let's see if they are any good
#linkage.markers <- file.path(base.dir, 'chip.exome.snps.csv')
linkage.markers <- read(linkage.markers.file)
colnames(linkage.markers) <- c('chr','snp','cm', 'pos')


# markers
# only keep linkage markers in chr
linkage.markers <- linkage.markers[which(linkage.markers$chr==chr),]
(nrow(map <- linkage.markers))
rownames(map) <- map[,'snp']
#x$cm <- x$cm-min(x$cm)
#print(nrow(x <- x[!duplicated(round(x$cm,0)),]))

# remove individuals which are trimmed/skipped from pedigree
pedigree <- read.table(file.path(base.dir,'pedigree/seq-ped.csv'),header=FALSE)
colnames(pedigree) <- c('family','ID','Father','Mother','Gender','Affection')
rownames(pedigree) <- pedigree$ID
# individuals to trim from pedigree
if (length(trim)>0) {
   remove.ids <- c()
   for (s in trim) {
   remove.ids <- unique(c(remove.ids, get.subfamily(s,pedigree)$ID))
   }
   pedigree <- pedigree[-match(remove.ids,pedigree$ID),]
}
# individual to skip from pedigreee
if (length(skip)>0) {
   for (s in skip) {
    pedigree <- pedigree.skip(s,pedigree)
   }
}

#
print('write plink extract file')
extract.file <- file.path(base.dir,'genotypes',sprintf('clean_%s.txt',chr))
print(extract.file)
write.table(linkage.markers$snp,file=extract.file,row.name=FALSE,quote=FALSE,col.name=FALSE)

print('extract markers for genotypes')
print(system(sprintf('plink  --file %s --extract %s --recode --out %s', file.path(base.dir,'genotypes','final'), extract.file, file.path(base.dir,'genotypes',sprintf('merlin_%s',chr)))))

# some markers might have been lost on extraction, need to rewrite file
print('map file')
map.file <- file.path(base.dir,'genotypes',sprintf('merlin_%s.map',chr))
print(map.file)
m <- read.table(map.file,header=FALSE)
colnames(m) <- c('chr','snp','pos','cm')
print(dim(map))
print(dim(map <- map[m$snp,]))
write.table(map[,c('chr','snp','cm')],file=map.file,row.name=FALSE,quote=FALSE,sep='\t',col.name=FALSE)

print('merlin dat file')
dat.file <- file.path(base.dir,'genotypes',sprintf('merlin_%s.dat',chr))
print(dat.file)
dat <- data.frame(A='M',Analysis=map$snp)
write.table(dat,file=dat.file,row.name=FALSE,quote=FALSE,sep='\t')

# read genotype ped file
ped.file <- file.path(base.dir, 'genotypes', sprintf('merlin_%s.ped',chr))
all.ped <- read.table(ped.file,header=FALSE)
# set colnames
colnames(all.ped)[1:6] <- c('fam','id','dad','mum','sex','affection')
for (i in seq(1,2*(nrow(map)),2)) { colnames(all.ped)[(i+6):(i+7)] <- paste(map[ceiling(i/2),2],1:2,sep='.') }
rownames(all.ped) <- all.ped$id

print('intersect')
print(pedigree$ID)
print(intersect(pedigree$ID,all.ped$id))
all.ped <- all.ped[intersect(pedigree$ID,all.ped$id),]

genome.samples <- c( 'JW11', 'JW15', 'JW3', 'JW5')
print(genome.samples <- intersect(genome.samples,pedigree$ID))
# genomes
for (sample in genome.samples) {
    print(sample)
    anno <- read(sprintf('/cluster/project8/vyp/pontikos/People/PetraLiskova/all/genomes/%s/%s_VEP_%s-annotations.csv',sample,sample,chr))
    #existing_variation <- strsplit( anno$Existing_variation, '&' )
    x <- as.data.frame(do.call('rbind', strsplit(anno$VARIANT_ID,'_')))
    colnames(x) <- c('chrom','pos','ref','alt')
    x$pos <- as.numeric(x$pos)
    anno <- cbind(anno,x)
    geno <- read(file.path(base.dir,sprintf('/all/genomes/%s/%s_VEP_%s-genotypes.csv',sample,sample,chr)))
    geno <- as.numeric(geno[,2])
    ped <- rep(0,ncol(all.ped))
    names(ped) <- colnames(all.ped)
    ped[2] <- sample
    for (i in seq(1,2*(nrow(map)),2)) {
        names(ped)[(i+6):(i+7)] <- paste(map[ceiling(i/2),2],1:2,sep='.')
        snp <- map[ceiling(i/2),2]
        j <- grep(sprintf("\\b%s\\b",snp),anno$Existing_variation)
        #print(snp)
        #print(length(j))
        if (length(j)>0) {
            if (is.na(geno[j])) {
                ped[6+i] <- 0
                ped[7+i] <- 0
            } else if (geno[j]==1) {
                ped[6+i] <- anno[j,'ref']
                ped[7+i] <- anno[j,'alt']
            } else  if (geno[j]==2) {
                ped[6+i] <- anno[j,'alt']
                ped[7+i] <- anno[j,'alt']
            } else  if (geno[j]==0) {
                ped[6+i] <- anno[j,'ref']
                ped[7+i] <- anno[j,'ref']
            }
         } else {
                ped[6+i] <- anno[i,'ref']
                ped[7+i] <- anno[i,'ref']
         }
    }
    all.ped <- rbind(all.ped,ped)
}

# now match up the alleles to the genotypes for consistency
genotype.samples <- grep('^P|^J',all.ped[,'id'],value=TRUE,invert=TRUE)
consistent.alleles <- function(all.ped, samples) {
    if (length(samples)==0) {
        return(all.ped)
    }
    for (j in 7:ncol(all.ped)) {
        old.alleles <- names(sort(table(all.ped[samples,j]),decreasing=TRUE))
        old.major <- old.alleles[[1]]
        if (length(old.alleles)==2) {
            old.minor <- old.alleles[[2]]
        } else {
            old.minor <- NULL
        }
        x <- all.ped[genotype.samples,j]
        x <- x[x!=0]
        new.alleles <- names(sort(table(x),decreasing=TRUE))
        new.major <- new.alleles[[1]]
        if (length(new.alleles)==2) {
            new.minor <- new.alleles[[2]]
        } else {
            new.minor <- NULL
        }
        if (!is.null(old.major) && !is.null(new.major) && old.major!=new.major && old.major!='0' && new.major!='0') {
            cat(colnames(all.ped)[j],'major :', old.major, '->', new.major, '\n')
            all.ped[which(all.ped[,j]==old.major),j] <- new.major
        }
        if (!is.null(old.minor) && !is.null(new.minor) && old.minor!=new.minor && old.minor!='0' && old.minor!='0') {
            cat(colnames(all.ped)[j],'minor :', old.minor, '->', new.minor, '\n')
            all.ped[which(all.ped[,j]==old.minor),j] <- new.minor
        }
    }
    return(all.ped)
}
all.ped <- consistent.alleles(all.ped,genome.samples)

# map file
print(map.file <- file.path(sprintf('merlin_%s.map',chr)))
write.table(map[,c('chr','snp','cm')],file=map.file,row.name=FALSE,quote=FALSE,sep='\t',col.name=FALSE)
# dat file
dat.file <- file.path(sprintf('merlin_%s.dat',chr))
write.table(dat,file=dat.file,row.name=FALSE,quote=FALSE,sep='\t')
# ped file
#final.ped.file <- file.path(base.dir,sprintf('merlin_%s.ped',chr))
final.ped.file <- file.path(sprintf('merlin_%s.ped',chr))
write.table(all.ped,file=final.ped.file,col.names=FALSE,row.names=FALSE,quote=FALSE)
write.csv(all.ped,file=gsub('.ped','.csv', final.ped.file),quote=FALSE,row.names=FALSE)

# exomes
exome.samples <- c("J2", "J3", "J5", "J6", "J8", "J9", "J10", "J11", "J12", "J13", "J14")
print(exome.samples <- intersect(exome.samples,pedigree$ID))
if (length(exome.samples)>0) {
anno <- read(file.path(base.dir,sprintf('exomes/VEP_%s-annotations.csv',chr)))
#existing_variation <- strsplit( anno$Existing_variation, '&' )
x <- as.data.frame(do.call('rbind', strsplit(anno$VARIANT_ID,'_')))
colnames(x) <- c('chrom','pos','ref','alt')
x$pos <- as.numeric(x$pos)
anno <- cbind(anno,x)
geno <- read(file.path(base.dir,sprintf('exomes/VEP_%s-genotypes.csv',chr)))
for (sample in exome.samples) {
    print(sample)
    ped <- rep(0,ncol(all.ped))
    names(ped) <- colnames(all.ped)
    ped[2] <- sample
    for (i in seq(1,2*(nrow(map)),2)) {
        names(ped)[(i+6):(i+7)] <- paste(map[ceiling(i/2),2],1:2,sep='.')
        snp <- map[ceiling(i/2),2]
        j <- grep(sprintf("\\b%s\\b",snp),anno$Existing_variation)
        if (length(j)>0) {
            if (is.na(geno[j,sample])) {
                ped[6+i] <- 0
                ped[7+i] <- 0
            } else if (geno[j,sample]==1) {
                ped[6+i] <- anno[j,'ref']
                ped[7+i] <- anno[j,'alt']
            } else  if (geno[j,sample]==2) {
                ped[6+i] <- anno[j,'alt']
                ped[7+i] <- anno[j,'alt']
            } else  if (geno[j,sample]==0) {
                ped[6+i] <- anno[j,'ref']
                ped[7+i] <- anno[j,'ref']
            }
         } else {
                print(snp)
                print(j)
                ped[6+i] <- anno[j,'ref']
                ped[7+i] <- anno[j,'ref']
         }
    }
    all.ped <- rbind(all.ped,ped)
}
all.ped <- consistent.alleles(all.ped,exome.samples)
}

if (FALSE) {
# merge samples analysed with both WES and WGS
wes.wgs <- data.frame( wes=c('J3','J5', 'J11'), wgs=c('JW3','JW5', 'JW11') )
for (i in 1:nrow(wes.wgs)) {
    wes.name <- wes.wgs[i,'wes']
    wgs.name <- wes.wgs[i,'wgs']
    wes <- all.ped[wes.name,]
    wgs <- all.ped[wgs.name,]
    for (i in 6:length(wgs)) {
        if ( wgs[i] == 0 || nchar(wgs[i])>1) {
            all.ped[wgs.name,i] <- wes[i]
        }
        if (wes[[i]]!=0 && wgs[i]!=wes[i]) {
            cat(names(wgs)[i],wgs[[i]],wes[[i]],'\n')
            # pick the most likely
            tab <- table(all.ped[,i])
            all.ped[wgs.name,i] <- names(which.max(tab[c(wgs[[i]],wes[[i]])]))
            #print(names(which.max(tab[c(wgs[[i]],wes[[i]])])))
        }
    }
}
all.ped <- consistent.alleles(all.ped,c(genome.samples,exome.samples))
}

# make sure that there is no more than 2 alleles per column
print(apply(all.ped[,7:ncol(all.ped)],2,function(x) sort(table(x))))
rownames(all.ped) <- all.ped$id

#
all.ped[,1:6] <- pedigree[all.ped$id,1:6]

# add connecting individuals missing from pedigree
print(missing.ids <- setdiff(pedigree$ID,all.ped$id))
print(dim(miss.ped <- matrix(0,nrow=length(missing.ids),ncol=ncol(all.ped))))
colnames(miss.ped) <- colnames(all.ped)
rownames(miss.ped) <- missing.ids
rownames(pedigree) <- pedigree$ID
print(dim(miss.ped))
print(pedigree[missing.ids,1:6])
print(miss.ped[missing.ids,1:6])
#miss.ped[missing.ids,1:6] <- pedigree[missing.ids,1:6]
miss.ped <- cbind(pedigree[missing.ids,1:6], miss.ped[missing.ids,-(1:6)])
colnames(miss.ped) <- colnames(all.ped)
print(miss.ped[missing.ids,1:6])
print(dim(all.ped))
print(dim(miss.ped))
all.ped <- rbind(all.ped,miss.ped)

# map file
print(map.file <- file.path(sprintf('merlin_%s.map',chr)))
write.table(map[,c('chr','snp','cm')],file=map.file,row.name=FALSE,quote=FALSE,sep='\t',col.name=FALSE)
# dat file
dat.file <- file.path(sprintf('merlin_%s.dat',chr))
write.table(dat,file=dat.file,row.name=FALSE,quote=FALSE,sep='\t')
# ped file
#final.ped.file <- file.path(base.dir,sprintf('merlin_%s.ped',chr))
final.ped.file <- file.path(sprintf('merlin_%s.ped',chr))
write.table(all.ped,file=final.ped.file,col.names=FALSE,row.names=FALSE,quote=FALSE)
write.csv(all.ped,file=gsub('.ped','.csv', final.ped.file),quote=FALSE,row.names=FALSE)

# run merlin
cmd <- "merlin -p %s -d  %s -m %s --npl --exp --bits 45 --megabytes 90000 --tabulate --markerNames > chr%s.results"
print(cmd <- sprintf(cmd,final.ped.file,dat.file,map.file,chr))
#print(system(cmd))
