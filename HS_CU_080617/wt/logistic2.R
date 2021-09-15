
library(GRanges)
library(segmentSeq)
chrlens <- c(Chr1=30427671,Chr2=19698289,Chr3=23459830,Chr4=18585056,Chr5=26975502)
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")

####################################################################################################################
# 481252 SNPs -> GRanges object 481257 length (+5=chrs?) -> defines intervals in which a crossover can de detected #
####################################################################################################################

snps <- read.table("BC.complete.tiger.txt")
snpGR <- do.call("c",lapply(1:5,function(chr){
    chrsnps <- snps[snps[,1]==chr,]
    GRanges(seqnames=names(chrlens)[chr],IRanges(start=c(1,chrsnps[,2]),end=c(chrsnps[,2],chrlens[chr])))
}))

##############################################################################################
# load 3,320 crossovers -> note Tom reduces to 3,115 ranges -> underestimates hot intervals? #
##############################################################################################

serra.wt <- read.table(file="serra.wt.coller.txt")
underwood.wt <- read.table(file="underwood.wt.coller.txt")
wt.cos <- rbind(serra.wt,underwood.wt)
n.wt <- 96+96+82+89+74
for(i in 1:5){
	wt.cos[which(wt.cos[,3]==i),3] <- chrs[i]
}
wt.cos.gr.coords <- GRanges(seqnames=wt.cos[,3],IRanges(start=wt.cos[,4],end=wt.cos[,5]))
cos.gr.coords <- wt.cos.gr.coords
cos.gr.coords <- sort(cos.gr.coords)
strand(cos.gr.coords) <- "*"
cos.gr.coords <- reduce(cos.gr.coords)

##############################################################################################################
# overlap to identify SNP windows which overlap a crossover, check type="within", overlap=3295 SNP intervals #
##############################################################################################################

coevents <- getOverlaps(snpGR,cos.gr.coords,overlapType="within",whichOverlaps=FALSE)

###########################
# regions? 480348 regions #
###########################

regions <- c(snpGR[!coevents],cos.gr.coords) 
regions <- sort(regions)

#############################################
# genes=27206, tes=29150 -> GRanges objects #
#############################################
###############################################################
# split TEs into families, define gene promoters, terminators #
###############################################################

genes <- read.table("all.genes.rpi.txt",sep="")
tes <- read.table("all.rpinuc.tes.txt",sep="")
genes <- GRanges(seqnames=paste("Chr",genes[,2],sep=""),IRanges(start=genes[,3],end=genes[,4]),strand=c("+", "-")[as.numeric(genes$strand=="REVERSE") + 1])
genes <- genes[order(as.character(seqnames(genes)),start(genes)),]
tes <- GRanges(seqnames=tes[,1],IRanges(start=tes[,2],end=tes[,3]),strand=tes$strand,sf=tes$class)
tes <- tes[order(as.character(seqnames(tes)),start(tes)),]
tesList <- lapply(levels(tes$sf),function(sf){
    tes[values(tes)$sf==sf,]
})
names(tesList) <- levels(tes$sf)

makePromoters <- function(genes){
    promoters <- genes
    end(promoters[strand(promoters)=="+"]) <- start(genes[strand(genes)=="+"])-1
    start(promoters[strand(promoters)=="+"]) <- start(genes[strand(genes)=="+"])-500
    start(promoters[strand(promoters)=="-"]) <- end(genes[strand(genes)=="-"])+1
    end(promoters[strand(promoters)=="-"]) <- end(genes[strand(genes)=="-"])+500
    start(promoters)[(start(promoters)<1)] <- 1
    promoters
}
promoters <- makePromoters(genes)

makeTerminators <- function(genes){
    terminators <- genes
    start(terminators[strand(terminators)=="+"]) <- end(genes[strand(genes)=="+"])-1
    end(terminators[strand(terminators)=="+"]) <- end(genes[strand(genes)=="+"])+500
    start(terminators[strand(terminators)=="-"]) <- start(genes[strand(genes)=="-"])-500
    end(terminators[strand(terminators)=="-"]) <- start(genes[strand(genes)=="-"])+1
    start(terminators)[(start(terminators)<1)] <- 1
    terminators
}
terminators <- makeTerminators(genes)

########################################################################
# overlap gene, TE, promoter annotations with regions (~SNP intervals) #
########################################################################

geneOverlap <- getOverlaps(regions,genes,whichOverlaps=FALSE)
promoterOverlap <- getOverlaps(regions,promoters,whichOverlaps=FALSE)
terminatorOverlap <- getOverlaps(regions,terminators,whichOverlaps=FALSE)
tesOverlap <- do.call("cbind",lapply(tesList,function(x){
    getOverlaps(regions,x,whichOverlaps=FALSE)
}))

##############################################################################################
# load SPO11-1 data, split according to region, calculate mean, max, min, sum in each region #
##############################################################################################
####################
# add in NUCS data #
####################

splitspo <- list()
for(spochr in dir(pattern="col_chr[1-5]_cov.norm.both.txt")){    
    spo11 <- read.table(spochr)[,1]
    chreg <- regions[seqnames(regions)==seqlevels(regions)[as.numeric(gsub("col_chr([1-5])_cov.norm.both.txt","\\1",spochr))]]
    splitspo <- c(splitspo,split(spo11,c(rep(1:length(chreg),width(chreg)-1),length(chreg))))
}
meanspo <- sapply(splitspo,mean)
maxspo <- sapply(splitspo,max)
minspo <- sapply(splitspo,min)
sumspo <- sapply(splitspo,sum)

splitnuc <- list()
for(i in 1:5){
    nuc <- read.table(file=paste("WT_nuc_chr",i,"_cov.norm.both.txt",sep=""))[,1]
    chreg <- regions[seqnames(regions)==chrs[i]]
    splitnuc <- c(splitnuc,split(nuc,c(rep(1:length(chreg),width(chreg)-1),length(chreg))))
}
meannuc <- sapply(splitnuc,mean)
maxnuc <- sapply(splitnuc,max)
minnuc <- sapply(splitnuc,min)
sumnuc <- sapply(splitnuc,sum)

###################################################################################################################
# define centromeres and banding regions - defines centromeres as 2Mb around CEN, bands are 1Mb intervals in arms #
###################################################################################################################

centromeres <- read.delim("centromere.txt",header=FALSE)
centGR <- GRanges(seqnames=paste("Chr",centromeres[,1],sep=""),IRanges(start=centromeres[,2]-2e6,end=centromeres[,3]+2e6))
centromere <- getOverlaps(regions,centGR,whichOverlaps=FALSE)
pericentGR <- sort(c(GRanges(seqnames(centGR),IRanges(pmax(1,start(centGR)-1-2e6),end=start(centGR)+1)),GRanges(seqnames(centGR),IRanges(end(centGR)+1,end=end(centGR)+1+2e6))))
pericent <- getOverlaps(regions,pericentGR,whichOverlaps=FALSE)
banding <- do.call("c",lapply(seqlevels(regions),function(chr){
    chrlen <- max(end(regions[seqnames(regions)==chr]))
    cent <- centGR[which(seqnames(centGR)==chr)]
    divL <- round((start(cent)-1)/2e6)
    divR <- round((chrlen-end(cent)+1)/2e6)
    ends <- c(seq.int(divL)*(start(cent)-1)/divL,end(cent),seq.int(divR)*(chrlen-end(cent)+1)/divR+end(cent))    
    GRanges(seqnames=chr,IRanges(start=c(1,ends[-length(ends)]+1),end=ends))
}))
bandRegions <- sapply(getOverlaps(regions,banding,whichOverlaps=TRUE), function(x) x[1])

##########################################################################################
# create data object for model -> contains CO overlap, SPO11-1 data, annotation overlaps #
##########################################################################################

dat <- cbind.data.frame(co = regions %in% cos.gr.coords,
                        meanspo=meanspo*1e4,maxspo=maxspo,minspo=minspo,sumspo=sumspo*1e4,
                        gene=geneOverlap,promoter=promoterOverlap,terminator=terminatorOverlap,tesOverlap,
                        band=as.factor(bandRegions),
                        width=width(regions),
                        centromere=centromere,
                        pericent=pericent)
colnames(dat) <- gsub("/", "_", colnames(dat))
colnames(dat) <- gsub("-", ".", colnames(dat))
dat$centromere[dat$centromere] <- as.character(seqnames(regions[centromere]))
dat$centromere[dat$centromere==FALSE] <- "Chr0"
dat$pericent[dat$pericent] <- as.character(seqnames(regions[pericent]))
dat$pericent[dat$pericent==FALSE] <- "Chr0"

#############################################################
# run binomial GLM model, link function is logit, plot data #
#############################################################
# predict function produces predicted values 
# Tom uses sumspo - why not meanspo - because width is a variable?

glmCO <- glm(co~band+sumspo+sumnuc*(gene+promoter+terminator+heli+ptmari+mudr+harbinger+hat+copia+enspm+sine+linel1+gypsy
+width),family=binomial(link="logit"),data=dat)
predict <- predict(glmCO,type="response")
sum.glm <- summary(glmCO)
coeff <- sum.glm$coefficients
write.csv(coeff,file="GLM.nucs.spo.summary.csv") 

glmCO <- glm(co~band+sumspo*(gene+promoter+terminator+heli+ptmari+mudr+harbinger+hat+copia+enspm+sine+linel1+gypsy
+width),family=binomial(link="logit"),data=dat)
predict <- predict(glmCO,type="response")
sum.glm <- summary(glmCO)
coeff <- sum.glm$coefficients
write.csv(coeff,file="GLM.summary.csv") 

###################################################
# plot predicted CO overlaps for annotation class #
###################################################

regid <- c(list(geneOverlap,promoterOverlap,terminatorOverlap),lapply(colnames(tesOverlap),function(tesn){tesOverlap[,tesn]}))
names(regid) <- c("gene","promoter","terminator",colnames(tesOverlap))
boxplot(lapply(regid,function(x) sapply(1:100,function(ii) sum(rbinom(sum(x),1,prob=predict[x]))/sum(dat$width[x])*1e6)),main="predicted CO/Mb by annotation")
tco <- sapply(regid,function(x) sum(dat$co[x])/sum(dat$width[x])*1e6)
points(x=1:length(regid),y=tco,col="red",pch = 19)

###################################################
# plot predicted CO for regions by SPO11-1 hexile #
###################################################

ssID <- split(1:nrow(dat),cut(dat$meanspo,breaks=sort(dat$meanspo)[round(0:7*(nrow(dat)/7))]))
boxplot(lapply(ssID, function(x) sapply(1:100, function(ii) sum(rbinom(length(x),1,prob=predict[x]))/sum(dat$width[x])*1e6)),main="predicted CO/Mb by SPO11-1 hexile")
tco <- sapply(ssID, function(x) sum(dat$co[x])/sum(dat$width[x])*1e6)
points(x=1:length(ssID),y=tco,col="red",pch=19)

##########################################
# plot predicted CO for chromosome bands #
##########################################

boxplot(lapply(split(1:nrow(dat),dat$band),function(x) sapply(1:100,function(ii) sum(rbinom(length(x),1,prob=predict[x]))/sum(dat$width[x])*1e6)),main="predicted CO/Mb by chromosome banding")
ssb <- split(1:nrow(dat),dat$band)
points(x=unique(dat$band),y=sapply(ssb,function(x) sum(dat$co[x])/sum(dat$width[x])*1e6),col="red",pch=19)

#########################################
# plot predicted CO for each SNP region #
#########################################

xy <- cbind(x=cumsum(width(regions)),y=(predict))
plot(xy[xy[,2]>0.01,],pch=".",main="Likelihood of CO at each SNP window");
abline(v=cumsum(chrlens),col="red",lty=3)
abline(v=c(0,cumsum(chrlens)[-5])+start(centGR),col="blue",lty=3)
abline(v=c(0,cumsum(chrlens)[-5])+end(centGR),col="blue",lty=3)
