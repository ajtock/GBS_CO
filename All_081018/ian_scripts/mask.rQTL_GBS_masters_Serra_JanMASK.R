
####################
# GBS masters data #
####################
###############################################################
# Underwood Col/Ler F2 2016 lane1 ind=96 cos=729 cos.ind=7.59 #
# Underwood Col/Ler F2 2016 lane2 ind=96 cos=751 cos.ind=7.82 #
# Serra Col/Ler F2 2016 lane1 ind=82 cos=611 cos.ind=7.45 #####
# Serra Col/Ler F2 2016 lane2 ind=89 cos=658 cos.ind=7.39 #
# Serra Col/Ler F2 2016 lane3 ind=74 cos=571 cos.ind=7.72 ####
# Serra SA Col/Ler F2 2016 lane1 ind=96 cos=714 cos.ind=7.43 #
# Serra SA Col/Ler F2 2016 lane2 ind=96 cos=723 cos.ind=7.53 #
# Serra JA Col/Ler F2 2016 lane1 ind=85 cos=601 cos.ind=7.07 #
# Serra JA Col/Ler F2 2016 lane2 ind=87 cos=653 cos.ind=7.51 ##########
# Underwood HEI10 Col/Ler F2 2016 lane1 ind=96 cos=1456 cos.ind=15.17 #
# Underwood HEI10 Col/Ler F2 2016 lane2 ind=96 cos=1472 cos.ind=15.33 #
# Underwood cmt3 Col/Ler F2 2016 lane1 ind=96 cos=708 cos.ind=7.34 ####
# Underwood cmt3 Col/Ler F2 2016 lane2 ind=96 cos=702 cos.ind=7.31 #
# Underwood cmt3 Col/Ler F2 2016 lane3 ind=96 cos=702 cos.ind=7.31 #
# Underwood cmt3 Col/Ler F2 2016 lane4 ind=96 cos=691 cos.ind=7.20 ########
# Serra recq4a recq4b Col/Ler F2 2017 lane1 ind=96 cos=2396 cos.ind=24.96 #
# Serra recq4a recq4b Col/Ler F2 2017 lane2 ind=95 cos=2387 cos.ind=25.13 #######
# Serra recq4a recq4b HEI10 Col/Ler F2 2017 lane1 ind=96 cos=3060 cos.ind=31.88 #
# Serra recq4a recq4b HEI10 Col/Ler F2 2017 lane2 ind=96 cos=2888 cos.ind=30.08 ##
# Lambing recq4a recq4b cmt3 Col/Ler F2 2017 lane1 ind=96 cos=2543 cos.ind=26.49 #
# Lambing recq4a recq4b cmt3 Col/Ler F2 2017 lane2 ind=96 cos=2524 cos.ind=26.29 #
# Lambing HEI10 cmt3 Col/Ler F2 2017 lane1 ind=96 cos=1476 cos.ind=15.38 #########
# Lambing HEI10 cmt3 Col/Ler F2 2017 lane2 ind=96 cos=1409 cos.ind=14.68 #
##########################################################################

#######################################################
# sorting for multiple hits over Serra PNAS genotypes #
#######################################################

hs.wt <- read.table(file="serra.wt.coller.txt")
hs.wt <- cbind(hs.wt,rep("hs.wt",times=length(hs.wt[,1])))
colnames(hs.wt) <- c("lane","lib","chr","start","end","cos","width","geno")
cu.hei10 <- read.table(file="underwood.HEI10.coller.txt")
cu.hei10 <- cbind(cu.hei10,rep("cu.hei10",times=length(cu.hei10[,1])))
colnames(cu.hei10) <- c("lane","lib","chr","start","end","cos","width","geno")
hs.recq <- read.table(file="serra.recq4a4b.coller.txt")
hs.recq <- cbind(hs.recq,rep("hs.recq",times=length(hs.recq[,1])))
colnames(hs.recq) <- c("lane","lib","chr","start","end","cos","width","geno")
hs.recqhei10 <- read.table(file="serra.recq4a4bHEI10.coller.txt")
hs.recqhei10 <- cbind(hs.recqhei10,rep("hs.recqhei10",times=length(hs.recqhei10[,1])))
colnames(hs.recqhei10) <- c("lane","lib","chr","start","end","cos","width","geno")
all.cos <- rbind(hs.wt,cu.hei10,hs.recq,hs.recqhei10)
all.chr.cos <- NULL
for(i in 1:5){
	print(i)
	chr.cos <- all.cos[which(all.cos[,3]==i),]	
	chr.cos <- chr.cos[order(chr.cos[,6]),]
	unique <- unique(chr.cos[,6])
	chr.anno <- NULL
	for(k in 1:length(unique)){
		print(k)
		match <- chr.cos[,6] %in% unique[k]
		anno <- rep(length(which(match==T)),times=length(which(match==T)))
		chr.anno <- c(chr.anno,anno)
		}
	chr.cos <- cbind(chr.cos,chr.anno)
	all.chr.cos <- rbind(all.chr.cos,chr.cos)
}
write.table(all.cos,file="all.pnascos.txt")
hist(all.chr.cos[,9],breaks=100,main="histogram of pnas crossover coords")

#######################################################################################################
# Raphael identified 74 apparant crossover likely caused by HEI10 translocation -> mask from analysis #
#######################################################################################################

74 crossovers
mask.coords <- c(20378158,20378485,20381452,160073,160663,161023,161262,161483)
match <- all.cos[,6] %in% mask.coords
mask.cos <- all.cos[-which(match==T),]
write.table(mask.cos,file="all.maskcos.txt")
all.chr.cos <- NULL
for(i in 1:5){
	print(i)
	chr.cos <- mask.cos[which(mask.cos[,3]==i),]	
	chr.cos <- chr.cos[order(chr.cos[,6]),]
	unique <- unique(chr.cos[,6])
	chr.anno <- NULL
	for(k in 1:length(unique)){
		print(k)
		match <- chr.cos[,6] %in% unique[k]
		anno <- rep(length(which(match==T)),times=length(which(match==T)))
		chr.anno <- c(chr.anno,anno)
		}
	chr.cos <- cbind(chr.cos,chr.anno)
	all.chr.cos <- rbind(all.chr.cos,chr.cos)
}
mask.cos <- all.chr.cos
levels <- unique(mask.cos[,9])
levels <- levels[order(levels)]
all.tal <- NULL
for(k in 1:length(levels)){
	print(k)
	tal <- length(unique(all.chr.cos[which(mask.cos[,9]==levels[k]),6]))
	all.tal <- c(all.tal,tal)
}
levels <- cbind(levels,all.tal)
write.csv(levels,file="pnascosmask.coord.tally.csv")

high.cos <- mask.cos[which(mask.cos[,9]>2),]
sort.cos <- high.cos[order(high.cos[,9]),]

##################
# summary tables #
##################

hswt.mask <- mask.cos[which(mask.cos[,8]=="hs.wt"),]
cuhei10.mask <- mask.cos[which(mask.cos[,8]=="cu.hei10"),]
hsrecq.mask <- mask.cos[which(mask.cos[,8]=="hs.recq"),]
hsrecqhei10.mask <- mask.cos[which(mask.cos[,8]=="hs.recqhei10"),]
hs.wt.n <- 245
cu.hei10.n <- 192
hs.recq.n <- 191
hs.recqhei10.n <- 192
samp <- c(hs.wt.n,cu.hei10.n,hs.recq.n,hs.recqhei10.n)
all.chr <- NULL
for(i in 1:5){
	print(i)
	all.chr <- rbind(all.chr,c(i,length(which(hswt.mask[,3]==i)),length(which(cuhei10.mask[,3]==i)),length(which(hsrecq.mask[3]==i)),length(which(hsrecqhei10.mask[,3]==i))))
}
all.chr <- rbind(all.chr,colSums(all.chr))
colnames(all.chr) <- c("chr","hs.wt","cu.hei10","hs.recq","hs.recqhei10")
all.ind <- NULL
for(j in 2:length(all.chr[1,])){
	all.ind <- cbind(all.ind,cos.ind <- all.chr[,j]/samp[j-1])
}
all.ind <- cbind(c(1,2,3,4,5,6),all.ind)
colnames(all.ind) <- c("chr","hs.wt","cu.hei10","hs.recq","hs.recqhei10")
all.chr <- rbind(all.chr,c(0,samp))
write.csv(all.chr,file="gbs.maskcos.summary.csv")
write.csv(all.ind,file="gbs.maskcos.ind.summary.csv")

########################
# chromosome idiograms #
########################
#####################
# histograms per F2 #
#####################

hs.wt.n <- 245
cu.hei10.n <- 192
hs.recq.n <- 191
hs.recqhei10.n <- 192

hswt.mask <- mask.cos[which(mask.cos[,8]=="hs.wt"),]
all.ind <- NULL
for(i in 1:3){
	print(i)
	lane <- hswt.mask[which(hswt.mask[,1]==i),]
	lane.libs <- unique(lane[,2])
	cos.ind <- NULL
	for(k in 1:length(lane.libs)){
		print(k)
		cos.ind <- c(cos.ind,length(which(lane[,2]==lane.libs[k])))
	}
	all.ind <- c(all.ind,cos.ind)
}
hswtmask.ind <- all.ind

cuhei10.mask <- mask.cos[which(mask.cos[,8]=="cu.hei10"),]
all.ind <- NULL
for(i in 1:2){
	print(i)
	lane <- cuhei10.mask[which(cuhei10.mask[,1]==i),]
	lane.libs <- unique(lane[,2])
	cos.ind <- NULL
	for(k in 1:length(lane.libs)){
		print(k)
		cos.ind <- c(cos.ind,length(which(lane[,2]==lane.libs[k])))
	}
	all.ind <- c(all.ind,cos.ind)
}
cuhei10mask.ind <- all.ind

hsrecq.mask <- mask.cos[which(mask.cos[,8]=="hs.recq"),]
all.ind <- NULL
for(i in 1:2){
	print(i)
	lane <- hsrecq.mask[which(hsrecq.mask[,1]==i),]
	lane.libs <- unique(lane[,2])
	cos.ind <- NULL
	for(k in 1:length(lane.libs)){
		print(k)
		cos.ind <- c(cos.ind,length(which(lane[,2]==lane.libs[k])))
	}
	all.ind <- c(all.ind,cos.ind)
}
hsrecqmask.ind <- all.ind

hsrecqhei10.mask <- mask.cos[which(mask.cos[,8]=="hs.recqhei10"),]
all.ind <- NULL
for(i in 1:2){
	print(i)
	lane <- hsrecqhei10.mask[which(hsrecqhei10.mask[,1]==i),]
	lane.libs <- unique(lane[,2])
	cos.ind <- NULL
	for(k in 1:length(lane.libs)){
		print(k)
		cos.ind <- c(cos.ind,length(which(lane[,2]==lane.libs[k])))
	}
	all.ind <- c(all.ind,cos.ind)
}
hsrecqhei10mask.ind <- all.ind

par(mfcol=c(4,1))
par(mar=c(1.8,1.8,1.8,1.8))
hist(hswtmask.ind,xlim=c(0,55),breaks=10000)
abline(v=mean(hswtmask.ind),col=2,lty=2)
hist(cuhei10mask.ind,xlim=c(0,55),breaks=10000)
abline(v=mean(cuhei10mask.ind),col=2,lty=2)
hist(hsrecqmask.ind,xlim=c(0,55),breaks=10000)
abline(v=mean(hsrecqmask.ind),col=2,lty=2)
hist(hsrecqhei10mask.ind,xlim=c(0,55),breaks=10000)
abline(v=mean(hsrecqhei10mask.ind),col=2,lty=2)

################################################################
# plotting crossovers per chromosome against chromosome length #
################################################################

wt.chrs <- c(1.8,1.3,1.4,1.3,1.7)
hei10.chrs <- c(4,2.7,2.8,2.1,3.5)
recq.chrs <- c(6.6,4.2,5.2,3,6.1)
recqhei10.chrs <- c(8.5,5.2,6.2,3.3,7.5)
chr.ends <- c(30427671,19698289,23459830,18585056,26975502)

ylim <- c(0,9.5)
xlim <- c(16000000,32000000)

plot(chr.ends,wt.chrs,col=1,ylim=ylim,xlim=xlim,type="p",pch=16)
par(new=T)
plot(chr.ends,hei10.chrs,col=4,ylim=ylim,xlim=xlim,type="p",pch=16)
par(new=T)
plot(chr.ends,recq.chrs,col="magenta",ylim=ylim,xlim=xlim,type="p",pch=16)
par(new=T)
plot(chr.ends,recqhei10.chrs,col=2,ylim=ylim,xlim=xlim,type="p",pch=16)

########################
# confidence intervals #
########################

t.test(hswtmask.ind)
95 percent confidence interval:
7.225235 7.795173
t.test(cuhei10mask.ind)
95 percent confidence interval:
 14.58816 15.55768
t.test(hsrecqmask.ind)
95 percent confidence interval:
 24.11092 25.97285
t.test(hsrecqhei10mask.ind)
95 percent confidence interval:
 29.80394 31.73772

dat <- hswtmask.ind
error <- qt(0.975,df=length(dat)-1)*sd(dat)/sqrt(length(dat))
lower <- mean(dat)-error
upper <- mean(dat)+error
print(lower)
7.225
print(mean(dat))
7.510
print(upper)
7.795

dat <- cuhei10mask.ind
error <- qt(0.975,df=length(dat)-1)*sd(dat)/sqrt(length(dat))
lower <- mean(dat)-error
upper <- mean(dat)+error
print(lower)
14.588
print(mean(dat))
15.073
print(upper)
15.558

dat <- hsrecqmask.ind
error <- qt(0.975,df=length(dat)-1)*sd(dat)/sqrt(length(dat))
lower <- mean(dat)-error
upper <- mean(dat)+error
print(lower)
24.042
print(mean(dat))
25.042
print(upper)
25.973

dat <- hsrecqhei10mask.ind
error <- qt(0.975,df=length(dat)-1)*sd(dat)/sqrt(length(dat))
lower <- mean(dat)-error
upper <- mean(dat)+error
print(lower)
29.804
print(mean(dat))
30.771
print(upper)
31.7

##########################################################
# parse GBS genotypes along chromosomes into rQTL format #
##########################################################

map.pwd <- "/Users/IanHenderson/Desktop/Lab_records/papers/Serra_PNAS/Fig_1_GBS/masters_mask/"
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chr.ends <- c(30427671,19698289,23459830,18585056,26975502)
all.chr.coords <- NULL
all.chr.lab <- NULL
for(i in 1:5){
	print(i)
	chr.coords <- c(seq(1,chr.ends[i],by=1000000),chr.ends[i])
	all.chr.coords <- c(all.chr.coords,chr.coords)
	chr.lab <- rep(i,length(chr.coords))
	all.chr.lab <- c(all.chr.lab,chr.lab)
}
one.samp <- seq(437,532,by=1)
two.samp <- seq(533,628,by=1)

pop.hei101.states <- NULL
pop.hei102.states <- NULL
pop.recq1.states <- NULL
pop.recqhei101.states <- NULL
pop.recqhei102.states <- NULL
for(k in 1:96){
	print(k)		
	hei101.dat <- read.table(file=paste(map.pwd,"n96.Hei10F2",k,".bchqsnvmask.smooth.co.txt",sep=""),header=F)
	hei102.dat <- read.table(file=paste(map.pwd,"n96.hei10lane2F2",k,".bchqsnvmask.smooth.co.txt",sep=""),header=F)
	recq1.dat <- read.table(file=paste(map.pwd,"n96.HS.Recq.series1",k,".bchqsnvmask.smooth.co.txt",sep=""),header=F)
	recqhei101.dat <- read.table(file=paste(map.pwd,"n96.HS.Recq.HEI10.series1",k,".bchqsnvmask.smooth.co.txt",sep=""),header=F)	
	recqhei102.dat <- read.table(file=paste(map.pwd,"n96.HS.Recq.HEI10.series2",k,".bchqsnvmask.smooth.co.txt",sep=""),header=F)
	print(k)
	all.recq1.states <- NULL
	all.hei101.states <- NULL
	all.hei102.states <- NULL
	all.recqhei101.states <- NULL	
	all.recqhei102.states <- NULL	
	for(i in 1:5){
		print(i)
		chr.coords <- c(seq(1,chr.ends[i],by=1000000),chr.ends[i])
		genorecq1.states <- rep(0,length(chr.coords))
		genohei101.states <- rep(0,length(chr.coords))
		genohei102.states <- rep(0,length(chr.coords))
		genorecqhei101.states <- rep(0,length(chr.coords))
		genorecqhei102.states <- rep(0,length(chr.coords))
		recq1.chr <- recq1.dat[which(recq1.dat[,2]==i),]
		for(m in 1:length(recq1.chr[,1])){
			print(m)
			genorecq1.states[which(chr.coords>=recq1.chr[m,3]&chr.coords<=recq1.chr[m,4])] <- recq1.chr[m,5]
		}		
		hei101.chr <- hei101.dat[which(hei101.dat[,2]==i),]
		for(m in 1:length(hei101.chr[,1])){
			print(m)
			genohei101.states[which(chr.coords>=hei101.chr[m,3]&chr.coords<=hei101.chr[m,4])] <- hei101.chr[m,5]
		}		
		hei102.chr <- hei102.dat[which(hei102.dat[,2]==i),]
		for(m in 1:length(hei102.chr[,1])){
			print(m)
			genohei102.states[which(chr.coords>=hei102.chr[m,3]&chr.coords<=hei102.chr[m,4])] <- hei102.chr[m,5]
		}
		recqhei101.chr <- recqhei101.dat[which(recqhei101.dat[,2]==i),]
		for(m in 1:length(recqhei101.chr[,1])){
			print(m)
			genorecqhei101.states[which(chr.coords>=recqhei101.chr[m,3]&chr.coords<=recqhei101.chr[m,4])] <- recqhei101.chr[m,5]
		}
		recqhei102.chr <- recqhei102.dat[which(recqhei102.dat[,2]==i),]
		for(m in 1:length(recqhei102.chr[,1])){
			print(m)
			genorecqhei102.states[which(chr.coords>=recqhei102.chr[m,3]&chr.coords<=recqhei102.chr[m,4])] <- recqhei102.chr[m,5]
		}
		genorecq1.states[which(genorecq1.states==1)] <- "A"
		genorecq1.states[which(genorecq1.states==2)] <- "H"
		genorecq1.states[which(genorecq1.states==3)] <- "B"
		all.recq1.states <- c(all.recq1.states,genorecq1.states)	
		genohei101.states[which(genohei101.states==1)] <- "A"
		genohei101.states[which(genohei101.states==2)] <- "H"
		genohei101.states[which(genohei101.states==3)] <- "B"
		all.hei101.states <- c(all.hei101.states,genohei101.states)
		genohei102.states[which(genohei102.states==1)] <- "A"
		genohei102.states[which(genohei102.states==2)] <- "H"
		genohei102.states[which(genohei102.states==3)] <- "B"
		all.hei102.states <- c(all.hei102.states,genohei102.states)
		genorecqhei101.states[which(genorecqhei101.states==1)] <- "A"
		genorecqhei101.states[which(genorecqhei101.states==2)] <- "H"
		genorecqhei101.states[which(genorecqhei101.states==3)] <- "B"
		all.recqhei101.states <- c(all.recqhei101.states,genorecqhei101.states)	
		genorecqhei102.states[which(genorecqhei102.states==1)] <- "A"
		genorecqhei102.states[which(genorecqhei102.states==2)] <- "H"
		genorecqhei102.states[which(genorecqhei102.states==3)] <- "B"
		all.recqhei102.states <- c(all.recqhei102.states,genorecqhei102.states)
		}
	pop.hei101.states <- rbind(pop.hei101.states,all.hei101.states)
	pop.hei102.states <- rbind(pop.hei102.states,all.hei102.states)
	pop.recq1.states <- rbind(pop.recq1.states,all.recq1.states)
	pop.recqhei101.states <- rbind(pop.recqhei101.states,all.recqhei101.states)
	pop.recqhei102.states <- rbind(pop.recqhei102.states,all.recqhei102.states)
}
lib.nums <- c(seq(1,33,by=1),seq(35,96,by=1))
pop.recq2.states <- NULL
for(k in 1:length(lib.nums)){
	print(k)
    recq2.dat <- read.table(file=paste(map.pwd,"n96_HS_Recq_series2",lib.nums[k],".bchqsnvmask.smooth.co.txt",sep=""),header=F) 
	all.recq2.states <- NULL
	for(i in 1:5){
		print(i)
		chr.coords <- c(seq(1,chr.ends[i],by=1000000),chr.ends[i])
		genorecq2.states <- rep(0,length(chr.coords))
		recq2.chr <- recq2.dat[which(recq2.dat[,2]==i),]
		for(m in 1:length(recq2.chr[,1])){
			print(m)
			genorecq2.states[which(chr.coords>=recq2.chr[m,3]&chr.coords<=recq2.chr[m,4])] <- recq2.chr[m,5]
		}		
		genorecq2.states[which(genorecq2.states==1)] <- "A"
		genorecq2.states[which(genorecq2.states==2)] <- "H"
		genorecq2.states[which(genorecq2.states==3)] <- "B"
		all.recq2.states <- c(all.recq2.states,genorecq2.states)
	}
	pop.recq2.states <- rbind(pop.recq2.states,all.recq2.states)
}
pop.hei10.dat <- rbind(all.chr.coords,all.chr.lab,pop.hei101.states,pop.hei102.states)
header.col <- c("id","",seq(1,length(pop.hei10.dat[,1])-2,by=1))
pop.hei10.dat <- cbind(header.col,pop.hei10.dat)
pop.hei10.dat <- as.data.frame(pop.hei10.dat,row.names=F,col.names=F)
write.table(pop.hei10.dat,file="hei10.f2mask.csv",row.names=F,col.names=F,sep=",")
pop.recq.dat <- rbind(all.chr.coords,all.chr.lab,pop.recq1.states,pop.recq2.states)
header.col <- c("id","",seq(1,length(pop.recq.dat[,1])-2,by=1))
pop.recq.dat <- cbind(header.col,pop.recq.dat)
pop.recq.dat <- as.data.frame(pop.recq.dat,row.names=F,col.names=F)
write.table(pop.recq.dat,file="recq.f2mask.csv",row.names=F,col.names=F,sep=",")
pop.recqhei10.dat <- rbind(all.chr.coords,all.chr.lab,pop.recqhei101.states,pop.recqhei102.states)
header.col <- c("id","",seq(1,length(pop.recqhei10.dat[,1])-2,by=1))
pop.recqhei10.dat <- cbind(header.col,pop.recqhei10.dat)
pop.recqhei10.dat <- as.data.frame(pop.recqhei10.dat,row.names=F,col.names=F)
write.table(pop.recqhei10.dat,file="recqhei10.f2mask.csv",row.names=F,col.names=F,sep=",")
lib.nums <- c(seq(1,7,by=1),seq(9,10,by=1),seq(13,17,by=1),seq(19,25,by=1),seq(29,30,by=1),seq(32,38,by=1),seq(40,50,by=1),seq(52,56,by=1),seq(58,73,by=1),seq(75,86,by=1),seq(88,92,by=1),seq(94,96,by=1))
pop.hswt1.states <- NULL
for(k in 1:length(lib.nums)){
	print(k)
	hswt1.dat <- read.table(file=paste(map.pwd,"n96.HS1.wtF2",lib.nums[k],".bchqsnvmask.smooth.co.txt",sep=""),header=F)
   	print(k)
	all.hswt1.states <- NULL
	for(i in 1:5){
		print(i)
		chr.coords <- c(seq(1,chr.ends[i],by=1000000),chr.ends[i])
		genohswt1.states <- rep(0,length(chr.coords))
		hswt1.chr <- hswt1.dat[which(hswt1.dat[,2]==i),]
		for(m in 1:length(hswt1.chr[,1])){
			print(m)
			genohswt1.states[which(chr.coords>=hswt1.chr[m,3]&chr.coords<=hswt1.chr[m,4])] <- hswt1.chr[m,5]
		}		
		genohswt1.states[which(genohswt1.states==1)] <- "A"
		genohswt1.states[which(genohswt1.states==2)] <- "H"
		genohswt1.states[which(genohswt1.states==3)] <- "B"
		all.hswt1.states <- c(all.hswt1.states,genohswt1.states)
	}
	pop.hswt1.states <- rbind(pop.hswt1.states,all.hswt1.states)
}
lib.nums <- c(seq(1,2,by=1),seq(4,6,by=1),seq(8,42,by=1),seq(44,55,by=1),seq(57,70,by=1),seq(72,75,by=1),seq(77,91,by=1),seq(93,96,by=1))
pop.hswt2.states <- NULL
for(k in 1:length(lib.nums)){
	print(k)
	hswt2.dat <- read.table(file=paste(map.pwd,"n96.HS2.wtF2",lib.nums[k],".bchqsnvmask.smooth.co.txt",sep=""),header=F)
	all.hswt2.states <- NULL
	for(i in 1:5){
		print(i)
		chr.coords <- c(seq(1,chr.ends[i],by=1000000),chr.ends[i])
		genohswt2.states <- rep(0,length(chr.coords))
		hswt2.chr <- hswt2.dat[which(hswt2.dat[,2]==i),]
		for(m in 1:length(hswt2.chr[,1])){
			print(m)
			genohswt2.states[which(chr.coords>=hswt2.chr[m,3]&chr.coords<=hswt2.chr[m,4])] <- hswt2.chr[m,5]
		}		
		genohswt2.states[which(genohswt2.states==1)] <- "A"
		genohswt2.states[which(genohswt2.states==2)] <- "H"
		genohswt2.states[which(genohswt2.states==3)] <- "B"
		all.hswt2.states <- c(all.hswt2.states,genohswt2.states)
	}
	pop.hswt2.states <- rbind(pop.hswt2.states,all.hswt2.states)
}
lib.nums <- c(seq(1,61,by=1),seq(68,80,by=1))
pop.hswt3.states <- NULL
for(k in 1:length(lib.nums)){
	print(k)
	hswt3.dat <- read.table(file=paste(map.pwd,"n96.HS.WT3F2",lib.nums[k],".bchqsnvmask.smooth.co.txt",sep=""),header=F)
	all.hswt3.states <- NULL
	for(i in 1:5){
		print(i)
		chr.coords <- c(seq(1,chr.ends[i],by=1000000),chr.ends[i])
		genohswt3.states <- rep(0,length(chr.coords))
		hswt3.chr <- hswt3.dat[which(hswt3.dat[,2]==i),]
		for(m in 1:length(hswt3.chr[,1])){
			print(m)
			genohswt3.states[which(chr.coords>=hswt3.chr[m,3]&chr.coords<=hswt3.chr[m,4])] <- hswt3.chr[m,5]
		}		
		genohswt3.states[which(genohswt3.states==1)] <- "A"
		genohswt3.states[which(genohswt3.states==2)] <- "H"
		genohswt3.states[which(genohswt3.states==3)] <- "B"
		all.hswt3.states <- c(all.hswt3.states,genohswt3.states)
	}
	pop.hswt3.states <- rbind(pop.hswt3.states,all.hswt3.states)
}
pop.hswt.dat <- rbind(all.chr.coords,all.chr.lab,pop.hswt1.states,pop.hswt2.states,pop.hswt3.states)
header.col <- c("id","",seq(1,length(pop.hswt.dat[,1])-2,by=1))
pop.hswt.dat <- cbind(header.col,pop.hswt.dat)
pop.hswt.dat <- as.data.frame(pop.hswt.dat,row.names=F,col.names=F)
write.table(pop.hswt.dat,file="hswt.f2mask.csv",row.names=F,col.names=F,sep=",")

###################################
# analyse genetic maps using rQTL #
###################################

library(qtl)
hswt.map <- read.cross(file="hswt.f2mask.csv",genotypes=c("A","H","B"))
hei10.map <- read.cross(file="hei10.f2mask.csv",genotypes=c("A","H","B"))
recq.map <- read.cross(file="recq.f2mask.csv",genotypes=c("A","H","B"))
recqhei10.map <- read.cross(file="recqhei10.f2mask.csv",genotypes=c("A","H","B"))

par(mfcol=c(2,2))
par(mar=c(1.8,1.8,1.8,1.8))
plot.rf(hswt.map,main="wt hs mask rf")
plot.rf(hei10.map,main="HEI10 mask rf")
plot.rf(recq.map,main="recq mask rf")
plot.rf(recqhei10.map,main="recq4a4b HEI10 mask rf")

hs.wt.n <- 245
cu.hei10.n <- 192
hs.recq.n <- 191
hs.recqhei10.n <- 192

hswt.haldane.map <- read.cross(file="hswt.f2mask.csv",genotypes=c("A","H","B"),map.function="haldane")
count.XO <- (sum(colSums(countXO(hswt.haldane.map,bychr=T))))/hs.wt.n
# 7.44
recq.haldane.map <- read.cross(file="recq.f2mask.csv",genotypes=c("A","H","B"),map.function="haldane")
count.XO <- (sum(colSums(countXO(recq.haldane.map,bychr=T))))/hs.recq.n
# 23.16
hei10.haldane.map <- read.cross(file="hei10.f2mask.csv",genotypes=c("A","H","B"),map.function="haldane")
count.XO <- (sum(colSums(countXO(hei10.haldane.map,bychr=T))))/cu.hei10.n
# 14.74
recqhei10.haldane.map <- read.cross(file="recqhei10.f2mask.csv",genotypes=c("A","H","B"),map.function="haldane")
count.XO <- (sum(colSums(countXO(recqhei10.haldane.map,bychr=T))))/hs.recqhei10.n
# 27.92

hswt.kosambi.map <- read.cross(file="hswt.f2mask.csv",genotypes=c("A","H","B"),map.function="kosambi")
count.XO <- (sum(colSums(countXO(hswt.kosambi.map,bychr=T))))/hs.wt.n
# 7.44
recq.kosambi.map <- read.cross(file="recq.f2mask.csv",genotypes=c("A","H","B"),map.function="kosambi")
count.XO <- (sum(colSums(countXO(recq.kosambi.map,bychr=T))))/hs.recq.n
# 23.16

hswt.cf.map <- read.cross(file="hswt.f2mask.csv",genotypes=c("A","H","B"),map.function="c-f")
(sum(colSums(countXO(hswt.cf.map,bychr=T))))/hs.wt.n
# 7.44
recq.cf.map <- read.cross(file="recq.f2mask.csv",genotypes=c("A","H","B"),map.function="c-f")
(sum(colSums(countXO(recq.cf.map,bychr=T))))/hs.recq.n
# 23.16

hswt.morgan.map <- read.cross(file="hswt.f2mask.csv",genotypes=c("A","H","B"),map.function="morgan")
(sum(colSums(countXO(hswt.morgan.map,bychr=T))))/hs.wt.n
# 7.44
recq.morgan.map <- read.cross(file="recq.f2mask.csv",genotypes=c("A","H","B"),map.function="morgan")
(sum(colSums(countXO(recq.morgan.map,bychr=T))))/hs.recq.n
# 23.16

par(mfcol=c(2,4))
par(mar=c(1.8,1.8,1.8,1.8))
plot.rf(hswt.haldane.map,main="hs.wt haldane")
plot.rf(recq.haldane.map,main="hs.recq haldane")
plot.rf(hswt.kosambi.map,main="hs.wt kosambi")
plot.rf(recq.kosambi.map,main="hs.recq kosambi")
plot.rf(hswt.cf.map,main="hs.wt cf")
plot.rf(recq.cf.map,main="hs.recq cf")
plot.rf(hswt.morgan.map,main="hs.wt morgan")
plot.rf(recq.morgan.map,main="hs.recq morgan")

par(mfcol=c(1,4))
par(mar=c(1.8,1.8,1.8,1.8))
plotMap(hswt.haldane.map,recq.haldane.map,main="wt vs recq haldane",ylim=c(400,0))
plotMap(hswt.kosambi.map,recq.kosambi.map,main="wt vs recq kosambi",ylim=c(400,0))
plotMap(hswt.cf.map,recq.cf.map,main="wt vs recq c-f",ylim=c(400,0))
plotMap(hswt.morgan.map,recq.morgan.map,main="wt vs recq morgan",ylim=c(400,0))

par(mfcol=c(1,3))
par(mar=c(1.8,1.8,1.8,1.8))
plotMap(hswt.haldane.map,hei10.haldane.map,main="wt vs recq haldane",ylim=c(550,0))
plotMap(hswt.haldane.map,recq.haldane.map,main="wt vs recq haldane",ylim=c(550,0))
plotMap(hswt.haldane.map,recqhei10.haldane.map,main="wt vs recq haldane",ylim=c(550,0))

wt.geno <- hswt.map$geno
wtchr.len <- NULL
for(i in 1:5){
	chr <- wt.geno[[i]]
	map <- chr$map
	wtchr.len <- c(wtchr.len,map[length(map)])
}
30427671 19698289 23459830 18585056 26975502 
93.73843 71.63953 74.72704 66.97380 90.52054 
wttot.len <- round(sum(wtchr.len))
398

hei10.geno <- hei10.map$geno
hei10chr.len <- NULL
for(i in 1:5){
	chr <- hei10.geno[[i]]
	map <- chr$map
	hei10chr.len <- c(hei10chr.len,map[length(map)])
}
30427671 19698289 23459830 18585056 26975502 
219.5239 155.8200 151.1703 113.9596 195.8517 
hei10tot.len <- round(sum(hei10chr.len))
836

recq.geno <- recq.map$geno
recqchr.len <- NULL
for(i in 1:5){
	chr <- recq.geno[[i]]
	map <- chr$map
	recqchr.len <- c(recqchr.len,map[length(map)])
}
30427671 19698289 23459830 18585056 26975502 
382.7160 240.3915 307.4288 164.5206 350.8401 
recqtot.len <- round(sum(recqchr.len))
1446

recqhei10.geno <- recqhei10.map$geno
recqhei10chr.len <- NULL
for(i in 1:5){
	chr <- recqhei10.geno[[i]]
	map <- chr$map
	recqhei10chr.len <- c(recqhei10chr.len,map[length(map)])
}
30427671 19698289 23459830 18585056 26975502 
517.2525 302.1269 381.2878 183.0497 466.6904 
recqhei10tot.len <- round(sum(recqhei10chr.len))
1850

################################################
# analysis variance in crossover values per F2 #
################################################

mean(hswtmask.ind)
7.51
sd(hswtmask.ind)
2.26

mean(cuhei10mask.ind)
15.07
sd(cuhei10mask.ind)
3.41

hsrecqmask.ind
mean(hsrecqmask.ind)
25.04
sd(hsrecqmask.ind)
6.52

hsrecqhei10mask.ind
mean(hsrecqhei10mask.ind)
30.77
sd(hsrecqhei10mask.ind)
6.79

library(vcd)
hs.wt.gf <- goodfit(hswtmask.ind,type="poisson",method="ML")
cu.hei10.gf <- goodfit(cuhei10mask.ind,type="poisson",method="ML")
hs.recq.gf <- goodfit(hsrecqmask.ind,type="poisson",method="ML")
hs.recqhei10.gf <- goodfit(hsrecqhei10mask.ind,type="poisson",method="ML")

plot(hs.wt.gf,main="hs.wt",xlim=c(0,55))
summary(hs.wt.gf)
                      X^2 df   P(> X^2)
Likelihood Ratio 25.5611 12 0.01237671

plot(cu.hei10.gf,main="cu.hei10",xlim=c(0,55))
summary(cu.hei10.gf)
  X^2 df    P(> X^2)
Likelihood Ratio 36.25779 17 0.004236534

plot(hs.recq.gf,main="hs.recq",xlim=c(0,55))
summary(hs.recq.gf)
   X^2 df     P(> X^2)
Likelihood Ratio 71.55466 29 1.852519e-05

plot(hs.recqhei10.gf,main="hs.recqhei10",xlim=c(0,55))
summary(hs.recqhei10.gf)
      X^2 df   P(> X^2)
Likelihood Ratio 49.8203 31 0.01746465

library(lawstat)

dat <- c(hswtmask.ind,cuhei10mask.ind)
group <- c(rep(1,length(hswtmask.ind)),rep(2,length(cuhei10mask.ind)))
levene.test(dat,group)
Test Statistic = 24.474, p-value = 1.078e-06

dat <- c(hswtmask.ind,hsrecqmask.ind)
group <- c(rep(1,length(hs.wt.ind)),rep(2,length(hs.recq.ind)))
levene.test(dat,group)
Test Statistic = 149.62, p-value < 2.2e-16

dat <- c(hswtmask.ind,hsrecqhei10mask.ind)
group <- c(rep(1,length(hswtmask.ind)),rep(2,length(hsrecqhei10mask.ind)))
levene.test(dat,group)
Test Statistic = 164.01, p-value < 2.2e-16

############################################################################################
# intercrossover distances - analyse the same number of random positions and compare ratio #
############################################################################################

hswt.mask <- mask.cos[which(mask.cos[,8]=="hs.wt"),]
cuhei10.mask <- mask.cos[which(mask.cos[,8]=="cu.hei10"),]
hsrecq.mask <- mask.cos[which(mask.cos[,8]=="hs.recq"),]
hsrecqhei10.mask <- mask.cos[which(mask.cos[,8]=="hs.recqhei10"),]

par(mfcol=c(2,2))
par(mar=c(1.8,1.8,1.8,1.8))

chr.ends <- c(30427671,19698289,23459830,18585056,26975502)

dat <- hswt.mask
lane.cos <- NULL
lane.ran <- NULL
for(j in 1:length(unique(dat[,1]))){
	lane <- dat[which(dat[,1]==j),]
	lib.cos <- NULL
	lib.ran <- NULL
	for(k in 1:length(unique(dat[,2]))){
		lib <- lane[which(lane[,2]==dat[k,2]),]
		chr.cos <- NULL
		chr.ran <- NULL
		for(i in 1:5){
			chr <- lib[which(lib[,3]==i),]
			if(length(chr[,1]>1)){
				ran <- round(runif(length(chr[,1]),min=1,max=chr.ends[i]))
				ran <- ran[order(ran)]
				chr.ran <- c(chr.ran,diff(ran))
				chr.cos <- c(chr.cos,diff(chr[,6]))
				print("TRUE")	
				}
			}
		lib.cos <- c(lib.cos,chr.cos)
		lib.ran <- c(lib.ran,chr.ran)
		}
	lane.cos <- c(lane.cos,lib.cos)
	lane.ran <- c(lane.ran,lib.ran)
}
plot(density(lane.cos),col=2,ylim=c(0,max(density(lane.ran)[[2]])),xlim=c(0,30000000),main="hs.wt.mask")
abline(v=mean(lane.cos),col=2,lty=2)
par(new=T)
plot(density(lane.ran),col=4,ylim=c(0,max(density(lane.ran)[[2]])),xlim=c(0,30000000),main="")
abline(v=mean(lane.ran),col=4,lty=2)
mean(lane.cos)/mean(lane.ran)
1.262246
code <- c(rep(1,length(lane.cos)),c(rep(2,length(lane.ran))))
dco <- c(lane.cos,lane.ran)
dat <- cbind(code,dco)
wilcox <- wilcox.test(dco ~ code, data=dat)	
wilcox$p.value
1.390334e-11
mean(lane.cos)
8718472
mean(lane.ran)
6907112

dat <- cuhei10.mask
lane.cos <- NULL
lane.ran <- NULL
for(j in 1:length(unique(dat[,1]))){
	lane <- dat[which(dat[,1]==j),]
	lib.cos <- NULL
	lib.ran <- NULL
	for(k in 1:length(unique(dat[,2]))){
		lib <- lane[which(lane[,2]==dat[k,2]),]
		chr.cos <- NULL
		chr.ran <- NULL
		for(i in 1:5){
			chr <- lib[which(lib[,3]==i),]
			if(length(chr[,1]>1)){
				ran <- round(runif(length(chr[,1]),min=1,max=chr.ends[i]))
				ran <- ran[order(ran)]
				chr.ran <- c(chr.ran,diff(ran))
				chr.cos <- c(chr.cos,diff(chr[,6]))
				print("TRUE")	
				}
			}
		lib.cos <- c(lib.cos,chr.cos)
		lib.ran <- c(lib.ran,chr.ran)
		}
	lane.cos <- c(lane.cos,lib.cos)
	lane.ran <- c(lane.ran,lib.ran)
}
plot(density(lane.cos),col=2,ylim=c(0,max(density(lane.ran)[[2]])),xlim=c(0,30000000),main="cu.hei10.mask")
abline(v=mean(lane.cos),col=2,lty=2)
par(new=T)
plot(density(lane.ran),col=4,ylim=c(0,max(density(lane.ran)[[2]])),xlim=c(0,30000000),main="")
abline(v=mean(lane.ran),col=4,lty=2)
mean(lane.cos)/mean(lane.ran)
1.176459
code <- c(rep(1,length(lane.cos)),c(rep(2,length(lane.ran))))
dco <- c(lane.cos,lane.ran)
dat <- cbind(code,dco)
wilcox <- wilcox.test(dco ~ code, data=dat)	
wilcox$p.value
3.880427e-12
mean(lane.cos)
6058884
mean(lane.ran)
5150103

dat <- hsrecq.mask
lane.cos <- NULL
lane.ran <- NULL
for(j in 1:length(unique(dat[,1]))){
	lane <- dat[which(dat[,1]==j),]
	lib.cos <- NULL
	lib.ran <- NULL
	for(k in 1:length(unique(dat[,2]))){
		lib <- lane[which(lane[,2]==dat[k,2]),]
		chr.cos <- NULL
		chr.ran <- NULL
		for(i in 1:5){
			chr <- lib[which(lib[,3]==i),]
			if(length(chr[,1]>1)){
				ran <- round(runif(length(chr[,1]),min=1,max=chr.ends[i]))
				ran <- ran[order(ran)]
				chr.ran <- c(chr.ran,diff(ran))
				chr.cos <- c(chr.cos,diff(chr[,6]))
				print("TRUE")	
				}
			}
		lib.cos <- c(lib.cos,chr.cos)
		lib.ran <- c(lib.ran,chr.ran)
		}
	lane.cos <- c(lane.cos,lib.cos)
	lane.ran <- c(lane.ran,lib.ran)
}
plot(density(lane.cos),col=2,ylim=c(0,max(density(lane.cos)[[2]])),xlim=c(0,30000000),main="hs.recq.mask")
abline(v=mean(lane.cos),col=2,lty=2)
par(new=T)
plot(density(lane.ran),col=4,ylim=c(0,max(density(lane.cos)[[2]])),xlim=c(0,30000000),main="")
abline(v=mean(lane.ran),col=4,lty=2)
mean(lane.cos)/mean(lane.ran)
1.114793
code <- c(rep(1,length(lane.cos)),c(rep(2,length(lane.ran))))
dco <- c(lane.cos,lane.ran)
dat <- cbind(code,dco)
wilcox <- wilcox.test(dco ~ code, data=dat)	
wilcox$p.value
0.07821
mean(lane.cos)
3843722
mean(lane.ran)
3447926

dat <- hsrecqhei10.mask
lane.cos <- NULL
lane.ran <- NULL
for(j in 1:length(unique(dat[,1]))){
	lane <- dat[which(dat[,1]==j),]
	lib.cos <- NULL
	lib.ran <- NULL
	for(k in 1:length(unique(dat[,2]))){
		lib <- lane[which(lane[,2]==dat[k,2]),]
		chr.cos <- NULL
		chr.ran <- NULL
		for(i in 1:5){
			chr <- lib[which(lib[,3]==i),]
			if(length(chr[,1]>1)){
				ran <- round(runif(length(chr[,1]),min=1,max=chr.ends[i]))
				ran <- ran[order(ran)]
				chr.ran <- c(chr.ran,diff(ran))
				chr.cos <- c(chr.cos,diff(chr[,6]))
				print("TRUE")	
				}
			}
		lib.cos <- c(lib.cos,chr.cos)
		lib.ran <- c(lib.ran,chr.ran)
		}
	lane.cos <- c(lane.cos,lib.cos)
	lane.ran <- c(lane.ran,lib.ran)
}
plot(density(lane.cos),col=2,ylim=c(0,max(density(lane.cos)[[2]])),xlim=c(0,30000000),main="hs.recqhei10.mask")
abline(v=mean(lane.cos),col=2,lty=2)
par(new=T)
plot(density(lane.ran),col=4,ylim=c(0,max(density(lane.cos)[[2]])),xlim=c(0,30000000),main="")
abline(v=mean(lane.ran),col=4,lty=2)
mean(lane.cos)/mean(lane.ran)
1.061456
code <- c(rep(1,length(lane.cos)),c(rep(2,length(lane.ran))))
dco <- c(lane.cos,lane.ran)
dat <- cbind(code,dco)
wilcox <- wilcox.test(dco ~ code, data=dat)	
wilcox$p.value
0.2571502
mean(lane.cos)
3300015
mean(lane.ran)
3089627

##################################
# plotting along the chromosomes #
##################################

hswt.mask <- mask.cos[which(mask.cos[,8]=="hs.wt"),]
cuhei10.mask <- mask.cos[which(mask.cos[,8]=="cu.hei10"),]
hsrecq.mask <- mask.cos[which(mask.cos[,8]=="hs.recq"),]
hsrecqhei10.mask <- mask.cos[which(mask.cos[,8]=="hs.recqhei10"),]

hs.wt.n <- 245
cu.hei10.n <- 192
hs.recq.n <- 191
hs.recqhei10.n <- 192
samp <- c(hs.wt.n,cu.hei10.n,hs.recq.n,hs.recqhei10.n)

library(GenomicRanges)
centromeres <- c(15086045,3607929,13587786,3956021,11725024)
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chr.ends <- c(30427671,19698289,23459830,18585056,26975502)
tha.cum <- cumsum(chr.ends)
tha.cum <- c(0,tha.cum)
tha.tot <- tha.cum[length(tha.cum)]
centromeres[2] <- centromeres[2]+tha.cum[2]
centromeres[3] <- centromeres[3]+tha.cum[3]
centromeres[4] <- centromeres[4]+tha.cum[4]
centromeres[5] <- centromeres[5]+tha.cum[5]

thal.r.genes <- read.csv(file="S1_Table.csv",header=T)
all.thalr.tss <- NULL
for(i in 1:5){
	chr.thalr.genes <- thal.r.genes[which(thal.r.genes[,3]==i),]
	thalr.tss <- chr.thalr.genes[,5]
	thalr.tss <- thalr.tss+tha.cum[i]
	all.thalr.tss <- c(all.thalr.tss,thalr.tss)
}
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
r.genes <- read.csv(file="S1_table.csv")
coller.snps <- read.table(file="Ler_hq.291K_50bp.txt")
coller.snvmnv <- read.table(file="SNV.MNV.allelic.txt")
coller.invs <- read.csv(file="ColLer_Inversions.csv")
yel.meth <- read.delim(file="Index1_CG_CHG_CHH_10kb_value_above_0.txt",header=T)

test <- seq(1,500,by=1)
j=8
ma <- rep(1,test[j])/test[j]
cos.matrix <- NULL
chr.invs <- NULL
for(i in 1:5){
	print(i)
	hswt.chr <- hswt.mask[which(hswt.mask[,3]==i),]
	cuhei10.chr <- cuhei10.mask[which(cuhei10.mask[,3]==i),]
	hsrecq.chr <- hsrecq.mask[which(hsrecq.mask[,3]==i),]
	hsrecqhei10.chr <- hsrecqhei10.mask[which(hsrecqhei10.mask[,3]==i),]
	print(i)
	windows <- seq(1,chr.ends[i],by=300000)
	windows <- c(windows,chr.ends[i])
	cum.windows <- windows+tha.cum[i]
	print(i)
	chr.coller.snps <- coller.snps[which(coller.snps[,2]==i),]
	snp.coords <- chr.coller.snps[,3]
	chr.coller.snvmnv <- coller.snvmnv[which(coller.snvmnv[,1]==chrs[i]),]
	snvmnv.coords <- chr.coller.snvmnv[,2]
	chr.coller.invs <- coller.invs[which(coller.invs[,1]==chrs[i]),]
	print(i)
	win.hswt <- NULL
	win.cuhei10 <- NULL
	win.hsrecq <- NULL
	win.hsrecqhei10 <- NULL
	win.snps <- NULL	
	win.snvmnv <- NULL
	for(j in 1:length(windows)-1){
		print(j)
		win.hswt <- c(win.hswt,length(which(hswt.chr[,6]>=windows[j]&hswt.chr[,6]<=windows[j+1])))
		win.cuhei10 <- c(win.cuhei10,length(which(cuhei10.chr[,6]>=windows[j]&cuhei10.chr[,6]<=windows[j+1])))
		win.hsrecq <- c(win.hsrecq,length(which(hsrecq.chr[,6]>=windows[j]&hsrecq.chr[,6]<=windows[j+1])))
		win.hsrecqhei10 <- c(win.hsrecqhei10,length(which(hsrecqhei10.chr[,6]>=windows[j]&hsrecqhei10.chr[,6]<=windows[j+1])))
		win.snps <- c(win.snps,length(which(snp.coords>=windows[j]&snp.coords<=windows[j+1])))
		win.snvmnv <- c(win.snvmnv,sum(chr.coller.snvmnv[which(snvmnv.coords>=windows[j]&snvmnv.coords<=windows[j+1]),9]))
	}
	win.hswt <- win.hswt/hs.wt.n
	win.cuhei10 <- win.cuhei10/cu.hei10.n
	win.hsrecq <- win.hsrecq/hs.recq.n
	win.hsrecqhei10 <- win.hsrecqhei10/hs.recqhei10.n
	filt.hsrecq <- filter(win.hsrecq,ma)
	which.na <- which(is.na(filt.hsrecq)==TRUE)
	left.na <- which.na[which(which.na<20)]
	left.val <- filt.hsrecq[left.na[length(left.na)]+1]
	filt.hsrecq[left.na] <- left.val
	right.na <- which.na[which(which.na>55)]
	right.val <- filt.hsrecq[right.na[1]-1]
	filt.hsrecq[right.na] <- right.val
	filt.hsrecqhei10 <- filter(win.hsrecqhei10,ma)
	which.na <- which(is.na(filt.hsrecqhei10)==TRUE)
	left.na <- which.na[which(which.na<20)]
	left.val <- filt.hsrecqhei10[left.na[length(left.na)]+1]
	filt.hsrecqhei10[left.na] <- left.val
	right.na <- which.na[which(which.na>55)]
	right.val <- filt.hsrecqhei10[right.na[1]-1]
	filt.hsrecqhei10[right.na] <- right.val
	filt.cuhei10 <- filter(win.cuhei10,ma)
	which.na <- which(is.na(filt.cuhei10)==TRUE)
	left.na <- which.na[which(which.na<20)]
	left.val <- filt.cuhei10[left.na[length(left.na)]+1]
	filt.cuhei10[left.na] <- left.val
	right.na <- which.na[which(which.na>55)]
	right.val <- filt.cuhei10[right.na[1]-1]
	filt.cuhei10[right.na] <- right.val
	filt.hswt <- filter(win.hswt,ma)
	which.na <- which(is.na(filt.hswt)==TRUE)
	left.na <- which.na[which(which.na<20)]
	left.val <- filt.hswt[left.na[length(left.na)]+1]
	filt.hswt[left.na] <- left.val
	right.na <- which.na[which(which.na>55)]
	right.val <- filt.hswt[right.na[1]-1]
	filt.hswt[right.na] <- right.val
	filt.snps <- filter(win.snps,ma)
	which.na <- which(is.na(filt.snps)==TRUE)
	left.na <- which.na[which(which.na<20)]
	left.val <- filt.snps[left.na[length(left.na)]+1]
	filt.snps[left.na] <- left.val
	right.na <- which.na[which(which.na>55)]
	right.val <- filt.snps[right.na[1]-1]
	filt.snps[right.na] <- right.val
	filt.snvs <- filter(win.snvmnv,ma)
	which.na <- which(is.na(filt.snvs)==TRUE)
	left.na <- which.na[which(which.na<20)]
	left.val <- filt.snvs[left.na[length(left.na)]+1]
	filt.snvs[left.na] <- left.val
	right.na <- which.na[which(which.na>55)]
	right.val <- filt.snvs[right.na[1]-1]
	filt.snvs[right.na] <- right.val
	chr.dat <- cbind(rep(i,length(windows)),cum.windows,filt.hswt,filt.cuhei10,filt.hsrecq,filt.hsrecqhei10,filt.snps,filt.snvs)
	cos.matrix <- rbind(cos.matrix,chr.dat)	
	cum.invs <- cbind(chr.coller.invs[,2]+tha.cum[i],chr.coller.invs[,3]+tha.cum[i])
	chr.invs <- rbind(chr.invs,cum.invs)
	print(i)
}

j=100
ma <- rep(1,test[j])/test[j]
chr.dat <- NULL
all.yel.dat <- NULL
for(i in 1:5){
	print(i)
	chip.dat <- read.table(file=paste(i,"_rec8.msh4.nkd.10kb.txt",sep=""))
	genetes.dat <- read.table(file=paste(i,"_chrcov_genetes.txt",sep="")) 
	nucs.dat <- read.table(file=paste(i,".WT_nuc.both.10kb.txt",sep=""))
	rpi.dat <- read.table(file=paste(i,".col.both.10kb.txt",sep=""))
	cum.windows <- chip.dat[,1]+tha.cum[i]
	dat <- cbind(cum.windows,chip.dat[,2:8],genetes.dat[,2:3],nucs.dat[,2],rpi.dat[,2])
	chr.dat <- rbind(chr.dat,dat)
	print(i)
	filt.rec8chip <- filter(chr.dat[,2],ma)
	filt.rec8in <- filter(chr.dat[,3],ma)
	filt.msh4chip <- filter(chr.dat[,4],ma)
	filt.msh4in <- filter(chr.dat[,5],ma)
	filt.nkd <- filter(chr.dat[,6],ma)
	filt.rec8norm <- filter(chr.dat[,7],ma)
	filt.msh4norm <- filter(chr.dat[,8],ma)
	filt.gene <- filter(chr.dat[,9],ma)
	filt.tes <- filter(chr.dat[,10],ma)
	filt.nucs <- filter(chr.dat[,11],ma)
	filt.rpi <- filter(chr.dat[,12],ma)
	print(i)
	yel.chr <- yel.meth[which(yel.meth[,2]==i),]
	yel.wins <- c(yel.chr[,3],chr.ends[i])
	ycum.wins <- yel.wins+tha.cum[i]
	yel.all <- c(yel.chr[,13],0)
	yel.cg <- c(yel.chr[,14],0)
	yel.chg <- c(yel.chr[,15],0)
	yel.chh <- c(yel.chr[,16],0)
	yel.dat <- cbind(ycum.wins,yel.all,yel.cg,yel.chg,yel.chh)
	print(i)
	filt.meth <- filter(yel.dat[,2],ma)
	which.na <- which(is.na(filt.meth)==TRUE)
	left.na <- which.na[which(which.na<100)]
	left.val <- filt.meth[left.na[length(left.na)]+1]
	filt.meth[left.na] <- left.val
	right.na <- which.na[which(which.na>100)]
	right.val <- filt.meth[right.na[1]-1]
	filt.meth[right.na] <- right.val
	print(i)
	yel.dat <- cbind(yel.dat,filt.meth)
	all.yel.dat <- rbind(all.yel.dat,yel.dat)
}

#####################################
# plotting for Serra et al Figure 3 #
#####################################

par(mfcol=c(3,1))
par(mar=c(1.8,1.8,1.8,1.8))
centromeres <- c(15086045,3607929,13587786,3956021,11725024)
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chr.ends <- c(30427671,19698289,23459830,18585056,26975502)
tha.cum <- cumsum(chr.ends)
tha.cum <- c(0,tha.cum)
tha.tot <- tha.cum[length(tha.cum)]
centromeres[2] <- centromeres[2]+tha.cum[2]
centromeres[3] <- centromeres[3]+tha.cum[3]
centromeres[4] <- centromeres[4]+tha.cum[4]
centromeres[5] <- centromeres[5]+tha.cum[5]

inversions <- read.csv(file="ColLer_inversions.csv",header=T)
knob <- inversions[34,]
knob.start <- knob[,2]+tha.cum[4]
knob.end <- knob[,3]+tha.cum[4]

xlim <- c(0,tha.tot)
ylim=c(0,0.14)
plot(cos.matrix[,2],cos.matrix[,3],type="l",xlim=xlim,ylim=ylim,main="hs.wt vs HEI10 vs recq4a4b vs HEI10recq4a4b",col=1)
par(new=T)
plot(cos.matrix[,2],cos.matrix[,4],type="l",xlim=xlim,ylim=ylim,col=4)
par(new=T)
plot(cos.matrix[,2],cos.matrix[,5],type="l",xlim=xlim,ylim=ylim,col="purple")
par(new=T)
plot(cos.matrix[,2],cos.matrix[,6],type="l",xlim=xlim,ylim=ylim,col=2)
abline(v=centromeres,col=1,lty=2,lwd=0.5)
abline(v=tha.cum,col=1,lty=1,lwd=0.3)
rug(all.thalr.tss,col=1,ticksize=0.05)

plot(cos.matrix[,2],cos.matrix[,3],type="l",xlim=xlim,ylim=c(0,0.045),main="hs.wt vs SNPs",col=1)
par(new=T)
plot(cos.matrix[,2],cos.matrix[,7],type="l",xlim=xlim,col="dark green",yaxt="n",lwd=1)
axis(side=4,col="dark green")
par(new=T)
plot(all.yel.dat[,1],all.yel.dat[,6],type="l",xlim=xlim,col=4,yaxt="n",lwd=1)
axis(side=4,col=4)
abline(v=centromeres,col=1,lty=2,lwd=0.5)
abline(v=tha.cum,col=1,lty=1,lwd=0.3)
rug(all.thalr.tss,col=1,ticksize=0.05)

ylim=c(0,0.15)
plot(cos.matrix[,2],cos.matrix[,3],type="l",xlim=xlim,ylim=ylim,main="hs.wt vs HEI10 recq4a4b vs SNPs")
par(new=T)
plot(cos.matrix[,2],cos.matrix[,6],type="l",xlim=xlim,ylim=ylim,col=2)
par(new=T)
plot(cos.matrix[,2],cos.matrix[,7],type="l",xlim=xlim,col="dark green",yaxt="n",lwd=1)
axis(side=4,col="dark green")
par(new=T)
plot(all.yel.dat[,1],all.yel.dat[,6],type="l",xlim=xlim,col=4,yaxt="n",lwd=1)
axis(side=4,col=4)
abline(v=centromeres,col=1,lty=2,lwd=0.5)
abline(v=tha.cum,col=1,lty=1,lwd=0.3)
rug(all.thalr.tss,col=1,ticksize=0.05)

##############################################
# correlating DNA methylation and crossovers #
##############################################

yel.coords <- all.yel.dat[,1]
cos.coords <- cos.matrix[,2]
yel.wins <- NULL
for(k in 1:length(cos.coords)-1){
	print(k)
	print(which(yel.coords>=cos.coords[k]&yel.coords<cos.coords[k+1]))
	yel.wins <- c(yel.wins,mean(all.yel.dat[which(yel.coords>=cos.coords[k]&yel.coords<cos.coords[k+1]),2]))
}
plot(cos.matrix[,2],cos.matrix[,3],type="l",col=2)
par(new=T)
plot(cos.matrix[,2],yel.wins,type="l",col=4)

cor.test(yel.wins,cos.matrix[,3])
# hs.wt -0.233 2.134e-06
cor.test(yel.wins,cos.matrix[,4])
# cu.hei10 -0.7398002  <2.2e-16
cor.test(yel.wins,cos.matrix[,5])
# hs.recq -0.828 <2.2e-16
cor.test(yel.wins,cos.matrix[,6])
# hs.recq -0.8101 <2.2e-16

#################################################################################
# seperate into pericentromeres/centromeres versus arms - measure and correlate #
#################################################################################
###################################################
# correlate at 10 kb, 50 kb, 100 kb, 500 kb, 1 Mb #
###################################################

hswt.mask <- mask.cos[which(mask.cos[,8]=="hs.wt"),]
cuhei10.mask <- mask.cos[which(mask.cos[,8]=="cu.hei10"),]
hsrecq.mask <- mask.cos[which(mask.cos[,8]=="hs.recq"),]
hsrecqhei10.mask <- mask.cos[which(mask.cos[,8]=="hs.recqhei10"),]
hs.wt.n <- 245
cu.hei10.n <- 192
hs.recq.n <- 191
hs.recqhei10.n <- 192

library(GenomicRanges)
centromeres <- c(15086045,3607929,13587786,3956021,11725024)
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chr.ends <- c(30427671,19698289,23459830,18585056,26975502)
tha.cum <- cumsum(chr.ends)
tha.cum <- c(0,tha.cum)
tha.tot <- tha.cum[length(tha.cum)]
centromeres[2] <- centromeres[2]+tha.cum[2]
centromeres[3] <- centromeres[3]+tha.cum[3]
centromeres[4] <- centromeres[4]+tha.cum[4]
centromeres[5] <- centromeres[5]+tha.cum[5]
thal.r.genes <- read.csv(file="S1_Table.csv",header=T)
all.thalr.tss <- NULL
for(i in 1:5){
	chr.thalr.genes <- thal.r.genes[which(thal.r.genes[,3]==i),]
	thalr.tss <- chr.thalr.genes[,5]
	thalr.tss <- thalr.tss+tha.cum[i]
	all.thalr.tss <- c(all.thalr.tss,thalr.tss)
}
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
coller.snps <- read.table(file="Ler_hq.291K_50bp.txt")
chr.ends <- c(30427671,19698289,23459830,18585056,26975502)
centromeres <- c(15086045,3607929,13587786,3956021,11725024)

win.size <- 500000
test <- seq(1,500,by=1)
k=30
ma <- rep(1,test[k])/test[k]
cos.matrix <- NULL
for(i in 1:5){
	print(i)
	hswt.chr <- hswt.mask[which(hswt.mask[,3]==i),]
	cuhei10.chr <- cuhei10.mask[which(cuhei10.mask[,3]==i),]
	hsrecq.chr <- hsrecq.mask[which(hsrecq.mask[,3]==i),]
	hsrecqhei10.chr <- hsrecqhei10.mask[which(hsrecqhei10.mask[,3]==i),]
	print(i)
	windows <- seq(1,chr.ends[i],by=win.size)
	windows <- c(windows,chr.ends[i])
	cum.windows <- windows+tha.cum[i]
	print(i)
	chr.coller.snps <- coller.snps[which(coller.snps[,2]==i),]
	snp.coords <- chr.coller.snps[,3]
	print(i)
	win.hswt <- NULL
	win.cuhei10 <- NULL
	win.hsrecq <- NULL
	win.hsrecqhei10 <- NULL
	win.snps <- NULL	
	for(j in 1:length(windows)-1){
		print(j)
		win.hswt <- c(win.hswt,length(which(hswt.chr[,6]>=windows[j]&hswt.chr[,6]<=windows[j+1])))
		win.cuhei10 <- c(win.cuhei10,length(which(cuhei10.chr[,6]>=windows[j]&cuhei10.chr[,6]<=windows[j+1])))
		win.hsrecq <- c(win.hsrecq,length(which(hsrecq.chr[,6]>=windows[j]&hsrecq.chr[,6]<=windows[j+1])))
		win.hsrecqhei10 <- c(win.hsrecqhei10,length(which(hsrecqhei10.chr[,6]>=windows[j]&hsrecqhei10.chr[,6]<=windows[j+1])))
		win.snps <- c(win.snps,length(which(snp.coords>=windows[j]&snp.coords<=windows[j+1])))
	}
	win.hswt <- win.hswt/hs.wt.n
	win.cuhei10 <- win.cuhei10/cu.hei10.n
	win.hsrecq <- win.hsrecq/hs.recq.n
	win.hsrecqhei10 <- win.hsrecqhei10/hs.recqhei10.n
	print(i)
	filt.lim <- 20
	print(i)
	filt.hswt <- filter(win.hswt,ma)
	which.na <- which(is.na(filt.hswt)==TRUE)
	left.na <- which.na[which(which.na<filt.lim)]
	left.val <- filt.hswt[left.na[length(left.na)]+1]
	filt.hswt[left.na] <- left.val
	right.na <- which.na[which(which.na>filt.lim)]
	right.val <- filt.hswt[right.na[1]-1]
	filt.hswt[right.na] <- right.val
	print(i)	
	filt.cuhei10 <- filter(win.cuhei10,ma)
	which.na <- which(is.na(filt.cuhei10)==TRUE)
	left.na <- which.na[which(which.na<filt.lim)]
	left.val <- filt.cuhei10[left.na[length(left.na)]+1]
	filt.cuhei10[left.na] <- left.val
	right.na <- which.na[which(which.na>filt.lim)]
	right.val <- filt.cuhei10[right.na[1]-1]
	filt.cuhei10[right.na] <- right.val
	print(i)	
	filt.hsrecq <- filter(win.hsrecq,ma)
	which.na <- which(is.na(filt.hsrecq)==TRUE)
	left.na <- which.na[which(which.na<filt.lim)]
	left.val <- filt.hsrecq[left.na[length(left.na)]+1]
	filt.hsrecq[left.na] <- left.val
	right.na <- which.na[which(which.na>filt.lim)]
	right.val <- filt.hsrecq[right.na[1]-1]
	filt.hsrecq[right.na] <- right.val
	print(i)
	filt.hsrecqhei10 <- filter(win.hsrecqhei10,ma)
	which.na <- which(is.na(filt.hsrecqhei10)==TRUE)
	left.na <- which.na[which(which.na<filt.lim)]
	left.val <- filt.hsrecqhei10[left.na[length(left.na)]+1]
	filt.hsrecqhei10[left.na] <- left.val
	right.na <- which.na[which(which.na>filt.lim)]
	right.val <- filt.hsrecqhei10[right.na[1]-1]
	filt.hsrecqhei10[right.na] <- right.val
	print(i)	
	filt.snps <- filter(win.snps,ma)
	which.na <- which(is.na(filt.snps)==TRUE)
	left.na <- which.na[which(which.na<filt.lim)]
	left.val <- filt.snps[left.na[length(left.na)]+1]
	filt.snps[left.na] <- left.val
	right.na <- which.na[which(which.na>filt.lim)]
	right.val <- filt.snps[right.na[1]-1]
	filt.snps[right.na] <- right.val
	print(i)
	chr.dat <- cbind(rep(i,length(windows)),cum.windows,win.hswt,win.cuhei10,win.hsrecq,win.hsrecqhei10,win.snps,filt.hswt,filt.cuhei10,filt.hsrecq,filt.hsrecqhei10,filt.snps)
	cos.matrix <- rbind(cos.matrix,chr.dat)	
	print(i)
}

# @ 50 kb

cor.test(cos.matrix[,8],cos.matrix[,12])
0.421786 <2.2e-16
cor.test(cos.matrix[,9],cos.matrix[,12])
-0.1634 <2.2e-16
cor.test(cos.matrix[,10],cos.matrix[,12])
-0.498 <2.2e-16
cor.test(cos.matrix[,11],cos.matrix[,12])
-0.494 <2.2e-16

# @ 100 kb

cor.test(cos.matrix[,8],cos.matrix[,12])
0.440 <2.2e-16
cor.test(cos.matrix[,9],cos.matrix[,12])
-0.300 <2.2e-16
cor.test(cos.matrix[,10],cos.matrix[,12])
-0.580 <2.2e-16
cor.test(cos.matrix[,11],cos.matrix[,12])
-0.584 <2.2e-16

# @ 200 kb

cor.test(cos.matrix[,8],cos.matrix[,12])
0.475 <2.2e-16
cor.test(cos.matrix[,9],cos.matrix[,12])
-0.523 <2.2e-16
cor.test(cos.matrix[,10],cos.matrix[,12])
-0.726 <2.2e-16
cor.test(cos.matrix[,11],cos.matrix[,12])
-0.741 <2.2e-16

# @ 300 kb

cor.test(cos.matrix[,8],cos.matrix[,12])
0.564 <2.2e-16
cor.test(cos.matrix[,9],cos.matrix[,12])
-0.640 <2.2e-16
cor.test(cos.matrix[,10],cos.matrix[,12])
-0.805 <2.2e-16
cor.test(cos.matrix[,11],cos.matrix[,12])
-0.810 <2.2e-16

# @ 400 kb

cor.test(cos.matrix[,8],cos.matrix[,12])
0.295 <2.2e-16
cor.test(cos.matrix[,9],cos.matrix[,12])
-0.629 <2.2e-16
cor.test(cos.matrix[,10],cos.matrix[,12])
-0.825 <2.2e-16
cor.test(cos.matrix[,11],cos.matrix[,12])
-0.818 <2.2e-16

cen.coords <- read.csv(file="cen.pericen.coords.csv")
cen.pericen <- read.csv(file="cen.pericen.coords.csv",header=T)
cen.pericen <- cen.pericen[,-1]
cen.pericen <- cen.pericen[,-1]
cen.pericen[2,] <- cen.pericen[2,]+tha.cum[2]
cen.pericen[3,] <- cen.pericen[3,]+tha.cum[3]
cen.pericen[4,] <- cen.pericen[4,]+tha.cum[4]
cen.pericen[5,] <- cen.pericen[5,]+tha.cum[5]

all.cen <- NULL
all.arm <- NULL
hswt.cen <- NULL
hswt.arm <- NULL
cuhei10.cen <- NULL
cuhei10.arm <- NULL
hsrecq.cen <- NULL
hsrecq.arm <- NULL
hsrecqhei10.cen <- NULL
hsrecqhei10.arm <- NULL
snps.cen <- NULL
snps.arm <- NULL
for(i in 1:5){
	print(i)
	chr.matrix <- cos.matrix[which(cos.matrix[,1]==i),]
	print(i)
	all.cen <- rbind(all.cen,chr.matrix[which(chr.matrix[,2]>cen.pericen[i,1]&chr.matrix[,2]<cen.pericen[i,4]),])
	all.arm <- rbind(all.arm,chr.matrix[-which(chr.matrix[,2]>cen.pericen[i,1]&chr.matrix[,2]<cen.pericen[i,4]),])
	snps.cen <- c(snps.cen,chr.matrix[which(chr.matrix[,2]>cen.pericen[i,1]&chr.matrix[,2]<cen.pericen[i,4]),7])
	snps.arm <- c(snps.arm,chr.matrix[-which(chr.matrix[,2]>cen.pericen[i,1]&chr.matrix[,2]<cen.pericen[i,4]),7])
	hswt.cen <- c(hswt.cen,chr.matrix[which(chr.matrix[,2]>cen.pericen[i,1]&chr.matrix[,2]<cen.pericen[i,4]),3])
	hswt.arm <- c(hswt.arm,chr.matrix[-which(chr.matrix[,2]>cen.pericen[i,1]&chr.matrix[,2]<cen.pericen[i,4]),3])
	cuhei10.cen <- c(cuhei10.cen,chr.matrix[which(chr.matrix[,2]>cen.pericen[i,1]&chr.matrix[,2]<cen.pericen[i,4]),4])
	cuhei10.arm <- c(cuhei10.arm,chr.matrix[-which(chr.matrix[,2]>cen.pericen[i,1]&chr.matrix[,2]<cen.pericen[i,4]),4])
	hsrecq.cen <- c(hsrecq.cen,chr.matrix[which(chr.matrix[,2]>cen.pericen[i,1]&chr.matrix[,2]<cen.pericen[i,4]),5])
	hsrecq.arm <- c(hsrecq.arm,chr.matrix[-which(chr.matrix[,2]>cen.pericen[i,1]&chr.matrix[,2]<cen.pericen[i,4]),5])
	hsrecqhei10.cen <- c(hsrecqhei10.cen,chr.matrix[which(chr.matrix[,2]>cen.pericen[i,1]&chr.matrix[,2]<cen.pericen[i,4]),6])
	hsrecqhei10.arm <- c(hsrecqhei10.arm,chr.matrix[-which(chr.matrix[,2]>cen.pericen[i,1]&chr.matrix[,2]<cen.pericen[i,4]),6])
}

############################################################################
# seperate into pericentromeres/centromeres versus arms - count crossovers #
############################################################################

hswt.mask <- mask.cos[which(mask.cos[,8]=="hs.wt"),]
cuhei10.mask <- mask.cos[which(mask.cos[,8]=="cu.hei10"),]
hsrecq.mask <- mask.cos[which(mask.cos[,8]=="hs.recq"),]
hsrecqhei10.mask <- mask.cos[which(mask.cos[,8]=="hs.recqhei10"),]
hs.wt.n <- 245
cu.hei10.n <- 192
hs.recq.n <- 191
hs.recqhei10.n <- 192

chr.ends <- c(30427671,19698289,23459830,18585056,26975502)
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
cen.coords <- read.csv(file="cen.pericen.coords.csv")

dat <- hswt.mask
dat.ann <- NULL
for(i in 1:5){
	print(i)
	chr.cen <- cen.coords[i,]
	chr.dat <- dat[which(dat[,3]==i),6]
	chr.dat <- chr.dat[order(chr.dat)]
	chr.dat <- cbind(chr.dat,rep(0,length(chr.dat)),rep(i,length(chr.dat)))
	dat.coords <- as.numeric(chr.dat[,1])
	chr.dat[which(dat.coords>=1&dat.coords<=chr.cen[1,3]),2] <- 1
	chr.dat[which(dat.coords>=chr.cen[1,3]&dat.coords<=chr.cen[1,4]),2] <- 2
	chr.dat[which(dat.coords>=chr.cen[1,4]&dat.coords<=chr.cen[1,5]),2] <- 3
	chr.dat[which(dat.coords>=chr.cen[1,5]&dat.coords<=chr.cen[1,6]),2] <- 2
	chr.dat[which(dat.coords>=chr.cen[1,6]&dat.coords<=chr.ends[i]),2] <- 1
	dat.ann <- rbind(dat.ann,chr.dat)
	print(dim(dat.ann))
}
dat.count <- NULL
for(i in 1:5){
	chr.dat <- dat.ann[which(dat.ann[,3]==i),]
	count <- c(length(chr.dat[which(chr.dat[,2]==1),2]),length(chr.dat[which(chr.dat[,2]==2),2]),length(chr.dat[which(chr.dat[,2]==3),2]))
	dat.count <- cbind(dat.count,count)
}
dat.count <- cbind(dat.count,rowSums(dat.count))
dat.count <- rbind(dat.count,colSums(dat.count))
write.csv(dat.count,file="mask.summary.hswt.csv")

dat <- cuhei10.mask
dat.ann <- NULL
for(i in 1:5){
	print(i)
	chr.cen <- cen.coords[i,]
	chr.dat <- dat[which(dat[,3]==i),6]
	chr.dat <- chr.dat[order(chr.dat)]
	chr.dat <- cbind(chr.dat,rep(0,length(chr.dat)),rep(i,length(chr.dat)))
	dat.coords <- as.numeric(chr.dat[,1])
	chr.dat[which(dat.coords>=1&dat.coords<=chr.cen[1,3]),2] <- 1
	chr.dat[which(dat.coords>=chr.cen[1,3]&dat.coords<=chr.cen[1,4]),2] <- 2
	chr.dat[which(dat.coords>=chr.cen[1,4]&dat.coords<=chr.cen[1,5]),2] <- 3
	chr.dat[which(dat.coords>=chr.cen[1,5]&dat.coords<=chr.cen[1,6]),2] <- 2
	chr.dat[which(dat.coords>=chr.cen[1,6]&dat.coords<=chr.ends[i]),2] <- 1
	dat.ann <- rbind(dat.ann,chr.dat)
	print(dim(dat.ann))
}
dat.count <- NULL
for(i in 1:5){
	chr.dat <- dat.ann[which(dat.ann[,3]==i),]
	count <- c(length(chr.dat[which(chr.dat[,2]==1),2]),length(chr.dat[which(chr.dat[,2]==2),2]),length(chr.dat[which(chr.dat[,2]==3),2]))
	dat.count <- cbind(dat.count,count)
}
dat.count <- cbind(dat.count,rowSums(dat.count))
dat.count <- rbind(dat.count,colSums(dat.count))
write.csv(dat.count,file="mask.summary.cuhei10.csv")

dat <- hsrecq.mask
dat.ann <- NULL
for(i in 1:5){
	print(i)
	chr.cen <- cen.coords[i,]
	chr.dat <- dat[which(dat[,3]==i),6]
	chr.dat <- chr.dat[order(chr.dat)]
	chr.dat <- cbind(chr.dat,rep(0,length(chr.dat)),rep(i,length(chr.dat)))
	dat.coords <- as.numeric(chr.dat[,1])
	chr.dat[which(dat.coords>=1&dat.coords<=chr.cen[1,3]),2] <- 1
	chr.dat[which(dat.coords>=chr.cen[1,3]&dat.coords<=chr.cen[1,4]),2] <- 2
	chr.dat[which(dat.coords>=chr.cen[1,4]&dat.coords<=chr.cen[1,5]),2] <- 3
	chr.dat[which(dat.coords>=chr.cen[1,5]&dat.coords<=chr.cen[1,6]),2] <- 2
	chr.dat[which(dat.coords>=chr.cen[1,6]&dat.coords<=chr.ends[i]),2] <- 1
	dat.ann <- rbind(dat.ann,chr.dat)
	print(dim(dat.ann))
}
dat.count <- NULL
for(i in 1:5){
	chr.dat <- dat.ann[which(dat.ann[,3]==i),]
	count <- c(length(chr.dat[which(chr.dat[,2]==1),2]),length(chr.dat[which(chr.dat[,2]==2),2]),length(chr.dat[which(chr.dat[,2]==3),2]))
	dat.count <- cbind(dat.count,count)
}
dat.count <- cbind(dat.count,rowSums(dat.count))
dat.count <- rbind(dat.count,colSums(dat.count))
write.csv(dat.count,file="mask.summary.hsrecq.csv")

dat <- hsrecqhei10.mask
dat.ann <- NULL
for(i in 1:5){
	print(i)
	chr.cen <- cen.coords[i,]
	chr.dat <- dat[which(dat[,3]==i),6]
	chr.dat <- chr.dat[order(chr.dat)]
	chr.dat <- cbind(chr.dat,rep(0,length(chr.dat)),rep(i,length(chr.dat)))
	dat.coords <- as.numeric(chr.dat[,1])
	chr.dat[which(dat.coords>=1&dat.coords<=chr.cen[1,3]),2] <- 1
	chr.dat[which(dat.coords>=chr.cen[1,3]&dat.coords<=chr.cen[1,4]),2] <- 2
	chr.dat[which(dat.coords>=chr.cen[1,4]&dat.coords<=chr.cen[1,5]),2] <- 3
	chr.dat[which(dat.coords>=chr.cen[1,5]&dat.coords<=chr.cen[1,6]),2] <- 2
	chr.dat[which(dat.coords>=chr.cen[1,6]&dat.coords<=chr.ends[i]),2] <- 1
	dat.ann <- rbind(dat.ann,chr.dat)
	print(dim(dat.ann))
}
dat.count <- NULL
for(i in 1:5){
	chr.dat <- dat.ann[which(dat.ann[,3]==i),]
	count <- c(length(chr.dat[which(chr.dat[,2]==1),2]),length(chr.dat[which(chr.dat[,2]==2),2]),length(chr.dat[which(chr.dat[,2]==3),2]))
	dat.count <- cbind(dat.count,count)
}
dat.count <- cbind(dat.count,rowSums(dat.count))
dat.count <- rbind(dat.count,colSums(dat.count))
write.csv(dat.count,file="mask.summary.hsrecqhei10.csv")

####################
# TEL-CEN analysis #
####################

chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
r.genes <- read.csv(file="S1_table.csv")
coller.snps <- read.table(file="Ler_hq.291K_50bp.txt")
coller.snvmnv <- read.table(file="SNV.MNV.allelic.txt")
coller.invs <- read.csv(file="ColLer_Inversions.csv")
yel.meth <- read.delim(file="Index1_CG_CHG_CHH_10kb_value_above_0.txt",header=T)

hswt.mask <- mask.cos[which(mask.cos[,8]=="hs.wt"),]
cuhei10.mask <- mask.cos[which(mask.cos[,8]=="cu.hei10"),]
hsrecq.mask <- mask.cos[which(mask.cos[,8]=="hs.recq"),]
hsrecqhei10.mask <- mask.cos[which(mask.cos[,8]=="hs.recqhei10"),]
hs.wt.n <- 245
cu.hei10.n <- 192
hs.recq.n <- 191
hs.recqhei10.n <- 192

wins <- seq(0,1,by=0.01)
chr.ends <- c(30427671,19698289,23459830,18585056,26975502)
centromeres <- c(15086045,3607929,13587786,3956021,11725024)
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")

snp.left <- NULL
snp.right <- NULL
for(i in 1:5){
	print(i)
	windows <- seq(1,chr.ends[i],by=100000)
	chr.snps <- coller.snps[which(coller.snps[,2]==i),]
	chr.snvs <- coller.snvmnv[which(coller.snvmnv[,2]==chrs[i]),]
	chr.meth <- yel.meth[which(yel.meth[,2]==i),]
	snp.coords <- chr.snps[,3]
	snv.coords <- chr.snvs[,2]
	win.snps <- NULL	
	win.snvmnv <- NULL
	win.meth <- NULL
	for(j in 1:length(windows)-1){
		win.snps <- c(win.snps,length(which(snp.coords>=windows[j]&snp.coords<=windows[j+1])))
		win.snvmnv <- c(win.snvmnv,sum(chr.snvs[which(snv.coords>=windows[j]&snv.coords<=windows[j+1]),9]))
		win.meth <- c(win.meth,mean(chr.meth[which(chr.meth[,3]>=windows[j]&chr.meth[,4]<=windows[j+1]),13]))
	}
	left.wins <- windows[which(windows<centromeres[i])]
	left.prop <- left.wins/centromeres[i]
	left.snps <- win.snps[which(windows<centromeres[i])]
	left.meth <- win.meth[which(windows<centromeres[i])]
	right.arm <- chr.ends[i]-centromeres[i]	
	right.wins <- windows[which(windows>centromeres[i])]
	right.snps <- win.snps[which(windows>centromeres[i])]
	right.snps <- rev(right.snps)
	right.meth <- win.meth[which(windows>centromeres[i])]
	right.meth <- rev(right.meth)	
	right.wins <- right.wins-centromeres[i]
	right.prop <- right.wins/right.arm
	right.prop <- 1-right.prop
	right.prop <- rev(right.prop)
	chr.left <- cbind(rep(i,length(left.prop)),left.prop,left.snps,left.meth)
	chr.right <- cbind(rep(i,length(right.prop)),right.prop,right.snps,right.meth)
	snp.left <- rbind(snp.left,chr.left)
	snp.right <- rbind(snp.right,chr.right)
}
snp.all <- rbind(snp.left,snp.right)
snp.wins <- NULL
meth.wins <- NULL
for(j in 1:length(wins)-1){
	snp.wins <- c(snp.wins,mean(snp.all[which(snp.all[,2]>=wins[j]&snp.all[,2]<wins[j+1]),3]))
	meth.wins <- c(meth.wins,mean(snp.all[which(snp.all[,2]>=wins[j]&snp.all[,2]<wins[j+1]),4]))
}

hswt.left <- NULL
hswt.right <- NULL
for(i in 1:5){
	print(i)
	hswt.coords <- hswt.mask[which(hswt.mask[,3]==i),6]	
	left.hswt <- hswt.coords[which(hswt.coords<centromeres[i])]
	left.prop <- left.hswt/centromeres[i]
	right.arm <- chr.ends[i]-centromeres[i]	
	right.hswt <- hswt.coords[which(hswt.coords>centromeres[i])]
	right.hswt <- right.hswt-centromeres[i]
	right.prop <- right.hswt/right.arm
	right.prop <- 1-right.prop
	right.prop <- rev(right.prop)
	chr.left <- cbind(rep(i,length(left.prop)),left.prop)
	chr.right <- cbind(rep(i,length(right.prop)),right.prop)
	hswt.left <- rbind(hswt.left,chr.left)
	hswt.right <- rbind(hswt.right,chr.right)
}
hswt.all <- c(hswt.left[,2],hswt.right[,2])
hswt.wins <- NULL
for(j in 1:length(wins)-1){
	win.hswt <- length(which(hswt.all>=wins[j]&hswt.all<wins[j+1]))
	hswt.wins <- c(hswt.wins,win.hswt)
}
hswt.norm <- hswt.wins/hs.wt.n

cuhei10.left <- NULL
cuhei10.right <- NULL
for(i in 1:5){
	print(i)
	cuhei10.coords <- cuhei10.mask[which(cuhei10.mask[,3]==i),6]	
	left.cuhei10 <- cuhei10.coords[which(cuhei10.coords<centromeres[i])]
	left.prop <- left.cuhei10/centromeres[i]
	right.arm <- chr.ends[i]-centromeres[i]	
	right.cuhei10 <- cuhei10.coords[which(cuhei10.coords>centromeres[i])]
	right.cuhei10 <- right.cuhei10-centromeres[i]
	right.prop <- right.cuhei10/right.arm
	right.prop <- 1-right.prop
	right.prop <- rev(right.prop)
	chr.left <- cbind(rep(i,length(left.prop)),left.prop)
	chr.right <- cbind(rep(i,length(right.prop)),right.prop)
	cuhei10.left <- rbind(cuhei10.left,chr.left)
	cuhei10.right <- rbind(cuhei10.right,chr.right)
}
cuhei10.all <- c(cuhei10.left[,2],cuhei10.right[,2])
cuhei10.wins <- NULL
for(j in 1:length(wins)-1){
	win.cuhei10 <- length(which(cuhei10.all>=wins[j]&cuhei10.all<wins[j+1]))
	cuhei10.wins <- c(cuhei10.wins,win.cuhei10)
}
cuhei10.norm <- cuhei10.wins/cu.hei10.n

hsrecq.left <- NULL
hsrecq.right <- NULL
for(i in 1:5){
	print(i)
	hsrecq.coords <- hsrecq.mask[which(hsrecq.mask[,3]==i),6]	
	left.hsrecq <- hsrecq.coords[which(hsrecq.coords<centromeres[i])]
	left.prop <- left.hsrecq/centromeres[i]
	right.arm <- chr.ends[i]-centromeres[i]	
	right.hsrecq <- hsrecq.coords[which(hsrecq.coords>centromeres[i])]
	right.hsrecq <- right.hsrecq-centromeres[i]
	right.prop <- right.hsrecq/right.arm
	right.prop <- 1-right.prop
	right.prop <- rev(right.prop)
	chr.left <- cbind(rep(i,length(left.prop)),left.prop)
	chr.right <- cbind(rep(i,length(right.prop)),right.prop)
	hsrecq.left <- rbind(hsrecq.left,chr.left)
	hsrecq.right <- rbind(hsrecq.right,chr.right)
}
hsrecq.all <- c(hsrecq.left[,2],hsrecq.right[,2])
hsrecq.wins <- NULL
for(j in 1:length(wins)-1){
	win.hsrecq <- length(which(hsrecq.all>=wins[j]&hsrecq.all<wins[j+1]))
	hsrecq.wins <- c(hsrecq.wins,win.hsrecq)
}
hsrecq.norm <- hsrecq.wins/hs.recq.n

hsrecqhei10.left <- NULL
hsrecqhei10.right <- NULL
for(i in 1:5){
	print(i)
	hsrecqhei10.coords <- hsrecqhei10.mask[which(hsrecqhei10.mask[,3]==i),6]	
	left.hsrecqhei10 <- hsrecqhei10.coords[which(hsrecqhei10.coords<centromeres[i])]
	left.prop <- left.hsrecqhei10/centromeres[i]
	right.arm <- chr.ends[i]-centromeres[i]	
	right.hsrecqhei10 <- hsrecqhei10.coords[which(hsrecqhei10.coords>centromeres[i])]
	right.hsrecqhei10 <- right.hsrecqhei10-centromeres[i]
	right.prop <- right.hsrecqhei10/right.arm
	right.prop <- 1-right.prop
	right.prop <- rev(right.prop)
	chr.left <- cbind(rep(i,length(left.prop)),left.prop)
	chr.right <- cbind(rep(i,length(right.prop)),right.prop)
	hsrecqhei10.left <- rbind(hsrecqhei10.left,chr.left)
	hsrecqhei10.right <- rbind(hsrecqhei10.right,chr.right)
}
hsrecqhei10.all <- c(hsrecqhei10.left[,2],hsrecqhei10.right[,2])
hsrecqhei10.wins <- NULL
for(j in 1:length(wins)-1){
	win.hsrecqhei10 <- length(which(hsrecqhei10.all>=wins[j]&hsrecqhei10.all<wins[j+1]))
	hsrecqhei10.wins <- c(hsrecqhei10.wins,win.hsrecqhei10)
}
hsrecqhei10.norm <- hsrecqhei10.wins/hs.recqhei10.n

####################
# TEL-CEN plotting #
####################

par(mfcol=c(2,3))
par(mar=c(1.8,1.8,1.8,1.8))
df <- 40

# all recombination
ylim=c(0,0.6)
plot(wins,hswt.norm,type="n",ylim=ylim)
lines(smooth.spline(wins,hswt.norm,df=df),col=1)
abline(h=mean(hswt.norm),col=1,lty=2,lwd=0.5)
lines(smooth.spline(wins,cuhei10.norm,df=df),col=4)
abline(h=mean(cuhei10.norm),col=4,lty=2,lwd=0.5)
lines(smooth.spline(wins,hsrecq.norm,df=df),col="purple")
abline(h=mean(hsrecq.norm),col="purple",lty=2,lwd=0.5)
lines(smooth.spline(wins,hsrecqhei10.norm,df=df),col=2)
abline(h=mean(hsrecqhei10.norm),col=2,lty=2,lwd=0.5)

# meth vs snps
meth.wins[1] <- 0
meth.wins[2] <- 0
plot(wins,meth.wins,type="n",col=4)
lines(smooth.spline(wins,meth.wins,df=df),col=4)
abline(h=mean(meth.wins),lty=2,lwd=0.5,col=4)
par(new=T)
ylim=c(100,400)
snp.wins[1] <- 0
plot(wins,snp.wins,type="n",ylim=ylim,yaxt="n",col="dark green")
lines(smooth.spline(wins,snp.wins,df=df),col="dark green")
axis(side=4)
abline(h=mean(snp.wins),lty=2,lwd=0.5,col="dark green")

# hs.wt vs snps
ylim=c(0,0.14)
plot(wins,hswt.norm,type="n",ylim=ylim)
lines(smooth.spline(wins,hswt.norm,df=df),col=1)
abline(h=mean(hswt.norm),lty=2,lwd=0.5,col=1)
par(new=T)
ylim=c(100,400)
snp.wins[1] <- 0
plot(wins,snp.wins,type="n",ylim=ylim,yaxt="n",col="dark green")
lines(smooth.spline(wins,snp.wins,df=df),col="dark green")
axis(side=4)
abline(h=mean(snp.wins),lty=2,lwd=0.5,col="dark green")

# hs.wt vs meth
ylim=c(0,0.14)
plot(wins,hswt.norm,type="n",ylim=ylim)
lines(smooth.spline(wins,hswt.norm,df=df),col=1)
abline(h=mean(hswt.norm),lty=2,lwd=0.5,col=1)
par(new=T)
meth.wins[1] <- 0
meth.wins[2] <- 0
plot(wins,meth.wins,type="n",col=4,yaxt="n")
lines(smooth.spline(wins,meth.wins,df=df),col=4)
axis(side=4)
abline(h=mean(meth.wins),lty=2,lwd=0.5,col=4)

# all.cos vs snps
ylim=c(0,0.6)
plot(wins,hswt.norm,type="n",ylim=ylim)
lines(smooth.spline(wins,hswt.norm,df=df),col=1)
abline(h=mean(hswt.norm),col=1,lty=2,lwd=0.5)
lines(smooth.spline(wins,cuhei10.norm,df=df),col=4)
abline(h=mean(cuhei10.norm),col=4,lty=2,lwd=0.5)
lines(smooth.spline(wins,hsrecq.norm,df=df),col="purple")
abline(h=mean(hsrecq.norm),col="purple",lty=2,lwd=0.5)
lines(smooth.spline(wins,hsrecqhei10.norm,df=df),col=2)
abline(h=mean(hsrecqhei10.norm),col=2,lty=2,lwd=0.5)
par(new=T)
ylim=c(100,400)
snp.wins[1] <- 0
plot(wins,snp.wins,type="n",ylim=ylim,yaxt="n",col="dark green")
lines(smooth.spline(wins,snp.wins,df=df),col="dark green")
axis(side=4)
abline(h=mean(snp.wins),lty=2,lwd=0.5,col="dark green")

# all.cos vs meth
ylim=c(0,0.6)
plot(wins,hswt.norm,type="n",ylim=ylim)
lines(smooth.spline(wins,hswt.norm,df=df),col=1)
abline(h=mean(hswt.norm),col=1,lty=2,lwd=0.5)
lines(smooth.spline(wins,cuhei10.norm,df=df),col=4)
abline(h=mean(cuhei10.norm),col=4,lty=2,lwd=0.5)
lines(smooth.spline(wins,hsrecq.norm,df=df),col="purple")
abline(h=mean(hsrecq.norm),col="purple",lty=2,lwd=0.5)
lines(smooth.spline(wins,hsrecqhei10.norm,df=df),col=2)
abline(h=mean(hsrecqhei10.norm),col=2,lty=2,lwd=0.5)
par(new=T)
meth.wins[1] <- 0
meth.wins[2] <- 0
plot(wins,meth.wins,type="n",col=4,yaxt="n")
lines(smooth.spline(wins,meth.wins,df=df),col=4)
axis(side=4)
abline(h=mean(meth.wins),lty=2,lwd=0.5,col=4)

######################################################
# testing for crossovers in the knob - all give zero #
######################################################

hswt.mask <- mask.cos[which(mask.cos[,8]=="hs.wt"),]
cuhei10.mask <- mask.cos[which(mask.cos[,8]=="cu.hei10"),]
hsrecq.mask <- mask.cos[which(mask.cos[,8]=="hs.recq"),]
hsrecqhei10.mask <- mask.cos[which(mask.cos[,8]=="hs.recqhei10"),]

# testing for crossovers in the knob - all give zero

chr.dat <- hswt.mask[which(hswt.mask[,3]==4),]
which(chr.dat[,4]>=1612606&chr.dat[,4]<=2782625)
which(chr.dat[,5]>=1612606&chr.dat[,5]<=2782625)

chr.dat <- cuhei10.mask[which(cuhei10.mask[,3]==4),]
which(chr.dat[,4]>=1612606&chr.dat[,4]<=2782625)
which(chr.dat[,5]>=1612606&chr.dat[,5]<=2782625)

chr.dat <- hsrecq.mask[which(hsrecq.mask[,3]==4),]
which(chr.dat[,4]>=1612606&chr.dat[,4]<=2782625)
which(chr.dat[,5]>=1612606&chr.dat[,5]<=2782625)

chr.dat <- hsrecqhei10.mask[which(hsrecqhei10.mask[,3]==4),]
which(chr.dat[,4]>=1612606&chr.dat[,4]<=2782625)
which(chr.dat[,5]>=1612606&chr.dat[,5]<=2782625)

#######
# all #
#######

coller.invs <- read.csv(file="ColLer_Inversions.csv")
tot.cos <- length(hswt.mask[,1])+length(cuhei10.mask[,1])+length(hsrecq.mask[,1])+length(hsrecqhei10.mask[,1])

library(segmentSeq)
library(GenomicAlignments)
chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chr.lens <- c(30427671,19698289,23459830,18585056,26975502)
all.cos <- NULL
all.ran <- NULL
for(i in 1:5){
	print(i)
	chr.invs <- coller.invs[which(coller.invs[,1]==chrs[i]),]
	cos1 <- hswt.mask[which(hswt.mask[,3]==i),]
	cos2 <- cuhei10.mask[which(cuhei10.mask[,3]==i),]
	cos3 <- hsrecq.mask[which(hsrecq.mask[,3]==i),]
	cos4 <- hsrecqhei10.mask[which(hsrecqhei10.mask[,3]==i),]
	chr.cos <- rbind(cos1,cos2,cos3,cos4)
	cos.gr.coords <- GRanges(seqnames=chrs[i],strand="+",ranges=IRanges(start=chr.cos$start,width=chr.cos$width+1))
	inv.gr.coords <- GRanges(seqnames=chrs[i],strand="+",ranges=IRanges(start=chr.invs$Start,width=chr.invs$Length+1))	
	ran.start <- round(runif(n=length(cos.gr.coords),min=1,max=chr.lens[i]-10000))	
	ran.end <- ran.start+width(cos.gr.coords)
	ran.gr.coords <- GRanges(seqnames=chrs[i],strand="+",ranges=IRanges(start=ran.start,end=ran.end))	
	cos.overlap <- findOverlaps(inv.gr.coords,cos.gr.coords,type="any")
	ran.overlap <- findOverlaps(inv.gr.coords,ran.gr.coords,type="any")
	all.cos <- c(all.cos,length(cos.overlap))
	all.ran <- c(all.ran,length(ran.overlap))	
}

all.cos <- sum(all.cos)
all.ran <- sum(all.ran)
cos.dat <- c(all.cos,tot.cos-all.cos)
ran.dat <- c(all.ran,tot.cos-all.ran)
all.dat <- cbind(cos.dat,ran.dat)
chisq.test(all.dat)
> <2.2e-16



