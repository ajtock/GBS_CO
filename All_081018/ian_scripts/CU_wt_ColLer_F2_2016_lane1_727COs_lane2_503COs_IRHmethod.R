
inDir <- "/projects/ajt200/GBS_CO/All_081018/cos_files/"

lib.nums <- seq(1,96,by=1)
CU_lane1.cos <- NULL
for(k in 1:length(lib.nums)){
        print(k)
        dat <- read.table(file=paste0(inDir,"n96.wtF2",k,".bchqsnvmask.smooth.co.txt",sep=""),header=F)
        all <- NULL
        for(i in 1:5){
                print(i)
                chr <- dat[which(dat[,2]==i),]
                if((length(chr[,1])>1)==T){
                        start <- chr[,4]
                        stop <- chr[,3]
                        stop <- stop[-1]
                        start <- start[-length(start)]
                        diff <- stop-start
                        mid <- round(diff/2)
                        cos <- start+mid
                        chrs <- rep(i,length(cos))
                        lib <- rep(paste("1.",lib.nums[k],sep=""),length(chrs))
                        width <- stop-start
                        bind <- cbind(lib,chrs,start,stop,cos,width)
                        all <- rbind(all,bind)
                        print(dim(all))
                        }
                }
        CU_lane1.cos <- rbind(CU_lane1.cos,all)
}

lib.nums <- seq(1,96,by=1)
CU_lane2.cos <- NULL
for(k in 1:length(lib.nums)){
        print(k)
        dat <- read.table(file=paste0(inDir,"n96.wt.lane2.",k,".bchqsnvmask.smooth.co.txt",sep=""),header=F)
        all <- NULL
        for(i in 1:5){
                print(i)
                chr <- dat[which(dat[,2]==i),]
                if((length(chr[,1])>1)==T){
                        start <- chr[,4]
                        stop <- chr[,3]
                        stop <- stop[-1]
                        start <- start[-length(start)]
                        diff <- stop-start
                        mid <- round(diff/2)
                        cos <- start+mid
                        chrs <- rep(i,length(cos))
                        lib <- rep(paste("2.",lib.nums[k],sep=""),length(chrs))
                        width <- stop-start
                        bind <- cbind(lib,chrs,start,stop,cos,width)
                        all <- rbind(all,bind)
                        print(dim(all))
                        }
                }
        CU_lane2.cos <- rbind(CU_lane2.cos,all)
}
