big_table<- "summary/bins_summary.txt"

data<-read.table(big_table,h=T)
data[,2]<-as.character(data[,2])
data[,3]<-as.character(data[,3])
#collect_pi
pigot<-c()
nsitesgot<-c()
piol<-c()
nsitesol<-c()

for (row_number in c(1:dim(data)[1])){
	line<-data[row_number,]
	#gotland pi
	filename<-paste("pi_bigpops/bins/gotland95/",as.character(line[3]),"/bin",as.character(line[1]),"gotland.bed",sep="")#gotland pi
	tempdata<-read.table(filename)
	print(c(row_number,filename))
	pigot[length(pigot)+1]<-sum(tempdata[,5],na.rm=T)/ sum(tempdata[,6],na.rm=T)
	nsitesgot[length(nsitesgot)+1]<-sum(tempdata[,6],na.rm=T)
	#oland
	filename<-paste("pi_bigpops/bins/oland95/",as.character(line[3]),"/bin",as.character(line[1]),"oland.bed",sep="")#gotland pi
	tempdata<-read.table(filename)
	print(c(row_number,filename))
	piol[length(piol)+1]<-sum(tempdata[,5],na.rm=T)/ sum(tempdata[,6],na.rm=T)
	nsitesol[length(nsitesol)+1]<-sum(tempdata[,6],na.rm=T)
}


mean_pi<-mean(cbind(piol,pigot))
mean_nsites<-mean(cbind(nsitesol,nsitesgot))
min_nsites<-min(cbind(nsitesol,nsitesgot))
data<-cbind(data,pigot,nsitesgot,piol,nsitesol,mean_pi,mean_nsites,min_nsites)

#real fst gotland 94

fst_real<-c()



for (row_number in c(1:dim(data)[1])){
	line<-data[row_number,]
	#gotland fst
	filename<-paste("/home/ludovic/nobackup_ludo/SBE/fsts/MF_fst/logsreal",as.character(line[3]),"/bin",as.character(line[1]),".bed",sep="")#gotland pi
	print (filename)
	system(paste("cat",filename,"| grep 'Weir and Cockerham mean Fst estimate' | cut -f 7 -d ' ' > temp"))
	tempdata<-read.table("temp")
	#print(c(row_number,filename))
	fst_real[length(fst_real)+1]<-as.numeric(tempdata)
}



#real fst gotland 94

fst_boot<-c()



for (row_number in c(1:dim(data)[1])){
	line<-data[row_number,]
	#gotland fst
	filename<-paste("/home/ludovic/nobackup_ludo/SBE/fsts/concatenated/logsboot",as.character(line[3]),"/bin",as.character(line[1]),".bed",sep="")#gotland pi
	print (filename)
	system(paste("cat",filename,"| grep 'Weir and Cockerham mean Fst estimate' | cut -f 7 -d ' ' > temp"))
}	
	tempdata<-read.table("temp") # should be 12000 values.. take 10000 non nas and then extract the percentile
	#print(c(row_number,filename))
	fst_real[length(fst_real)+1]<-as.numeric(tempdata)
}



