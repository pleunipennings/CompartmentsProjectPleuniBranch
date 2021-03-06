#compile and move program
system("./make")
system("mv main ../main")

#set parameters
R = 4
c1 = 0.99; c2 = 0.99; s=0.01; dy=1; dx=0.1; 
m="2E-03";
#m=0;
uf1 = "1E-05"; ub1 = "1E-04"; uf2 = "0" ; ub2 = "0";
x0 = "1E+04"; x1 = "1E+06"; x2 = 0; x3 = 0; #the number of uninfected cells in the compartments
seed = 1
nr = 1; #number of runs
tcf = 10; #threshold of when to stop simulations
numruns=500
rununtilfailureOrngens=0; #if 0 run until failure if >0 run this num of generations
outputmutantinfo=0; #if 1 output info on mutants for debugging, if 0 output only migtimes

migtimes=vector()

#remove old output files
if (TRUE){
if (file.exists("../firstfilled.txt")) system ("rm ../firstfilled.txt")
if (file.exists("../texit.txt")) system ("rm ../texit.txt")
if (file.exists("../path_sucmut.txt")) system ("rm ../path_sucmut.txt")
if (file.exists("../tfailure.txt")) system ("rm ../tfailure.txt")
if (file.exists("../vload.txt")) system ("rm ../vload.txt")
if (file.exists("../path_sdc.txt")) system ("rm ../path_sdc.txt")
if (file.exists("../mmt.txt")) system ("rm ../mmt.txt")
if (file.exists("../output.txt")) system ("rm ../output.txt")
}
	
#run simulation
WD <- getwd()
setwd("/Users/pleunipennings/Documents/Research/HIV/CompartmentsProject/CodePleuniBranch")
#pdf("fig.pdf")

for (seed in 1:numruns){
	system (paste("rm ./output",seed,".txt",sep=""))
	
#create parameter file
	filetowrite="./parameters.txt"
	write(paste("R:",R),file=filetowrite)
	write(paste("c1:",c1),file=filetowrite,append=TRUE)
	write(paste("c2:",c2),file=filetowrite,append=TRUE)
	write(paste("s:",s),file=filetowrite,append=TRUE)
	write(paste("dy:",dy),file=filetowrite,append=TRUE)
	write(paste("dx:",dx),file=filetowrite,append=TRUE)
	write(paste("m:",m),file=filetowrite,append=TRUE)
	write(paste("uf1:",uf1),file=filetowrite,append=TRUE)
	write(paste("ub1:",ub1),file=filetowrite,append=TRUE)
	write(paste("uf2:",uf2),file=filetowrite,append=TRUE)
	write(paste("ub2:",ub2),file=filetowrite,append=TRUE)
	write(paste("x0:",x0),file=filetowrite,append=TRUE)
	write(paste("x1:",x1),file=filetowrite,append=TRUE)
	write(paste("x2:",x2),file=filetowrite,append=TRUE)
	write(paste("x3:",x3),file=filetowrite,append=TRUE)
	write(paste("nr:",nr),file=filetowrite,append=TRUE)
	write(paste("tcf:",tcf),file=filetowrite,append=TRUE)
	write(paste("seed:",seed),file=filetowrite,append=TRUE)
	write(paste("ruf:",rununtilfailureOrngens),file=filetowrite,append=TRUE)
	write(paste("omi:",outputmutantinfo),file=filetowrite,append=TRUE)
	
#	print(system.time(system (paste("./main >> ./output",seed,".txt",sep=""))))
	print(paste("sim",seed,"out of",numruns))
	system (paste("./main >> ./output",seed,".txt",sep=""))
#read output files
#scan("tfailure.txt")->FailureTimes

#Data<-data.frame("Run" = 1:length(FailureTimes),"FailureTimes" = FailureTimes)
if(FALSE){
	png("FailureTimes.png")

library(MASS)
x.est <- fitdistr(Data$FailureTimes, "exponential")$estimate

hist(Data$FailureTimes,ylim=c(0,0.001),freq=FALSE, breaks=seq(0,90000,by=50),xlim=c(0,10000))
curve(dexp(x, rate = x.est), add = TRUE, col = "red", lwd = 2)
curve(dexp(x, rate = 0.001095), add = TRUE, col = "green", lwd = 2)

dev.off()
}

#migtimes<-c(migtimes,scan(paste("./output",seed,".txt",sep="")))
#	hist(x,main=seed)

}	
#dev.off()

#read.table("output.txt",sep="\t")->x
#x[,(length(x[1,])+1)]<-0;
#x[,(length(x[1,])+1)]<-0;
#names(x)<-c("i","j","k","t","WT","mutsize","deltat","mutsizedeltat")
#x$deltat[2:length(x[,1])]<-x$t[2:length(x[,1])]-x$t[1:(length(x[,1])-1)]
#get average size of mutant pop
#x$mutsizedeltat=x$mutsize*x$deltat
#print("mean number of mutants")
#print(sum(x$mutsizedeltat)/max(x$t))

if (!is.null(WD)) setwd(WD)

if (FALSE){
pdf("MutSelBal.pdf")
read.table("../output1.txt",sep="\t")->x1
names(x1)<-c("i","j","k","t","wt","mut")
#plot(x1$t,x1$wt,ylim=c(0,max(x1$wt)*1.2))
plot(x1$t,x1$mut,ylim=c(0,max(c(x1$mut,x1$mut))),t="l",lwd=2)
abline(v=x1$t[which(x1$k==2)],lty=2,lwd=0.5)

migevents<-vector()
for (i in 2:numruns){
	read.table(paste("../output",i,".txt",sep=""))->x
	names(x)<-c("i","j","k","t","wt","mut")
	points(x$t,x$mut,col=i,t="l",lwd=2)
	abline(v=x$t[which(x$k==2)],col=i,lty=2,lwd=0.5)
	migevents<-c(migevents,x$t[which(x$k==2)])
}
	
	dev.off()
	s=0
	for (i in 1:numruns){
		read.table(paste("../output",i,".txt",sep=""))->x
		names(x)<-c("i","j","k","t","wt","mut")
		mean1=mean(x$mut[1:floor(length(x$mut)/2)])
		mean2=mean(x$mut[ceiling(length(x$mut)/2):length(x$mut)])
		print(paste(i,round(mean1,2),round(mean2,2),length(which(x$k==2))))
		if (mean1>mean2)s=s+1
	}
	print(s)
	
}

if (FALSE){
pdf("MigrationTimesIn4000runs.pdf")
plot(sort(migtimes))
abline(v=length(migtimes)/2)
abline(v=length(migtimes)/4)
abline(v=length(migtimes)*3/4)
abline(h=c(250,500,750))
dev.off()
}

for (i in 1:numruns){
system (paste("rm ../output",i,".txt",sep=""))
}

scan("../tfailure.txt")->FailureTimes
hist(FailureTimes)
print(paste("mean Failure time",round(mean(FailureTimes))))

library(MASS)
x.est <- fitdistr(FailureTimes, "exponential")$estimate

predictedrate=as.numeric(x0)*0.075*(as.numeric(uf1)/s)*as.numeric(m)*(R-1)/R
predictedrate=predictedrate/10

if (FALSE){
pdf(paste("mu",uf1,"s",s,"R_1.081.pdf",sep=""))
hist(FailureTimes,ylim=c(0,0.001),freq=FALSE, breaks=seq(0,90000,by=50),xlim=c(0,20000))
abline(v=mean(FailureTimes),col=2)
abline(v=1/predictedrate,col=3)
curve(dexp(x, rate = x.est), add = TRUE, col = "red", lwd = 2)
curve(dexp(x, rate = predictedrate), add = TRUE, col = "green", lwd = 2)
dev.off()
}