#compile and move program
system("./make")
system("mv main ../main")

#set parameters
R = 4
c1 = 0.99; c2 = 0.99; s=0.05; dy=1; dx=0.1; m="1E-03";
uf1 = "1E-04"; ub1 = "0"; uf2 = "0" ; ub2 = "0";
x0 = "1E+03"; x1 = "1E+06"; x2 = 0; x3 = 0; #the number of uninfected cells in the compartments
seed = 1
nr = 1; #number of runs
tcf = 10; #threshold of when to stop simulations

#create parameter file
filetowrite="../parameters.txt"
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

for (i in 1:2){
print(system.time(system ("../main >> ../output.txt")))

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

	scan("../output.txt")->x
	hist(x)

#read.table("output.txt",sep="\t")->x
#x[,(length(x[1,])+1)]<-0;
#x[,(length(x[1,])+1)]<-0;
#names(x)<-c("i","j","k","t","WT","mutsize","deltat","mutsizedeltat")
#x$deltat[2:length(x[,1])]<-x$t[2:length(x[,1])]-x$t[1:(length(x[,1])-1)]
#get average size of mutant pop
#x$mutsizedeltat=x$mutsize*x$deltat
#print("mean number of mutants")
#print(sum(x$mutsizedeltat)/max(x$t))
}

