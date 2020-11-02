source("~/Dropbox/HingeModels/molt.simulator.R")


require(moult)
library(dplyr)
require(chngpt)
require(reshape)


# import simulated data

simsAll<-read.csv("~/Dropbox/HingeModels/simulatedData.csv")

#I'll want three fit functions . fit.DH, fit.UZ, and fit.Pimm

#then I'll want three evaluate functions, eval.DH, eval.UZ, and eval.Pimm
#this takes 100 simuations of various n's


fit.UZ<-function(sims,samples){
	for(i in 1:100){
		sims[[i]]$pscore<-sims[[i]]$pscore/100
		sims[[i]]$V3<-factor(sims[[i]]$V3,levels=samples)
			}
		fits<-list()
		fit<-list()

		for(j in 1:100){
		for(i in factor(samples)){
			fit[[i]]<-NA
			tryCatch({
		fit[[i]]<-moult(pscore~julian,data=sims[[j]][sims[[j]]$V3==i,])	
			fits[[j]]<-fit
			}, error=function(e){})
			}
		}
	ci<-list()
	ci_end<-list()
	cis<-list()
	cis_end<-list()
	for(j in 1:length(fits)){
    	for(i in 1:length(samples)){
    		tryCatch({
        		fit=fits[[j]][[i]]
        		se=sqrt(diag(vcov(fit)))
        		ci[[i]]<-cbind(coef(fit)-1.96*se, coef(fit)+1.96*se)
        		comb=c(1,1,0)
        		est =  comb%*%coef(fit)
        		se_end=sqrt(comb%*%vcov(fit)%*%comb)     
       	 ci_end[[i]]<-cbind(as.numeric(est-1.96*se_end), as.numeric(est+1.96*se_end))
             			}, error=function(e){})
    				}
		names(ci)<-names(fits[[1]])
    	cis[[j]]<-ci
    	names(ci_end)<-names(fits[[1]])
    	cis_end[[j]]<-ci_end
						}
	UZ.mat<-list()
	for (i in 1:100){
		UZ.mat[[i]]<-matrix(nrow=length(samples),ncol=10)
		rownames(UZ.mat[[i]])<-as.character(samples)
		colnames(UZ.mat[[i]])<-c("duration_min","duration_max","start_min","start_max","SD_min","SD_max","term_min","term_max","init","term")
					}
					newData<-data.frame(pscore=c(0,100))
	for(j in 1:100){
	for(i in 1:length(samples)){
				tryCatch({
    		 	UZ.mat[[j]][i,1]<-cis[[j]][[i]][,1][1]
      			UZ.mat[[j]][i,2]<-cis[[j]][[i]][,2][1]
     			UZ.mat[[j]][i,3]<-cis[[j]][[i]][,1][2]
     			UZ.mat[[j]][i,4]<-cis[[j]][[i]][,2][2]
     			UZ.mat[[j]][i,5]<-cis[[j]][[i]][,1][3]
     			UZ.mat[[j]][i,6]<-cis[[j]][[i]][,2][3]
     			UZ.mat[[j]][i,7]<-cis_end[[j]][[i]][,1]
     			UZ.mat[[j]][i,8]<-cis_end[[j]][[i]][,2]
     			UZ.mat[[j]][i,9]<-unlist(fits[[j]][[i]])$`optim.par.(Intercept)`
     			UZ.mat[[j]][i,10]<-unlist(fits[[j]][[i]])$`optim.par.(Intercept)`+unlist(fits[[j]][[i]])$optim.par1
     			        		}, error=function(e){})
					}
		UZ.mat[[j]]<-data.frame(UZ.mat[[j]])
		UZ.mat[[j]]$n<-rownames(UZ.mat[[j]])	
						}

		UZ.results<-bind_rows(UZ.mat)

	return(UZ.results)
}



fit.DH<-function(sims,samples){

	hinge<-list()
	hinges<-list()
	for(i in 1:100){
		hinges[[i]]<-list()
					}
	for(j in 1:100){
		cat("\n")
		cat("Sim",j,"/100","\n")
	for(i in unique(sims[[j]]$V3)){
		tryCatch({
			cat(i," ")
     		 hinges[[j]][i]<-summary(double.hinge(x=sims[[j]][sims[[j]]$V3==i,]$julian,y=sims[[j]][sims[[j]]$V3==i,]$pscore, lower.y=0, upper.y=100, var.type="bootstrap"))
  						})
					}
				}
			hinge.mats<-list()

		for (i in 1:length(hinges)){
			hinge.mats[[i]]<-matrix(nrow=length(unique(sims[[1]]$V3)),ncol=7)
			colnames(hinge.mats[[i]])<-c("n","startDate","endDate","startmin","startmax","endmin","endmax")
				for (j in 1: nrow(hinge.mats[[i]])){
					hinge.mats[[i]][j,2]<-unlist(hinges[[i]][j])[[1]]
					hinge.mats[[i]][j,3]<-unlist(hinges[[i]][j])[[2]]
					hinge.mats[[i]][j,4]<-unlist(hinges[[i]][j])[[7]]
					hinge.mats[[i]][j,5]<-unlist(hinges[[i]][j])[[10]]
					hinge.mats[[i]][j,6]<-unlist(hinges[[i]][j])[[8]]
					hinge.mats[[i]][j,7]<-unlist(hinges[[i]][j])[[11]]
			}
			hinge.mats[[i]]<-data.frame(hinge.mats[[i]])
			hinge.mats[[i]]$n<-samples
			hinge.mats[[i]]$correct<-NA
			hinge.mats[[i]]$termCorrect<-NA

		}

		hinge.mat<-bind_rows(hinge.mats)
		return(hinge.mat)
	}




fit.Pimm<-function(sims,samples){
	for(i in 1:length(sims)){
		sims[[i]]<-sims[[i]][!sims[[i]]$pscore==0,]
		sims[[i]]<-sims[[i]][!sims[[i]]$pscore==100,]
	}
	
	pimm.mat<-matrix(nrow=length(samples),ncol=4)
	rownames(pimm.mat)<-unique(samples)
	colnames(pimm.mat)<-c("end_min","end_max","start_min","start_max")
	pimm.mat<-data.frame(pimm.mat)
	pimm.mat$n<-samples
	fit<-list()
	fits<-list()
	for(k in factor(samples)){
	fit[[k]]<-NA}
	for(i in 1:100){
		{
		fits[[i]]<-fit
	}
}
	for(j in 1:100){
		cat("\n")
		cat("Sim",j,"/100","\n")
	for(i in unique(sims[[j]]$V3)){
		tryCatch({
			cat(i," ")
     		 fits[[j]][[i]]<-lm(julian~pscore,data=sims[[j]][sims[[j]]$V3==i,])
     		 fits[[j]][[i]]$terms
  						})
					}
				}
	ci<-list()
	cis<-list()
	for(j in 1:length(fits)){
    	for(i in 1:length(samples)){
    		newData<-data.frame(pscore=c(0,100))
tryCatch({
    		ci[[i]]<-predict(fits[[j]][[i]], newdata = newData, interval = "confidence")
    		             			}, error=function(e){})
    			}
    	cis[[j]]<-ci
						}

	Pimm.mat<-list()
	for (i in 1:100){
		Pimm.mat[[i]]<-matrix(nrow=length(samples),ncol=6)
		rownames(Pimm.mat[[i]])<-as.character(samples)
		colnames(Pimm.mat[[i]])<-c("end_min","end_max","start_min","start_max","start","term")
					}
	for(j in 1:100){
	for(i in 1:length(samples)){
				tryCatch({
    		 	Pimm.mat[[j]][i,1]<-cis[[j]][[i]][,2][2]
      			Pimm.mat[[j]][i,2]<-cis[[j]][[i]][,3][2]
     			Pimm.mat[[j]][i,3]<-cis[[j]][[i]][,2][1]
     			Pimm.mat[[j]][i,4]<-cis[[j]][[i]][,3][1]
     			Pimm.mat[[j]][i,5]<-cis[[j]][[i]][,1][1]
     			Pimm.mat[[j]][i,6]<-cis[[j]][[i]][,1][2]

     			        		}, error=function(e){})
					}
		Pimm.mat[[j]]<-data.frame(Pimm.mat[[j]])
		Pimm.mat[[j]]$n<-rownames(Pimm.mat[[j]])	
						}

		Pimm.results<-bind_rows(Pimm.mat)

	return(Pimm.results)

	}



evaluate.UZ<-function(UZ.results,init,duration,samples){

		UZ.results<-bind_rows(UZ.results)

		if(length(duration)==1){

		UZ.results$correct<-UZ.results$start_min<init+1&UZ.results$start_max>init-1
		UZ.results$termCorrect<-UZ.results$term_min<(init+duration)-1&UZ.results$term_max>(init+duration)+1
		} else {
		for (i in duration){
		UZ.results$correct<-UZ.results$start_min<init+1&UZ.results$start_max>init-1
		UZ.results$termCorrect[UZ.results$n==i]<-UZ.results$term_min[UZ.results$n==i]<(init+as.numeric(i))-1&UZ.results$term_max[UZ.results$n==i]>(init+as.numeric(i))+1
			}
		}

		#count up corrects by n

		Corrects<-data.frame(matrix(nrow=length(unique(UZ.results$n)),ncol=6))
		colnames(Corrects)<-c("Correct","Incorrect","termCorrect","termIncorrect","model","n")
		Corrects$n<-samples
		Corrects$model<-"UZ"

				for (i in unique(UZ.results$n)){
				tryCatch({
					if(is.na(table(UZ.results$correct[UZ.results$n==i])[2])){
					Corrects$Correct[Corrects$n==i]<-100
					Corrects$Incorrect<-0}else{	
			Corrects$Correct[Corrects$n==i]<-table(UZ.results$correct[UZ.results$n==i])[2]
			Corrects$Incorrect[Corrects$n==i]<-table(UZ.results$correct[UZ.results$n==i])[1]
		}
		if(is.na(table(UZ.results$termCorrect[UZ.results$n==i])[2])){
			Corrects$termCorrect[Corrects$n==i]<-100
			Corrects$termIncorrect[Corrects$n==i]<-0}else{

			Corrects$termCorrect[Corrects$n==i]<-table(UZ.results$termCorrect[UZ.results$n==i])[2]
			Corrects$termIncorrect[Corrects$n==i]<-table(UZ.results$termCorrect[UZ.results$n==i])[1]
		}
					}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
		}





		Corrects[is.na(Corrects)]<-0
		Corrects$noModel<-100-(as.numeric(Corrects$Correct)+as.numeric(Corrects$Incorrect))
		Corrects$TermnoModel<-100-(as.numeric(Corrects$termCorrect)+as.numeric(Corrects$termIncorrect))
		Corrects[,1]<-as.numeric(Corrects[,1])
		Corrects[,2]<-as.numeric(Corrects[,2])
		Corrects[,3]<-as.numeric(Corrects[,3])
		Corrects[,4]<-as.numeric(Corrects[,4])
		Corrects[,7]<-as.numeric(Corrects[,7])




		Corrects.melt<-melt(Corrects,id=c("n","model"))
		Corrects.melt$value<-as.numeric(as.character(Corrects.melt$value))/100
		Corrects.melt$param<-"Initiation"
		Corrects.melt$param[25:48]<-"Termination"
		Corrects.melt$param[61:72]<-"Termination"
		Corrects.melt$param<-factor(Corrects.melt$param)
		Corrects.melt$variable[25:36]<-"Correct"
		Corrects.melt$variable[37:48]<-"Incorrect"
		Corrects.melt$variable[61:72]<-"noModel"

#Maybe also return distance between known and modelled parameters


return(Corrects.melt)
}




evaluate.DH<-function(DH.results,init,duration,samples){

		hinge.mat<-DH.results

if(length(duration)==1){



					hinge.mat$correct<-hinge.mat$startmin<init+1&hinge.mat$startmax>init-1
					hinge.mat$termCorrect<-hinge.mat$endmin<(init+duration+1)&hinge.mat$endmax>(init+duration-1)
					}else{

						for(i in duration){
					hinge.mat$correct<-hinge.mat$startmin<init+1&hinge.mat$startmax>init-1
					hinge.mat$termCorrect[hinge.mat$n==i]<-hinge.mat$endmin[hinge.mat$n==i]<(init+as.numeric(i)+1)&hinge.mat$endmax[hinge.mat$n==i]>(init+as.numeric(i)-1)
					}}

				

		hinges.all<-hinge.mat

		Corrects.hinges<-data.frame(matrix(nrow=length(unique(hinges.all$n)),ncol=6))
		colnames(Corrects.hinges)<-c("Correct","Incorrect","model","n","termCorrect","termIncorrect")
		Corrects.hinges$n<-samples



		Corrects.hinges$model<-"hinge"

		Corrects.hinges$model<-"hinge"

		for (i in unique(hinges.all$n)){
				tryCatch({
			num<-as.numeric(i)
			Corrects.hinges$Correct[Corrects.hinges$n==i]<-sum(hinges.all$correct[hinges.all$n==i]==TRUE)
			Corrects.hinges$Incorrect[Corrects.hinges$n==i]<-sum(hinges.all$correct[hinges.all$n==i]==FALSE)
			Corrects.hinges$termCorrect[Corrects.hinges$n==i]<-sum(hinges.all$termCorrect[hinges.all$n==i]==TRUE)
			Corrects.hinges$termIncorrect[Corrects.hinges$n==i]<-sum(hinges.all$termCorrect[hinges.all$n==i]==FALSE)

					}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
		}

		Corrects.hinges$"NA"<-"0"
		Corrects.hinges$noModel<-"0"
		Corrects.hinges$TermnoModel<-"0"

		CA.n<-Corrects.hinges
		dur.melt<-melt(CA.n,id=c("n","model"))


		dur.melt$value<-as.numeric(as.character(dur.melt$value))/100
		dur.melt$param<-"Initiation"


		dur.melt$variable[25:36]<-"Correct"
		dur.melt$variable[37:48]<-"Incorrect"
		dur.melt$variable[49:84]<-"noModel"
		dur.melt$param[25:60]<-"Termination"

return(dur.melt)


		}


evaluate.Pimm<-function(Pimm.results,init,duration,samples){



if(length(duration)==1){

			Pimm.results$correct<-Pimm.results$start_min<init+1&Pimm.results$start_max>init-1
		Pimm.results$termCorrect<-Pimm.results$end_min<(init+duration)+1&Pimm.results$end_max>(init+duration)-1

		}else{
			for(i in duration){
					Pimm.results$correct<-Pimm.results$start_min<init+1&Pimm.results$start_max>init-1
					Pimm.results$termCorrect[Pimm.results$n==i]<-Pimm.results$end_min[Pimm.results$n==i]<(init+as.numeric(i))+1&Pimm.results$end_max[Pimm.results$n==i]>(init+as.numeric(i)-1)
				}
		}




		#count up corrects by n

		Corrects<-data.frame(matrix(nrow=length(unique(Pimm.results$n)),ncol=6))
		colnames(Corrects)<-c("Correct","Incorrect","termCorrect","termIncorrect","model","n")
		Corrects$n<-samples
		Corrects$model<-"Pimm"

		for (i in unique(Pimm.results$n)){
				tryCatch({
					if(is.na(table(Pimm.results$correct[Pimm.results$n==i])[2])){
					Corrects$Correct[Corrects$n==i]<-table(Pimm.results$correct[Pimm.results$n==i])[1]
					Corrects$Incorrect<-0}else{	
			Corrects$Correct[Corrects$n==i]<-table(Pimm.results$correct[Pimm.results$n==i])[2]
			Corrects$Incorrect[Corrects$n==i]<-table(Pimm.results$correct[Pimm.results$n==i])[1]
		}
		if(is.na(table(Pimm.results$termCorrect[Pimm.results$n==i])[2])){
			Corrects$termCorrect[Corrects$n==i]<-table(Pimm.results$termCorrect[Pimm.results$n==i])[1]
			Corrects$termIncorrect[Corrects$n==i]<-0}else{

			Corrects$termCorrect[Corrects$n==i]<-table(Pimm.results$termCorrect[Pimm.results$n==i])[2]
			Corrects$termIncorrect[Corrects$n==i]<-table(Pimm.results$termCorrect[Pimm.results$n==i])[1]
		}
					}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
		}


		Corrects[is.na(Corrects)]<-0
		Corrects$noModel<-100-(as.numeric(Corrects$Correct)+as.numeric(Corrects$Incorrect))
		Corrects$TermnoModel<-100-(as.numeric(Corrects$termCorrect)+as.numeric(Corrects$termIncorrect))
		Corrects[,1]<-as.numeric(Corrects[,1])
		Corrects[,2]<-as.numeric(Corrects[,2])
		Corrects[,3]<-as.numeric(Corrects[,3])
		Corrects[,4]<-as.numeric(Corrects[,4])
		Corrects[,7]<-as.numeric(Corrects[,7])




		Corrects.melt<-melt(Corrects,id=c("n","model"))
		Corrects.melt$value<-as.numeric(as.character(Corrects.melt$value))/100
		Corrects.melt$param<-"Initiation"
		Corrects.melt$param[25:48]<-"Termination"
		Corrects.melt$param[61:72]<-"Termination"
		Corrects.melt$param<-factor(Corrects.melt$param)
		Corrects.melt$variable[25:36]<-"Correct"
		Corrects.melt$variable[37:48]<-"Incorrect"
		Corrects.melt$variable[61:72]<-"noModel"

return(Corrects.melt)
}

