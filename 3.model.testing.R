#Model testing

library(ggplot2)
library(dplyr)
library(moult)
library(chngpt)
library(reshape)
library(gridExtra)
library(grid)
library(cowplot)
library(ggformula)

source("~/Dropbox/HingeModels/molt.simulator.R")
source("~/Dropbox/HingeModels/modelFitFunctions.R")

#######################
#Sample Size variation#
#######################


#simulate data
samples<-c(5,10,15,20,30,40,50,75,100,150,200,250)

sims<-list()

for(j in 1:100){
	sims[[as.character(j)]]<-molt.oneSim(samples=samples,init=100,duration=100,var=10)
				}

for(i in 1:100){
	sims[[i]]<-bind_rows(sims[[i]])
}				

#fit to UZ
UZ.sample.size<-fit.UZ(sims=sims,samples=samples)
write.csv(UZ.sample.size,file="~/Dropbox/HingeModels/UZ_sample_size.csv")
UZ.eval.sample.size<-evaluate.UZ(UZ.sample.size,init=100,duration=100,samples=samples)
write.csv(UZ.eval.sample.size,file="~/Dropbox/HingeModels/UZ_sample_size_eval.csv")

#fit to DH

DH.sample.size<-fit.DH(sims=sims,samples=samples)
write.csv(DH.sample.size,file="~/Dropbox/HingeModels/DH_sample_size.csv")

DH.eval.sample.size<-evaluate.DH(DH.sample.size,init=100,duration=100,samples=samples)
write.csv(DH.eval.sample.size,file="~/Dropbox/HingeModels/DH_sample_size_eval.csv")


#fit to Pimm

Pimm.sample.size<-fit.Pimm(sims=sims,samples=samples)
write.csv(Pimm.sample.size,file="~/Dropbox/HingeModels/Pimm_sample_size.csv")
Pimm.eval.sample.size<-evaluate.Pimm(Pimm.sample.size,init=100,duration=100,samples=samples)
write.csv(Pimm.eval.sample.size,file="~/Dropbox/HingeModels/Pimm_sample_size_eval.csv")






#######################
##start variation###
#######################


#simulate data
v.samples<-c(0.5,1,2,3,5,8,13,21,34,55,89,144)

v.sims<-list()

for(j in 1:100){
	v.sims[[as.character(j)]]<-v.molt.oneSim(samples=50,init=100,duration=100,var=v.samples)
				}

for(i in 1:100){
	v.sims[[i]]<-bind_rows(v.sims[[i]])
}	


#fit to UZ
UZ.var<-fit.UZ(sims=v.sims,samples=v.samples)
write.csv(UZ.var,file="~/Dropbox/HingeModels/UZ_var.csv")
UZ.eval.var<-evaluate.UZ(UZ.var,init=100,duration=100,samples=v.samples)
write.csv(UZ.eval.var,file="~/Dropbox/HingeModels/UZ_var_eval.csv")

#fit to DH

DH.var<-fit.DH(sims=v.sims,samples=v.samples)
write.csv(DH.var,file="~/Dropbox/HingeModels/DH_Var.csv")
DH.eval.var<-evaluate.DH(DH.var,init=100,duration=100,samples=v.samples)
write.csv(DH.eval.var,file="~/Dropbox/HingeModels/DH_Var_eval.csv")


#fit to Pimm

Pimm.var<-fit.Pimm(sims=v.sims,samples=v.samples)
write.csv(Pimm.var,file="~/Dropbox/HingeModels/Pimm_var.csv")
Pimm.eval.var<-evaluate.Pimm(Pimm.var,init=100,duration=100,samples=v.samples)
write.csv(Pimm.eval.var,file="~/Dropbox/HingeModels/Pimm_var_eval.csv")




#######################
##Duration variation###
#######################

#simulate data
d.samples<-c(1,2,5,10,20,30,40,50,100,150,200,250)


d.sims<-list()

for(j in 1:100){
	d.sims[[as.character(j)]]<-d.molt.oneSim(samples=50,init=100,duration=d.samples,var=10)
				}

for(i in 1:100){
	d.sims[[i]]<-bind_rows(d.sims[[i]])
}	



#fit to UZ
UZ.dur<-fit.UZ(sims=d.sims,samples=d.samples)
write.csv(UZ.dur,file="~/Dropbox/HingeModels/UZ_dur.csv")
UZ.eval.dur<-evaluate.UZ(UZ.dur,init=100,duration=d.samples,samples=d.samples)
write.csv(UZ.eval.dur,file="~/Dropbox/HingeModels/UZ_dur_eval.csv")

#fit to DH

DH.dur<-fit.DH(sims=d.sims,samples=d.samples)
write.csv(DH.dur,file="~/Dropbox/HingeModels/DH_dur.csv")
DH.eval.dur<-evaluate.DH(DH.dur,init=100,duration=d.samples,samples=d.samples)
write.csv(DH.eval.dur,file="~/Dropbox/HingeModels/DH_dur_eval.csv")


#fit to Pimm

Pimm.dur<-fit.Pimm(sims=d.sims,samples=d.samples)
write.csv(Pimm.dur,file="~/Dropbox/HingeModels/Pimm_dur.csv")
Pimm.eval.dur<-evaluate.Pimm(Pimm.dur,init=100,duration=d.samples,samples=d.samples)
write.csv(Pimm.eval.dur,file="~/Dropbox/HingeModels/Pimm_dur_eval.csv")


#######################
##Hidden Populations###
#######################

h.samples<-c(-50,-30,-20,-10,-5,0,5,10,20,30,50,60)




h.sims<-list()

for(j in 1:100){
h.sims[[as.character(j)]]<-h.molt.oneSim()
cat(j," ")
}


for(i in 1:100){
	h.sims[[i]]<-bind_rows(h.sims[[i]])
}


#fit to UZ
UZ.hidden<-fit.UZ(sims=h.sims,samples=h.samples)
write.csv(UZ.hidden,file="~/Dropbox/HingeModels/UZ_hidden.csv")
UZ.eval.hidden<-evaluate.UZ(UZ.hidden,init=100,duration=100,samples=h.samples)
write.csv(UZ.eval.hidden,file="~/Dropbox/HingeModels/UZ_hidden_eval.csv")

#fit to DH

DH.hidden<-fit.DH(sims=h.sims,samples=h.samples)
write.csv(DH.hidden,file="~/Dropbox/HingeModels/DH_hidden.csv")
DH.eval.hidden<-evaluate.DH(DH.hidden,init=100,duration=100,samples=h.samples)
write.csv(DH.eval.hidden,file="~/Dropbox/HingeModels/DH_hidden_eval.csv")


#fit to Pimm

Pimm.hidden<-fit.Pimm(sims=h.sims,samples=h.samples)
write.csv(Pimm.hidden,file="~/Dropbox/HingeModels/Pimm_hidden.csv")
Pimm.eval.hidden<-evaluate.Pimm(Pimm.hidden,init=100,duration=100,samples=h.samples)
write.csv(Pimm.eval.hidden,file="~/Dropbox/HingeModels/Pimm_hidden_eval.csv")



#######################
##Rate Variation###
#######################

rv.samples<-c(0,.05,.1,.2,.5,.75,1,1.5,2,3,4,5)




rv.sims<-list()

for(j in 1:100){
rv.sims[[as.character(j)]]<-rv.molt.oneSim()
cat(j," ")
}


for(i in 1:100){
	rv.sims[[i]]<-bind_rows(rv.sims[[i]])
}



#fit to UZ
UZ.rv<-fit.UZ(sims=rv.sims,samples=rv.samples)
write.csv(UZ.rv,file="~/Dropbox/HingeModels/UZ_rv.csv")
UZ.eval.rv<-evaluate.UZ(UZ.rv,init=100,duration=100,samples=rv.samples)
write.csv(UZ.eval.rv,file="~/Dropbox/HingeModels/UZ_rv_eval.csv")

#fit to DH


DH.rv<-fit.DH(sims=rv.sims,samples=rv.samples)
write.csv(DH.rv,file="~/Dropbox/HingeModels/DH_rv.csv")
DH.eval.rv<-evaluate.DH(DH.rv,init=100,duration=100,samples=rv.samples)
write.csv(DH.eval.rv,file="~/Dropbox/HingeModels/DH_rv_eval.csv")


#fit to Pimm

Pimm.rv<-fit.Pimm(sims=rv.sims,samples=rv.samples)
write.csv(Pimm.rv,file="~/Dropbox/HingeModels/Pimm_rv.csv")
Pimm.eval.rv<-evaluate.Pimm(Pimm.rv,init=100,duration=100,samples=rv.samples)
write.csv(Pimm.eval.rv,file="~/Dropbox/HingeModels/Pimm_rv_eval.csv")



