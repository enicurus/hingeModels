###Simulate molt data

#Author: Ryan S. Terrill
#ornithoterrill at gmail dot com
#Last update: 17 September 2019


#n=numer of samples
#init = julian day of molt initation
#duration = molt duration
#var = SD of molt initiation
#start= julian day to start sampling
#end = julian day to end sampling

molt.oneSim<-function(samples,init,duration,var){
	simData<-list()
for (i in samples){
	simData[[as.character(i)]]<-molt.simulate(n=i,init=init,duration=duration,var=var)
	}
for (i in samples){
	simData[[as.character(i)]][,3]<-as.character(i)
	}
	return(simData)
}

v.molt.oneSim<-function(samples,init,duration,var){
	v.simData<-list()
for (i in v.samples){
	v.simData[[as.character(i)]]<-molt.simulate(n=samples,init=init,duration=duration,var=i)
	}
for (i in v.samples){
	v.simData[[as.character(i)]][,3]<-as.character(i)
	}
	return(v.simData)
}


d.molt.oneSim<-function(samples,init,duration,var){
	d.simData<-list()
for (i in duration){
	d.simData[[as.character(i)]]<-molt.simulate(n=samples,init=init,duration=i,var=var)
	}
for (i in duration){
	d.simData[[as.character(i)]][,3]<-as.character(i)
	}
	return(d.simData)
}


h.molt.oneSim<-function(){
	h.simData<-list()
for (i in h.samples){
	h.simData[[as.character(i)]]<-molt.simulate.hidden(n=100,init=100,duration=100,var=15,offset=i)
	}
for (i in h.samples){
	h.simData[[as.character(i)]][,3]<-as.character(i)
	}
	return(h.simData)
}

rv.molt.oneSim<-function(){
	rv.simData<-list()
for (i in rv.samples){
	rv.simData[[as.character(i)]]<-molt.simulate.slopevar(n=100,init=100,duration=100,var=15,durVar=i)
	}
for (i in rv.samples){
	rv.simData[[as.character(i)]][,3]<-as.character(i)
	}
	return(rv.simData)
}



molt.simulate<-function(n,init,duration,var,start=1,end=365){
	num<-sample(start:end,n,replace=TRUE)
	out<-matrix(ncol=2,nrow=n)
	if(length(duration)==1){
	term=init+duration
	for(i in 1:n){
		initiation<-rnorm(n,mean=init,sd=var)
			if(num[i]<initiation[i]){
				out[i,]<-c(num[i],0)
			}
		else {
			out[i,]<-c(num[i],(((100/(term-init))*num[i]))-initiation[i]*(100/(term-init)))
			}
						termination<-list()
		termination[i]<-term-(init-initiation[i])	
		if(num[i]>termination[i]){
			out[i,]<-c(num[i],100)
	}
	}
}

else{
		term<-list()
		for (i in 1:length(duration)){
			term[i]<-init+duration[i]
		}

			for(i in 1:n){
		initiation<-rnorm(n,mean=init,sd=var)
			if(num[i]<initiation[i]){
				out[i,]<-c(num[i],0)
			}
		else {
			out[i,]<-c(num[i],(((100/(term[i]-init))*num[i]))-initiation[i]*(100/(term[i]-init)))
			}
						termination<-list()
		termination[i]<-term[i]-(init-initiation[i])	
		if(num[i]>termination[i]){
			out[i,]<-c(num[i],100)
	}
	}
}
	out<-data.frame(out);colnames(out)<-c("julian","pscore")
return(out)
}



# include slope variance

molt.simulate.slopevar<-function(n,init,duration,var,durVar,start=1,end=365){
	num<-sample(start:end,n,replace=TRUE)
	out<-matrix(ncol=2,nrow=n)
	for(i in 1:n){
			term=rnorm(n,mean=duration,sd=durVar)+init
		initiation<-rnorm(n,mean=init,sd=var)
			if(num[i]<initiation[i]){
				out[i,]<-c(num[i],0)
			}
		else {
			out[i,]<-c(num[i],(((100/(term[i]-init))*num[i]))-initiation[i]*(100/(term[i]-init)))
			}
						termination<-list()
		termination[i]<-term[i]-(init-initiation[i])	
		if(num[i]>termination[i]){
			out[i,]<-c(num[i],100)
	}
	}
	out<-data.frame(out);colnames(out)<-c("julian","pscore")
return(out)
}

test=molt.simulate.slopevar(n=200,init=50,duration=150,var=10,durVar=100)

#plot(test)


#### simulate sigmoidal data

# here's how to make some sigmoidal data -- adapt this to the sims
#sigmoid<-function(x){
#	1/(1+exp(-x))
#}

#x <- seq(-5, 5, 0.1)

#plot(x+10, sigmoid(x), col='blue')

## Skewed start and end - add a "hidden" population that starts early/late, end earllate
## the hidden population can be set with hiddenPopSize


molt.simulate.hidden<-function(n,init,duration,var,durVar=0,offset,hiddenPopSize=n/5,start=1,end=365){
	num<-sample(start:end,n,replace=TRUE)
	out<-matrix(ncol=2,nrow=n)
	for(i in 1:n){
			term=rnorm(n,mean=duration,sd=durVar)+init
		initiation<-rnorm(n,mean=init,sd=var)
			if(num[i]<initiation[i]){
				out[i,]<-c(num[i],0)
			}
		else {
			out[i,]<-c(num[i],(((100/(term[i]-init))*num[i]))-initiation[i]*(100/(term[i]-init)))
			}
						termination<-list()
		termination[i]<-term[i]-(init-initiation[i])	
		if(num[i]>termination[i]){
			out[i,]<-c(num[i],100)
	}
	}
	out<-data.frame(out);colnames(out)<-c("julian","pscore")
	num2<-sample(start:end,hiddenPopSize,replace=TRUE)
	out2<-matrix(ncol=2,nrow=hiddenPopSize)
	for(i in 1:hiddenPopSize){
			term=rnorm(hiddenPopSize,mean=duration,sd=durVar)+init
		initiation2<-rnorm(hiddenPopSize,mean=init+offset,sd=var)
			if(num2[i]<initiation2[i]){
				out2[i,]<-c(num2[i],0)
			}
		else {
			out2[i,]<-c(num2[i],(((100/(term[i]-init))*num2[i]))-initiation2[i]*(100/(term[i]-init)))
			}
						termination2<-list()
		termination2[i]<-term[i]-(init-initiation2[i])	
		if(num2[i]>termination2[i]){
			out2[i,]<-c(num2[i],100)
	}
	}
	
	out<-data.frame(out);colnames(out)<-c("julian","pscore")
	out2<-data.frame(out2);colnames(out2)<-c("julian","pscore")
		out<-rbind(out,out2)
return(out)
}

