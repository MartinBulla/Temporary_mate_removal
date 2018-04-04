# script used to impute an escape distance for a bird from nest s628
	
{# run first	
	wd = "C:/Users/mbulla/Documents/Dropbox/Science/Projects/MS/Removal/Analyses/"
	require(Amelia)
	require(plyr)
	require(lattice)
	
		Mode <- function(x) {
			ux <- unique(x)
			ux[which.max(tabulate(match(x, ux)))]
			}
			
	u=read.csv(paste(wd,'Data/experiment_metadata.csv', sep=""),stringsAsFactors=FALSE)
		load(file=paste(wd,'Data/experimental.Rdata', sep=""))
		
		b=b[-which(b$nest%in%c('s410','s514','s524','s624')),] # excludes nest where wrong bird captures, nests with desertion and predation prior or during experiment
		b$inc_eff[b$nest=='s806']=b$inc_eff_2[b$nest=='s806'] # incubation efficency calculated only from rfid readings, as no temperature data obtained, in all others it is only temperature based
		b$bird_ID=paste(b$nest,b$sex,sep="")
		b$lat=u$lat[match(b$nest,u$nest)]
		b$lon=u$lon[match(b$nest,u$nest)]		
		# pre-define variables
			b$bout_start=as.POSIXct(b$bout_start)
			b$exper=as.factor(b$exper)
			b$sex=as.factor(b$sex)
			b$exper_sex=as.factor(ifelse(b$exper=='c', 'c', ifelse(b$sex=='f', 't_female', 't_male')))
			b$exper_sex_2=relevel(b$exper_sex,ref="t_female")
			b$exper_sex_int=as.factor(ifelse(b$exper=='c' & b$sex=='f', 'c_female',ifelse(b$exper=='c' & b$sex=='m', 'c_male', ifelse(b$sex=='f', 't_female', 't_male'))))
			
			# time as midpoint of bout
				b$midbout=b$bout_start+60*60*b$bout_length/2
				b$time = as.numeric(difftime(b$midbout,trunc(b$midbout,"day"),   #DO NOT USE rfid_T$day here!! (the time will be wrong)
								units = "hours"))
				b$rad=b$time*pi/12
				b$sin_=sin(b$rad)
				b$cos_=cos(b$rad)
		# create separate dataset for control and treated experimental period
			#bb=ddply(b,.(nest), summarise, n=length(nest))
			bb=b[b$exper=='c',]
			bt=b[b$exper=='t',]
		p_ff = read.csv(file=paste(wd,'Data/prop_inc.csv', sep=""), sep=",",stringsAsFactors = FALSE)	
			bt$p_f=p_ff$med_[match(bt$nest,p_ff$nest)]
			bt$prop=ifelse(bt$sex=='f', bt$p_f, 1-bt$p_f)
		
		v=read.csv(paste(wd,'Data/escape.csv', sep=""), stringsAsFactors=FALSE)
				v = v[-which(v$nest%in%c('s410','s514','s524','s624')),] # limit to experimental individulas
				es=ddply(v,.(nest,sex),summarise, esc=median(distm))
				bt$esc=es$esc[match(paste(bt$nest, bt$sex),paste(es$nest,es$sex))]	
}		
{# impute		
	data1 = bt[,c('bird_ID_filled', 'sex', 'esc', 'inc_eff', 'disturb', 't_ambient_med','prop')]
						idvars1 = c('bird_ID_filled')
						noms1=c('sex')
						priors=matrix(c(1:25,rep(3,25),rep(0,25),rep(80,25), rep(.99,25)),25,5) # first is number of rows; second is column number of imputed variable (2); third is minimum (0), fourth is maximum (80m); fifthe is confidence (.99)
						l=list()
						for(i in 1:1000){
								am = amelia(x=data1, idvars = idvars1,noms=noms1,priors=priors, m = 1)
								#am = amelia(x=data1, idvars = idvars1, noms=noms1,m = 1)
								l[[i]]=am$imputations[[1]]
								}
						u=do.call(rbind,l)
						boxplot(u$esc[u$bird_ID==229194794])
						densityplot(~u$esc[u$bird_ID==229194794])
						median(u$esc[u$bird_ID==229194794]) # 34.05681
						mean(u$esc[u$bird_ID==229194794])
						Mode(u$esc[u$bird_ID==229194794])
}