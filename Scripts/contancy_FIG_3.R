{# tools
{# SETTINGS & DATA
	{# do you want plots within R (PNG=FALSE) or as PNG (PNG = TRUE)?
		PNG=TRUE
		#PNG = TRUE
	}
	{# define working directories
	     wd="C:/Users/mbulla/Documents/Dropbox/Science/Projects/MS/Removal/Analyses/"	
		 outdir="C:/Users/mbulla/Documents/Dropbox/Science/Projects/MS/Removal/Outputs/"
		#out3="C:/Users/mbulla/Documents/Dropbox/Science/Projects/MS/Removal/Tabels/"
		#out2="C:/Users/mbulla/Documents/Dropbox/Science/Projects/MS/Removal/Figs/"	
	}
	{# load packages, constants and data
		source(paste(wd, 'Scripts/Constants&Functions.R',sep=""))
		#source(paste(wd, 'Scripts/Prepare_Data.R',sep=""))
	}
}
	
	{# run first
		#u=read.csv(paste(wd,'Data/freeze/experiment_metadata (3).csv', sep=""),stringsAsFactors=FALSE)
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
			bt$sex=as.factor(bt$sex)
			bb$inc_eff_treated=bt$inc_eff[match(bb$nest, bt$nest)]
			bb$col_=ifelse(bb$sex=="f", "#FCB42C","#535F7C")
			
	}
				
			{# prepare bout lengths -for definition of control treatment
				load(file=paste(wd,'bout.Rdata', sep=""))
				d=b[-which(b$nest%in%c('s410','s514','s524','s624')),]
				d$bout_start=as.POSIXct(d$bout_start)
				d=d[-which(d$nest=='s402' & d$bout_start>as.POSIXct("2013-06-19 12:00:00")& d$bout_start<as.POSIXct("2013-06-21 02:00:00")),] 
					#d[which(d$nest=='s404' & d$exper=='t'),3:9]
					#d[which(d$nest=='s404'),c(4:9,12)]
					#min(d$bout_ID[which(!is.na(d$exper) & d$exper=='t' & d$nest=='s404')])
				dd=d[which(d$bout_ID<min(d$bout_ID[which(!is.na(d$exper) & d$exper=='t')])&  d$bout_type=='inc' & d$bout_ID>1),]
						d[which(d$nest=='s404' & d$bout_type=='inc' & d$bout_ID>1 & d$bout_ID<min(d$bout_ID[which(!is.na(d$exper) & d$exper=='t')])),]
						#d[which(d$nest=='s808' & d$bout_start<d$bout_start[which(!is.na(d$exper) & d$exper=='t')][1]),]
						#d$bout_start[which(d$nest=='s808' & !is.na(d$exper) & d$exper=='t')][1]
				dsplit=split(d,d$nest)
					foo=lapply(dsplit,function(x) {
												#x=dsplit$"s404"
													x=x[which(x$bout_ID<min(x$bout_ID[which(!is.na(x$exper) & x$exper=='t')])&  x$bout_type=='inc' & x$bout_ID>1),]# limit to before treatment data
													tst=x$sex[nrow(x)]
													x=x[which(x$sex!=tst),]
													
													if(nrow(x)>3){x=x[(nrow(x)-2):nrow(x),]} # keep only three bouts
																				
												return(x)
													}
													)
														
					dd=do.call(rbind, foo)
				ggplot(dd,aes(x=nest,y=bout_length/60,fill=nest))+geom_boxplot()+geom_point()
				#dd[which(dd$nest=='s402'),]
				#dd[which(dd$nest=='s808'),]
				dd=ddply(dd,.(nest,sex),summarise, bout=median(bout_length))
												# test if sex correctly assigned -YES
													g=read.csv(paste(wd,'experiment_metadata.csv', sep=""), stringsAsFactors=FALSE)
													g=g[g$type=='exp',]
													dd$taken=g$taken[match(dd$nest,g$nest)]
													dd[dd$sex==dd$taken,]
			}									
			{# define nests - order according to constancy in treated
				load(file=paste(wd,'experimental.Rdata', sep=""))
				b=b[-which(b$nest%in%c('s410','s514','s524','s624')),]
				bt=b[b$exper=='t',]
				bt$inc_eff[bt$nest=='s806']=bt$inc_eff_2[bt$nest=='s806']
				
				bb=b[b$exper=='c',]
				bb$inc_eff[bt$nest=='s806']=bb$inc_eff_2[bb$nest=='s806']
				
					#plot(bt$inc_eff, bt$inc_eff_3)
					#abline(a=0,b=1)
				#bt[order(bt$inc_eff,decreasing = TRUE),c('nest','inc_eff']
				nests_=bt$nest[order(bt$inc_eff,decreasing = FALSE)]
				# create fake dummy for plotting
				nn=data.frame(nest=nests_, n=c(1:length(nests_)), stringsAsFactors=FALSE)
				nn$sex=bt$sex[match(nn$nest,bt$nest)]
				nn$col_=ifelse(nn$sex=="f", "#FCB42C","#535F7C")
				nn$sex_symbol=ifelse(nn$sex=="f", "\u2640","\u2642")
				nn$bt=round(bt$inc_eff[match(nn$nest,bt$nest)],2)
				nn$bb=round(bb$inc_eff[match(nn$nest,bb$nest)],2)
				nn$bt=ifelse(nchar(nn$bt)==3, paste(nn$bt,"0",sep=""), nn$bt)
				nn$bb=ifelse(nchar(nn$bb)==3, paste(nn$bb,"0",sep=""), nn$bb)
				
			}	
			{# prepare main metadata with the above and return of the captive bird and for labeling
				load(file=paste(wd,'bout.Rdata', sep=""))
				d=b[-which(b$nest%in%c('s410','s514','s524','s624')),]
				d=d[-which(d$bout_type=='gap'),]
				d$bout_start=as.POSIXct(d$bout_start)
				m=read.csv(paste(wd,'experiment_metadata.csv', sep=""), stringsAsFactors=FALSE)
				m=m[-which(m$after==""),]
				d$after=m$after[match(d$nest, m$nest)]
				
				d[d$nest=='s406',c('exper','sex','bout_ID','bout_ID_tr', 'bout_start')]
				
				# create start of experimental and return of the captive bird for those where bird returned
				dsplit=split(d[d$after%in%c("gap_return", "return"),],d$nest[d$after%in%c("gap_return", "return")])
				foo=lapply(dsplit,function(x) {
												#x=dsplit$"s805"
													#print(x$nest[1])
													x=x[!is.na(x$sex),]
													tst=unique(x$sex[!is.na(x$exper) & x$exper=='t'])
													#print(length(tst))
													tst2=max(x$bout_start[!is.na(x$exper) & x$exper=='t'])
													if(tst=='f'){x=data.frame(nest=x$nest[1],start_=min(x$bout_start[!is.na(x$exper) & x$exper=='t']), end_=min(x$bout_start[x$bout_start>tst2 & x$sex=='m']))}
													if(tst=='m'){x=data.frame(nest=x$nest[1], start_=min(x$bout_start[!is.na(x$exper) & x$exper=='t']),end_=min(x$bout_start[x$bout_start>tst2 & x$sex=='f']))}
													return(x)
													}
													)
				e=do.call(rbind,foo)
				e$release=as.POSIXct(m$release[match(e$nest,m$nest)])
				e$bout_a=as.numeric(difftime(e$end_,e$release, units='hours'))
				e$after='return'
				
				
				# median bout to estimate period (when bird could have returned) for birds with partner never returning
					#summary(e)
					bout_med=max(e$bout_a)
				
					d2=d[which(!d$after%in%c("gap_return", "return")),]
					e2=ddply(d2,.(nest), summarise, after=unique(after), start_=min(bout_start[!is.na(exper) & exper=='t']))
					e2$release=as.POSIXct(m$release[match(e2$nest,m$nest)])
					e2$end_=e2$release+bout_med*60*60
					e2$bout_a=bout_med
					
				# bring two datasets together
					e=rbind(e,e2)
				
				# add bout length for definition of control
					e$c_bout=dd$bout[match(e$nest,dd$nest)]
					e$end_c=e$start_+e$c_bout*60
					
					e$r_time=as.numeric(difftime(e$release,e$end_c,units='hours'))
					e$r_time[e$an!=""]=NA
					e$e_time=as.numeric(difftime(e$end_,e$end_c,units='hours'))
				
					e$n=nn$n[match(e$nest,nn$nest)]
					e$an=ifelse(e$after=='desertion', 'deserted',ifelse(e$after=='continues', 'continues',""))
							f=data.frame(nest=c('s409','s510','s516','s711','s806'),days=c(5,10,9,4,3))
							f$days=paste(f$days, 'days')
					e$an[e$an=='continues']=f$days[match(e$nest[e$an=='continues'], f$nest)]
					
					e$sex=nn$sex[match(e$nest,nn$nest)]
					ann_text <- data.frame(datetime_c=31, roll_con=0.5, lab=e$an,
						   n = e$n, sex=e$sex, nest=e$nest, col_line=NA)

			}
	}			
{# prepare nest datasets 	- add this part once sqlite is ready
				db_=FALSE # in preparation and usable once database with raw data is ready
				if(db_ == TRUE){
				conMy=dbConnect(MySQL(),user='root',host='127.0.0.1', password='',dbname='')
				l=list()
				for (j in 1:length(nests_)){
						print(nests_[j])
						pp=paste("extract_2013.mb_inc_2013_",nests_[j],sep="")	
						start_=e$start_[e$nest==nests_[j]]
						end_=e$end_[e$nest==nests_[j]]
						release_=e$release[e$nest==nests_[j]]
						b=wet(conMy, paste("SELECT pk, datetime_, inc_t, incubation FROM  ", pp," where datetime_>='",start_,"' and datetime_<='",end_,"'"))
						b$datetime_= as.POSIXct(b$datetime_)
						b$datetime_c=as.numeric(difftime(b$datetime_,e$end_c[e$nest==nests_[j]], units='hours'))
						
						if(nests_[j]=='s806'){b$inc_t=b$incubation
											  b$inc_t[is.na(b$inc_t)]=0}
						
						b$roll_con=rollmean(b$inc_t, 12*60, align="center", fill="extend") # 12=1min
						b$nest=nests_[j]
						b$col_line=ifelse(b$datetime_<=release_,'black','lightgrey')
						l[[j]]=b
						}
				u=do.call(rbind,l)
				# add order for plotting
				u$n=nn$n[match(u$nest,nn$nest)]
				u$sex=nn$sex[match(u$nest,nn$nest)]
				u$col_=nn$col_[match(u$nest,nn$nest)]
				
				save(u, file = paste(wd, 'constancy_for_Fig_3.Rdata'))
				}else{load(paste(wd, 'constancy_for_Fig_3.Rdata'))}				
				}
				