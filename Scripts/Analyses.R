# R Script generates statistical outputs for 

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

 # INTRODUCTION
	{# Figure 1 - predictions
		if(PNG == TRUE) {
					png(paste(out2,"Figure_1", ".png", sep = ""), width=0.875*2, height=0.875*2, units = "in", res = 600)	 
					}else{
					dev.new(width=0.875*2, height=0.875*2) #dev.new(width=3.5, height=1.97)
					}	
			
		#par(mfrow=c(2,2),mar=c(2,2,0.7,0.7), ps=12, mgp=c(1.2,0.35,0), las=1, cex=1,tcl=-0.15,bty="o", oma = c(2, 2, 0.5, 0.5)) # cex is multiplied as in multiplot figure the cex is reduced to 0.66	
		par(mfrow=c(2,2),mar=c(0,0,0,0), ps=12, mgp=c(1.2,0.35,0), las=1, cex=1,tcl=-0.15,bty="o", oma = c(2, 2, 0.5, 0.5),col.axis="grey30",col.lab="grey30", col.main="grey30", fg="grey70") 
		
		{#(a)	
							#par(mar=c(0,2,2.7,0), ps=12, mgp=c(0.25,0,0), las=1, cex.lab=0.7, cex.axis=0.6, tcl=-0.15,bty="o")
							#par(mar=c(0,0,1,0))
							plot(c(0,1), c(0,1),xlab=NA, ylab=NA, xlim=c(-0.125,1.125),ylim=c(-0.05,1.1),type="n",xaxt='n',yaxt='n')#,yaxt='n')#, )#, 	
							axis(2, at = seq(0,1, 0.5),cex.axis=0.5)	
							mtext(expression(bold("a")),side=3,line=-0.65, cex=0.7, las=1,adj=0.04, col='black')
							polygon(x=c(c(0.4, 0.8),rev(c(0.4, 0.8))), y=c(c(-0.096,-0.096),c(1.155,1.155)), col='grey95', border=FALSE)
							lines(c(0,0.4,0.4,0.8,0.8,1),c(0.9,0.9,0,0,0.9,0.9), lwd=1.5, col='black')
							#abline(v=0.4, col='red', lty=3, lwd=0.5)
							box(col = "grey70")
		}
		{#(b) 
							#par(mar=c(0,0,2.7,2), ps=12, mgp=c(0.25,0,0), las=1, cex.lab=0.7, cex.axis=0.6, tcl=-0.15,bty="o")
							plot(c(0,1), c(0,1),xlab=NA, ylab=NA,xlim=c(-0.125,1.125),ylim=c(-0.05,1.1), type="n",xaxt='n',yaxt='n')#,yaxt='n')#, )#, 	
							mtext(expression(bold("b")),side=3,line=-0.65, cex=0.6, las=1,adj=0.04, col='black')
							polygon(x=c(c(0.4, 0.8),rev(c(0.4, 0.8))), y=c(c(-0.096,-0.096),c(1.155,1.155)), col='grey95', border=FALSE)
							lines(c(0,1),c(0.9,0.9), lwd=1.5, col='black')
							#abline(v=0.4, col='red', lty=3, lwd=0.5)
							box(col = "grey70")
		}
		{#(c) 
							#par(mar=c(2.7,2,0,0), ps=12, mgp=c(0.25,0,0), las=1, cex.lab=0.7, cex.axis=0.6, tcl=-0.15,bty="o")
							#par(mgp=c(0.2,0.15,0))
							plot(c(0,1), c(0,1),xlab=NA, ylab=NA, xlim=c(-0.125,1.125),ylim=c(-0.05,1.1),type="n",xaxt='n',yaxt='n')#,yaxt='n')#, )#, 	
							axis(2,  at = seq(0,1, 0.5),cex.axis=0.5)
							#axis(1, col = "grey40", col.axis = "grey20", at = seq(0,1, 0.5),labels=FALSE)# seq(0,24,12),cex.axis=0.6)
							axis(1, at = seq(0,1, 0.5),labels=seq(0,24,12),cex.axis=0.5,mgp=c(0,-0.15,0))
							#text(seq(0,1, 0.5),  par("usr")[3], labels = seq(0,24,12),  xpd = TRUE, cex=0.6) 
							#mtext(text = seq(0,24,12),seq(0,1, 0.5),  par("usr")[3],   outer = TRUE, cex=0.6) 
							polygon(x=c(c(0.4, 0.8),rev(c(0.4, 0.8))), y=c(c(-0.096,-0.096),c(1.155,1.155)), col='grey95', border=FALSE)
							mtext(expression(bold("c")),side=3,line=-0.65, cex=0.7, las=1,adj=0.04, col='black')
							lines(c(0,0.6,0.6,1),c(0.9,0.9,0,0), lwd=1.5, col='black')
							#abline(v=0.4, col='red', lty=3, lwd=0.5)
							box(col = "grey70")
		}
		{#(d) 
						#par(mar=c(2.7,0,0,2), ps=12, mgp=c(0.5,0.35,0), las=1, cex.lab=0.7, cex.axis=0.6, tcl=-0.15,bty="o")
						plot(c(0,1), c(0,1),xlab=NA, ylab=NA,xlim=c(-0.125,1.125),ylim=c(-0.05,1.1), type="n",xaxt='n',yaxt='n')#,yaxt='n')#, )#, 
						mtext(expression(bold("d")),side=3,line=-0.65, cex=0.6, las=1,adj=0.04, col='black')
						axis(1, at = seq(0,1, 0.5),labels=seq(0,24,12),cex.axis=0.5,mgp=c(0,-0.15,0))
						polygon(x=c(c(0.4, 0.8),rev(c(0.4, 0.8))), y=c(c(-0.096,-0.096),c(1.155,1.155)), col='grey95', border=FALSE)
						lines(c(0,0.6,0.6,1),c(0.9,0.9,0.55,0.55), lwd=1.5, col='black')
						#abline(v=0.4, col='red', lty=3, lwd=0.5)
						box(col = "grey70")				
				}
				
		mtext('Time [h]',side=1,line=0.6, cex=0.6, las=1, outer=TRUE, col="grey30")
		mtext('Nest attendance',side=2,line=1.3, cex=0.6, las=0,outer=TRUE,col="grey30")
			
		if(PNG == TRUE) {dev.off()}
	}
				

# METHODS
	{# Recordind incubation and escape distance
	   {# escape distance - (note: sometimes we had two estimates (our and gps) for a single event)
		load("C:\\Users\\mbulla\\Documents\\Dropbox\\Science\\Projects\\MS\\Removal\\Analyses\\Data\\escape.Rdata")
		write.csv(vv,paste(wd,'Data/escape.csv', sep=""),row.names=FALSE)
		str(vv)
		table(vv$nest_sex)
		ggplot(vv,aes(x = nest_sex, y = distm)) + geom_boxplot() + geom_dotplot(aes(fill = nest_sex),binaxis = 'y', stackdir = 'center',  position = position_dodge())
		densityplot(~vv$distm)
		densityplot(~log(vv$distm+0.5))
		
		vg = vv[vv$type == 'gps',]	
		ve = vv[!paste(vv$nest_sex, vv$datetime_)%in%paste(vg$nest_sex, vg$datetime_),]
		v = rbind(vg,ve)
		x = data.frame(table(v$nest_sex))
		summary(factor(x$Freq))
		summary(x$Freq)
		nrow(x)+1 #1 bird had no escape distance estimate)
		
		m = lmer(distm ~ 1 + (1|nest_sex), v)
		summary(m)
		
		xx = x$Var1[x$Freq>1]
		m = lmer(distm ~ 1 + (1|nest_sex), v[v$nest_sex%in%xx,])
			summary(m)
		table(v$nest_sex[v$nest_sex%in%xx])
		
		vv = v
		vv$time = as.numeric(difftime(vv$datetime_,trunc(vv$datetime_,"day"), units = "hours"))
		vv$day_j=as.numeric(format(as.Date(trunc(vv$datetime_, "day")),"%j")) - as.numeric(format(as.Date(trunc(min(vv$datetime_), "day")),"%j"))+1
		ggplot(vv,aes(x = day_j, y = distm)) + geom_point() +stat_smooth(method = 'lm')
			 load("C:\\Users\\mbulla\\Documents\\Dropbox\\Science\\Projects\\MS\\Removal\\Analyses\\Data\\inc_start_estimation_2013_new.Rdata")
			 e$nest = tolower(e$nest)
			vv$inc_start = e$inc_start[match(vv$nest,e$nest)] 	
			vv$day_inc=as.numeric(difftime(vv$datetime_,vv$inc_start, units = "days"))	
		ggplot(vv,aes(x = day_inc, y = distm)) + geom_point() +stat_smooth(method = 'lm')
		ggplot(vv,aes(x = day_inc, y = distm)) + geom_point() +stat_smooth()
		
		m = lmer(distm ~ 1 + (1|nest_sex), vv)
		m = lmer(distm ~ day_j + (day_j|nest_sex), vv)
		m = lmer(distm ~ day_inc + (day_inc|nest_sex), vv)
		m = lmer(distm ~ 1 + (day_inc|nest_sex), vv)
		75.56/(75.56+0.08+269)
		pp = profile(m)
		confint(pp)
		m = lmer(distm ~ day_inc + (day_inc|nest_sex), vv[vv$nest_sex%in%xx,])
		plot(allEffects(m))
		summary(glht(m))
		summary(m)
		
		v = vv[paste(vv$datetime_,vv$nest_sex)%in%paste(vv$datetime_,vv$nest_sex)[duplicated(paste(vv$datetime_,vv$nest_sex))],]
		v = v[order(v$nest_sex, v$datetime_),]
		table(paste(v$nest_sex, v$datetime))
		pd <- position_dodge(0.4)
		ggplot(v,aes(x=type, y = distm, col = factor(paste(nest_sex, datetime_))))+geom_point(position = pd) +geom_line(position = pd, aes(group = factor(paste(nest_sex, datetime_)), col = factor(paste(nest_sex, datetime_)))) +
				theme(legend.position = "none") + ylab('escape distance [m]') + xlab('type of estimation')
		

			
		
}
{# escape distance 2011-2012
 {# tools	
	wet_=dbGetQuery
	
	
	
	db_barr12="M:\\Science\\Projects\\PHD\\FIELD_METHODS\\Barrow\\2012\\DATA\\Database\\Barrow_2012.sqlite"
	con_barr12 = dbConnect(dbDriver("SQLite"),dbname = db_barr12)
	
	db_barr11="M:\\Science\\Projects\\PHD\\FIELD_METHODS\\Barrow\\2011\\DATA\\Database\\Barrow_2011.sqlite"
	con_barr11 = dbConnect(dbDriver("SQLite"),dbname = db_barr11)
}	
 {# Barrow 2011
	v=dbGetQuery(con_barr11,"SELECT date_time as datetime_, nest, dist_left_nest dist_est FROM NESTS where dist_left_nest is not null and species='sesa'")
		v[nchar(v$datetime_)!=nchar('2012-06-12 20:36:00'),]
		v = v[nchar(v$datetime_)>nchar('2012-06-12 20'),]
		v$datetime_ = as.POSIXct(v$datetime_,format="%Y-%m-%d %H:%M:%S", tz = 'America/Anchorage')		
	
	load("M:\\Science\\Projects\\PHD\\STATS\\DATASETS\\weather_SESA_2011_MS_BeEc.Rdata")
		d$nest = tolower(d$nest)
	v = v[v$nest %in% unique(d$nest),]	
			l = list()
			for(i in 1:nrow(v)){
							fi=v[i,]
							di=d[which(d$nest==fi$nest & fi$datetime_>=d$bout_start & fi$datetime_<=d$bout_end+60),] # adding 1 min to the bout end to make the selection more robust
							if(nrow(di)==0){#fi$sex=NA
											 print(fi) #all these visits where before we had rfid on and hence sex cannot be determined
											}else{	fi$sex=di$sex
													fi$bird_ID_filled = di$bird_ID_filled
													l[[i]]=fi
													print(i)
												 }
							
							
							}
	v=do.call(rbind,l)		
	v$year = 2011	
	v11 =v
	}
 {# Barrow 2012
	require(sp)
	v=dbGetQuery(con_barr12,"SELECT visits_pk pk,datetime_, nest, gps_left, wp_left, dist_left_nest dist_est FROM VISITS where gps_left is not null")
	g=dbGetQuery(con_barr12,"SELECT*FROM GPS")
	n=dbGetQuery(conMy,"SELECT nest, species, latit, longit FROM avesatbarrow.nests where year_=2012 and species = 'sesa'")
	v$lat_g=g$latit[match(paste(v$gps_left,v$wp_left),paste(g$gps,g$wp))]
	v$lon_g=g$longit[match(paste(v$gps_left,v$wp_left),paste(g$gps,g$wp))]
	v=v[!is.na(v$lat_g),]
	#v=v[-which(v$pk%in%c(45,259,607,1385)),]
	v$lat=n$latit[match(v$nest,tolower(n$nest))]
	v$lon=n$longit[match(v$nest,tolower(n$nest))]
	v$sp=n$species[match(v$nest,tolower(n$nest))]
		for(i in 1:nrow(v) ) { 
			if(is.na(v$lat[i])){v$distm[i]==NA} else {
					v[i,'distm'] = spDists( as.matrix(v[i, c("lon_g", "lat_g")]), v[i, c("lon", "lat")], longlat=TRUE)*1000
						}
						}
	v=v[v$distm<1250,]					
	v$issues=0				
	v1=v
	# add nests where dist estimated
	v=dbGetQuery(con_barr12,"SELECT visits_pk pk,datetime_, nest, gps_left, wp_left, dist_left_nest dist_est FROM VISITS where dist_left_nest is not null and wp_left is NULL")
	g=dbGetQuery(con_barr12,"SELECT*FROM GPS")
	n=dbGetQuery(conMy,"SELECT nest, species, latit, longit FROM avesatbarrow.nests where year_=2012")
	v$lat_g=g$latit[match(paste(v$gps_left,v$wp_left),paste(g$gps,g$wp))]
	v$lon_g=g$longit[match(paste(v$gps_left,v$wp_left),paste(g$gps,g$wp))]
	
	v$lat=n$latit[match(v$nest,tolower(n$nest))]
	v$lon=n$longit[match(v$nest,tolower(n$nest))]
	v$sp=n$species[match(v$nest,tolower(n$nest))]
	
	v$distm=NA
	v$issues=0
	v2=v
	
		v=rbind(v1,v2)
		v$issues[!is.na(v$dist_est) & !is.na(v$distm)]=1
		v_=v[v$issues==1,]
		ggplot(v_,aes(y=v_$dist_est,x=v_$distm))+geom_point()+stat_smooth(method="lm")
	v2$distm=v2$dist_est
	v1$type="gps"
	v2$type="est"
	v=rbind(v1,v2)
	v = v[!is.na(v$datetime_),]
	v[nchar(v$datetime_)!=nchar('2012-06-12 20:36:00'),]
	v$datetime_ = as.POSIXct(v$datetime_, tz = 'America/Anchorage')		
	
	load("M:\\Science\\Projects\\PHD\\STATS\\DATASETS\\C_SESA_inc_start.Rdata")
		d1 = d
	load("M:\\Science\\Projects\\PHD\\STATS\\DATASETS\\C_SESA_short_inc_start.Rdata")
	d = rbind(d1,d)
	d$nest= tolower(d$nest)	
	v = v[v$nest %in% unique(d$nest),]	
			l = list()
			for(i in 1:nrow(v)){
							fi=v[i,]
							di=d[which(d$nest==fi$nest & fi$datetime_>=d$bout_start & fi$datetime_<=d$bout_end+60),] # adding 1 min to the bout end to make the selection more robust
							if(nrow(di)==0){#fi$sex=NA
											 print(fi) #all these visits where before we had rfid on and hence sex cannot be determined
											}else{	fi$sex=di$sex
													fi$bird_ID_filled = di$bird_ID_filled
													l[[i]]=fi
													print(i)
												 }
							
							
							}
	v=do.call(rbind,l)		
	v$year = 2012	
	v12 =v
}
 {# combine
	names(v11)[names(v11) == 'dist_est'] = 'distm'	
	v11$type = 'est'
	v12 = v12[v12$sp == 'SESA',]
	v = rbind(v12[c('year','datetime_','nest','sex', 'bird_ID_filled','distm','type')], v11)
	nrow(v)
	save(v, file = paste(wd,'escpate_2011-2012.Rdata', sep=''))
	}

 {# analyses
	load(file = paste(wd,'escpate_2011-2012.Rdata', sep=''))
	nrow(v)
	#nrow(v[is.na(v$bird_ID_filled),])
	#v[v$nest == 's510',]
	v$bird_ID_filled[is.na(v$bird_ID_filled)] = 'f'
	v$nest_sex = paste(v$nest, v$sex)
	nrow(v[duplicated(paste(v$nest_sex,v$datetime_)),])
	#v[paste(v$nest_sex,v$datetime_) %in% paste(v$nest_sex,v$datetime_)[duplicated(paste(v$nest_sex,v$datetime_))],]
	#v = v[!duplicated(paste(v$nest_sex,v$datetime_)),]
	#vg = v[v$type == 'gps',]	
	#ve = v[!paste(v$nest_sex, v$datetime_)%in%paste(vg$nest_sex, vg$datetime_),]
	#v = rbind(vg,ve)
	#nrow(v)
	
	x = data.frame(table(v$bird_ID_filled))
	summary(factor(x$Freq))
	summary(x$Freq)
	
	densityplot(~v$distm)
	
	v = v[v$distm<200,]
	
	m = lmer(distm ~ 1 + (1|bird_ID_filled), v)
	m = lmer(distm ~ 1 + (1|bird_ID_filled), v[v$type=='est',])
	summary(m)
	
	xx = x$Var1[x$Freq>1]
	m = lmer(distm ~ 1 + (1|bird_ID_filled), v[v$bird_ID_filled%in%xx,])
		summary(m)
	table(v$nest_sex[v$nest_sex%in%xx])
	
	vv = v
	vv$time = as.numeric(difftime(vv$datetime_,trunc(vv$datetime_,"day"), units = "hours"))
	vv$day_j=as.numeric(format(as.Date(trunc(vv$datetime_, "day")),"%j")) - as.numeric(format(as.Date(trunc(min(vv$datetime_), "day")),"%j"))+1
	ggplot(vv,aes(x = day_j, y = distm)) + geom_point() +stat_smooth(method = 'lm')
		 load("M:\\Science\\Projects\\PHD\\STATS\\DATASETS\\inc_start_est_2011_SESA_NON-SESA.Rdata")
		 e$year = 2011
		 e1 = e
		 
		 load("M:\\Science\\Projects\\PHD\\STATS\\DATASETS\\inc_start_est_2012_SESA_NON-SESA.Rdata")
		  e$year = 2012
		 e2 = e
		 e = rbind(e1,e2)
		 e = e[!is.na(e$inc_est),]
		 e$inc_est = as.POSIXct(e$inc_est,format="%Y-%m-%d %H:%M:%S", tz = 'America/Anchorage')		
		  e$nest = tolower(e$nest)
		vv$inc_start = e$inc_est[match(paste(vv$year, vv$nest),paste(e$year,e$nest))] 	
		vv$day_inc=as.numeric(difftime(vv$datetime_,vv$inc_start, units = "days"))	
	ggplot(vv,aes(x = day_inc, y = distm)) + geom_point() +stat_smooth(method = 'lm')
	ggplot(vv,aes(x = day_inc, y = distm)) + geom_point() +stat_smooth()
	
	m = lmer(distm ~ 1 + (1|bird_ID_filled), vv)
	m = lmer(distm ~ day_j + (day_j|bird_ID_filled), vv)
	m = lmer(distm ~ day_inc + (day_inc|bird_ID_filled), vv,na.action='na.omit')
	m = lmer(distm ~ 1 + (day_inc|bird_ID_filled), vv,na.action='na.omit')
	m = lmer(distm ~ day_inc + (day_inc|bird_ID_filled), vv[vv$bird_ID_filled%in%xx,])
	m = lmer(distm ~ 1 + (day_inc|bird_ID_filled), vv[vv$bird_ID_filled%in%xx,])
	plot(allEffects(m))
	summary(glht(m))
	summary(m)
	pp = profile(m)
	confint(pp)
	library(MCMCglmm)
		mm <- MCMCglmm(distm ~ 1,random = ~us(day_inc):bird_ID_filled, data=vv[vv$bird_ID_filled%in%xx & !is.na(vv$day_inc),])
summary(mm)
}
}	

	}
	{# Experimental procedure 
		{# period between capture of the 'removed' parent and return of the 'focal' parent
			e=read.csv(file=paste(wd,'experiment_times.csv', sep=""), sep=",",stringsAsFactors =FALSE)
			e$taken = as.POSIXct(e$taken)
			e$start_ = as.POSIXct(e$start_)
			summary(as.numeric(difftime(e$start_,e$taken, units = 'hours')))
		}
		{# correlation between temperature based and rfid based constancy
			load(file=paste(wd,'bout.Rdata', sep=""))
			bb=b[-which(is.na(b$inc_eff)|is.na(b$inc_eff_2)),] # use only bouts with both constancies
			cor(bb$inc_eff,bb$inc_eff_2, method='spearman')
			cor(bb$inc_eff,bb$inc_eff_2, method='pearson')
			nrow(bb)
			
			dd_=dd[-which(is.na(dd$inc_eff)|is.na(dd$inc_eff_2)),] # use only bouts with both constancies
			cor(dd_$inc_eff,dd_$inc_eff_2, method='spearman')
			cor(dd_$inc_eff,dd_$inc_eff_2, method='pearson')
		}
	}
	
	{# Explaining the diversity in compensation - Correlation nest attendance vs compensation
				load(file=paste(wd,'experimental.Rdata', sep=""))
				b=b[-which(b$nest%in%c('s410','s514','s524','s624')),] 
				b$inc_eff[b$nest=='s806']=b$inc_eff_2[b$nest=='s806'] 
				bb=b[b$exper=='c',]
				bt=b[b$exper=='t',]
				bb$inc_eff_treated=bt$inc_eff[match(bb$nest, bt$nest)]
				bb$compensation=bb$inc_eff_treated/bb$inc_eff
				
				cor(bb$compensation,bb$inc_eff_treated,method='pearson')
				cor(bb$compensation,bb$inc_eff_treated,method='spearman')
				plot(bb$inc_eff_treated~bb$compensation)
			}
	{# Mass of removed individuals
		d =read.csv(paste(wd,'Data/captivity.csv', sep=""),stringsAsFactors=FALSE)
		d = d[d$phase == 'start' & !is.na(d$mass),] # removing NAs where fat was scored by two observes
		summary(d$mass)
		nrow(d)
	}
	{# relative mass loss in removed parents (without two dead females)
		d =read.csv(paste(wd,'captivity.csv', sep=""),stringsAsFactors=FALSE)
		dd = d[!is.na(d$mass) & d$phase %in%c('start','end'),] # removes the dead females 257188570 257188566
		length(unique(dd$ring_num))
		v = ddply(dd,.(ring_num), summarise, rel_mass = (mass[phase=='end'] - mass[phase=='start'])/mass[phase=='start'], abs_mass = (mass[phase=='end'] - mass[phase=='start']))
		v[duplicated(v$ring_num),]
		length(unique(v$ring_num))
		summary(v)
		#unique(dd$ring_num)[!unique(dd$ring_num)%in%v$ring_num]
		cor.test(v$rel_mass, v$abs_mass)
		ggplot(v,aes(x=rel_mass, y = abs_mass)) + geom_point() + stat_smooth()
		ggplot(v,aes(x=rel_mass, y = abs_mass)) + geom_point() + stat_smooth(method = 'lm')
		
		#u[!u$ID_taken%in%v$ring_num,] #257188570, 257188566
		#v[!v$ring_num%in%u$ID_taken,] 
		}

	{# semipalmated sandpiper energetic requirements
		# based on equation in Norton 1973 page 64 - given in cubic centimeter per gram and hour
					# cc = 0.001 liters (thus to get the results in literes we divide them by by 1000)
					# O2 in liters multiply by 20.10 to get KJ
			#27g SESA and 6.2C - median tundra temperature in Barrow
			27*((10.69+(-0.203)*6.2)/1000)*20.1 # kj/27g/h # 4.699782	
			24*27*((10.69+(-0.203)*6.2)/1000)*20.1 # kj/27g/d # 4.699782	
			
		# based on Ashkenazei and Safriel 1979
			# 19-59 duing incubation
	}
	{# Ethics
		# starved birds (s402, s510 deaath
			p = data.frame(bird_ID =c(257188558,257188571,257188564,255174350,257188570,257188566,255174353,257188572) , nest = c('s502','s702','s509','s711','s402','s510','s516','s404'), stringsAsFactors=FALSE)
		# number of mealworms at the end
			d =read.csv(paste(wd,'captivity.csv', sep=""),stringsAsFactors=FALSE)
			d = d[d$phase == 'end',]
			dd = d[!is.na(d$worms_num),]
			summary(dd$worms_num)
			summary(factor(dd$worms_num))
		# number of returned birds according to sex and whether they were starved
			d =read.csv(paste(wd,'experiment_metadata.csv', sep=""),stringsAsFactors=FALSE)
			d = d[d$type == 'exp',]
			d$sex=ifelse(d$ID_taken==d$IDfemale, 'f','m')
			d$starved = ifelse(d$nest%in%p$nest, 'yes','no')
			d$starved[d$comments=='died'] = 'died'
			summary(factor(d$sex))
			summary(factor(d$back))
			table(paste(d$sex,d$back),d$starved)
		# return times of the starved/non-starved parents	
	}
# RESULTS
	{# run first
		u=read.csv(paste(wd,'experiment_metadata.csv', sep=""),stringsAsFactors=FALSE)
		load(file=paste(wd,'experimental.Rdata', sep=""))
		
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

	{# Compensation for absence of parental care 
		{# Supplementary Table 1
			{# 1 - simple model 
			m=lmer(inc_eff~exper+(1|bird_ID),b, REML=FALSE) # best
			 
			pred=c('Intercept (control)','Period (treated)')
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
				# Fixed effects
					v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))
					ci=apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
					oi=data.frame(model='1',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
					rownames(oi) = NULL
						oi$estimate_r=round(oi$estimate,3)
						oi$lwr_r=round(oi$lwr,3)
						oi$upr_r=round(oi$upr,3)
						#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
					oii=oi[c('model','type',"effect", "estimate_r","lwr_r",'upr_r')]	
				# Random effects
					l=data.frame(summary(m)$varcor)
						ri=data.frame(model='1', type='random (var)',effect=l$grp, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA)
					o1=rbind(oii,ri)
			
			}
			{# 2 - model controlled for bout length
			 m2=lmer(inc_eff~exper+scale(bout_length)+(1|bird_ID),b, REML=FALSE)
			 pred=c('Intercept (control)','Period (treated)','Bout (scaled)')
						nsim <- 2000
						bsim <- sim(m2, n.sim=nsim)  
				# Fixed effects
					v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))
					ci=apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
					oi=data.frame(model='2',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
					rownames(oi) = NULL
						oi$estimate_r=round(oi$estimate,3)
						oi$lwr_r=round(oi$lwr,3)
						oi$upr_r=round(oi$upr,3)
						#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
					oii=oi[c('model','type',"effect", "estimate_r","lwr_r",'upr_r')]	
				# Random effects
					l=data.frame(summary(m)$varcor)
						ri=data.frame(model='2',type='random (var)',effect=l$grp, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA)
					o2=rbind(oii,ri)
			
						
			}	
			{# 3 - model with sex 
			 m3=lmer(inc_eff~exper*sex+(1|bird_ID),b, REML=FALSE)
			 pred=c('Intercept (control)','Period (treated)','Sex(male)','Period:sex')
						nsim <- 2000
						bsim <- sim(m3, n.sim=nsim)  
				# Fixed effects
					v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))
					ci=apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
					oi=data.frame(model='3',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
					rownames(oi) = NULL
						oi$estimate_r=round(oi$estimate,3)
						oi$lwr_r=round(oi$lwr,3)
						oi$upr_r=round(oi$upr,3)
						#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
					oii=oi[c('model','type',"effect", "estimate_r","lwr_r",'upr_r')]	
				# Random effects
					l=data.frame(summary(m)$varcor)
						ri=data.frame(model='3',type='random (var)',effect=l$grp, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA)
					o3=rbind(oii,ri)

			}				
			
	
			{# combine and export to excel table
						sname = tempfile(fileext='.xls')
						wb = loadWorkbook(sname,create = TRUE)	
						createSheet(wb, name = "output")
						writeWorksheet(wb, rbind(o1,o2,o3), sheet = "output")
						createSheet(wb, name = "output_AIC")
						writeWorksheet(wb, rbind(o), sheet = "output_AIC")
						saveWorkbook(wb)
						shell(sname)
			}
			
			   densityplot(b$inc_eff)
			{# model assumptions simple
						dev.new(width=6,height=9)
						
						par(mfrow=c(4,3))
						m=lmer(inc_eff~exper+(1|bird_ID),b, REML=FALSE)
						
						scatter.smooth(fitted(m),resid(m),col='red');abline(h=0, lty=2)
									 
						scatter.smooth(fitted(m),sqrt(abs(resid(m))), col='red')
						
							
						qqnorm(resid(m), main=list("Normal Q-Q Plot: residuals", cex=0.8),col='red') 
						qqline(resid(m))
			 
						qqnorm(unlist(ranef(m)$bird_ID[1]), main = " Bird")
						qqline(unlist(ranef(m)$bird_ID[1]))
						
						
						scatter.smooth(resid(m)~b$exper);abline(h=0, lty=2, col='red')
						plot(resid(m)~b$exper);abline(h=0, lty=2, col='red')
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=b$lon, y=b$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
				
				}
			{# model assumptions 2
						dev.new(width=6,height=9)
						
						par(mfrow=c(4,3))
						
						m=lmer(inc_eff~exper+scale(bout_length)+(1|bird_ID),b, REML=FALSE)
						scatter.smooth(fitted(m),resid(m),col='red');abline(h=0, lty=2)
						scatter.smooth(fitted(m),sqrt(abs(resid(m))), col='red')
						
							
						qqnorm(resid(m), main=list("Normal Q-Q Plot: residuals", cex=0.8),col='red') 
						qqline(resid(m))
			 
						qqnorm(unlist(ranef(m)$bird_ID[1]), main = " Bird")
						qqline(unlist(ranef(m)$bird_ID[1]))
						
						
						scatter.smooth(resid(m)~b$exper);abline(h=0, lty=2, col='red')
						plot(resid(m)~b$exper);abline(h=0, lty=2, col='red')
						scatter.smooth(resid(m)~scale(b$bout_length));abline(h=0, lty=2, col='red')
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=b$lon, y=b$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
				
				}
			{# model assumptions 3
						dev.new(width=6,height=9)
						par(mfrow=c(4,3))
						m=lmer(inc_eff~exper*sex+(1|bird_ID),b, REML=FALSE)
						
						scatter.smooth(fitted(m),resid(m),col='red');abline(h=0, lty=2)
						scatter.smooth(fitted(m),sqrt(abs(resid(m))), col='red')
						
							
						qqnorm(resid(m), main=list("Normal Q-Q Plot: residuals", cex=0.8),col='red') 
						qqline(resid(m))
			 
						qqnorm(unlist(ranef(m)$bird_ID[1]), main = " Bird")
						qqline(unlist(ranef(m)$bird_ID[1]))
						
						
						scatter.smooth(resid(m)~b$exper);abline(h=0, lty=2, col='red')
						plot(resid(m)~b$exper);abline(h=0, lty=2, col='red')
						plot(resid(m)~b$sex);abline(h=0, lty=2, col='red')
						plot(resid(m)~factor(paste(b$exper,b$sex,sep=" ")));abline(h=0, lty=2, col='red')
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=b$lon, y=b$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
				
				}
		}
		{# % compensation and probability that compensation is higher than 0.5
				  m=lmer(inc_eff~exper+(1|nest),b, REML=FALSE) 
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
					# % compensation
					quantile(((bsim@fixef[,1]+bsim@fixef[,2])/bsim@fixef[,1]),probs=c(0.025,0.5,0.975))
								mean((bsim@fixef[,1]+bsim@fixef[,2])/bsim@fixef[,1])
					# probability			
					sum(((bsim@fixef[,1]+bsim@fixef[,2])/bsim@fixef[,1])>0.5)/nsim
					p=((bsim@fixef[,1]+bsim@fixef[,2])/bsim@fixef[,1])
					sum(length(p[p>0.55 & p<0.65]))/nsim
				
						
		}
		{# distribution of compensation
			summary(100*bt$inc_eff/bb$inc_eff)
			dev.new(width=3.5/2,height=1.85)
			par(mar=c(2.2,2.1,0.5,0.1), ps=12, mgp=c(1.2,0.35,0), las=1, cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE, col.axis="grey30",font.main = 1, col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt
			dd <- density(100*bt$inc_eff/bb$inc_eff) # returns the density data
			plot(dd) # plots the results 
			
			densityplot(~100*bt$inc_eff/bb$inc_eff)
			densityplot(~100*bt$inc_eff/bb$inc_eff, groups=as.character(bt$sex))
			histogram(~100*bt$inc_eff/bb$inc_eff, breaks = 9)
			histogram(~100*bt$inc_eff/bb$inc_eff, type = 'count', breaks = 9)
			hist(bt$inc_eff/bb$inc_eff, breaks=10, xlim=c(0,1.1))
			plot(density(100*bt$inc_eff/bb$inc_eff), xlim=c(0,102))
		}	
	
		{# Figure 2
			{# run first
				k=0.1 # distance
				b$exper_sex=ifelse(b$exper=='c',ifelse(b$sex=='f', 1-k,k+1),ifelse(b$sex=='f', 4-k,k+4))
							{# model predictions
				 m=lmer(inc_eff~exper+(1|nest),b, REML=FALSE) # best
						pred=c('Intercept (control)','Period (treated)')
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
				
				# coefficients
					v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
				
				# values to predict for		
					newD=data.frame(exper=c('c','t'))
						
				# exactly the model which was used has to be specified here 
					X <- model.matrix(~ exper,data=newD)	
								
				# calculate predicted values and creditability intervals
					newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
							predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
							for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@fixef[i,]
							newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
							newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
					pp=newD	
					pt_=pp[pp$exper=='t',]
					pc_=pp[pp$exper=='c',]
				}
			
			}	
			{# plot dist, boxplot - no y-axis labels
			dev.new(width=3.5,height=1.85)
			#png(paste(out2,"Figure_2.png", sep=""), width=3.5,height=1.85,units="in",res=600)
			#png(paste(out2,"Figure_2_dist-box_in_sex.png", sep=""), width=3.5,height=1.85,units="in",res=600)
			#jpeg(paste(out2,"Figure_2_dist-box_smaller_font.jpeg", sep=""), width=3.5,height=1.85,units="in",res=100)
			{# plot distribution
				par(mfrow=c(1,2),mar=c(2.2,2.1,0.5,0.1), ps=12, mgp=c(1.2,0.35,0), las=1, cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE, col.axis="grey30",font.main = 1, col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt
			
				plot(bb$inc_eff_treated~bb$inc_eff, col=bb$col_,bg=adjustcolor(bb$col_, alpha.f = 0.4), pch=21,xlim=c(0,1), ylim=c(0,1),  ylab=NA,xlab=NA, xaxt='n',yaxt='n', type='n')
					lines(c(0,1),c(0,1), lty=3, col="lightgrey")
					axis(1, at=seq(0,1,by=0.25), labels=FALSE)
					axis(2, at=seq(0,1,by=0.25), labels=FALSE)
					axis(2, at=seq(0,1,by=0.5), labels=c('0.0','0.5','1.0'))
					text(seq(0,1,by=0.5), par("usr")[3]-0.08, labels = c('0.0','0.5','1.0'),  xpd = TRUE, cex=0.5, col="grey30") #labels = z_g$genus
					
					#mtext('Control period',side=1,line=0.6, cex=0.6, las=1, col='grey30')
					#mtext('Incubation constancy',side=1,line=0.6, cex=0.6, las=1, col='grey30')
					mtext('Treated period',side=2,line=-0.7, cex=0.5, las=3, col='grey45')
					mtext('Nest attendance',side=2,line=1, cex=0.6, las=3, col='grey30')
					
					mtext('Control period',side=1,line=-1, cex=0.5, las=1, col='grey45')
					mtext('Nest attendance',side=1,line=0.4, cex=0.6, las=1, col='grey30')
					#text(c(0.5), par("usr")[3]-0.18, labels = c('Control period'),  xpd = TRUE, cex=0.5, col="grey30")
					#text(c(0.5), par("usr")[3]-0.28, labels = c('Incubation constancy'),  xpd = TRUE, cex=0.6, col="grey30")
					
					text(x=0.05,y=0.95, labels='\u2640', col='#FCB42C', cex=0.6)
					text(x=0.12,y=0.98, labels='\u2642', col='#535F7C', cex=0.6)
					
					mtext(expression(bold("a")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
					# data
					points(bb$inc_eff_treated~bb$inc_eff, col=bb$col_,bg=adjustcolor( bb$col_, alpha.f = 0.4), pch=21, cex=0.5)
				
				# predictions	
					# 95%CI for control
					arrows(x0=pc_$lwr, y0=pt_$pred,x1=pc_$upr, y1=pt_$pred,
						code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)
					# 95%CI for treated
					arrows(x0=pc_$pred, y0=pt_$lwr,x1=pc_$pred, y1=pt_$upr,
						code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)	
					
					#points(y=(0.939-0.384),x=0.939, pch=21, cex=0.8,col="white",bg="white")	
					
					points(y=pt_$pred,x=pc_$pred, pch=20, cex=0.9,col="red")
					#points(y=(0.939-0.384),x=0.939, pch=21, cex=0.8,col="red",bg=adjustcolor("red", alpha.f = 0.4))	
				
				#legend("topleft", legend=c('\u2640','\u2642'), pch=21, col=c( '#FCB42C','#535F7C'),bg=adjustcolor(c( '#FCB42C','#535F7C')), cex=0.6, bty='n')
				}
			{# boxplot - no points 
				par(mar=c(2.2,1,0.5,1.2))
					
					boxplot(inc_eff ~ exper_sex, data = b,
										ylab = NULL,xaxt='n', yaxt='n',
										#at=c(1,2,3.5,4.5),
										at=c(1,2,3.7,4.7), 
										type='n',
										outcex=0.5,outpch=20,boxwex=0.5,whisklty=1,staplelty=0,#medlwd=1, 
										lwd = 1,
										border=c('#FCB42C','#535F7C','#FCB42C','#535F7C'),
										col = adjustcolor("white", alpha.f = 0), # trick for PNGs, to show what is underneath the boxplot else can be taken out
										#outcol="darkgrey",boxcol='darkgrey',whiskcol='darkgrey',staplecol='darkgrey',medcol='darkgrey', 
										par(bty='l')
										)					
					
					axis(1, at=c(1.5,4.2), labels=FALSE)
					
					text(c(1,2), par("usr")[3]+0.07, labels = c('\u2640','\u2642'), font=4, xpd = TRUE, cex=0.6, col=c('#FCB42C','#535F7C'))#col="grey30") #labels 
					text(c(1.5,4.2), par("usr")[3]-0.08, labels = c('Control','Treated'),  xpd = TRUE, cex=0.5, col="grey30")
					text(c(2.85), par("usr")[3]-0.18, labels = c('Period'),  xpd = TRUE, cex=0.6, col="grey30")
					axis(2, at=seq(0,1,by=0.25), labels=FALSE)
					
					mtext(expression(bold("b")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
					
					# predictions	
						points(y=pp$pred,x=c(1.5,4.2), pch=20, cex=0.9,col="red")
					# 95%CI 
						arrows(x0=c(1.5,4.2), y0=pp$lwr,x1=c(1.5,4.2), y1=pp$upr,
						code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)
								
					
					}


			dev.off()
			}
			
			{# not used
						{# boxplot - points grey, box color
				par(mar=c(2.2,1,0.5,1.2))
					
				boxplot(inc_eff ~ exper_sex, data = b, 
										 ylab =NULL, xaxt='n',yaxt='n',
										par(bty='n'),
										#at=c(1,2,3.5,4.5),
										at=c(1,2,3.7,4.7), type='n',
										outcex=0.5, outpch=20,boxwex=0.5,whisklty=1,staplelty=0,#medlwd=1, 
										lwd = 0.25, ylim=c(0,1),
										outcol="darkgrey",boxcol='darkgrey',whiskcol='darkgrey',staplecol='darkgrey',medcol='darkgrey'
										) # col=z_g$cols, border=z_g$cols
										
					stripchart(inc_eff ~ factor(exper_sex), vertical = TRUE, data = b, method = "jitter", add = TRUE, 
										at=c(1,2,3.7,4.7),
										pch = 21,cex=0.5, 
										col="gray63",
										bg=adjustcolor("gray63", alpha.f = 0.4)
										)
										
					boxplot(inc_eff ~ exper_sex, data = b,
										ylab = NULL,xaxt='n', yaxt='n',
										#at=c(1,2,3.5,4.5),
										at=c(1,2,3.7,4.7), 
										type='n',
										outcex=0.5,outpch=20,boxwex=0.5,whisklty=1,staplelty=0,#medlwd=1, 
										lwd = 1,
										border=c('#FCB42C','#535F7C','#FCB42C','#535F7C'),
										col = adjustcolor("white", alpha.f = 0), # trick for PNGs, to show what is underneath the boxplot else can be taken out
										#outcol="darkgrey",boxcol='darkgrey',whiskcol='darkgrey',staplecol='darkgrey',medcol='darkgrey', 
										par(bty='l'),
										add=TRUE
										)					
					
					axis(1, at=c(1.5,4.2), labels=FALSE)
					
					text(c(1,2), par("usr")[3]+0.07, labels = c('\u2640','\u2642'), font=4, xpd = TRUE, cex=0.6, col=c('#FCB42C','#535F7C'))#col="grey30") #labels 
					text(c(1.5,4.2), par("usr")[3]-0.08, labels = c('Control','Treated'),  xpd = TRUE, cex=0.5, col="grey30")
					text(c(2.85), par("usr")[3]-0.18, labels = c('Period'),  xpd = TRUE, cex=0.6, col="grey30")
					axis(2, at=seq(0,1,by=0.25), labels=FALSE)
					
					mtext(expression(bold("b")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
					
					# predictions	
						points(y=pp$pred,x=c(1.5,4.2), pch=20, cex=0.9,col="red")
					# 95%CI 
						arrows(x0=c(1.5,4.2), y0=pp$lwr,x1=c(1.5,4.2), y1=pp$upr,
						code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)
								
					
					}

			{# plot dist, boxplot - no y-axis labels, 
			dev.new(width=3.5,height=1.85)
			#png(paste(out2,"Figure_2_dist-box_smaller_font.png", sep=""), width=3.5,height=1.85,units="in",res=600)
			{# plot distribution
				par(mfrow=c(1,2),mar=c(2.2,2.1,0.5,0.1), ps=12, mgp=c(1.2,0.35,0), las=1, cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE, col.axis="grey30",font.main = 1, col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt
			
				plot(bb$inc_eff_treated~bb$inc_eff, col=bb$col_,bg=adjustcolor(bb$col_, alpha.f = 0.4), pch=21,xlim=c(0,1), ylim=c(0,1),  ylab=NA,xlab=NA, xaxt='n',yaxt='n', type='n')
					lines(c(0,1),c(0,1), lty=3, col="lightgrey")
					axis(1, at=seq(0,1,by=0.25), labels=FALSE)
					axis(2, at=seq(0,1,by=0.25), labels=FALSE)
					axis(2, at=seq(0,1,by=0.5), labels=c('0.0','0.5','1.0'))
					text(seq(0,1,by=0.5), par("usr")[3]-0.08, labels = c('0.0','0.5','1.0'),  xpd = TRUE, cex=0.5, col="grey30") #labels = z_g$genus
					text(c(0.5), par("usr")[3]-0.18, labels = c('Control period'),  xpd = TRUE, cex=0.5, col="grey30")
					text(c(0.5), par("usr")[3]-0.28, labels = c('Nest attendance'),  xpd = TRUE, cex=0.6, col="grey30")
					#mtext('Control period',side=1,line=0.6, cex=0.6, las=1, col='grey30')
					#mtext('Incubation constancy',side=1,line=0.6, cex=0.6, las=1, col='grey30')
					mtext('Treated period',side=2,line=1, cex=0.5, las=3, col='grey30')
					mtext('Nest attendance',side=2,line=1.6, cex=0.6, las=3, col='grey30')
					
					text(x=0.05,y=0.95, labels='\u2640', col='#FCB42C', cex=0.7)
					text(x=0.12,y=0.98, labels='\u2642', col='#535F7C', cex=0.7)
					
					mtext(expression(bold("a")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
					# data
					points(bb$inc_eff_treated~bb$inc_eff, col=bb$col_,bg=adjustcolor( bb$col_, alpha.f = 0.4), pch=21, cex=0.5)
				
				# predictions	
					# 95%CI for control
					arrows(x0=pc_$lwr, y0=pt_$pred,x1=pc_$upr, y1=pt_$pred,
						code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)
					# 95%CI for treated
					arrows(x0=pc_$pred, y0=pt_$lwr,x1=pc_$pred, y1=pt_$upr,
						code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)	
					
					#points(y=(0.939-0.384),x=0.939, pch=21, cex=0.8,col="white",bg="white")	
					
					points(y=pt_$pred,x=pc_$pred, pch=20, cex=0.9,col="red")
					#points(y=(0.939-0.384),x=0.939, pch=21, cex=0.8,col="red",bg=adjustcolor("red", alpha.f = 0.4))	
				
				#legend("topleft", legend=c('\u2640','\u2642'), pch=21, col=c( '#FCB42C','#535F7C'),bg=adjustcolor(c( '#FCB42C','#535F7C')), cex=0.6, bty='n')
				}
			{# boxplot - points grey, box color
				par(mar=c(2.2,1.5,0.5,0.7))
					
				boxplot(inc_eff ~ exper_sex, data = b, 
										 ylab =NULL, xaxt='n',yaxt='n',
										par(bty='n'),
										#at=c(1,2,3.5,4.5),
										at=c(1,2,3.7,4.7), type='n',
										outcex=0.5, outpch=20,boxwex=0.5,whisklty=1,staplelty=0,#medlwd=1, 
										lwd = 0.25, ylim=c(0,1),
										outcol="darkgrey",boxcol='darkgrey',whiskcol='darkgrey',staplecol='darkgrey',medcol='darkgrey'
										) # col=z_g$cols, border=z_g$cols
										
					stripchart(inc_eff ~ factor(exper_sex), vertical = TRUE, data = b, method = "jitter", add = TRUE, 
										at=c(1,2,3.7,4.7),
										pch = 21,cex=0.5, 
										col="gray63",
										bg=adjustcolor("gray63", alpha.f = 0.4)
										)
										
					boxplot(inc_eff ~ exper_sex, data = b,
										ylab = NULL,xaxt='n', yaxt='n',
										#at=c(1,2,3.5,4.5),
										at=c(1,2,3.7,4.7), 
										type='n',
										outcex=0.5,outpch=20,boxwex=0.5,whisklty=1,staplelty=0,#medlwd=1, 
										lwd = 1,
										border=c('#FCB42C','#535F7C','#FCB42C','#535F7C'),
										col = adjustcolor("white", alpha.f = 0), # trick for PNGs, to show what is underneath the boxplot else can be taken out
										#outcol="darkgrey",boxcol='darkgrey',whiskcol='darkgrey',staplecol='darkgrey',medcol='darkgrey', 
										par(bty='l'),
										add=TRUE
										)					
					
					axis(1, at=c(1,2,3.7,4.7), labels=FALSE)
					text(c(1,2,3.7,4.7), par("usr")[3]-0.08, labels = c('\u2640','\u2642','\u2640','\u2642'), font=4, xpd = TRUE, cex=0.6, col=c('#FCB42C','#535F7C','#FCB42C','#535F7C'))#col="grey30") #labels 
					text(c(1.5,4.2), par("usr")[3]-0.18, labels = c('Control','Treated'),  xpd = TRUE, cex=0.5, col="grey30")
					text(c(2.85), par("usr")[3]-0.28, labels = c('Period'),  xpd = TRUE, cex=0.6, col="grey30")
					
					axis(2, at=seq(0,1,by=0.25), labels=FALSE)
									
					#mtext(expression(bold("b")),side=3,line=-0.3, cex=0.7, las=1,adj=0.02, col="grey30")
					#mtext(expression(bold("a")),side=1,line=-1.1, cex=0.7, las=1,adj=0.03, col="grey30")
					mtext(expression(bold("b")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
					# predictions	
						points(y=pp$pred,x=c(1.5,4.2), pch=20, cex=0.9,col="red")
					# 95%CI 
						arrows(x0=c(1.5,4.2), y0=pp$lwr,x1=c(1.5,4.2), y1=pp$upr,
						code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)
								
					
					}

			dev.off()
			}
						
			{# plot boxplot, dist - no y-axis labels
			dev.new(width=3.5,height=1.85)
			#png(paste(out2,"Figure_2_smaller_font.png", sep=""), width=3.5,height=1.85,units="in",res=600)
			par(mfrow=c(1,2),mar=c(2.2,2.2,0.5,0), ps=12, mgp=c(1.2,0.35,0), las=1, cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE, col.axis="grey30",font.main = 1, col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt
			{# boxplot - points grey, box color
				
					
					boxplot(inc_eff ~ exper_sex, data = b, 
										 ylab =NULL, xaxt='n',yaxt='n',
										par(bty='n'),
										#at=c(1,2,3.5,4.5),
										at=c(1,2,3.7,4.7), type='n',
										outcex=0.5, outpch=20,boxwex=0.5,whisklty=1,staplelty=0,#medlwd=1, 
										lwd = 0.25, ylim=c(0,1),
										outcol="darkgrey",boxcol='darkgrey',whiskcol='darkgrey',staplecol='darkgrey',medcol='darkgrey'
										) # col=z_g$cols, border=z_g$cols
										
					stripchart(inc_eff ~ factor(exper_sex), vertical = TRUE, data = b, method = "jitter", add = TRUE, 
										at=c(1,2,3.7,4.7),
										pch = 21,cex=0.5, 
										col="gray63",
										bg=adjustcolor("gray63", alpha.f = 0.4)
										)
										
					boxplot(inc_eff ~ exper_sex, data = b,
										ylab = 'Nest attendance',xaxt='n', 
										#at=c(1,2,3.5,4.5),
										at=c(1,2,3.7,4.7), 
										type='n',
										outcex=0.5,outpch=20,boxwex=0.5,whisklty=1,staplelty=0,#medlwd=1, 
										lwd = 1,
										border=c('#FCB42C','#535F7C','#FCB42C','#535F7C'),
										col = adjustcolor("white", alpha.f = 0), # trick for PNGs, to show what is underneath the boxplot else can be taken out
										#outcol="darkgrey",boxcol='darkgrey',whiskcol='darkgrey',staplecol='darkgrey',medcol='darkgrey', 
										par(bty='l'),
										add=TRUE
										)					
					
					axis(1, at=c(1,2,3.7,4.7), labels=FALSE)
					text(c(1,2,3.7,4.7), par("usr")[3]-0.08, labels = c('\u2640','\u2642','\u2640','\u2642'), font=2, xpd = TRUE, cex=0.6, col=c('#FCB42C','#535F7C','#FCB42C','#535F7C'))#col="grey30") #labels 
					text(c(1.5,4.2), par("usr")[3]-0.18, labels = c('Control','Treated'),  xpd = TRUE, cex=0.5, col="grey30")
					text(c(2.85), par("usr")[3]-0.28, labels = c('Period'),  xpd = TRUE, cex=0.6, col="grey30")
					
					#mtext(expression(bold("b")),side=3,line=-0.3, cex=0.7, las=1,adj=0.02, col="grey30")
					#mtext(expression(bold("a")),side=1,line=-1.1, cex=0.7, las=1,adj=0.03, col="grey30")
					mtext(expression(bold("a")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
					# predictions	
						points(y=pp$pred,x=c(1.5,4.2), pch=20, cex=0.9,col="red")
					# 95%CI 
						arrows(x0=c(1.5,4.2), y0=pp$lwr,x1=c(1.5,4.2), y1=pp$upr,
						code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)
								
					
					}
			{# plot distribution
				par(mar=c(2.2,1.5,0.5,0.7))
				plot(bb$inc_eff_treated~bb$inc_eff, col=bb$col_,bg=adjustcolor(bb$col_, alpha.f = 0.4), pch=21,xlim=c(0,1), ylim=c(0,1),  ylab=NA,xlab=NA, xaxt='n',yaxt='n', type='n')
					lines(c(0,1),c(0,1), lty=3, col="lightgrey")
					axis(1, at=seq(0,1,by=0.2), labels=FALSE)
					axis(2, at=seq(0,1,by=0.2), labels=FALSE)
					text(seq(0,1,by=0.2), par("usr")[3]-0.08, labels = seq(0,1,by=0.2),  xpd = TRUE, cex=0.5, col="grey30") #labels = z_g$genus
					text(c(0.5), par("usr")[3]-0.18, labels = c('Control period'),  xpd = TRUE, cex=0.5, col="grey30")
					text(c(0.5), par("usr")[3]-0.28, labels = c('Nest attendance'),  xpd = TRUE, cex=0.6, col="grey30")
					#mtext('Control period',side=1,line=0.6, cex=0.6, las=1, col='grey30')
					#mtext('Incubation constancy',side=1,line=0.6, cex=0.6, las=1, col='grey30')
					mtext('Treated period',side=2,line=0.3, cex=0.5, las=3, col='grey30')
					
					text(x=0.05,y=0.95, labels='\u2640', col='#FCB42C', cex=0.7)
					text(x=0.11,y=0.98, labels='\u2642', col='#535F7C', cex=0.7)
					
					mtext(expression(bold("b")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
					# data
					points(bb$inc_eff_treated~bb$inc_eff, col=bb$col_,bg=adjustcolor( bb$col_, alpha.f = 0.4), pch=21, cex=0.5)
				
				# predictions	
					# 95%CI for control
					arrows(x0=pc_$lwr, y0=pt_$pred,x1=pc_$upr, y1=pt_$pred,
						code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)
					# 95%CI for treated
					arrows(x0=pc_$pred, y0=pt_$lwr,x1=pc_$pred, y1=pt_$upr,
						code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)	
					
					#points(y=(0.939-0.384),x=0.939, pch=21, cex=0.8,col="white",bg="white")	
					
					points(y=pt_$pred,x=pc_$pred, pch=20, cex=0.9,col="red")
					#points(y=(0.939-0.384),x=0.939, pch=21, cex=0.8,col="red",bg=adjustcolor("red", alpha.f = 0.4))	
				
				#legend("topleft", legend=c('\u2640','\u2642'), pch=21, col=c( '#FCB42C','#535F7C'),bg=adjustcolor(c( '#FCB42C','#535F7C')), cex=0.6, bty='n')
				}
			dev.off()
			}
					
			{# plot boxplot, dist
			dev.new(width=3.5,height=1.85)
			png(paste(out2,"Figure_2.png", sep=""), width=3.5,height=1.85,units="in",res=600)
			par(mfrow=c(1,2),mar=c(2.2,2.2,0.5,0), ps=12, mgp=c(1.2,0.35,0), las=1, cex.lab=0.7,cex.main=0.7, cex.axis=0.6, tcl=-0.15,bty="l",xpd=TRUE, col.axis="grey30",font.main = 1, col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt
			{# boxplot - points grey, box color
				k=0.1 # distance
				b$exper_sex=ifelse(b$exper=='c',ifelse(b$sex=='f', 1-k,k+1),ifelse(b$sex=='f', 4-k,k+4))
					
					boxplot(inc_eff ~ exper_sex, data = b, 
										 ylab =NULL, xaxt='n',yaxt='n',
										par(bty='n'),
										#at=c(1,2,3.5,4.5),
										at=c(1,2,3.7,4.7), type='n',
										outcex=0.5, outpch=20,boxwex=0.5,whisklty=1,staplelty=0,#medlwd=1, 
										lwd = 0.25, ylim=c(0,1),
										outcol="darkgrey",boxcol='darkgrey',whiskcol='darkgrey',staplecol='darkgrey',medcol='darkgrey'
										) # col=z_g$cols, border=z_g$cols
										
					stripchart(inc_eff ~ factor(exper_sex), vertical = TRUE, data = b, method = "jitter", add = TRUE, 
										at=c(1,2,3.7,4.7),
										pch = 21,cex=0.5, 
										col="gray63",
										bg=adjustcolor("gray63", alpha.f = 0.4)
										)
										
					boxplot(inc_eff ~ exper_sex, data = b,
										ylab = 'Nest attendance',xaxt='n', 
										#at=c(1,2,3.5,4.5),
										at=c(1,2,3.7,4.7), 
										type='n',
										outcex=0.5,outpch=20,boxwex=0.5,whisklty=1,staplelty=0,#medlwd=1, 
										lwd = 1,
										border=c('#FCB42C','#535F7C','#FCB42C','#535F7C'),
										col = adjustcolor("white", alpha.f = 0), # trick for PNGs, to show what is underneath the boxplot else can be taken out
										#outcol="darkgrey",boxcol='darkgrey',whiskcol='darkgrey',staplecol='darkgrey',medcol='darkgrey', 
										par(bty='l'),
										add=TRUE
										)					
					
					axis(1, at=c(1,2,3.7,4.7), labels=FALSE)
					text(c(1,2,3.7,4.7), par("usr")[3]-0.08, labels = c('\u2640','\u2642','\u2640','\u2642'), font=2, xpd = TRUE, cex=0.6, col=c('#FCB42C','#535F7C','#FCB42C','#535F7C'))#col="grey30") #labels 
					text(c(1.5,4.2), par("usr")[3]-0.18, labels = c('Control','Treated'),  xpd = TRUE, cex=0.6, col="grey30")
					text(c(2.85), par("usr")[3]-0.28, labels = c('Period'),  xpd = TRUE, cex=0.7, col="grey30")
					# predictions	
						points(y=pp$pred,x=c(1.5,4.2), pch=20, cex=0.9,col="red")
					# 95%CI 
						arrows(x0=c(1.5,4.2), y0=pp$lwr,x1=c(1.5,4.2), y1=pp$upr,
						code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)
								
					
					}
						
			{# plot distribution
				plot(bb$inc_eff_treated~bb$inc_eff, col=bb$col_,bg=adjustcolor(bb$col_, alpha.f = 0.4), pch=21,xlim=c(0,1), ylim=c(0,1),  ylab=NA,xlab=NA, xaxt='n', type='n')
					lines(c(0,1),c(0,1), lty=3, col="lightgrey")
					axis(1, at=seq(0,1,by=0.2), labels=FALSE)
					text(seq(0,1,by=0.2), par("usr")[3]-0.08, labels = seq(0,1,by=0.2),  xpd = TRUE, cex=0.6, col="grey30") #labels = z_g$genus
					text(c(0.5), par("usr")[3]-0.18, labels = c('Control period'),  xpd = TRUE, cex=0.6, col="grey30")
					text(c(0.5), par("usr")[3]-0.28, labels = c('Nest attendance'),  xpd = TRUE, cex=0.7, col="grey30")
					#mtext('Control period',side=1,line=0.6, cex=0.6, las=1, col='grey30')
					#mtext('Nest attendance',side=1,line=0.6, cex=0.6, las=1, col='grey30')
					mtext('Treated period',side=2,line=1.2, cex=0.6, las=3, col='grey30')
					text(x=0.05,y=0.95, labels='\u2640', col='#FCB42C', cex=0.6)
					text(x=0.11,y=0.98, labels='\u2642', col='#535F7C', cex=0.6)
					# data
					points(bb$inc_eff_treated~bb$inc_eff, col=bb$col_,bg=adjustcolor( bb$col_, alpha.f = 0.4), pch=21, cex=0.5)
				
				# predictions	
					# 95%CI for control
					arrows(x0=pc_$lwr, y0=pt_$pred,x1=pc_$upr, y1=pt_$pred,
						code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)
					# 95%CI for treated
					arrows(x0=pc_$pred, y0=pt_$lwr,x1=pc_$pred, y1=pt_$upr,
						code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)	
					
					#points(y=(0.939-0.384),x=0.939, pch=21, cex=0.8,col="white",bg="white")	
					
					points(y=pt_$pred,x=pc_$pred, pch=20, cex=0.9,col="red")
					#points(y=(0.939-0.384),x=0.939, pch=21, cex=0.8,col="red",bg=adjustcolor("red", alpha.f = 0.4))	
				
				#legend("topleft", legend=c('\u2640','\u2642'), pch=21, col=c( '#FCB42C','#535F7C'),bg=adjustcolor(c( '#FCB42C','#535F7C')), cex=0.6, bty='n')
				}
			dev.off()
			}
			
			{# plot dist, boxplot
			dev.new(width=3.5,height=1.75)
			#png(paste(out2,"constancy_treated_control_relationship.png", sep=""), width=2.5,height=2.5,units="in",res=600)
			par(mfrow=c(1,2),mar=c(1.5,2.2,0.7,0), ps=12, mgp=c(1.2,0.35,0), las=1, cex.lab=0.7,cex.main=0.7, cex.axis=0.6, tcl=-0.15,bty="l",xpd=TRUE, col.axis="grey30",font.main = 1, col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt
			
			{# plot distribution
			plot(bb$inc_eff_treated~bb$inc_eff, col=bb$col_,bg=adjustcolor( bb$col_, alpha.f = 0.4), pch=21,xlim=c(0,1), ylim=c(0,1), ylab='Treated period', xlab=NA, xaxt='n', type='n')
				title("Nest attendance", line = 0.25)
				lines(c(0,1),c(0,1), lty=3, col="lightgrey")
				axis(1, at=seq(0,1,by=0.2), labels=FALSE)
				text(seq(0,1,by=0.2), par("usr")[3]-0.08, labels = seq(0,1,by=0.2),  xpd = TRUE, cex=0.6, col="grey30") #labels = z_g$genus
				mtext('Control period',side=1,line=0.6, cex=0.7, las=1, col='grey30')
				text(x=0.05,y=0.95, labels='\u2640', col='#FCB42C', cex=0.6)
				text(x=0.11,y=0.98, labels='\u2642', col='#535F7C', cex=0.6)
				# data
					points(bb$inc_eff_treated~bb$inc_eff, col=bb$col_,bg=adjustcolor( bb$col_, alpha.f = 0.4), pch=21, cex=0.5)
				
				# predictions	
					# 95%CI for control
					arrows(x0=pc_$lwr, y0=pt_$pred,x1=pc_$upr, y1=pt_$pred,
						code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)
					# 95%CI for treated
					arrows(x0=pc_$pred, y0=pt_$lwr,x1=pc_$pred, y1=pt_$upr,
						code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)	
					
					#points(y=(0.939-0.384),x=0.939, pch=21, cex=0.8,col="white",bg="white")	
					
					points(y=pt_$pred,x=pc_$pred, pch=20, cex=0.9,col="red")
					#points(y=(0.939-0.384),x=0.939, pch=21, cex=0.8,col="red",bg=adjustcolor("red", alpha.f = 0.4))	
				
				#legend("topleft", legend=c('\u2640','\u2642'), pch=21, col=c( '#FCB42C','#535F7C'),bg=adjustcolor(c( '#FCB42C','#535F7C')), cex=0.6, bty='n')
				}
			{# boxplot - points grey, box color
				k=0.1 # distance
				b$exper_sex=ifelse(b$exper=='c',ifelse(b$sex=='f', 1-k,k+1),ifelse(b$sex=='f', 4-k,k+4))
					
					boxplot(inc_eff ~ exper_sex, data = b, 
										 ylab =NULL, xaxt='n',yaxt='n',
										par(bty='n'),
										#at=c(1,2,3.5,4.5),
										at=c(1,2,3.7,4.7), type='n',
										outcex=0.5, outpch=20,boxwex=0.5,whisklty=1,staplelty=0,#medlwd=1, 
										lwd = 0.25,
										outcol="darkgrey",boxcol='darkgrey',whiskcol='darkgrey',staplecol='darkgrey',medcol='darkgrey'
										) # col=z_g$cols, border=z_g$cols
										
					stripchart(inc_eff ~ factor(exper_sex), vertical = TRUE, data = b, method = "jitter", add = TRUE, 
										at=c(1,2,3.7,4.7),
										pch = 21,cex=0.5, 
										col="gray63",
										bg=adjustcolor("gray63", alpha.f = 0.4)
										)
										
					boxplot(inc_eff ~ exper_sex, data = b,
										ylab = 'Nest attendance',xaxt='n', 
										#at=c(1,2,3.5,4.5),
										at=c(1,2,3.7,4.7), 
										type='n',
										outcex=0.5, outpch=20,boxwex=0.5,whisklty=1,staplelty=0,#medlwd=1, 
										lwd = 0.25,
										border=c('#FCB42C','#535F7C','#FCB42C','#535F7C'),
										#outcol="darkgrey",boxcol='darkgrey',whiskcol='darkgrey',staplecol='darkgrey',medcol='darkgrey', 
										par(bty='l'),
										add=TRUE
										)					
					
					axis(1, at=c(1,2,3.7,4.7), labels=FALSE)
					text(c(1,2,3.7,4.7), par("usr")[3]-0.08, labels = c('\u2640','\u2642','\u2640','\u2642'),  xpd = TRUE, cex=0.6, col=c('#FCB42C','#535F7C','#FCB42C','#535F7C'))#col="grey30") #labels 
					text(c(1.5,4.2), par("usr")[3]-0.16, labels = c('Control','Treated'),  xpd = TRUE, cex=0.6, col="grey30")
					
					# predictions	
						points(y=pp$pred,x=c(1.5,4.2), pch=20, cex=0.9,col="red")
					# 95%CI 
						arrows(x0=c(1.5,4.2), y0=pp$lwr,x1=c(1.5,4.2), y1=pp$upr,
						code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)
								
					
					}
			
			{# boxplot - points color, box black
				k=0.1 # distance
				b$exper_sex=ifelse(b$exper=='c',ifelse(b$sex=='f', 1-k,k+1),ifelse(b$sex=='f', 4-k,k+4))
					
					boxplot(inc_eff ~ exper_sex, data = b, 
										 ylab =NULL, xaxt='n',yaxt='n',
										par(bty='n'),
										#at=c(1,2,3.5,4.5),
										at=c(1,2,3.7,4.7), type='n',
										outcex=0.5, outpch=20,boxwex=0.5,whisklty=1,staplelty=0,#medlwd=1, 
										lwd = 0.25,
										outcol="darkgrey",boxcol='darkgrey',whiskcol='darkgrey',staplecol='darkgrey',medcol='darkgrey'
										) # col=z_g$cols, border=z_g$cols
										
					stripchart(inc_eff ~ factor(exper_sex), vertical = TRUE, data = b, method = "jitter", add = TRUE, 
										at=c(1,2,3.7,4.7),
										pch = 20, 
										col = adjustcolor(c('#FCB42C','#535F7C','#FCB42C','#535F7C'),alpha.f=0.4)
										)
					boxplot(inc_eff ~ exper_sex, data = b,
										ylab = 'Nest attendance',xaxt='n', 
										#at=c(1,2,3.5,4.5),
										at=c(1,2,3.7,4.7), 
										type='n',
										outcex=0.5, outpch=20,boxwex=0.5,whisklty=1,staplelty=0,#medlwd=1, 
										lwd = 0.25,
										border='black',#outcol="darkgrey",boxcol='darkgrey',whiskcol='darkgrey',staplecol='darkgrey',medcol='darkgrey', 
										par(bty='l'),
										add=TRUE
										)					
					
					axis(1, at=c(1,2,3.7,4.7), labels=FALSE)
					text(c(1,2,3.7,4.7), par("usr")[3]-0.08, labels = c('\u2640','\u2642','\u2640','\u2642'),  xpd = TRUE, cex=0.6, col=c('#FCB42C','#535F7C','#FCB42C','#535F7C'))#col="grey30") #labels 
					text(c(1.5,4.2), par("usr")[3]-0.16, labels = c('Control','Treated'),  xpd = TRUE, cex=0.6, col="grey30")
					}
			}						
			}
		}		
		{# Figure 3
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

			{# all in plot with facets
				{# prepare nest datasets 	
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
				{# plot 7x4 with labels -with constancy lables and sex label, black and grey line
					nn$label=paste(nn$bb,nn$bt,"   ",nn$sex_symbol,sep="        ")
					nn$label=paste("     ",nn$bb,"        ",nn$bt,"                    ",nn$sex_symbol,sep="")
					
					# label constancy
					nnn=nn$label
					names(nnn)=nn$n
					
					# label nest
					nnnest=nn$nest
					names(nnnest)=nn$n
				
				
				#dev.new(width=7.4, height=4.35) # for 4*7
				dev.new(width=7, height=4.35) # for 4*7
				  ggplot(u,aes(x=datetime_c,y=roll_con, fill=col_line))+
						#stat_smooth(method="lm", fill="grey95",alpha=1, aes(colour=genus_))+
						#scale_colour_manual(values=go$cols, name="Genus")+
						#facet_wrap(~n, nrow=5, labeller=as_labeller(nnnest))+ #label is nests
						facet_wrap(~n, nrow=5, labeller=as_labeller(nnn))+ #label is constancy
						geom_point(aes(color=col_line), shape=20,size=0.01)+#colour="grey38",aes(fill=genus_))+
						#geom_line(aes(color=col_line), size=0.05)+#colour="grey38",aes(fill=genus_))+
						scale_color_manual(values=c("black","gray90"), guide=FALSE)+
						#scale_fill_manual(values=go$cols, name="Genus")+ 
						labs(x="Time [h]", y="Nest attendance")+
						scale_x_continuous(limits=c(-16,33), breaks=seq(-15,30,by=15), labels=c(seq(-15,30,by=15)))+
						geom_vline(xintercept=0, lty=3,lwd=0.5, col='red')+
						#geom_vline(data=e,aes(xintercept=r_time), lty=3,lwd=0.5, col='red')+
						geom_text(data = ann_text,aes(label =lab),size=1.9, angle = 90, colour='grey50',hjust='middle',vjust='middle')+
						#scale_x_continuous(breaks=15, labels=15)+
						scale_y_continuous(limits=c(-0.05,1.05),breaks=seq(0,1,by=0.25), labels=c('0.0','','0.5','','1.0'))+
						#coord_cartesian(ylim=c(0,1), xlim=c(-15,15))+ #coord_cartesian(ylim=c(-5,61))+ # needs to be 10 if it is not to touch x-axis
						theme_light()+
						theme(	axis.line=element_line(colour="grey70", size=0.25),
								#panel.border=element_rect(colour="white"),
								panel.border=element_rect(colour="grey70", size=0.25),
								panel.grid = element_blank(),
								
								axis.title=element_text(size=7, colour="grey30"),
								axis.title.y = element_text(vjust=1),
								axis.title.x = element_text(vjust=0.2),
								axis.text=element_text(size=6),# margin=units(0.5,"mm")),
								axis.ticks.length=unit(0.5,"mm"),
								#axis.ticks.margin,
								
								strip.text.x =element_text(size = 6, color="grey30",  margin=margin(1,1,1,1,"mm")), #grey50
								strip.background=element_rect(fill="grey99",colour="grey70", size=0.25),
								#strip.background = element_blank(), 
								#strip.text = element_blank(),
								panel.margin = unit(0, "mm"),
								legend.position="none"
								#legend.background=element_rect(colour="white"),
								#legend.key=element_rect(fill="grey99", colour="white"),
								#legend.text=element_text(size=7, colour="grey30"),
								#legend.title=element_text(size=8, colour="grey30")
								)
						##### despite 4.35 width it saved it as 4.34
						ggsave(paste(out2,"constancy_per_nest_with_constancy_sex_labels_after_period_in_grey_smaller_na.png", sep=""),width=7.4, height=4.35, units='in',dpi=600)
						#ggsave(paste(out2,"constancy_per_nest_with_nest_labels_sex_after_period_3.png", sep=""))
				}
				
				{# not used				
				{# 7x4 with labels -with constancy lables and color sex
				nn$label=paste(nn$bb,nn$bt,"              ",sep="       ")
				# label constancy
				nnn=nn$label
				names(nnn)=nn$n
				
				# label nest
				nnnest=nn$nest
				names(nnnest)=nn$n
				
				dev.new(width=7.4, height=4.35) # for 4*7
				  ggplot(u,aes(x=datetime_c,y=roll_con, fill=sex))+
						#stat_smooth(method="lm", fill="grey95",alpha=1, aes(colour=genus_))+
						#scale_colour_manual(values=go$cols, name="Genus")+
						facet_wrap(~n, nrow=5, labeller=as_labeller(nnnest))+ #label is nests
						#facet_wrap(~n, nrow=5, labeller=as_labeller(nnn))+ #label is constancy
						geom_point(aes(color=sex), shape=20,size=0.25)+#colour="grey38",aes(fill=genus_))+
						scale_color_manual(values=c("#FCB42C","#535F7C"), guide=FALSE)+
						#scale_fill_manual(values=go$cols, name="Genus")+ 
						labs(x="Time in experimental period [h]", y="Nest attendance")+
						scale_x_continuous(limits=c(-16,33), breaks=seq(-15,30,by=15), labels=c(seq(-15,30,by=15)))+
						geom_vline(xintercept=0, lty=3,lwd=0.5, col='red')+
						geom_vline(data=e,aes(xintercept=r_time), lty=3,lwd=0.5, col='red')+
						geom_text(data = ann_text,aes(label =lab),size=1.9, angle = 90, colour='grey50',hjust='middle',vjust='middle')+
						#scale_x_continuous(breaks=15, labels=15)+
						scale_y_continuous(limits=c(-0.05,1.05),breaks=seq(0,1,by=0.25), labels=seq(0,1,by=0.25))+
						#coord_cartesian(ylim=c(0,1), xlim=c(-15,15))+ #coord_cartesian(ylim=c(-5,61))+ # needs to be 10 if it is not to touch x-axis
						theme_light()+
						theme(	axis.line=element_line(colour="grey70", size=0.25),
								#panel.border=element_rect(colour="white"),
								panel.border=element_rect(colour="grey70", size=0.25),
								panel.grid = element_blank(),
								
								axis.title=element_text(size=7, colour="grey50"),
								axis.title.y = element_text(vjust=1),
								axis.title.x = element_text(vjust=0.2),
								axis.text=element_text(size=6),# margin=units(0.5,"mm")),
								axis.ticks.length=unit(0.5,"mm"),
								#axis.ticks.margin,
								
								strip.text.x =element_text(size = 6, color="grey30",  margin=margin(1,1,1,1,"mm")), #grey50
								strip.background=element_rect(fill="grey99",colour="grey70", size=0.25),
								#strip.background = element_blank(), 
								#strip.text = element_blank(),
								panel.margin = unit(0, "mm"),
								legend.position="none"
								#legend.background=element_rect(colour="white"),
								#legend.key=element_rect(fill="grey99", colour="white"),
								#legend.text=element_text(size=7, colour="grey30"),
								#legend.title=element_text(size=8, colour="grey30")
								)
						
						ggsave(paste(out2,"constancy_per_nest_with_constancy_labels_sex_after_period_3.png", sep=""))
						ggsave(paste(out2,"constancy_per_nest_with_nest_labels_sex_after_period_3.png", sep=""))
				}
			
				{# 7x4 with nest labels
				dev.new(width=7.4, height=4.35) # for 4*7
				ggplot(u,aes(x=datetime_c,y=roll_con))+
						#stat_smooth(method="lm", fill="grey95",alpha=1, aes(colour=genus_))+
						#scale_colour_manual(values=go$cols, name="Genus")+
						facet_wrap(~nest, nrow=4)+
						geom_point(shape=20,size=0.25)+#colour="grey38",aes(fill=genus_))+
						#scale_fill_manual(values=go$cols, name="Genus")+ 
						labs(x="Time in experimental period [h]", y="Nest attendance")+
						scale_x_continuous(limits=c(-16,16), breaks=seq(-15,15,by=5), labels=c(seq(-15,15,by=5)))+
						geom_vline(xintercept=0, lty=3,lwd=0.5, col='red')+
						#scale_x_continuous(breaks=15, labels=15)+
						scale_y_continuous(limits=c(-0.05,1.05),breaks=seq(0,1,by=0.25), labels=seq(0,1,by=0.25))+
						#coord_cartesian(ylim=c(0,1), xlim=c(-15,15))+ #coord_cartesian(ylim=c(-5,61))+ # needs to be 10 if it is not to touch x-axis
						theme_light()+
						theme(	axis.line=element_line(colour="grey70", size=0.25),
								#panel.border=element_rect(colour="white"),
								panel.border=element_rect(colour="grey70", size=0.25),
								panel.grid = element_blank(),
								
								axis.title=element_text(size=7, colour="grey50"),
								axis.title.y = element_text(vjust=1),
								axis.title.x = element_text(vjust=0.2),
								axis.text=element_text(size=6),# margin=units(0.5,"mm")),
								axis.ticks.length=unit(0.5,"mm"),
								#axis.ticks.margin,
								
								strip.text.x =element_text(size = 6, color="grey30",  margin=margin(1,1,1,1,"mm")), #grey50
								strip.background=element_rect(fill="grey99",colour="grey70", size=0.25),
								#strip.background = element_blank(), 
								#strip.text = element_blank(),
								panel.margin = unit(0, "mm"),
								
								#legend.background=element_rect(colour="white"),
								legend.key=element_rect(fill="grey99", colour="white"),
								legend.text=element_text(size=7, colour="grey30"),
								legend.title=element_text(size=8, colour="grey30")
								)
						
						ggsave(paste(out2,"constancy_per_nest_with_nest_labels_tick.png", sep=""))
				}
				}
			}
		}
		}
		
	{# Explaining the diversity in compensation 
		{# run first - prepare data
			{# define colors
				col_t='#FCB42C'
				col_t_light='#fef0d4'
				col_c='#535F7C'
				col_c_light='#cbcfd7'
				col_2011=rgb(0,100,0,255, maxColorValue=255)
				col_2011_light=rgb(102,162,102,255, maxColorValue=255)
			}
			{# add proprotion of incubation for focal bird from three days with available data prior to treatment
			l=list()
			for(i in 1:nrow(bb)){
								x=bb[bb$nest==bb$nest[i],]
								day=as.Date(trunc(x$bout_start, "day"))
								pp=paste("extract_2013.mb_inc_2013_",bb$nest[i],sep="")	
								if('s402'==b$nest[i]){
								y=wet(conMy, paste("SELECT pk, datetime_, inc_t, incubation, sex_inc_filled FROM  ", pp," where exchange_gap is NULL and datetime_<'",day,"' and datetime_>='",day-4,"'", sep=""))
								y$datetime_= as.POSIXct(y$datetime_)
								y$day=as.Date(trunc(y$datetime_, "day"))
								y=y[y$day!=as.Date("2013-06-20"),]
								}else{
									y=wet(conMy, paste("SELECT pk, datetime_, inc_t, incubation, sex_inc_filled FROM  ", pp," where exchange_gap is NULL and datetime_<'",day,"' and datetime_>='",day-3,"'", sep=""))
									y$datetime_= as.POSIXct(y$datetime_)
									y$day=as.Date(trunc(y$datetime_, "day"))
									}
								l[[i]]=ddply(y,.(day), summarise, nest=x$nest, p_f=length(sex_inc_filled[sex_inc_filled=='f'])/length(sex_inc_filled))
									#length(y$sex_inc_filled[y$sex_inc_filled=='f' & y$day=='2013-06-19'])/length(y$sex_inc_filled[y$day=='2013-06-19'])
								
								print(bb$nest[i])
								}
			p_f=do.call(rbind,l)
			p_ff=ddply(p_f,.(nest),summarise, mean_=mean(p_f), med_=median(p_f))
			bt$p_f=p_ff$med_[match(bt$nest,p_ff$nest)]
			bt$prop=ifelse(bt$sex=='f', bt$p_f, 1-bt$p_f)
			}
			{# add escape distance - imputed for a bird from nest s628
				load(file=paste(wd,'escape.Rdata',sep="")) #vv - all estimates, es - aggregated per nest and sex
				bt$esc=es$esc[match(paste(bt$nest, bt$sex),paste(es$nest,es$sex))]
				bt$esc[bt$nest=='s628']=34.05681 # value based on imputation
			}
			{# add incubation start
				load("C:\\Users\\mbulla\\Documents\\Dropbox\\Science\\Projects\\MS\\Removal\\Analyses\\Data\\inc_start_estimation_2013_new.Rdata")
				 e$nest = tolower(e$nest)
				 e$inc_start = as.character(e$inc_start)
				 e$inc_start = as.POSIXct(e$inc_start, tz = 'UTC')
				 bt$inc_start = e$inc_start[match(bt$nest,e$nest)] 	
				 bt$day_inc=as.numeric(difftime(bt$bout_start,bt$inc_start, units = "days"))	
			}
			{# load 2011 data
				g=read.csv(paste(wd,'B-bout_length_constancy.csv', sep=""), stringsAsFactors=FALSE)
						g=g[-which(is.na(g$inc_eff)),]
						h=wet(conMy, "SELECT nest, latit as lat, longit as lon FROM avesatbarrow.nests where species='sesa' and year_=2011")
						g$lat=h$lat[match(g$nest,h$nest)]
						g$lon=h$lon[match(g$nest,h$nest)]		
						g$bout_start=as.POSIXct(g$bout_start)
						g$midbout=g$bout_start+60*g$bout_length/2
						g$time = as.numeric(difftime(g$midbout,trunc(g$midbout,"day"),   #DO NOT USE rfid_T$day here!! (the time will be wrong)
												units = "hours"))
						g$rad=g$time*pi/12
						g$sin_=sin(g$rad)
						g$cos_=cos(g$rad)
						
						g$bird_ID=paste(g$nest,g$sex,sep=" ")
				
		}
		}			
					{# script used to impute an escape distance for a bird from nest s628
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
					
		{# explore relationship between time of day and temperature
			#dev.new(width = 1.9685, height = 1.9685)
			ggplot(bt, aes(x = time, y = t_ambient_med)) + geom_point(size=1.15)+ stat_smooth(size=0.8, color='orange',fill='orange') + 
			xlab('Time of day [h]') + ylab('Median tundra temperature [°C]') + annotate(geom="text", x=24, y=27, label="a",fontface =2, size = 8*0.352777778)+ scale_x_continuous(limits = c(0,24), breaks = c(0,6,12,18,24)) + scale_y_continuous(limits = c(-5,27), breaks = c(-5,0,5,10,15,20,25))+ 
			theme_light()+
			theme(	axis.ticks.length=unit(0.5,"mm"),
					
					axis.line.y = element_line(color="grey40", size = 0.25),
					axis.text.y=element_text(size=6, color='grey20'),# margin=units(0.5,"mm")),
					axis.title.y = element_text(size=7, color='grey20', margin = margin( r = 6)),		
					
					#axis.line.y = element_line(color="grey70", size = 0.25),
					axis.text.x = element_text(size=6, colour = "white"),
					axis.title.x = element_text(size=7, colour = "white"),
					axis.ticks.x = element_blank(),
							
					panel.background=element_blank(),
					panel.border=element_blank(),
					panel.grid.major=element_blank(),
					panel.grid.minor=element_blank(),
					plot.background=element_blank()
					)
			ggsave(paste(out2,"Supplementary_Figure_Xa_new.png", sep=""),width = 5, height = 5, units = "cm",  dpi = 300)
			
			ggplot(bt, aes(x = time, y = t_amb_avg)) + geom_point(size=1.15)+ stat_smooth(size=0.8, color='orange',fill='orange') +
			xlab('Time of day [h]') + ylab('Mean tundra temperature [°C]') + annotate(geom="text", x=24, y=27, label="b",fontface =2,size = 8*0.352777778)+ scale_x_continuous(limits = c(0,24), breaks = c(0,6,12,18,24)) + scale_y_continuous(limits = c(-5,27), breaks = c(-5,0,5,10,15,20,25))+ 
			theme_light()+
			theme(	axis.ticks.length=unit(0.5,"mm"),
					
					axis.line.y = element_line(color="grey40", size = 0.25),
					axis.text.y=element_text(size=6, color='grey20'),# margin=units(0.5,"mm")),
					axis.title.y = element_text(size=7, color='grey20', margin = margin( r = 6)),		
					
					#axis.line.y = element_line(color="grey70", size = 0.25),
					axis.text.x = element_text(size=6, colour = "white"),
					axis.title.x = element_text(size=7, colour = "white"),
					axis.ticks.x = element_blank(),
							
					panel.background=element_blank(),
					panel.border=element_blank(),
					panel.grid.major=element_blank(),
					panel.grid.minor=element_blank(),
					plot.background=element_blank()
					)
			ggsave(paste(out2,"Supplementary_Figure_Xb_new.png", sep=""),width = 5, height = 5, units = "cm",  dpi = 300)
			
			ggplot(bt, aes(x = time, y = t_station_med)) + geom_point(size=1.15)+ stat_smooth(size=0.8, color='orange',fill='orange') +
			xlab('Time of day [h]') + ylab('Median ambient temperature [°C]') + annotate(geom="text", x=24, y=27, label="c",fontface =2,size = 8*0.352777778)+ scale_x_continuous(limits = c(0,24), breaks = c(0,6,12,18,24)) + scale_y_continuous(limits = c(-5,27), breaks = c(-5,0,5,10,15,20,25))+ 
			theme_light()+
			theme(	axis.ticks.length=unit(0.5,"mm"),
					
					axis.line.y = element_line(color="grey40", size = 0.25),
					axis.text.y=element_text(size=6, color='grey20'),# margin=units(0.5,"mm")),
					axis.title.y = element_text(size=7, color='grey20', margin = margin( r = 6)),		
					
					axis.line.x = element_line(color="grey40", size = 0.25),
					axis.text.x = element_text(size=6, colour = "grey20"),
					axis.title.x = element_text(size=7, colour = "grey20"),
					axis.ticks.x = element_blank(),
							
					panel.background=element_blank(),
					panel.border=element_blank(),
					panel.grid.major=element_blank(),
					panel.grid.minor=element_blank(),
					plot.background=element_blank()
					)
			ggsave(paste(out2,"Supplementary_Figure_Xc_new.png", sep=""),width = 5, height = 5, units = "cm",  dpi = 300)
			# not used
			ggplot(bt, aes(x = t_ambient_med, y = esc)) + geom_point()+ stat_smooth()
			ggplot(bt, aes(x = t_ambient_med, y = inc_eff)) + geom_point()+ stat_smooth(method='lm')
			ggplot(bt, aes(x = t_amb_avg, y = inc_eff)) + geom_point()+ stat_smooth(method='lm')
			ggplot(bt, aes(x = t_station_med, y = inc_eff)) + geom_point()+ stat_smooth(method='lm')
			m=lm(inc_eff~t_ambient_med+esc+prop, bt)
		}
		{# explore relationship between esc and day_inc and inc_start
			ggplot(bt, aes(x = inc_start, y = day_inc)) + geom_point()+ stat_smooth(method='lm')
			ggplot(bt, aes(x = inc_start, y = esc)) + geom_point()+ stat_smooth(method='lm')
			m=lm(inc_eff~day_inc+prop, bt)
			AICc(m,nobs = 25)
		}
		{# explore predictors with sex
			ggplot(bt, aes(x=sex, y = time)) + geom_boxplot() + 
				geom_dotplot(aes(fill = sex),binaxis = 'y', stackdir = 'center',  position = position_dodge())+
				ylab('Time of day [h]') + xlab('Sex')+annotate(geom="text", x=2.5, y=20, label="a",fontface =2)+ 
				scale_y_continuous(limits = c(0,20), breaks = c(0,5,10,15,20))+scale_x_discrete(labels=c("f" = "\u2640", "m" = "\u2642"))+ 	theme( 
				legend.position="none",
				#axis.text.x = element_text(colour = "white"), 
				#axis.title.x = element_text(colour = "white"),
				#axis.ticks.x = element_blank(),
				axis.line.x = element_line(color="black", size = 0.25),
				axis.title.y = element_text(margin = margin( r = 9)),
				axis.line.y = element_line(color="black", size = 0.25),
				panel.background=element_blank(),
				panel.border=element_blank(),
				panel.grid.major=element_blank(),
				panel.grid.minor=element_blank(),
				plot.background=element_blank()
				) 
			ggsave(paste(out2,"Supplementary_Figure_XXa.png", sep=""),width = 5, height = 6, units = "cm",  dpi = 300)
			
			ggplot(bt, aes(x=sex, y = t_ambient_med)) + geom_boxplot() + geom_dotplot(aes(fill = sex),binaxis = 'y', stackdir = 'center',  position = position_dodge())+ylab('Tundra temperature [°C]')+ xlab('Sex')+annotate(geom="text", x=2.5, y=30, label="b",fontface =2)+ scale_y_continuous(limits = c(0,30), breaks = c(0,5,10,15,20,25,30))+scale_x_discrete(labels=c("f" = "\u2640", "m" = "\u2642"))+ 	theme( 
				legend.position="none",
				#axis.text.x = element_text(colour = "white"), 
				#axis.title.x = element_text(colour = "white"),
				#axis.ticks.x = element_blank(),
				axis.line.x = element_line(color="black", size = 0.25),
				axis.title.y = element_text(margin = margin( r = 9)),
				axis.line.y = element_line(color="black", size = 0.25),
				panel.background=element_blank(),
				panel.border=element_blank(),
				panel.grid.major=element_blank(),
				panel.grid.minor=element_blank(),
				plot.background=element_blank()
				) 
			ggsave(paste(out2,"Supplementary_Figure_XXb.png", sep=""),width = 5, height = 6, units = "cm",  dpi = 300)
			
			ggplot(bt, aes(x=sex, y = prop)) + geom_boxplot() + geom_dotplot(aes(fill = sex),binaxis = 'y', stackdir = 'center',  position = position_dodge())+ylab('Proportion of incubation')+ xlab('Sex')+annotate(geom="text", x=2.5, y=0.56, label="c",fontface =2)+ scale_y_continuous(limits = c(0.40,0.56))+scale_x_discrete(labels=c("f" = "\u2640", "m" = "\u2642"))+ 	theme( 
				legend.position="none",
				#axis.text.x = element_text(colour = "white"), 
				#axis.title.x = element_text(colour = "white"),
				#axis.ticks.x = element_blank(),
				axis.line.x = element_line(color="black", size = 0.25),
				axis.title.y = element_text(margin = margin( r = 9)),
				axis.line.y = element_line(color="black", size = 0.25),
				panel.background=element_blank(),
				panel.border=element_blank(),
				panel.grid.major=element_blank(),
				panel.grid.minor=element_blank(),
				plot.background=element_blank()
				) 
			ggsave(paste(out2,"Supplementary_Figure_XXc.png", sep=""),width = 5, height = 6, units = "cm",  dpi = 300)
			ggplot(bt, aes(x=sex, y = esc)) + geom_boxplot() + geom_dotplot(aes(fill = sex),binaxis = 'y', stackdir = 'center',  position = position_dodge())+ylab('Escape distance [m]')+ xlab('Sex')+annotate(geom="text", x=2.5, y=65, label="d",fontface =2)+ scale_y_continuous(limits = c(0,65))+scale_x_discrete(labels=c("f" = "\u2640", "m" = "\u2642"))+ 	theme( 
				legend.position="none",
				#axis.text.x = element_text(colour = "white"), 
				#axis.title.x = element_text(colour = "white"),
				#axis.ticks.x = element_blank(),
				axis.line.x = element_line(color="black", size = 0.25),
				axis.title.y = element_text(margin = margin( r = 9)),
				axis.line.y = element_line(color="black", size = 0.25),
				panel.background=element_blank(),
				panel.border=element_blank(),
				panel.grid.major=element_blank(),
				panel.grid.minor=element_blank(),
				plot.background=element_blank()
				) 
			ggsave(paste(out2,"Supplementary_Figure_XXd.png", sep=""),width = 5, height = 6, units = "cm",  dpi = 300)
			
			
			ggplot(bt, aes(x=sex, y = log(esc))) + geom_boxplot() + geom_dotplot(aes(fill = sex),binaxis = 'y', stackdir = 'center',  position = position_dodge())
			
			
			
		}

		{# Supplementary Table 2
			{# 1 - Temperation, proportion, escape
				m=lm(inc_eff~scale(t_ambient_med) + scale(prop)+scale(esc), bt)
						pred=c('Intercept','T','Proportion','Escape')
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
					v <- apply(bsim@coef, 2,quantile, prob=c(0.5))
					ci=apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
					oi=data.frame(model='1',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
					rownames(oi) = NULL
						oi$estimate_r=round(oi$estimate,3)
						oi$lwr_r=round(oi$lwr,3)
						oi$upr_r=round(oi$upr,3)
						#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
					o1=oi[c('model','type',"effect", "estimate_r","lwr_r",'upr_r')]	
					
			}
			{# 2 - Temperature model - treated
					m=lm(inc_eff~scale(t_ambient_med), bt)
						pred=c('Intercept','T')
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
					v <- apply(bsim@coef, 2,quantile, prob=c(0.5))
					ci=apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
					oi=data.frame(model='2',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
					rownames(oi) = NULL
						oi$estimate_r=round(oi$estimate,3)
						oi$lwr_r=round(oi$lwr,3)
						oi$upr_r=round(oi$upr,3)
						#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
					o2=oi[c('model','type',"effect", "estimate_r","lwr_r",'upr_r')]	
					
			}
			{# 3 - Time model - treated
					m=lm(inc_eff~sin_+cos_, bt)
						pred=c('Intercept','Sin','Cos')
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
					v <- apply(bsim@coef, 2,quantile, prob=c(0.5))
					ci=apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
					oi=data.frame(model='3',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
					rownames(oi) = NULL
						oi$estimate_r=round(oi$estimate,3)
						oi$lwr_r=round(oi$lwr,3)
						oi$upr_r=round(oi$upr,3)
						#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
					o3=oi[c('model','type',"effect", "estimate_r","lwr_r",'upr_r')]	
				}
			{# 4 - Proportion 
					m=lm(inc_eff~scale(prop), bt)
						pred=c('Intercept','Proportion')
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
					v <- apply(bsim@coef, 2,quantile, prob=c(0.5))
					ci=apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
					oi=data.frame(model='4',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
					rownames(oi) = NULL
						oi$estimate_r=round(oi$estimate,3)
						oi$lwr_r=round(oi$lwr,3)
						oi$upr_r=round(oi$upr,3)
						#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
					o4=oi[c('model','type',"effect", "estimate_r","lwr_r",'upr_r')]	
					
				}
			{# 5 - Escape distance
					m=lm(inc_eff~scale(esc), bt)
						pred=c('Intercept','escape distance')
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
					v <- apply(bsim@coef, 2,quantile, prob=c(0.5))
					ci=apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
					oi=data.frame(model='5',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
					rownames(oi) = NULL
						oi$estimate_r=round(oi$estimate,3)
						oi$lwr_r=round(oi$lwr,3)
						oi$upr_r=round(oi$upr,3)
						#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
					o5=oi[c('model','type',"effect", "estimate_r","lwr_r",'upr_r')]	
					
					}				
			{# AICc comparison - find out which model fits better
			  mf=lm(inc_eff~prop+esc+t_ambient_med, bt)
			  m0=lm(inc_eff~t_amb_avg, bt)
			  m1=lm(inc_eff~sin_+cos_, bt)
			  m2=lm(inc_eff~prop, bt)
			  m3=lm(inc_eff~esc, bt)
				o=data.frame(model=c('mf','m0','m1','m2','m3'), AIC=c(AICc(mf, nobs=25),AICc(m0, nobs=25),AICc(m1, nobs=25),AICc(m2,nobs=25),AICc(m3,nobs=25)))
						o$delta=o$AIC-min(o$AIC)
						o$prob=exp(-0.5*o$delta)/sum(exp(-0.5*o$delta))
						o$ER=max(o$prob)/o$prob
						#o[order(o$delta),]
						o$AIC=round(o$AIC,2)
						o$delta=round(o$delta,2)
						o$prob=round(o$prob,3)
						o$ER=round(o$ER,2)
						
		}
		
			{# combine and export to excel table
						sname = tempfile(fileext='.xls')
						wb = loadWorkbook(sname,create = TRUE)	
						createSheet(wb, name = "output")
						writeWorksheet(wb, rbind(o1,o2,o3,o4,o5), sheet = "output")
						createSheet(wb, name = "AIC")
						writeWorksheet(wb, o, sheet = "AIC")
						saveWorkbook(wb)
						shell(sname)
			}
			
			
			{# model assumptions temperature - OK
						dev.new(width=6,height=9)
						par(mfrow=c(4,3))
						m=lm(inc_eff~scale(t_ambient_med) + scale(prop)+scale(esc), bt)
									
						plot(m)
						
						scatter.smooth(resid(m)~bt$t_ambient_med);abline(h=0, lty=2)
						scatter.smooth(resid(m)~bt$prop);abline(h=0, lty=2)
						scatter.smooth(resid(m)~bt$esc);abline(h=0, lty=2)
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=bt$lon, y=bt$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
				
				}
			{# model assumptions temperature - OK
						dev.new(width=6,height=9)
						par(mfrow=c(4,3))
						m=lm(inc_eff~t_ambient_med, bt)
									
						plot(m)
						
						scatter.smooth(resid(m)~bt$sin_);abline(h=0, lty=2)
						scatter.smooth(resid(m)~bt$cos_);abline(h=0, lty=2)
						
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=bt$lon, y=bt$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
				
				}
			{# model assumptions time - OK
						dev.new(width=6,height=9)
						par(mfrow=c(4,3))
						m=lm(inc_eff~sin_+cos_, bt)
									
						plot(m)
						
						scatter.smooth(resid(m)~bt$sin_);abline(h=0, lty=2)
						scatter.smooth(resid(m)~bt$cos_);abline(h=0, lty=2)
						
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=bt$lon, y=bt$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
				
				}
			{# model assumptions proportion - OK
						
						dev.new(width=6,height=9)
						par(mfrow=c(4,3))
						
						m=lm(inc_eff~prop, bt)
						
						plot(m)
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=bt$lon, y=bt$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
				
				}
			{# model assumptions escape - OK
						dev.new(width=6,height=9)
						par(mfrow=c(4,3))
						m=lm(inc_eff~esc, bt)
														
						plot(m)
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=bt$lon, y=bt$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
				
				}
				
		}
		{# Supplementary Table 2 - with sex not used
			{# 1 - Temperation, proportion, escape
				m=lm(inc_eff~prop+esc+t_ambient_med, bt)
						pred=c('Intercept','Proportion','Escape','T')
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
					v <- apply(bsim@coef, 2,quantile, prob=c(0.5))
					ci=apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
					oi=data.frame(model='1',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
					rownames(oi) = NULL
						oi$estimate_r=round(oi$estimate,3)
						oi$lwr_r=round(oi$lwr,3)
						oi$upr_r=round(oi$upr,3)
						#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
					o1=oi[c('model','type',"effect", "estimate_r","lwr_r",'upr_r')]	
					
					m=lm(inc_eff~sin_*sex+cos_*sex, bt)
						pred=c('Intercept','Sin','Cos','Sex', 'sin:sex','cos:sex')
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
					v <- apply(bsim@coef, 2,quantile, prob=c(0.5))
					ci=apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
					oi=data.frame(model='1',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
					rownames(oi) = NULL
						oi$estimate_r=round(oi$estimate,3)
						oi$lwr_r=round(oi$lwr,3)
						oi$upr_r=round(oi$upr,3)
						#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
					o1s=oi[c('model','type',"effect", "estimate_r","lwr_r",'upr_r')]	
					o1=rbind(o1,o1s)
			}
			{# 2 - Temperature model - treated
					m=lm(inc_eff~scale(t_ambient_med), bt)
						pred=c('Intercept','Sin','Cos')
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
					v <- apply(bsim@coef, 2,quantile, prob=c(0.5))
					ci=apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
					oi=data.frame(model='1',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
					rownames(oi) = NULL
						oi$estimate_r=round(oi$estimate,3)
						oi$lwr_r=round(oi$lwr,3)
						oi$upr_r=round(oi$upr,3)
						#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
					o2=oi[c('model','type',"effect", "estimate_r","lwr_r",'upr_r')]	
					
					m=lm(inc_eff~sin_*sex+cos_*sex, bt)
						pred=c('Intercept','Sin','Cos','Sex', 'sin:sex','cos:sex')
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
					v <- apply(bsim@coef, 2,quantile, prob=c(0.5))
					ci=apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
					oi=data.frame(model='1',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
					rownames(oi) = NULL
						oi$estimate_r=round(oi$estimate,3)
						oi$lwr_r=round(oi$lwr,3)
						oi$upr_r=round(oi$upr,3)
						#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
					o2s=oi[c('model','type',"effect", "estimate_r","lwr_r",'upr_r')]	
					o2=rbind(o2,o2s)
			}
			{# 3 - Time model - treated
					m=lm(inc_eff~sin_+cos_, bt)
						pred=c('Intercept','Sin','Cos')
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
					v <- apply(bsim@coef, 2,quantile, prob=c(0.5))
					ci=apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
					oi=data.frame(model='1',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
					rownames(oi) = NULL
						oi$estimate_r=round(oi$estimate,3)
						oi$lwr_r=round(oi$lwr,3)
						oi$upr_r=round(oi$upr,3)
						#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
					o3=oi[c('model','type',"effect", "estimate_r","lwr_r",'upr_r')]	
					
					m=lm(inc_eff~sin_*sex+cos_*sex, bt)
						pred=c('Intercept','Sin','Cos','Sex', 'sin:sex','cos:sex')
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
					v <- apply(bsim@coef, 2,quantile, prob=c(0.5))
					ci=apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
					oi=data.frame(model='1',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
					rownames(oi) = NULL
						oi$estimate_r=round(oi$estimate,3)
						oi$lwr_r=round(oi$lwr,3)
						oi$upr_r=round(oi$upr,3)
						#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
					o3s=oi[c('model','type',"effect", "estimate_r","lwr_r",'upr_r')]	
					o3=rbind(o3,o3s)
			}
			{# 4 - Proportion 
					m=lm(inc_eff~prop, bt)
						pred=c('Intercept','Proportion')
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
					v <- apply(bsim@coef, 2,quantile, prob=c(0.5))
					ci=apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
					oi=data.frame(model='2',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
					rownames(oi) = NULL
						oi$estimate_r=round(oi$estimate,3)
						oi$lwr_r=round(oi$lwr,3)
						oi$upr_r=round(oi$upr,3)
						#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
					o4=oi[c('model','type',"effect", "estimate_r","lwr_r",'upr_r')]	
					
					m=lm(inc_eff~prop*sex, bt)
						pred=c('Intercept','Proportion','Sex', 'Proportion:sex')
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
					v <- apply(bsim@coef, 2,quantile, prob=c(0.5))
					ci=apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
					oi=data.frame(model='1',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
					rownames(oi) = NULL
						oi$estimate_r=round(oi$estimate,3)
						oi$lwr_r=round(oi$lwr,3)
						oi$upr_r=round(oi$upr,3)
						#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
					o4s=oi[c('model','type',"effect", "estimate_r","lwr_r",'upr_r')]	
					o4=rbind(o4,o4s)
						
			}	
			{# 5 - Escape distance
					m=lm(inc_eff~esc, bt)
						pred=c('Intercept','escape distance')
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
					v <- apply(bsim@coef, 2,quantile, prob=c(0.5))
					ci=apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
					oi=data.frame(model='3',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
					rownames(oi) = NULL
						oi$estimate_r=round(oi$estimate,3)
						oi$lwr_r=round(oi$lwr,3)
						oi$upr_r=round(oi$upr,3)
						#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
					o5=oi[c('model','type',"effect", "estimate_r","lwr_r",'upr_r')]	
					
					m=lm(inc_eff~esc*sex, bt)
						pred=c('Intercept','escape distance','Sex', 'escape distance:sex')
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
					v <- apply(bsim@coef, 2,quantile, prob=c(0.5))
					ci=apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
					oi=data.frame(model='1',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
					rownames(oi) = NULL
						oi$estimate_r=round(oi$estimate,3)
						oi$lwr_r=round(oi$lwr,3)
						oi$upr_r=round(oi$upr,3)
						#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
					o5s=oi[c('model','type',"effect", "estimate_r","lwr_r",'upr_r')]	
					o5=rbind(o5,o5s)	
				}				
			{# AICc comparison - find out which model fits better
			  mf=lm(inc_eff~prop+esc+t_ambient_med, bt)
			  m0=lm(inc_eff~t_amb_avg, bt)
			  m0s=lm(inc_eff~t_ambient_med*sex, bt)
			  m1=lm(inc_eff~sin_+cos_, bt)
			  m1s=lm(inc_eff~sin_*sex+cos_*sex, bt)
			  m2=lm(inc_eff~prop, bt)
			  m2s=lm(inc_eff~prop*sex, bt)
			  m3=lm(inc_eff~esc, bt)
			  m3s=lm(inc_eff~esc*sex, bt)
				o=data.frame(model=c('mf','m0','m0s','m1','m1s','m2','m2s','m3','m3s'), AIC=c(AICc(mf, nobs=25),AICc(m0, nobs=25),AICc(m0s, nobs=25),AICc(m1, nobs=25),AICc(m1s, nobs=25),AICc(m2,nobs=25),AICc(m2s,nobs=25),AICc(m3,nobs=25),AICc(m3s,nobs=25)))
						o$delta=o$AIC-min(o$AIC)
						o$prob=exp(-0.5*o$delta)/sum(exp(-0.5*o$delta))
						o$ER=max(o$prob)/o$prob
						#o[order(o$delta),]
						o$AIC=round(o$AIC,2)
						o$delta=round(o$delta,2)
						o$prob=round(o$prob,3)
						o$ER=round(o$ER,2)
						
		}
		
			{# combine and export to excel table
						sname = tempfile(fileext='.xls')
						wb = loadWorkbook(sname,create = TRUE)	
						createSheet(wb, name = "output")
						writeWorksheet(wb, rbind(o1,o2,o3), sheet = "output")
						createSheet(wb, name = "AIC")
						writeWorksheet(wb, o, sheet = "AIC")
						saveWorkbook(wb)
						shell(sname)
			}
			
			{# model assumptions temperature
						dev.new(width=6,height=9)
						par(mfrow=c(4,3))
						m=lm(inc_eff~t_ambient_med, bt)
									
						plot(m)
						
						scatter.smooth(resid(m)~bt$sin_);abline(h=0, lty=2)
						scatter.smooth(resid(m)~bt$cos_);abline(h=0, lty=2)
						
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=bt$lon, y=bt$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
				
				}
			
			{# model assumptions time - OK
						dev.new(width=6,height=9)
						par(mfrow=c(4,3))
						m=lm(inc_eff~sin_+cos_, bt)
									
						plot(m)
						
						scatter.smooth(resid(m)~bt$sin_);abline(h=0, lty=2)
						scatter.smooth(resid(m)~bt$cos_);abline(h=0, lty=2)
						
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=bt$lon, y=bt$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
				
				}
			{# model assumptions proportion - OK
						
						dev.new(width=6,height=9)
						par(mfrow=c(4,3))
						
						m=lm(inc_eff~prop, bt)
						
						plot(m)
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=bt$lon, y=bt$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
				
				}
			{# model assumptions escape - OK
						dev.new(width=6,height=9)
						par(mfrow=c(4,3))
						m=lm(inc_eff~esc, bt)
														
						plot(m)
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=bt$lon, y=bt$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
				
				}
				
		}
		{# Supplementary Table 3
				{# 4 - Time model - control
					m=lm(inc_eff~sin_+cos_, b[b$exper=='c',])
						pred=c('Intercept','Sin','Cos')
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
					v <- apply(bsim@coef, 2,quantile, prob=c(0.5))
					ci=apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
					oi=data.frame(model='time_c',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
					rownames(oi) = NULL
						oi$estimate_r=round(oi$estimate,3)
						oi$lwr_r=round(oi$lwr,3)
						oi$upr_r=round(oi$upr,3)
						#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
					o4=oi[c('model','type',"effect", "estimate_r","lwr_r",'upr_r')]	
			}
				{# 5 - Time model - 2011
					m=lmer(inc_eff~sin_+cos_+(sin_+cos_|bird_ID), g)
					# simulation		
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
						apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
						ci=apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
						oi=data.frame(model='time_2011',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
						rownames(oi) = NULL
						oi$estimate_r=round(oi$estimate,3)
						oi$lwr_r=round(oi$lwr,3)
						oi$upr_r=round(oi$upr,3)
						#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
						o5=oi[c('model','type',"effect", "estimate_r","lwr_r",'upr_r')]	
					
						# Random effects
						l=data.frame(summary(m)$varcor)
						l=l[is.na(l$var2),]
						l$var2=NULL
						ri=data.frame(model='time_2011', type='random (var)',effect=l$var, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA)
					o5=rbind(o5,ri)	
					
			}
				{# combine and export to excel table
						sname = tempfile(fileext='.xls')
						wb = loadWorkbook(sname,create = TRUE)	
						createSheet(wb, name = "output")
						writeWorksheet(wb, rbind(o4,o5), sheet = "output")
						saveWorkbook(wb)
						shell(sname)
			}
				
				{# model assumptions control
						dev.new(width=6,height=9)
						par(mfrow=c(4,3))
						bc=b[b$exper=='c',]
						m=lm(inc_eff~sin_+cos_, b[b$exper=='c',])
									
						plot(m)
						
						scatter.smooth(resid(m)~bc$sin_);abline(h=0, lty=2)
						scatter.smooth(resid(m)~bc$cos_);abline(h=0, lty=2)
						
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=bc$lon, y=bc$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
				
				}
				{# model assumptions 2011 - auto-cor
						dev.new(width=6,height=9)
									
						m=lmer(inc_eff~sin_+cos_+(sin_+cos_|nest), g)
						g=g[-which(is.na(g$bout_lag)),]	
						m=lmer(inc_eff~bout_lag+sin_+cos_+(sin_+cos_|nest), g)
						par(mfrow=c(4,3))
						
						scatter.smooth(fitted(m),resid(m),col='red');abline(h=0, lty=2)
						scatter.smooth(fitted(m),sqrt(abs(resid(m))), col='red')
													
						qqnorm(resid(m), main=list("Normal Q-Q Plot: residuals", cex=0.8),col='red') 
						qqline(resid(m))
			 
						qqnorm(unlist(ranef(m)$nest[1]), main = " Bird")
						qqline(unlist(ranef(m)$nest[1]))
						
						qqnorm(unlist(ranef(m)$nest[2]), main = " Slope(sin)")
						qqline(unlist(ranef(m)$nest[2]))
						
						qqnorm(unlist(ranef(m)$nest[3]), main = " Slope(cos)")
						qqline(unlist(ranef(m)$nest[3]))
						
						scatter.smooth(resid(m)~g$sin_);abline(h=0, lty=2, col='red')
						scatter.smooth(resid(m)~g$cos_);abline(h=0, lty=2, col='red')
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=g$lon, y=g$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
				
				}
		
		}
		
			xx=data.frame(	nest = c("s402","s409","s510","s516","s520","s711","s806"), 
							uni = c(0,0,1,0,1,0,0), stringsAsFactors=FALSE)
							xx$esc=bt$esc[match(xx$nest, bt$nest)]
			ggplot(xx,aes(x=factor(uni), y=esc))+geom_boxplot()+geom_point()
			m=glm(uni~esc, xx, family='binomial')
				summary(m)
				plot(allEffects(m))
						
		
		{# Figure 4
			{# run first - prepare predictions
				{# effect of time
					{# treated
						#m=lm(inc_eff~sin(rad)+cos(rad), b[b$exper=='t',])
						m=lm(inc_eff~sin_+cos_, b[b$exper=='t',])
						
						# simulation		
							nsim <- 2000
							bsim <- sim(m, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
						# coefficients
							v <- apply(bsim@coef, 2, quantile, prob=c(0.5))
						# predicted values		
							newD=data.frame(time_=seq(0,24,0.25))
									newD$rad=2*pi*newD$time_ / 24
									newD$sin_=sin(newD$rad)
									newD$cos_=cos(newD$rad)
						# exactly the model which was used has to be specified here 
							X <- model.matrix(~ sin_ + cos_,data=newD)	
										
						# calculate predicted values and creditability intervals
							newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
									predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
									for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@coef[i,]
									newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
									newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
									#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
									#newD=newD[order(newD$t_tundra),]
								ptt=newD
					}			
					{# control
						#m=lm(inc_eff~sin(rad)+cos(rad), b[b$exper=='t',])
						m=lm(inc_eff~sin_+cos_, b[b$exper=='c',])
						summary(m)
						plot(allEffects(m))
					
						# simulation		
							nsim <- 2000
							bsim <- sim(m, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
						# coefficients
							v <- apply(bsim@coef, 2, quantile, prob=c(0.5))
						# predicted values		
							newD=data.frame(time_=seq(0,24,0.25))
									newD$rad=2*pi*newD$time_ / 24
									newD$sin_=sin(newD$rad)
									newD$cos_=cos(newD$rad)
						# exactly the model which was used has to be specified here 
							X <- model.matrix(~ sin_ + cos_,data=newD)	
										
						# calculate predicted values and creditability intervals
							newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
									predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
									for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@coef[i,]
									newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
									newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
									#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
									#newD=newD[order(newD$t_tundra),]
							ptc=newD
					}			
					{# constancy from 2011
						
						m=lmer(inc_eff~sin_+cos_+(sin_+cos_|bird_ID), g)
							#summary(m)
							#plot(allEffects(m))
						
							# simulation		
								nsim <- 2000
								bsim <- sim(m, n.sim=nsim)  
								apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							
							# coefficients
								v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))	
							# predicted values		
								newD=data.frame(time_=seq(0,24,0.25))
										newD$rad=2*pi*newD$time_ / 24
										newD$sin_=sin(newD$rad)
										newD$cos_=cos(newD$rad)
								
							# exactly the model which was used has to be specified here 
								X <- model.matrix(~ sin_ + cos_,data=newD)	
											
							# calculate predicted values and creditability intervals
								newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
										predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
										for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@fixef[i,]
										newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
										newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
										#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
										#newD=newD[order(newD$t_tundra),]
								pt11=newD
		}
				}	
				{# temperature
					m=lm(inc_eff~prop+esc+t_ambient_med, bt)
						  
						# simulation		
							nsim <- 2000
							bsim <- sim(m, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
						# coefficients
							v <-apply(bsim@coef, 2, quantile, prob=c(0.5))
						# predicted values		
							newD=data.frame(prop=mean(bt$prop),
											esc = mean(bt$esc),
											t_ambient_med = seq(min(bt$t_ambient_med), max(bt$t_ambient_med), length.out=200)
											)
								
						# exactly the model which was used has to be specified here 
							X <- model.matrix(~ prop+esc+t_ambient_med,data=newD)	
										
						# calculate predicted values and creditability intervals
							newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
									predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
									for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@coef[i,]
									newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
									newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
									#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
									#newD=newD[order(newD$t_tundra),]
								pte=newD
				}
				{# proportion of incubation priot to removal
					m=lm(inc_eff~prop+esc+t_ambient_med, bt)
						  
						# simulation		
							nsim <- 2000
							bsim <- sim(m, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
						# coefficients
							v <-apply(bsim@coef, 2, quantile, prob=c(0.5))
						# predicted values		
							newD=data.frame(t_ambient_med=mean(bt$t_ambient_med),
											esc = mean(bt$esc),
											prop = seq(min(bt$prop), max(bt$prop), length.out=200)
											)
								
						# exactly the model which was used has to be specified here 
							X <- model.matrix(~ prop+esc+t_ambient_med,data=newD)	
										
						# calculate predicted values and creditability intervals
							newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
									predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
									for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@coef[i,]
									newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
									newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
									#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
									#newD=newD[order(newD$t_tundra),]
								pp_=newD
				}
				{# escape distance
					m=lm(inc_eff~prop+esc+t_ambient_med, bt)
						  
						# simulation		
							nsim <- 2000
							bsim <- sim(m, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
						# coefficients
							v <-apply(bsim@coef, 2, quantile, prob=c(0.5))
						# predicted values		
							newD=data.frame(t_ambient_med=mean(bt$t_ambient_med),
											prop = mean(bt$prop),
											esc = seq(min(bt$esc), max(bt$esc), length.out=200)
											)
								
						# exactly the model which was used has to be specified here 
							X <- model.matrix(~ prop+esc+t_ambient_med,data=newD)	
										
						# calculate predicted values and creditability intervals
							newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
									predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
									for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@coef[i,]
									newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
									newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
									#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
									#newD=newD[order(newD$t_tundra),]
								pe=newD
				}
			}
			{# plot
				{# with polygons
				  #dev.new(width=3.5,height=1.85)
				  #dev.new(width=4.5,height=1.85)
				  png(paste(out2,"Figure_4_polygons_a-d.png", sep=""), width=4.5,height=1.85,units="in",res=600)
				  par(mfrow=c(1,4),mar=c(0.0,0,0,0.4),oma = c(1.9, 1.8, 0.2, 0.5),ps=12, mgp=c(1.2,0.35,0), las=1, cex=1, col.axis="grey30",font.main = 1, col.lab="grey30", col.main="grey30", fg="grey70", bty='n') # 0.6 makes font 7pt, 0.7 8pt
				{# time
				#par(mar=c(2.2,2.1,0.5,0.1), ps=12,	cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
				
				par(ps=12,	cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="n",xpd=TRUE)
				plot(ptt$pred~ptt$time_, pch=19,xlim=c(0,24), ylim=c(0,1), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
									
					axis(1, at=seq(0,24,by=6),labels=c(0,'',12,'',24),cex.axis=0.5,mgp=c(0,-0.15,0))
						mtext('Time of day\n[h]',side=1,line=1, cex=0.6, las=1, col='grey30')
						
					axis(2, at=seq(0,1,by=0.25), labels=c('0.0','','0.5','','1.0'))
						mtext('Nest attendance',side=2,line=1, cex=0.6, las=3, col='grey30')
					
					mtext(expression(bold("a")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
					
					# treated
							polygon(c(ptt$time_, rev(ptt$time_)), c(ptt$lwr, 
								rev(ptt$upr)), border=NA, col=adjustcolor(col_t ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(ptt$time_, ptt$pred, col=col_t,lwd=1)
							
					
					# control
							polygon(c(ptc$time_, rev(ptc$time_)), c(ptc$lwr, 
								rev(ptc$upr)), border=NA, col=adjustcolor(col_c ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(ptc$time_, ptc$pred, col=col_c,lwd=1)
					# treated control points		
							points(bt$time, bt$inc_eff, col=col_t,bg=adjustcolor(col_t ,alpha.f = 0.4), pch=21,cex=0.5)
							points(bb$time, bb$inc_eff, col=col_c,bg=adjustcolor(col_c ,alpha.f = 0.4), pch=21,cex=0.5)
					
					# 2011
						points(g$time[g$inc_eff<0.8], g$inc_eff[g$inc_eff<0.8], col=adjustcolor(col_2011 ,alpha.f = 0.5), pch=20,cex=0.1)		
					
					text(x=0,y=0.85, labels='control', col=col_c, cex=0.5, adj=0)
					text(x=0,y=0.75, labels='treated', col=col_t, cex=0.5, adj=0)		
					text(x=0,y=0.65, labels='2011', col=col_2011, cex=0.4, adj=0)	
					
					
				}
				{# temperature
					par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="n",xpd=TRUE)
					plot(pte$pred~pte$t_ambient_med, pch=19,xlim=c(0,28), ylim=c(0,1), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
									
					axis(1, at=seq(0,28,by=7),labels=seq(0,28,by=7),cex.axis=0.5,mgp=c(0,-0.15,0))
						mtext('Tundra temperature\n[°C]',side=1,line=1, cex=0.6, las=1, col='grey30')
						#mtext('Tundra temperature [°C]',side=1,line=0.4, cex=0.6, las=1, col='grey30')
						
					#axis(2, at=seq(0,1,by=0.25), labels=FALSE)
						#mtext('Nest attendance',side=2,line=1, cex=0.6, las=3, col='grey30')
					
					mtext(expression(bold("b")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
					
					# treated
							polygon(c(pte$t_ambient_med, rev(pte$t_ambient_med)), c(pte$lwr, 
								rev(pte$upr)), border=NA, col=adjustcolor(col_t ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pte$t_ambient_med, pte$pred, col=col_t,lwd=1)
							points(bt$t_ambient_med, bt$inc_eff, col=col_t,bg=adjustcolor(col_t ,alpha.f = 0.4), pch=21,cex=0.5)
				}
				{# proportion
					par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="n",xpd=TRUE)
					plot(pp_$pred~pp_$prop, pch=19,xlim=c(40,55), ylim=c(0,1), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
									
					axis(1, at=seq(40,55,by=5),labels=seq(40,55,by=5),cex.axis=0.5,mgp=c(0,-0.15,0))
						mtext('Incubation share\n[%]',side=1,line=1, cex=0.6, las=1, col='grey30')
						
					#axis(2, at=seq(0,1,by=0.25), labels=FALSE)
						#mtext('Nest attendance',side=2,line=1, cex=0.6, las=3, col='grey30')
					
					mtext(expression(bold("c")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
					
					# treated
							polygon(c(pp_$prop*100, rev(pp_$prop*100)), c(pp_$lwr, 
								rev(pp_$upr)), border=NA, col=adjustcolor(col_t ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pp_$prop*100, pp_$pred, col=col_t,lwd=1)
							points(bt$prop*100, bt$inc_eff, col=col_t,bg=adjustcolor(col_t ,alpha.f = 0.4), pch=21,cex=0.5)
				}
				{# escape
					par(ps=12,	cex =1, cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="n",xpd=TRUE)
					plot(pe$pred~pe$esc, pch=19,xlim=c(0,65), ylim=c(0,1), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
									
					axis(1, at=seq(0,65,by=15),labels=c(0,'',30,'',60),cex.axis=0.5,mgp=c(0,-0.15,0))
						mtext('Escape distance\n[m]',side=1,line=1, cex=0.6, las=1, col='grey30')
						
					#axis(2, at=seq(0,1,by=0.25), labels=FALSE)
						#mtext('Nest attendance',side=2,line=1, cex=0.6, las=3, col='grey30')
					
					mtext(expression(bold("d")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
					
					# treated
							polygon(c(pe$esc, rev(pe$esc)), c(pe$lwr, 
								rev(pe$upr)), border=NA, col=adjustcolor(col_t ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pe$esc, pe$pred, col=col_t,lwd=1)
							points(bt$esc, bt$inc_eff, col=col_t,bg=adjustcolor(col_t ,alpha.f = 0.4), pch=21,cex=0.5)
				}
				   dev.off()
				}										
				{# 2011 add small panel on top
					#dev.new(width=3.5,height=1.85/5)
					png(paste(out2,"Figure_4_2011_abcd.png", sep=""), width=4.5,height=1.85/5,units="in",res=600)
						  par(mfrow=c(1,4),mar=c(0.0,0,0,0.4),oma = c(0.01, 1.8, 0.2, 0.5),ps=12, mgp=c(1.2,0.35,0), las=1, cex=1, col.axis="grey30",font.main = 1, col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt
					
						par(ps=12,	cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="n",xpd=TRUE)
						plot(pt11$pred~pt11$time_, pch=19,xlim=c(0,24), ylim=c(0.75,1), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
											
							axis(2, at=seq(0.8,1,by=0.2), labels=c('0.8','1.0'))
							
							
							points(g$time[g$inc_eff>=0.8], g$inc_eff[g$inc_eff>=0.8], col=adjustcolor(col_2011 ,alpha.f = 0.5), pch=20,cex=0.1)		
							polygon(c(pt11$time_, rev(pt11$time_)), c(pt11$lwr, 
								rev(pt11$upr)), border=NA, col=adjustcolor(col_c ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pt11$time_, pt11$pred, col=col_2011,lwd=1)
								
							#mtext(expression(bold("a")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
					dev.off()
			}		
		
		
		{# run first - prepare predictions
				{# effect of time
					{# treated
						#m=lm(inc_eff~sin(rad)+cos(rad), b[b$exper=='t',])
						m=lm(inc_eff~sin_+cos_, b[b$exper=='t',])
						
						# simulation		
							nsim <- 2000
							bsim <- sim(m, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
						# coefficients
							v <- apply(bsim@coef, 2, quantile, prob=c(0.5))
						# predicted values		
							newD=data.frame(time_=seq(0,24,0.25))
									newD$rad=2*pi*newD$time_ / 24
									newD$sin_=sin(newD$rad)
									newD$cos_=cos(newD$rad)
						# exactly the model which was used has to be specified here 
							X <- model.matrix(~ sin_ + cos_,data=newD)	
										
						# calculate predicted values and creditability intervals
							newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
									predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
									for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@coef[i,]
									newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
									newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
									#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
									#newD=newD[order(newD$t_tundra),]
								ptt=newD
					}			
					{# control
						#m=lm(inc_eff~sin(rad)+cos(rad), b[b$exper=='t',])
						m=lm(inc_eff~sin_+cos_, b[b$exper=='c',])
						summary(m)
						plot(allEffects(m))
					
						# simulation		
							nsim <- 2000
							bsim <- sim(m, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
						# coefficients
							v <- apply(bsim@coef, 2, quantile, prob=c(0.5))
						# predicted values		
							newD=data.frame(time_=seq(0,24,0.25))
									newD$rad=2*pi*newD$time_ / 24
									newD$sin_=sin(newD$rad)
									newD$cos_=cos(newD$rad)
						# exactly the model which was used has to be specified here 
							X <- model.matrix(~ sin_ + cos_,data=newD)	
										
						# calculate predicted values and creditability intervals
							newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
									predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
									for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@coef[i,]
									newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
									newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
									#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
									#newD=newD[order(newD$t_tundra),]
							ptc=newD
					}			
					{# constancy from 2011
						
						m=lmer(inc_eff~sin_+cos_+(sin_+cos_|nest), g)
							#summary(m)
							#plot(allEffects(m))
						
							# simulation		
								nsim <- 2000
								bsim <- sim(m, n.sim=nsim)  
								apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							
							# coefficients
								v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))	
							# predicted values		
								newD=data.frame(time_=seq(0,24,0.25))
										newD$rad=2*pi*newD$time_ / 24
										newD$sin_=sin(newD$rad)
										newD$cos_=cos(newD$rad)
								
							# exactly the model which was used has to be specified here 
								X <- model.matrix(~ sin_ + cos_,data=newD)	
											
							# calculate predicted values and creditability intervals
								newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
										predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
										for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@fixef[i,]
										newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
										newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
										#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
										#newD=newD[order(newD$t_tundra),]
								pt11=newD
		}
				}	
				{# proportion of incubation priot to removal
						m=lm(inc_eff~prop,bt)
						  
						# simulation		
							nsim <- 2000
							bsim <- sim(m, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
						# coefficients
							v <-apply(bsim@coef, 2, quantile, prob=c(0.5))
						# predicted values		
							newD=data.frame(prop=seq(0.40,0.55,0.001))
								
						# exactly the model which was used has to be specified here 
							X <- model.matrix(~ prop,data=newD)	
										
						# calculate predicted values and creditability intervals
							newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
									predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
									for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@coef[i,]
									newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
									newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
									#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
									#newD=newD[order(newD$t_tundra),]
								pp_=newD
				}
				{# escape distance
						m=lm(inc_eff~esc,bt)
						  
						# simulation		
							nsim <- 2000
							bsim <- sim(m, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
						# coefficients
							v <-apply(bsim@coef, 2, quantile, prob=c(0.5))
						# predicted values		
							newD=data.frame(esc=seq(0,65,0.2))
									
						# exactly the model which was used has to be specified here 
							X <- model.matrix(~ esc,data=newD)	
										
						# calculate predicted values and creditability intervals
							newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
									predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
									for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@coef[i,]
									newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
									newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
									#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
									#newD=newD[order(newD$t_tundra),]
								pe=newD
				}
			}
		 
			}
					{# NOT USED with lines
					dev.new(width=3.5,height=1.85)
					#png(paste(out2,"Figure_4_lines.png", sep=""), width=3.5,height=1.85,units="in",res=600)
					
					  par(mfrow=c(1,3),mar=c(0.0,0,0,0.4),oma = c(1.8, 1.8, 0.5, 0.5),ps=12, mgp=c(1.2,0.35,0), las=1, cex=1, col.axis="grey30",font.main = 1, col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt
					{# time
					par(ps=12,	cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
					plot(ptt$pred~ptt$time_, pch=19,xlim=c(0,24), ylim=c(0,1), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
										
						axis(1, at=seq(0,24,by=6),labels=c(0,'',12,'',24),cex.axis=0.5,mgp=c(0,-0.15,0))
							mtext('Time [h]',side=1,line=0.4, cex=0.6, las=1, col='grey30')
							
						axis(2, at=seq(0,1,by=0.25), labels=c('0.0','','0.5','','1.0'))
							mtext('Nest attendance',side=2,line=1, cex=0.6, las=3, col='grey30')
						
						mtext(expression(bold("a")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
						
						# treated
								
								lines(ptt$time_, ptt$upr, col=col_t,lwd=1, lty=3)
								lines(ptt$time_, ptt$lwr, col=col_t,lwd=1, lty=3)
								lines(ptt$time_, ptt$pred, col=col_t,lwd=1)
								points(bt$time, bt$inc_eff, col=col_t,bg=adjustcolor(col_t ,alpha.f = 0.4), pch=21,cex=0.5)
						
						# control
								lines(ptc$time_, ptc$upr, col=col_c,lwd=1, lty=3)
								lines(ptc$time_, ptc$lwr, col=col_c,lwd=1, lty=3)
								lines(ptc$time_, ptc$pred, col=col_c,lwd=1)
								points(bb$time, bb$inc_eff, col=col_c,bg=adjustcolor(col_c ,alpha.f = 0.4), pch=21,cex=0.5)
								
						text(x=0,y=0.85, labels='control', col=col_c, cex=0.5, adj=0)
						text(x=0,y=0.75, labels='treated', col=col_t, cex=0.5, adj=0)		
						
					}
					{# proportion
						par(ps=12,	cex =1 cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
						plot(pp_$pred~pp_$prop, pch=19,xlim=c(40,55), ylim=c(0,1), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
										
						axis(1, at=seq(40,55,by=5),labels=seq(40,55,by=5),cex.axis=0.5,mgp=c(0,-0.15,0))
							mtext('Incubation share [%]',side=1,line=0.4, cex=0.6, las=1, col='grey30')
							
						axis(2, at=seq(0,1,by=0.25), labels=FALSE)
							#mtext('Nest attendance',side=2,line=1, cex=0.6, las=3, col='grey30')
						
						mtext(expression(bold("b")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
						
						# treated
								lines(pp_$prop*100, pp_$upr, col=col_t,lwd=1, lty=3)
								lines(pp_$prop*100, pp_$lwr, col=col_t,lwd=1, lty=3)
								lines(pp_$prop*100, pp_$pred, col=col_t,lwd=1)
								points(bt$prop*100, bt$inc_eff, col=col_t,bg=adjustcolor(col_t ,alpha.f = 0.4), pch=21,cex=0.5)
					}
					{# escape
						par(ps=12,	cex =1, cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
						plot(pe$pred~pe$esc, pch=19,xlim=c(0,65), ylim=c(0,1), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
										
						axis(1, at=seq(0,65,by=15),labels=c(0,'',30,'',60),cex.axis=0.5,mgp=c(0,-0.15,0))
							mtext('Escape distance [m]',side=1,line=0.4, cex=0.6, las=1, col='grey30')
							
						axis(2, at=seq(0,1,by=0.25), labels=FALSE)
							#mtext('Nest attendance',side=2,line=1, cex=0.6, las=3, col='grey30')
						
						mtext(expression(bold("c")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
						
						# treated
								lines(pe$esc, pe$upr, col=col_t,lwd=1, lty=3)
								lines(pe$esc, pe$lwr, col=col_t,lwd=1, lty=3)
								lines(pe$esc, pe$pred, col=col_t,lwd=1)
								points(bt$esc, bt$inc_eff, col=col_t,bg=adjustcolor(col_t ,alpha.f = 0.4), pch=21,cex=0.5)
					}
					   dev.off()
					}						
					{# based on a single model - not used
				{# proportion of incubation priot to removal
						m=lm(inc_eff~prop,bt)
						  
						# simulation		
							nsim <- 2000
							bsim <- sim(m, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
						# coefficients
							v <-apply(bsim@coef, 2, quantile, prob=c(0.5))
						# predicted values		
							newD=data.frame(prop=seq(0.40,0.55,0.001))
								
						# exactly the model which was used has to be specified here 
							X <- model.matrix(~ prop,data=newD)	
										
						# calculate predicted values and creditability intervals
							newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
									predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
									for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@coef[i,]
									newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
									newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
									#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
									#newD=newD[order(newD$t_tundra),]
								pp_=newD
				}
				{# escape distance
						m=lm(inc_eff~esc,bt)
						  
						# simulation		
							nsim <- 2000
							bsim <- sim(m, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
						# coefficients
							v <-apply(bsim@coef, 2, quantile, prob=c(0.5))
						# predicted values		
							newD=data.frame(esc=seq(0,65,0.2))
									
						# exactly the model which was used has to be specified here 
							X <- model.matrix(~ esc,data=newD)	
										
						# calculate predicted values and creditability intervals
							newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
									predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
									for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@coef[i,]
									newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
									newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
									#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
									#newD=newD[order(newD$t_tundra),]
								pe=newD
				}
				}
	
		}
	} 
	
	{# Post treatment effects -  only for nests where somebody returned
		
		{# number of days the widowed birds incubate after the treatment was over
			# s402 1.5 (died female)
			# s409 5
			# s510 10 (died female)
			# s516 9
			# s526 0
			# s711 4
			# s806 3
			summary(data.frame(nest=c('s402','s409','s510','s516','s526','s711','s806'),days=c(1.5,5,10,9,0,4,3)))
}
		
		{# prepare data
			u=read.csv(paste(wd,'experiment_metadata.csv', sep=""), stringsAsFactors=FALSE)	
			e=read.csv(file=paste(wd,'experiment_times.csv', sep=""), sep=",",stringsAsFactors =FALSE)
			e$taken=as.POSIXct(e$taken)
			e$release=as.POSIXct(e$release)
			load(file=paste(wd,'bout.Rdata', sep=""))
				#ggplot(b[!b$nest=='s806' & b$bout_type=='inc',], aes(x=inc_eff_3, y=inc_eff))+geom_point()+stat_smooth()
				#ggplot(b[!b$nest=='s806' & b$bout_type=='inc',], aes(x=inc_eff_3, y=inc_eff))+geom_point()+stat_smooth()+coord_cartesian(ylim = c(0.8,1), xlim = c(0.8,1))
			b$inc_eff[b$nest=='s702' & b$bout_ID==12]=b$inc_eff_2[b$nest=='s702' & b$bout_ID==12] # bird came for short time perioe and conctancy based on temperature was zero, this makes analyses dificult, hence we use constancy baased on tag readings
			d=b[-which(b$nest%in%c('s410','s514','s524','s624')),]
			d$bout_start=as.POSIXct(d$bout_start)
			d$bout_end=as.POSIXct(d$bout_end)
			# limit to nests where a removed bird returned
				d=d[-which(d$nest%in%c('s402','s409','s510','s516','s526','s711','s806')),] 
				e=e[-which(e$nest%in%c('s402','s409','s510','s516',"s526",'s711','s806')),] 
			
			# create zero exchange gaps where no exchange gap detected
			d=d[order(d$nest, d$bout_start),]
			d=ddply(d,.(nest), transform, before=as.character((c(NA,bout_type[-length(bout_type)]))))
			d$pk=c(1:nrow(d))
			l=list()
			for(i in 1:nrow(d)){
								di=d[i,]
								if(di$bout_type=='inc' & !is.na(di$before) & di$before=='inc'){
										di$pk=di$pk-0.5
										di$bout_ID=di$bout_ID-0.5
										di$bout_type='gap'
										di$bout_end=di$bout_start
										di$bout_length=di$inc_eff=di$inc_eff_2=di$inc_eff_3=di$disturb=di$disturb_log=di$dist_visit=di$dist_v_log=di$capture=0
										l[[i]]=di
										}
								}		
				d2=do.call(rbind, l)
			d=rbind(d,d2)
			d$before=NULL
			d=d[order(d$nest, d$bout_ID, d$bout_start),]
				# test correct assignment - YES
					#d_=ddply(d,.(nest), transform, before=as.character((c(NA,bout_type[-length(bout_type)]))))
					#d_=d_[which(!is.na(d_$before)),]
					#d_[d_$before==d_$bout_type,]
			
			li=list()
			lg=list()
			for(i in 1:nrow(e)){
							print(e$nest[i])
							ei=e[e$nest==e$nest[i],]
							di=d[d$nest==e$nest[i],]
							di=di[1:(nrow(di)-1),]
							{# treated bird 
								# before 
									db=di[di$sex==ei$sex & di$bout_end<(ei$taken-30*60),]
									db$pk=c(1:nrow(db))
									
								  # incubation bout
									dbt=db[db$bout_type=='inc' ,] # limit to bouts ending minimum 30 minutes before treatment capture
									dbt=dbt[(nrow(dbt)-2):nrow(dbt),] # take 3 bouts prior to treatment
									dbt$exper='b'
									dbt$treated_bird='y'
								  # gaps
											btg=db[db$bout_type=='gap' ,]
											btg=btg[(nrow(btg)-2):nrow(btg),] 
											btg$exper='b'
											btg$treated_bird='y'
								# after # if long gap and then still treated incubates, it is considered as his new bout
									da=di[di$sex==ei$sex & di$bout_start>(ei$release),]									
									da$pk=c(1:nrow(da))
								 # incubation bout
									dat=da[da$bout_type=='inc' ,] 
									dat=dat[1:3,] # take 3 bouts after to treatment
									dat$exper='a'
									dat$treated_bird='y'
								  # gaps
										atg=da[da$bout_type=='gap' ,] 
										atg=atg[1:3,]
										atg$exper='a'
										atg$treated_bird='y'
								inc_t=rbind(dbt, dat)
								gap_t=rbind(btg, atg)
								}
							{# control bird 
								# before 
									db=di[di$sex!=ei$sex & di$bout_end<(ei$taken-30*60),]
									db$pk=c(1:nrow(db))
								  # incubation bout
									dbt=db[db$bout_type=='inc' ,] # limit to bouts ending minimum 30 minutes before treatment capture
									dbt=dbt[(nrow(dbt)-2):nrow(dbt),] # take 3 bouts prior to treatment
									dbt$exper='b'
									dbt$treated_bird='n'
								 # gaps
											btg=db[db$bout_type=='gap' ,]
											btg=btg[(nrow(btg)-2):nrow(btg),] 
											btg$exper='b'
											btg$treated_bird='n'
								# after
									da=di[di$sex!=ei$sex & di$bout_start>(ei$release),]
									da$pk=c(1:nrow(da))
									
								  # incubation bout
									dat=da[da$bout_type=='inc' ,] # limit to bouts ending minimum 30 minutes before treatment capture
									dat=dat[1:3,] # take 3 bouts after to treatment
									dat$exper='a'
									dat$treated_bird='n'
								  # gaps
										atg=da[da$bout_type=='gap' ,] 
										atg=atg[1:3,]
										atg$exper='a'
										atg$treated_bird='n'
										
								inc_c=rbind(dbt, dat)
								gap_c=rbind(btg, atg)
								}
							li[[i]]=rbind(inc_t,inc_c)
							lg[[i]]=rbind(gap_t,gap_c)
							}
			inc=do.call(rbind,li)
			inc=inc[which(!is.na(inc$bout_length)),]
			
			# remove bouts where predation/desertion occured
				#d=d[-which(d$nest=='s519' & d$bout_start>'2013-06-30 08:00:00' | d$nest=='s712' & d$bout_end>'2013-06-28 01:00:00'),] 
				inc=inc[-which(inc$nest=='s519' & inc$bout_start>'2013-06-30 08:00:00' | inc$nest=='s712' & inc$bout_end>'2013-06-28 01:00:00'),] 
				inc=inc[order(inc$nest, inc$bout_start),]	
			# check which nest has less bouts 
				#s712 - nest was depradated and hence we have only 2 after bouts for each parent
			incsplit=split(inc, paste(inc$nest,inc$exper,inc$treated_bird, sep=" "))
			foo=lapply(incsplit,function(x) {
											#x=incsplit$"s808 b y"
											if(nrow(x)!=3){print(paste(x$nest[1],x$exper[1],x$treated_bird[1], sep=" "))}
											})
			gap=do.call(rbind,lg)
			gap=gap[which(!is.na(gap$bout_length)),]
			
			# remove gaps where predation/desertion occured
				#d=d[-which(d$nest=='s519' & d$bout_start>'2013-06-30 08:00:00' | d$nest=='s712' & d$bout_end>'2013-06-28 01:00:00'),] 
				gap=gap[!(gap$nest=='s519' & gap$bout_start>'2013-06-30 08:00:00' | gap$nest=='s712' & gap$bout_end>'2013-06-28 01:00:00'),] 
				#s712 - nest was depradated and hence we have only 2 after gaps for each parent
			gapsplit=split(gap, paste(gap$nest,gap$exper,gap$treated_bird, sep=" "))
			foo=lapply(gapsplit,function(x) {
											#x=gapsplit$"s808 b y"
											if(nrow(x)!=3){print(paste(x$nest[1],x$exper[1],x$treated_bird[1], sep=" "))}
											})
			gap=gap[order(gap$nest, gap$bout_start),]	
			load(file = paste(wd,"inc_start_estimation_2013_new.Rdata",sep=""))	
			e$nest=tolower(e$nest)
			
			inc$inc_start=e$inc_start[match(inc$nest,e$nest)]
			inc$bout_start=as.POSIXct(inc$bout_start)
			inc$inc_start=as.POSIXct(inc$inc_start)
			inc$bout_length=inc$bout_length/60
			inc$sex=as.factor(inc$sex)
			inc$exper=as.factor(inc$exper)
			inc$exper=relevel(inc$exper,ref="b")
			inc$treated_bird=as.factor(inc$treated_bird)
			inc$bird_ID=paste(inc$nest,inc$sex)
			
			# start of bout within incubation
				inc$bout_start_j  = as.numeric(format(inc$bout_start ,"%j")) - as.numeric(format(inc$inc_start,"%j")) +1
			# start of bout within season	
				inc$b_st_season_j = as.numeric(format(inc$bout_start,"%j")) - min(as.numeric(format(inc$inc_start,"%j")))+1 
			# start of incubation within season
				inc$inc_start_j = as.numeric(format(inc$inc_start,"%j")) - min(as.numeric(format(inc$inc_start,"%j")))+1
			
			gap$present=ifelse(gap$bout_length==0,0,1)
			gap$inc_start=e$inc_start[match(gap$nest,e$nest)]
			gap$bout_start=as.POSIXct(gap$bout_start)
			gap$inc_start=as.POSIXct(gap$inc_start)
			gap$sex=as.factor(gap$sex)
			gap$sex=as.factor(gap$sex)
			gap$exper=as.factor(gap$exper)
			gap$exper=relevel(gap$exper,ref="b")
			gap$treated_bird=as.factor(gap$treated_bird)
			gap$bird_ID=paste(gap$nest,gap$sex)
			# start of bout within oncpubation
				gap$bout_start_j  = as.numeric(format(gap$bout_start ,"%j")) - as.numeric(format(gap$inc_start,"%j")) +1
			# start of bout within season	
				gap$b_st_season_j = as.numeric(format(gap$bout_start,"%j")) - min(as.numeric(format(gap$inc_start,"%j")))+1 
			# start of gapubation within season
				gap$inc_start_j = as.numeric(format(gap$inc_start,"%j")) - min(as.numeric(format(gap$inc_start,"%j")))+1
				
			# bout_start mean centered within nest
						incsplit=split(inc,paste(inc$nest))
									foo=lapply(incsplit,function(x) {
																	#x=incsplit$"s404"
																	#x$t_a=c(x$treat[-1], NA) 	
																	x$bout_start_j_c=x$bout_start_j - mean(x$bout_start_j)
																	return(x)
																	}
																	)
								
								inc=do.call(rbind, foo)	
						gapsplit=split(gap,paste(gap$nest))
									foo=lapply(gapsplit,function(x) {
																	#x=gapsplit$"s404"
																	#x$t_a=c(x$treat[-1], NA) 	
																	x$bout_start_j_c=x$bout_start_j - mean(x$bout_start_j)
																	return(x)
																	}
																	)
								
								gap=do.call(rbind, foo)			
			gapp=gap[gap$bout_length>0,] # data with detected exchange gap
			
			# scale within experimental period
					mean_b=mean(inc$bout_length[inc$exper=='b'])
					mean_a=mean(inc$bout_length[inc$exper=='a'])
					sd_b=sd(inc$bout_length[inc$exper=='b'])
					sd_a=sd(inc$bout_length[inc$exper=='a'])
					inc$bout_length_c=ifelse(inc$exper=='b',(inc$bout_length-mean_b)/sd_b,(inc$bout_length-mean_a)/sd_a)
					
					mean_b=mean(inc$inc_eff[inc$exper=='b'])
					mean_a=mean(inc$inc_eff[inc$exper=='a'])
					sd_b=sd(inc$inc_eff[inc$exper=='b'])
					sd_a=sd(inc$inc_eff[inc$exper=='a'])
					inc$inc_eff_c=ifelse(inc$exper=='b',(inc$inc_eff-mean_b)/sd_b,(inc$inc_eff-mean_a)/sd_a)
				
					mean_b=mean(gapp$bout_length[gapp$exper=='b'])
					mean_a=mean(gapp$bout_length[gapp$exper=='a'])
					sd_b=sd(gapp$bout_length[gapp$exper=='b'])
					sd_a=sd(gapp$bout_length[gapp$exper=='a'])
					gapp$bout_length_c=ifelse(gapp$exper=='b',(gapp$bout_length-mean_b)/sd_b,(gapp$bout_length-mean_a)/sd_a)
						
						inc$lat=u$lat[match(inc$nest,u$nest)]
						inc$lon=u$lon[match(inc$nest,u$nest)]		
						gap$lat=u$lat[match(gap$nest,u$nest)]
						gap$lon=u$lon[match(gap$nest,u$nest)]						
						gapp$lat=u$lat[match(gapp$nest,u$nest)]
						gapp$lon=u$lon[match(gapp$nest,u$nest)]		
		}
		{# define colors
				col_t='#FCB42C'
				col_t_light='#fef0d4'
				col_c='#535F7C'
				col_c_light='#cbcfd7'
				col_2011=rgb(0,100,0,255, maxColorValue=255)
				col_2011_light=rgb(102,162,102,255, maxColorValue=255)
				inc$col_=ifelse(inc$exper=='b', col_t,col_c)
				gapp$col_=ifelse(gapp$exper=='b', col_t,col_c)
			}
		
		# starved birds (s402, s510 deaath
			p = data.frame(bird_ID =c(257188558,257188571,257188564,255174350,257188570,257188566,255174353,257188572) , nest = c('s502','s702','s509','s711','s402','s510','s516','s404'), stringsAsFactors=FALSE)
		{# how long it took birds to return - also relative mass loss (N = 25 nests minus two were female died)
			ee=read.csv(file=paste(wd,'experiment_times.csv', sep=""), sep=",",stringsAsFactors =FALSE)
			ee$release=as.POSIXct(ee$release)
			ee$starved = ifelse(ee$nest%in%p$nest, 'yes','no')
			ee[ee$starved =='yes',]
			nrow(ee[ee$starved =='no',])
			ee_ = ee[ee$sex=='m',]
			table(ee_$after, ee_$starved)
			
			inc$release=ee$release[match(inc$nest,ee$nest)]
			inc_=ddply(inc[inc$exper=='a' & inc$treated_bird=='n',],.(nest), transform, st_=min(bout_start))
			inc_=inc_[inc_$bout_start==inc_$st_,]
			
			nrow(inc_)
			inc_[,c('nest','release','bout_start')]
			summary(as.numeric(difftime(inc_$bout_start,inc_$release,unit='hours')))
			inc_$return_time = as.numeric(difftime(inc_$bout_start,inc_$release,unit='hours'))
			inc_$starved = ifelse(inc_$nest%in%p$nest, 'yes','no')
			summary(factor(inc_$starved))
			
			ee$return_time = inc_$return_time[match(ee$nest, inc_$nest)]
			ggplot(inc_,aes(x = starved, y = return_time)) + geom_boxplot() +geom_point(aes(fill = starved),shape = 21,  position = position_jitterdodge())
			ggplot(ee[ee$after=='return',],aes(x = starved, y = return_time)) + geom_boxplot() +geom_point(aes(fill = starved),shape = 21,  position = position_jitterdodge())
			
			ggplot(ee[ee$after=='return',],aes(x = starved, y = return_time)) + geom_boxplot() +  geom_dotplot(aes(fill = starved),binaxis = 'y', stackdir = 'center',  position = position_dodge())+
			ylab('Return time [h]')+xlab('Starved')+
			scale_y_continuous(limits = c(0,20))+ 	
			theme( 
				legend.position="none",
				#axis.text.x = element_text(colour = "white"), 
				#axis.title.x = element_text(colour = "white"),
				#axis.ticks.x = element_blank(),
				axis.line.x = element_line(color="black", size = 0.25),
				axis.title.y = element_text(margin = margin( r = 9)),
				axis.line.y = element_line(color="black", size = 0.25),
				panel.background=element_blank(),
				panel.border=element_blank(),
				panel.grid.major=element_blank(),
				panel.grid.minor=element_blank(),
				plot.background=element_blank()
				) 
			ggsave(paste(out2,"Supplementary_Figure_2.png", sep=""),width = 5, height = 6, units = "cm",  dpi = 300)
			
			
			u=read.csv(paste(wd,'experiment_metadata.csv', sep=""), stringsAsFactors=FALSE)	
			d =read.csv(paste(wd,'captivity.csv', sep=""),stringsAsFactors=FALSE)
				dd = d[!is.na(d$mass) & d$phase %in%c('start','end'),]
				v = ddply(dd,.(ring_num), summarise, rel_mass = (mass[phase=='end'] - mass[phase=='start'])/mass[phase=='start'], abs_mass = (mass[phase=='end'] - mass[phase=='start']))
			u$rel_mass = v$rel_mass[match(u$ID_taken, v$ring_num)]
			ee$mass_loss = 	-u$mass_loss[match(ee$nest, u$nest)]	#uu = u[9:nrow(u),]
			ee$rel_mass_loss = 	u$rel_mass[match(ee$nest, u$nest)]	#uu = u[9:nrow(u),]
			ee$returned = ifelse(ee$after == 'return','yes','no')
			ee$sex = ifelse(ee$sex =='f','m','f')
			ggplot(ee[ee$starved=='no',],aes(x = returned, y = mass_loss)) + geom_boxplot() +  geom_dotplot(aes(fill = sex),binaxis = 'y', stackdir = 'center',  position = position_dodge())+ylab('mass change [g]')
			ggplot(ee[ee$starved=='no',],aes(x = returned, y = mass_loss, fill=sex)) + geom_boxplot() +  geom_dotplot(aes(fill = sex),binaxis = 'y', stackdir = 'center',  position = position_dodge())+ylab('mass change [g]')
			
			
			med = ddply(ee,.(sex,starved,returned),summarise,mass_loss = median(mass_loss, na.rm =TRUE))
			ggplot(ee,aes(x = returned, y = mass_loss, fill=sex)) + geom_boxplot() + geom_dotplot(aes(fill = sex, shape = starved),binaxis = 'y', stackdir = 'center',  position = position_dodge())+ylab('mass change [g]')
			ggplot(ee,aes(x = returned, y = mass_loss, fill=sex)) + geom_boxplot() +  geom_dotplot(aes(fill = paste(sex, starved)),binaxis = 'y', stackdir = 'center',  position = position_dodge())+ylab('mass change [g]')
			
			ggplot(ee,aes(x = returned, y = mass_loss, fill=paste(sex, starved))) + geom_boxplot() + geom_dotplot(aes(fill = paste(sex, starved)),binaxis = 'y', stackdir = 'center',  position = position_dodge(),dotsize = 1.5)+xlab('Returned')+ylab('Mass change [g]') 
			dev.new(width = 2.3622, height =1.9685)
			ggplot(ee,aes(x = returned, y = mass_loss, fill=paste(sex, starved))) +  geom_dotplot(aes(fill = paste(sex, starved)),binaxis = 'y', stackdir = 'center',  position = position_dodge(),dotsize = 1.5)+xlab('Returned')+ylab('Mass change [g]') +geom_point(data = med, col = 'black',shape= 95, size =5 ) + #guides(shape = guide_legend(override.aes = list(size = 0.2)))+
			theme( 
				axis.ticks.length=unit(0.5,"mm"),
				#legend.position="none",
				#axis.text.x = element_text(colour = "white"), 
				#axis.title.x = element_text(colour = "white"),
				#axis.ticks.x = element_blank(),
				axis.line.x = element_line(color="grey40", size = 0.25),
				axis.text.x = element_text(size=6, colour = "grey20"),
				axis.title.x = element_text(size=7, colour = "grey20"),
				
				axis.title.y = element_text(size=7, color='grey20', margin = margin( r = 6)),
				axis.line.y = element_line(color="black", size = 0.25),
				axis.text.y=element_text(size=6, color='grey20'),# margin=units(0.5,"mm")),
				
				legend.key = element_blank(),
				legend.text = element_text(size=6, colour = "grey20"),
				legend.title = element_blank(),#element_text(size=7, colour = "grey20"),
				legend.key.size = unit(0.1,"points"),
				panel.background=element_blank(),
				panel.border=element_blank(),
				panel.grid.major=element_blank(),
				panel.grid.minor=element_blank(),
				plot.background=element_blank()
				) 
			ggsave(paste(out2,"Supplementary_Figure_5new.png", sep=""),width = 6, height = 5, units = "cm",  dpi = 300)
			
		
			ggplot(ee[ee$starved=='no',],aes(x = returned, y = rel_mass_loss)) + geom_boxplot() +  geom_dotplot(aes(fill = sex),binaxis = 'y', stackdir = 'center',  position = position_dodge())+ylab('relative mass change')
			ggplot(ee[ee$starved=='no',],aes(x = returned, y = rel_mass_loss, fill=sex)) + geom_boxplot() +  geom_dotplot(aes(fill = sex),binaxis = 'y', stackdir = 'center',  position = position_dodge())+ylab('relative mass change')
		
			ggplot(ee[ee$starved=='no',],aes(x = mass_loss, y = return_time, col=sex)) + geom_point() +  stat_smooth(method = 'lm') + ylab('return time [h]')+ xlab('mass change [g]')
			ggplot(ee[ee$starved=='no',],aes(x = mass_loss, y = return_time)) + geom_point() +  stat_smooth(method = 'lm') + ylab('return time [h]')+ xlab('mass change [g]')
			ggplot(ee,aes(x = mass_loss, y = return_time, col=sex)) + geom_point() +  stat_smooth(method = 'lm') + ylab('return time [h]')+ xlab('mass change [g]')
			
			densityplot(ee$mass_loss, group = ee$sex)
			densityplot(ee$return_time, group = ee$sex)
			
			ee$sex = as.factor(ee$sex)
			m = lm(return_time~mass_loss*sex,ee[ee$starved=='no',])
			m0 = lm(return_time~mass_loss,ee[ee$starved=='no',])
			AICc(m)
			AICc(m0)
				summary(glht(m))
				summary(m)
				plot(allEffects(m))
				
			m = lm(return_time~mass_loss*sex,ee)
			m0 = lm(return_time~mass_loss,ee)
			AICc(m)
			AICc(m0)
				summary(glht(m))
				summary(m0)
				plot(allEffects(m))	
			ggplot(ee[ee$starved=='no',],aes(x = rel_mass_loss, y = return_time, col=sex)) + geom_point() +  stat_smooth(method = 'lm') + ylab('return time [h]')+ xlab('relative mass change')
			ggplot(ee[ee$starved=='no',],aes(x = rel_mass_loss, y = return_time)) + geom_point() +  stat_smooth(method = 'lm') + ylab('return time [h]')+ xlab('relative mass change')
			ggplot(ee,aes(x = rel_mass_loss, y = return_time, col=sex)) + geom_point() +  stat_smooth(method = 'lm') + ylab('return time [h]')+ xlab('relative mass change')
			{# Supplementary Table Z
				nrow(ee[!(is.na(ee$return_time)|is.na(ee$mass_loss)),])
				m = lm(return_time~scale(mass_loss)*sex,ee)	
					pred=c('Intercept (female)','Mass_loss','Sex (male)','Mass:Sex')
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
				# Fixed effects
					v <- apply(bsim@coef, 2, quantile, prob=c(0.5))
					ci=apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
						oi=data.frame(model='1',response='return time',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
					rownames(oi) = NULL
						oi$estimate_r=round(oi$estimate,2)
						oi$lwr_r=round(oi$lwr,2)
						oi$upr_r=round(oi$upr,2)
						#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
					oii=oi[c('model','response','type',"effect", "estimate_r","lwr_r",'upr_r')]	
							sname = tempfile(fileext='.xls')
							wb = loadWorkbook(sname,create = TRUE)	
							createSheet(wb, name = "output")
							writeWorksheet(wb, oii, sheet = "output")
							saveWorkbook(wb)
							shell(sname)
			
			}
			}
		{# effect of cage on before treatment
			cc =  read.csv(paste(wd,'cages.csv', sep=""), stringsAsFactors = FALSE)
			cc$on = as.POSIXct(cc$on)
			cc$off = as.POSIXct(cc$off)
			incb = inc[inc$exper=='b',]
			incb$cage  = 'n'
			for(i in 1:nrow(incb)){
				n = incb$nest[i]
				ci = cc[cc$nest == n,]
				if(nrow(ci)!=0){
				incb$cage[i] = ifelse(incb$bout_start[i]>ci$on & incb$bout_start[i]<ci$off, 'y','n')
				}
				print(i)
			}
			summary(factor(incb$cage))
			table(incb$nest, incb$cage)
			ggplot(incb,aes(x = cage, y = bout_length)) + geom_boxplot() + geom_dotplot(aes(fill = cage),binaxis = 'y', stackdir = 'center',  position = position_dodge())
			ggplot(incb,aes(x = cage, y = inc_eff)) + geom_boxplot() + geom_dotplot(aes(fill = cage),binaxis = 'y', stackdir = 'center',  position = position_dodge())
			
		     incb$cage = as.factor(incb$cage)
			 m1=lmer(bout_length~cage+(1|nest)+(1|bird_ID),incb, REML=FALSE)
			 m1=lmer(inc_eff~cage+(1|nest)+(1|bird_ID),incb, weights=sqrt(bout_length),REML=FALSE)
			 plot(allEffects(m1))
			 summary(glht(m1))
			 m1=lmer(bout_length~cage+(1|bird_ID),incb, REML=FALSE)
		}
		{# Supplementary Table 4
			{# constancy
				{# 1 - model of interest
				 m1=lmer(inc_eff~exper*bout_start_j_c+(bout_start_j_c|bird_ID),inc, weights=sqrt(bout_length),REML=FALSE)
				 #mv=lme(inc_eff~exper*bout_start_j_c, data=inc, random=~1|bird_ID, weights=varIdent(form=~1|exper), method="ML") # best
				 pred=c('Intercept (before)','Period (after)','Day (centered)','Period:Day')
							nsim <- 2000
							bsim <- sim(m1, n.sim=nsim)  
					# Fixed effects
						v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))
						ci=apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
						oi=data.frame(model='1',response='constancy',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,], stringsAsFactors=FALSE)
						rownames(oi) = NULL
							oi$estimate_r=round(oi$estimate,3)
							oi$lwr_r=round(oi$lwr,3)
							oi$upr_r=round(oi$upr,3)
							#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
						oii=oi[c('model','response','type',"effect", "estimate_r","lwr_r",'upr_r')]	
					# Random effects
						l=data.frame(summary(m1)$varcor)
						l=l[is.na(l$var2),]
						l$var2=NULL
						ri=data.frame(model='1',response='constancy',type='random (var)',effect=l$var, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA, stringsAsFactors=FALSE)
							ri$effect[nrow(ri)]='residual'
						o1=rbind(oii,ri)
				
							
				}	
				{# 2 - simple model 
					 m=lmer(inc_eff~exper+(1|bird_ID),inc, REML=FALSE)
					 pred=c('Intercept (before)','Period (after)')
								nsim <- 2000
								bsim <- sim(m, n.sim=nsim)  
						# Fixed effects
							v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))
							ci=apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							oi=data.frame(model='2',response='constancy',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
							rownames(oi) = NULL
								oi$estimate_r=round(oi$estimate,3)
								oi$lwr_r=round(oi$lwr,3)
								oi$upr_r=round(oi$upr,3)
								#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
							oii=oi[c('model','response','type',"effect", "estimate_r","lwr_r",'upr_r')]	
						# Random effects
							l=data.frame(summary(m)$varcor)
								ri=data.frame(model='2',response='constancy', type='random (var)',effect=l$grp, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA)
							o2=rbind(oii,ri)
				}
				{# 3 - model with bird (focal/removed)
				 m3=lmer(inc_eff~exper*treated_bird*bout_start_j_c+(bout_start_j_c|bird_ID),inc, REML=FALSE)
				 pred=c('Intercept (before,removed)','Period (after)','Bird(focal)','Day','Period:bird', 'Period:Day','bird:Day','period:bird:day')
							nsim <- 2000
							bsim <- sim(m3, n.sim=nsim)  
					# Fixed effects
						v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))
						ci=apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
						oi=data.frame(model='3',response='constancy',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,], stringsAsFactors=FALSE)
						rownames(oi) = NULL
							oi$estimate_r=round(oi$estimate,3)
							oi$lwr_r=round(oi$lwr,3)
							oi$upr_r=round(oi$upr,3)
							#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
						oii=oi[c('model','response','type',"effect", "estimate_r","lwr_r",'upr_r')]	
					# Random effects
						l=data.frame(summary(m3)$varcor)
						l=l[is.na(l$var2),]
						l$var2=NULL
						ri=data.frame(model='1',response='constancy',type='random (var)',effect=l$var, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA, stringsAsFactors=FALSE)
							ri$effect[nrow(ri)]='residual'
						o3=rbind(oii,ri)

				}				
				  o_con=(rbind(o1,o2,o3))
				{# AICc # number of observations = number of birds*2 (2 = before and after period) length(unique(paste(inc$bird_ID,inc$exper)))
						o=data.frame(response='constancy',model=c('period:day','period','period:bird:day'), AIC=c(AICc(m1, nobs=36*2),AICc(m,nobs=36*2),AICc(m3,nobs=36*2)))
						o$delta=o$AIC-min(o$AIC)
						o$prob=exp(-0.5*o$delta)/sum(exp(-0.5*o$delta))
						o$ER=max(o$prob)/o$prob
						o$AIC=round(o$AIC,1)
						o$delta=round(o$delta,2)
						o$prob=round(o$prob,3)
						o$ER=round(o$ER,2)
						o_1=o
						#o[order(o$delta),]
						
						
						
				}
			}	
			{# bout_length
				{# 1 - model of interest
				 m1=lmer(bout_length~exper*bout_start_j_c+(bout_start_j_c|bird_ID),inc, REML=FALSE)
				 #mv=lme(inc_eff~exper*bout_start_j_c, data=inc, random=~bout_start_j_c|bird_ID, weights=varIdent(form=~1|exper), method="ML") # best
				 pred=c('Intercept (before)','Period (after)','Day (centered)','Period:Day')
							nsim <- 2000
							bsim <- sim(m1, n.sim=nsim)  
					# Fixed effects
						v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))
						ci=apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
						oi=data.frame(model='1',response='bout_length',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,], stringsAsFactors=FALSE)
						rownames(oi) = NULL
							oi$estimate_r=round(oi$estimate,2)
							oi$lwr_r=round(oi$lwr,2)
							oi$upr_r=round(oi$upr,2)
							#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
						oii=oi[c('model','response','type',"effect", "estimate_r","lwr_r",'upr_r')]	
					# Random effects
						l=data.frame(summary(m1)$varcor)
						l=l[is.na(l$var2),]
						l$var2=NULL
						ri=data.frame(model='1',response='bout_length',type='random (var)',effect=l$var, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA, stringsAsFactors=FALSE)
							ri$effect[nrow(ri)]='residual'
						o1=rbind(oii,ri)
				
							
				}	
				{# 2 - simple model 
					 m=lmer(bout_length~exper+(1|bird_ID),inc, REML=FALSE)
					 pred=c('Intercept (before)','Period (after)')
								nsim <- 2000
								bsim <- sim(m, n.sim=nsim)  
						# Fixed effects
							v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))
							ci=apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							oi=data.frame(model='2',response='bout_length',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
							rownames(oi) = NULL
								oi$estimate_r=round(oi$estimate,2)
								oi$lwr_r=round(oi$lwr,2)
								oi$upr_r=round(oi$upr,2)
								#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
							oii=oi[c('model','response','type',"effect", "estimate_r","lwr_r",'upr_r')]	
						# Random effects
							l=data.frame(summary(m)$varcor)
								ri=data.frame(model='2',response='bout_length', type='random (var)',effect=l$grp, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA)
							o2=rbind(oii,ri)
				}
				{# 3 - model with bird (focal/removed)
				 m3=lmer(bout_length~exper*treated_bird*bout_start_j_c+(bout_start_j_c|bird_ID),inc, REML=FALSE)
				 pred=c('Intercept (before,removed)','Period (after)','Bird(focal)','Day','Period:bird', 'Period:Day','bird:Day','period:bird:day')
							nsim <- 2000
							bsim <- sim(m3, n.sim=nsim)  
					# Fixed effects
						v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))
						ci=apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
						oi=data.frame(model='3',response='bout_length',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,], stringsAsFactors=FALSE)
						rownames(oi) = NULL
								oi$estimate_r=round(oi$estimate,2)
								oi$lwr_r=round(oi$lwr,2)
								oi$upr_r=round(oi$upr,2)
							#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
						oii=oi[c('model','response','type',"effect", "estimate_r","lwr_r",'upr_r')]	
					# Random effects
						l=data.frame(summary(m3)$varcor)
						l=l[is.na(l$var2),]
						l$var2=NULL
						ri=data.frame(model='1',response='bout_length',type='random (var)',effect=l$var, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA, stringsAsFactors=FALSE)
							ri$effect[nrow(ri)]='residual'
						o3=rbind(oii,ri)

				}				
				   o_bout=(rbind(o1,o2,o3))
				{# AICc # number of observations = number of birds*2 (2 = before and after period) length(unique(paste(inc$bird_ID,inc$exper)))
						o=data.frame(response='bout_length',model=c('period:day','period','period:bird:day'), AIC=c(AICc(m1, nobs=36*2),AICc(m,nobs=36*2),AICc(m3,nobs=36*2)))
						o$delta=o$AIC-min(o$AIC)
						o$prob=round(exp(-0.5*o$delta)/sum(exp(-0.5*o$delta)),4)
						o$ER=max(o$prob)/o$prob
						o$AIC=round(o$AIC,1)
						o$delta=round(o$delta,2)
						o$prob=round(o$prob,3)
						o$ER=round(o$ER,2)
						o_2=o #o[order(o$delta),]
						
						
				}
			}				
			{# combine and export to excel table
							sname = tempfile(fileext='.xls')
							wb = loadWorkbook(sname,create = TRUE)	
							createSheet(wb, name = "output")
							writeWorksheet(wb, rbind(o_con,o_bout), sheet = "output")
							createSheet(wb, name = "output_AIC")
							writeWorksheet(wb, rbind(o_1,o_2), sheet = "output_AIC")
							saveWorkbook(wb)
							shell(sname)
									
			}
						
			{# model assumptions
				{# constancy
					{# 1
						dev.new(width=6,height=9)
									
						m=lmer(inc_eff~exper*bout_start_j_c+(bout_start_j_c|bird_ID),inc, REML=FALSE)
						par(mfrow=c(4,3))
						
						scatter.smooth(fitted(m),resid(m),col='red');abline(h=0, lty=2)
						scatter.smooth(fitted(m),sqrt(abs(resid(m))), col='red')
													
						qqnorm(resid(m), main=list("Normal Q-Q Plot: residuals", cex=0.8),col='red') 
						qqline(resid(m))
			 
						qqnorm(unlist(ranef(m)$bird_ID[1]), main = " Bird")
						qqline(unlist(ranef(m)$bird_ID[1]))
						
						qqnorm(unlist(ranef(m)$bird_ID[2]), main = " Slope")
						qqline(unlist(ranef(m)$bird_ID[2]))
						
						scatter.smooth(resid(m)~inc$bout_start_j_c);abline(h=0, lty=2, col='red')
						scatter.smooth(resid(m)~inc$exper);abline(h=0, lty=2, col='red')
						plot(resid(m)~factor(inc$exper)); abline(h=0, lty=2, col='red')
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=inc$lon, y=inc$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
					}
					{# 2
						dev.new(width=6,height=9)
									
						m=lmer(inc_eff~exper+(1|bird_ID),inc, REML=FALSE)
						par(mfrow=c(4,3))
						
						scatter.smooth(fitted(m),resid(m),col='red');abline(h=0, lty=2)
						scatter.smooth(fitted(m),sqrt(abs(resid(m))), col='red')
													
						qqnorm(resid(m), main=list("Normal Q-Q Plot: residuals", cex=0.8),col='red') 
						qqline(resid(m))
			 
						qqnorm(unlist(ranef(m)$bird_ID[1]), main = " Bird")
						qqline(unlist(ranef(m)$bird_ID[1]))
						
						scatter.smooth(resid(m)~inc$exper);abline(h=0, lty=2, col='red')
						plot(resid(m)~factor(inc$exper)); abline(h=0, lty=2, col='red')
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=inc$lon, y=inc$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
					}
					{# 3
						dev.new(width=6,height=9)
									
						m=lmer(inc_eff~exper*treated_bird*bout_start_j_c+(bout_start_j_c|bird_ID),inc, REML=FALSE)
						par(mfrow=c(4,4))
						
						scatter.smooth(fitted(m),resid(m),col='red');abline(h=0, lty=2)
						scatter.smooth(fitted(m),sqrt(abs(resid(m))), col='red')
													
						qqnorm(resid(m), main=list("Normal Q-Q Plot: residuals", cex=0.8),col='red') 
						qqline(resid(m))
			 
						qqnorm(unlist(ranef(m)$bird_ID[1]), main = " Bird")
						qqline(unlist(ranef(m)$bird_ID[1]))
						
						qqnorm(unlist(ranef(m)$bird_ID[2]), main = " Slope")
						qqline(unlist(ranef(m)$bird_ID[2]))
						
						scatter.smooth(resid(m)~inc$bout_start_j_c);abline(h=0, lty=2, col='red')
						scatter.smooth(resid(m)~inc$exper);abline(h=0, lty=2, col='red')
						plot(resid(m)~factor(inc$exper)); abline(h=0, lty=2, col='red')
						plot(resid(m)~factor(inc$treated_bird)); abline(h=0, lty=2, col='red')
						plot(resid(m)~interaction(inc$exper,inc$treated_bird)); abline(h=0, lty=2, col='red')
						
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=inc$lon, y=inc$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
					}
				}	
				{# bout - auto-cor
					{# 1
						dev.new(width=6,height=9)
									
						m=lmer(bout_length~exper*bout_start_j_c+(bout_start_j_c|bird_ID),inc, REML=FALSE)
						par(mfrow=c(4,3))
						
						scatter.smooth(fitted(m),resid(m),col='red');abline(h=0, lty=2)
						scatter.smooth(fitted(m),sqrt(abs(resid(m))), col='red')
													
						qqnorm(resid(m), main=list("Normal Q-Q Plot: residuals", cex=0.8),col='red') 
						qqline(resid(m))
			 
						qqnorm(unlist(ranef(m)$bird_ID[1]), main = " Bird")
						qqline(unlist(ranef(m)$bird_ID[1]))
						
						qqnorm(unlist(ranef(m)$bird_ID[2]), main = " Slope")
						qqline(unlist(ranef(m)$bird_ID[2]))
						
						scatter.smooth(resid(m)~inc$bout_start_j_c);abline(h=0, lty=2, col='red')
						scatter.smooth(resid(m)~inc$exper);abline(h=0, lty=2, col='red')
						plot(resid(m)~factor(inc$exper)); abline(h=0, lty=2, col='red')
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=inc$lon, y=inc$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
					}
					{# 2
						dev.new(width=6,height=9)
									
						m==lmer(bout_length~exper+(1|bird_ID),inc, REML=FALSE)
						par(mfrow=c(4,3))
						
						scatter.smooth(fitted(m),resid(m),col='red');abline(h=0, lty=2)
						scatter.smooth(fitted(m),sqrt(abs(resid(m))), col='red')
													
						qqnorm(resid(m), main=list("Normal Q-Q Plot: residuals", cex=0.8),col='red') 
						qqline(resid(m))
			 
						qqnorm(unlist(ranef(m)$bird_ID[1]), main = " Bird")
						qqline(unlist(ranef(m)$bird_ID[1]))
						
						scatter.smooth(resid(m)~inc$exper);abline(h=0, lty=2, col='red')
						plot(resid(m)~factor(inc$exper)); abline(h=0, lty=2, col='red')
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=inc$lon, y=inc$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
					}
					{# 3
						dev.new(width=6,height=9)
									
						m=lmer(bout_length~exper*treated_bird*bout_start_j_c+(bout_start_j_c|bird_ID),inc, REML=FALSE)
						par(mfrow=c(4,4))
						
						scatter.smooth(fitted(m),resid(m),col='red');abline(h=0, lty=2)
						scatter.smooth(fitted(m),sqrt(abs(resid(m))), col='red')
													
						qqnorm(resid(m), main=list("Normal Q-Q Plot: residuals", cex=0.8),col='red') 
						qqline(resid(m))
			 
						qqnorm(unlist(ranef(m)$bird_ID[1]), main = " Bird")
						qqline(unlist(ranef(m)$bird_ID[1]))
						
						qqnorm(unlist(ranef(m)$bird_ID[2]), main = " Slope")
						qqline(unlist(ranef(m)$bird_ID[2]))
						
						scatter.smooth(resid(m)~inc$bout_start_j_c);abline(h=0, lty=2, col='red')
						scatter.smooth(resid(m)~inc$exper);abline(h=0, lty=2, col='red')
						plot(resid(m)~factor(inc$exper)); abline(h=0, lty=2, col='red')
						plot(resid(m)~factor(inc$treated_bird)); abline(h=0, lty=2, col='red')
						plot(resid(m)~interaction(inc$exper,inc$treated_bird)); abline(h=0, lty=2, col='red')
						
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=inc$lon, y=inc$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
					}
				}
			}					
		}
		{# Supplementary Table 5			
			{# gap presence
				table(gap$bird_ID_filled, gap$exper)
				length(unique(gap$bird_ID_filled))
				{# 1 - model of interest
				 m1=glmer(present~exper*bout_start_j_c+(bout_start_j_c|bird_ID),gap, family='binomial',control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) # uses different optimizer (results same as without it, but without it model does not converge)
				 pred=c('Intercept (before)','Period (after)','Day (centered)','Period:Day')
							nsim <- 2000
							bsim <- sim(m1, n.sim=nsim)  
					# Fixed effects
						v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))
						ci=apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
						oi=data.frame(model='1',response='gap_y_n',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,], stringsAsFactors=FALSE)
						rownames(oi) = NULL
							oi$estimate_r=round(oi$estimate,2)
							oi$lwr_r=round(oi$lwr,2)
							oi$upr_r=round(oi$upr,2)
							#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
						oii=oi[c('model','response','type',"effect", "estimate_r","lwr_r",'upr_r')]	
					# Random effects
						l=data.frame(summary(m1)$varcor)
						l=l[is.na(l$var2),]
						l$var2=NULL
						ri=data.frame(model='1',response='gap_y_n',type='random (var)',effect=l$var, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA, stringsAsFactors=FALSE)
							ri$effect[nrow(ri)]='residual'
						o1=rbind(oii,ri)
				
							
				}	
				{# 2 - simple model 
					 m=glmer(present~exper+(1|bird_ID),gap, family='binomial')
					 pred=c('Intercept (before)','Period (after)')
								nsim <- 2000
								bsim <- sim(m, n.sim=nsim)  
						# Fixed effects
							v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))
							ci=apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							oi=data.frame(model='2',response='gap_y_n',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
							rownames(oi) = NULL
								oi$estimate_r=round(oi$estimate,2)
								oi$lwr_r=round(oi$lwr,2)
								oi$upr_r=round(oi$upr,2)
								#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
							oii=oi[c('model','response','type',"effect", "estimate_r","lwr_r",'upr_r')]	
						# Random effects
							l=data.frame(summary(m)$varcor)
								ri=data.frame(model='2',response='gap_y_n', type='random (var)',effect=l$grp, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA)
							o2=rbind(oii,ri)
				}
				{# 3 - model with bird (focal/removed)
				 m3=glmer(present~exper*treated_bird*bout_start_j_c+(bout_start_j_c|bird_ID),gap, family='binomial',control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
				 pred=c('Intercept (before,removed)','Period (after)','Bird(focal)','Day','Period:bird', 'Period:Day','bird:Day','period:bird:day')
							nsim <- 2000
							bsim <- sim(m3, n.sim=nsim)  
					# Fixed effects
						v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))
						ci=apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
						oi=data.frame(model='3',response='gap_y_n',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,], stringsAsFactors=FALSE)
						rownames(oi) = NULL
								oi$estimate_r=round(oi$estimate,2)
								oi$lwr_r=round(oi$lwr,2)
								oi$upr_r=round(oi$upr,2)
							#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
						oii=oi[c('model','response','type',"effect", "estimate_r","lwr_r",'upr_r')]	
					# Random effects
						l=data.frame(summary(m3)$varcor)
						l=l[is.na(l$var2),]
						l$var2=NULL
						ri=data.frame(model='3',response='gap_y_n',type='random (var)',effect=l$var, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA, stringsAsFactors=FALSE)
							ri$effect[nrow(ri)]='residual'
						o3=rbind(oii,ri)

				}				
				   o_gapyn=(rbind(o1,o2,o3))
				{# AICc # number of observations = number of birds*2 (2 = before and after period) length(unique(paste(gap$bird_ID,gap$exper)))
						o=data.frame(response='gapyn',model=c('period:day','period','period:bird:day'), AIC=c(AICc(m1, nobs=36*2),AICc(m,nobs=36*2),AICc(m3,nobs=36*2)))
						o$delta=o$AIC-min(o$AIC)
						o$prob=round(exp(-0.5*o$delta)/sum(exp(-0.5*o$delta)),4)
						o$ER=max(o$prob)/o$prob
						o$AIC=round(o$AIC,1)
						o$delta=round(o$delta,2)
						o$prob=round(o$prob,3)
						o$ER=round(o$ER,2)
						o$ER_pd=NA
						o_1=o
						#o[order(o$delta),]
				}
			}				
			{# gap
				table(gapp$bird_ID_filled, gapp$exper)
				length(unique(gapp$bird_ID_filled))
				length(unique(gapp$bird_ID))
				{# 1 - model of interest
				 m1=lmer(bout_length~exper*bout_start_j_c+(bout_start_j_c|bird_ID),gapp, REML=FALSE)
				 pred=c('Intercept (before)','Period (after)','Day (centered)','Period:Day')
							nsim <- 2000
							bsim <- sim(m1, n.sim=nsim)  
					# Fixed effects
						v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))
						ci=apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
						oi=data.frame(model='1',response='gap_min',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,], stringsAsFactors=FALSE)
						rownames(oi) = NULL
								oi$estimate_r=round(oi$estimate,2)
								oi$lwr_r=round(oi$lwr,2)
								oi$upr_r=round(oi$upr,2)
							#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
						oii=oi[c('model','response','type',"effect", "estimate_r","lwr_r",'upr_r')]	
					# Random effects
						l=data.frame(summary(m1)$varcor)
						l=l[is.na(l$var2),]
						l$var2=NULL
						ri=data.frame(model='1',response='gap_min',type='random (var)',effect=l$var, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA, stringsAsFactors=FALSE)
							ri$effect[nrow(ri)]='residual'
						o1=rbind(oii,ri)
				
							
				}	
				{# 2 - simple model 
					 m=lmer(bout_length~exper+(1|bird_ID),gapp, REML=FALSE)
					 pred=c('Intercept (before)','Period (after)')
								nsim <- 2000
								bsim <- sim(m, n.sim=nsim)  
						# Fixed effects
							v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))
							ci=apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							oi=data.frame(model='2',response='gap_min',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
							rownames(oi) = NULL
								oi$estimate_r=round(oi$estimate,2)
								oi$lwr_r=round(oi$lwr,2)
								oi$upr_r=round(oi$upr,2)
								#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
							oii=oi[c('model','response','type',"effect", "estimate_r","lwr_r",'upr_r')]	
						# Random effects
							l=data.frame(summary(m)$varcor)
								ri=data.frame(model='2',response='gap_min', type='random (var)',effect=l$grp, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA)
							o2=rbind(oii,ri)
				}
				{# 3 - model with bird (focal/removed)
				 m3=lmer(bout_length~exper*treated_bird*bout_start_j_c+(bout_start_j_c|bird_ID),gapp, REML=FALSE)
				 pred=c('Intercept (before,removed)','Period (after)','Bird(focal)','Day','Period:bird', 'Period:Day','bird:Day','period:bird:day')
							nsim <- 2000
							bsim <- sim(m3, n.sim=nsim)  
					# Fixed effects
						v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))
						ci=apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
						oi=data.frame(model='3',response='gap_min',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,], stringsAsFactors=FALSE)
						rownames(oi) = NULL
								oi$estimate_r=round(oi$estimate,2)
								oi$lwr_r=round(oi$lwr,2)
								oi$upr_r=round(oi$upr,2)
							#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
						oii=oi[c('model','response','type',"effect", "estimate_r","lwr_r",'upr_r')]	
					# Random effects
						l=data.frame(summary(m3)$varcor)
						l=l[is.na(l$var2),]
						l$var2=NULL
						ri=data.frame(model='1',response='gap_min',type='random (var)',effect=l$var, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA, stringsAsFactors=FALSE)
							ri$effect[nrow(ri)]='residual'
						o3=rbind(oii,ri)

				}				
				   o_gap=(rbind(o1,o2,o3))
				{# AICc # number of observations = 68 number of birds*2 (2 = before and after period); not all birds had gap in both length(unique(paste(gapp$bird_ID,gapp$exper)))
							# length(unique(paste(gapp$bird_ID[gapp$exper=='b'])))
							# length(unique(paste(gapp$bird_ID[gapp$exper=='a'])))
							b_i=unique(gapp$bird_ID[gapp$exper=='b'])
							a_i=unique(gapp$bird_ID[gapp$exper=='a'])
							b_i[!b_i%in%a_i] # s519 f not in after
							a_i[!a_i%in%b_i] # s622 f not in before
							length(unique(c(b_i,a_i)))
							
						o=data.frame(response='gap',model=c('period:day','period','period:bird:day'), AIC=c(AICc(m1, nobs=68),AICc(m,nobs=68),AICc(m3,nobs=68)))
						o$delta=o$AIC-min(o$AIC)
						o$prob=exp(-0.5*o$delta)/sum(exp(-0.5*o$delta))
						o$ER=max(o$prob)/o$prob
						o$ER_pd=o$prob[o$model=='period:day']/o$prob # comparison of period and period:day model
						o$AIC=round(o$AIC,1)
						o$delta=round(o$delta,2)
						o$prob=round(o$prob,3)
						o$ER=round(o$ER,2)
						o_2=o #o[order(o$delta),]
				}
			}				
			{# combine and export to excel table
							sname = tempfile(fileext='.xls')
							wb = loadWorkbook(sname,create = TRUE)	
							createSheet(wb, name = "output")
							writeWorksheet(wb, rbind(o_gapyn,o_gap), sheet = "output")
							createSheet(wb, name = "output_AIC")
							writeWorksheet(wb, rbind(o_1,o_2), sheet = "output_AIC")
							saveWorkbook(wb)
							shell(sname)
							
				}
			
			{# model assumptions
				{# gap presence
					{# 1
						dev.new(width=6,height=9)
									
						m=glmer(present~exper*bout_start_j_c+(bout_start_j_c|bird_ID),gap, family='binomial',control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
						
						par(mfrow=c(4,3))
						
						scatter.smooth(fitted(m),resid(m),col='red');abline(h=0, lty=2)
						scatter.smooth(fitted(m),sqrt(abs(resid(m))), col='red')
													
						qqnorm(resid(m), main=list("Normal Q-Q Plot: residuals", cex=0.8),col='red') 
						qqline(resid(m))
			 
						qqnorm(unlist(ranef(m)$bird_ID[1]), main = " Bird")
						qqline(unlist(ranef(m)$bird_ID[1]))
						
						qqnorm(unlist(ranef(m)$bird_ID[2]), main = " Slope")
						qqline(unlist(ranef(m)$bird_ID[2]))
						
						scatter.smooth(resid(m)~inc$bout_start_j_c);abline(h=0, lty=2, col='red')
						scatter.smooth(resid(m)~inc$exper);abline(h=0, lty=2, col='red')
						plot(resid(m)~factor(inc$exper)); abline(h=0, lty=2, col='red')
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=inc$lon, y=inc$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
					}
					{# 2
						dev.new(width=6,height=9)
									
						m=glmer(present~exper+(1|bird_ID),gap, family='binomial')
						par(mfrow=c(4,3))
						
						scatter.smooth(fitted(m),resid(m),col='red');abline(h=0, lty=2)
						scatter.smooth(fitted(m),sqrt(abs(resid(m))), col='red')
													
						qqnorm(resid(m), main=list("Normal Q-Q Plot: residuals", cex=0.8),col='red') 
						qqline(resid(m))
			 
						qqnorm(unlist(ranef(m)$bird_ID[1]), main = " Bird")
						qqline(unlist(ranef(m)$bird_ID[1]))
						
						scatter.smooth(resid(m)~inc$exper);abline(h=0, lty=2, col='red')
						plot(resid(m)~factor(inc$exper)); abline(h=0, lty=2, col='red')
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=inc$lon, y=inc$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
					}
					{# 3
						dev.new(width=6,height=9)
									
						m=glmer(present~exper*treated_bird*bout_start_j_c+(bout_start_j_c|bird_ID),gap, family='binomial',control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
						par(mfrow=c(4,4))
						
						scatter.smooth(fitted(m),resid(m),col='red');abline(h=0, lty=2)
						scatter.smooth(fitted(m),sqrt(abs(resid(m))), col='red')
													
						qqnorm(resid(m), main=list("Normal Q-Q Plot: residuals", cex=0.8),col='red') 
						qqline(resid(m))
			 
						qqnorm(unlist(ranef(m)$bird_ID[1]), main = " Bird")
						qqline(unlist(ranef(m)$bird_ID[1]))
						
						qqnorm(unlist(ranef(m)$bird_ID[2]), main = " Slope")
						qqline(unlist(ranef(m)$bird_ID[2]))
						
						scatter.smooth(resid(m)~inc$bout_start_j_c);abline(h=0, lty=2, col='red')
						scatter.smooth(resid(m)~inc$exper);abline(h=0, lty=2, col='red')
						plot(resid(m)~factor(inc$exper)); abline(h=0, lty=2, col='red')
						plot(resid(m)~factor(inc$treated_bird)); abline(h=0, lty=2, col='red')
						plot(resid(m)~interaction(inc$exper,inc$treated_bird)); abline(h=0, lty=2, col='red')
						
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=inc$lon, y=inc$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
					}
				}	
				{# gap
					{# 1
						dev.new(width=6,height=9)
									
						m=lmer(bout_length~exper*bout_start_j_c+(bout_start_j_c|bird_ID),gap, REML=FALSE)
						par(mfrow=c(4,3))
						
						scatter.smooth(fitted(m),resid(m),col='red');abline(h=0, lty=2)
						scatter.smooth(fitted(m),sqrt(abs(resid(m))), col='red')
													
						qqnorm(resid(m), main=list("Normal Q-Q Plot: residuals", cex=0.8),col='red') 
						qqline(resid(m))
			 
						qqnorm(unlist(ranef(m)$bird_ID[1]), main = " Bird")
						qqline(unlist(ranef(m)$bird_ID[1]))
						
						qqnorm(unlist(ranef(m)$bird_ID[2]), main = " Slope")
						qqline(unlist(ranef(m)$bird_ID[2]))
						
						scatter.smooth(resid(m)~inc$bout_start_j_c);abline(h=0, lty=2, col='red')
						scatter.smooth(resid(m)~inc$exper);abline(h=0, lty=2, col='red')
						plot(resid(m)~factor(inc$exper)); abline(h=0, lty=2, col='red')
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=inc$lon, y=inc$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
					}
					{# 2
						dev.new(width=6,height=9)
									
						m=lmer(bout_length~exper+(1|bird_ID),gap, REML=FALSE)
						par(mfrow=c(4,3))
						
						scatter.smooth(fitted(m),resid(m),col='red');abline(h=0, lty=2)
						scatter.smooth(fitted(m),sqrt(abs(resid(m))), col='red')
													
						qqnorm(resid(m), main=list("Normal Q-Q Plot: residuals", cex=0.8),col='red') 
						qqline(resid(m))
			 
						qqnorm(unlist(ranef(m)$bird_ID[1]), main = " Bird")
						qqline(unlist(ranef(m)$bird_ID[1]))
						
						scatter.smooth(resid(m)~inc$exper);abline(h=0, lty=2, col='red')
						plot(resid(m)~factor(inc$exper)); abline(h=0, lty=2, col='red')
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=inc$lon, y=inc$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
					}
					{# 3
						dev.new(width=6,height=9)
									
						m=lmer(bout_length~exper*treated_bird*bout_start_j_c+(bout_start_j_c|bird_ID),gap, REML=FALSE)
						par(mfrow=c(4,4))
						
						scatter.smooth(fitted(m),resid(m),col='red');abline(h=0, lty=2)
						scatter.smooth(fitted(m),sqrt(abs(resid(m))), col='red')
													
						qqnorm(resid(m), main=list("Normal Q-Q Plot: residuals", cex=0.8),col='red') 
						qqline(resid(m))
			 
						qqnorm(unlist(ranef(m)$bird_ID[1]), main = " Bird")
						qqline(unlist(ranef(m)$bird_ID[1]))
						
						qqnorm(unlist(ranef(m)$bird_ID[2]), main = " Slope")
						qqline(unlist(ranef(m)$bird_ID[2]))
						
						scatter.smooth(resid(m)~inc$bout_start_j_c);abline(h=0, lty=2, col='red')
						scatter.smooth(resid(m)~inc$exper);abline(h=0, lty=2, col='red')
						plot(resid(m)~factor(inc$exper)); abline(h=0, lty=2, col='red')
						plot(resid(m)~factor(inc$treated_bird)); abline(h=0, lty=2, col='red')
						plot(resid(m)~interaction(inc$exper,inc$treated_bird)); abline(h=0, lty=2, col='red')
						
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=inc$lon, y=inc$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
					}
				}
			}					
		}	
		
		{# Figure 5
			{# run first - prepare predictions
				{# constancy
					{# simple model (prediction for day 	
							m=lmer(inc_eff~exper+(1|bird_ID),inc, REML=FALSE)
							#summary(m)
							#plot(allEffects(m))
						
							# simulation		
								nsim <- 2000
								bsim <- sim(m, n.sim=nsim)  
								apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							
							# coefficients
								v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))	
							# predicted values		
								newD=data.frame(exper=c('a','b'))
								
							# exactly the model which was used has to be specified here 
								X <- model.matrix(~ exper,data=newD)	
											
							# calculate predicted values and creditability intervals
								newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
										predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
										for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@fixef[i,]
										newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
										newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
										#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
										#newD=newD[order(newD$t_tundra),]
								pcs=newD
								pcsB=pcs[pcs$exper=='b',]
								pcsA=pcs[pcs$exper=='a',]
					}			
					{# day model	
						m=lmer(inc_eff~exper*bout_start_j_c+(bout_start_j_c|bird_ID),inc, REML=FALSE)
							#summary(m)
							#plot(allEffects(m))
						
							# simulation		
								nsim <- 2000
								bsim <- sim(m, n.sim=nsim)  
								apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							
							# coefficients
								v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))	
							# predicted values		
								newDb=data.frame(bout_start_j_c=seq(min(inc$bout_start_j_c),0,length.out=100),
												exper='b')
								newDa=data.frame(bout_start_j_c=seq(0,max(inc$bout_start_j_c),length.out=100),
												exper='a')				
								newD=rbind(newDb,newDa)			
							# exactly the model which was used has to be specified here 
								X <- model.matrix(~ exper*bout_start_j_c,data=newD)	
											
							# calculate predicted values and creditability intervals
								newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
										predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
										for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@fixef[i,]
										newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
										newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
										#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
										#newD=newD[order(newD$t_tundra),]
								pcon=newD
								pconB=pcon[pcon$exper=='b',]
								pconA=pcon[pcon$exper=='a',]
					}			
				}
				{# bout
					{# simple model
						m=lmer(bout_length~exper+(1|bird_ID),inc, REML=FALSE)
							#summary(m)
							#plot(allEffects(m))
						
							# simulation		
								nsim <- 2000
								bsim <- sim(m, n.sim=nsim)  
								apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							
							# coefficients
								v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))	
							# predicted values		
								newD=data.frame(exper=c('a','b'))
								
							# exactly the model which was used has to be specified here 
								X <- model.matrix(~ exper,data=newD)	
											
							# calculate predicted values and creditability intervals
								newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
										predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
										for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@fixef[i,]
										newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
										newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
										#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
										#newD=newD[order(newD$t_tundra),]
								pbs=newD
								pbsB=pbs[pbs$exper=='b',]
								pbsA=pbs[pbs$exper=='a',]
							
					}
					{# day model
						m=lmer(bout_length~exper*bout_start_j_c+(bout_start_j_c|bird_ID),inc, REML=FALSE)
							#summary(m)
							#plot(allEffects(m))
						
							# simulation		
								nsim <- 2000
								bsim <- sim(m, n.sim=nsim)  
								apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							
							# coefficients
								v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))	
							# predicted values		
								newDb=data.frame(bout_start_j_c=seq(min(inc$bout_start_j_c),0,length.out=100),
												exper='b')
								newDa=data.frame(bout_start_j_c=seq(0,max(inc$bout_start_j_c),length.out=100),
												exper='a')				
								newD=rbind(newDb,newDa)			
							# exactly the model which was used has to be specified here 
								X <- model.matrix(~ exper*bout_start_j_c,data=newD)	
											
							# calculate predicted values and creditability intervals
								newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
										predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
										for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@fixef[i,]
										newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
										newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
										#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
										#newD=newD[order(newD$t_tundra),]
								pbout=newD
								pboutB=pbout[pbout$exper=='b',]
								pboutA=pbout[pbout$exper=='a',]
					}
				
				}
				{# gap length
					{# simple model	
						m=lmer(log(bout_length)~exper+(1|bird_ID),gapp, REML=FALSE)
							#summary(m)
							#plot(allEffects(m))
						
							# simulation		
								nsim <- 2000
								bsim <- sim(m, n.sim=nsim)  
								apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							
							# coefficients
								v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))	
							# predicted values		
								newD=data.frame(exper=c('a','b'))
								
							# exactly the model which was used has to be specified here 
								X <- model.matrix(~ exper,data=newD)	
											
							# calculate predicted values and creditability intervals
								newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
										predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
										for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@fixef[i,]
										newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
										newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
										#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
										#newD=newD[order(newD$t_tundra),]
								pgs=newD
								pgsB=pgs[pgs$exper=='b',]
								pgsA=pgs[pgs$exper=='a',]
					}			
					{# day model	
						m=lmer(log(bout_length)~exper*bout_start_j_c+(bout_start_j_c|bird_ID),gapp, REML=FALSE)
							#summary(m)
							#plot(allEffects(m))
						
							# simulation		
								nsim <- 2000
								bsim <- sim(m, n.sim=nsim)  
								apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							
							# coefficients
								v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))	
							# predicted values		
								newDb=data.frame(bout_start_j_c=seq(min(gapp$bout_start_j_c),0,length.out=100),
												exper='b')
								newDa=data.frame(bout_start_j_c=seq(0,max(gapp$bout_start_j_c),length.out=100),
												exper='a')				
								newD=rbind(newDb,newDa)			
							# exactly the model which was used has to be specified here 
								X <- model.matrix(~ exper*bout_start_j_c,data=newD)	
											
							# calculate predicted values and creditability intervals
								newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
										predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
										for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@fixef[i,]
										newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
										newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
										#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
										#newD=newD[order(newD$t_tundra),]
								pgap=newD
								pgapB=pgap[pgap$exper=='b',]
								#pgapB$pred=ifelse(pgapB$pred<=0,log(0.015),log(as.numeric(pgapB$pred)))
								#pgapB$lwr=ifelse(pgapB$lwr>0,log(pgapB$lwr),log(0.015))
								#pgapB$upr=ifelse(pgapB$upr>0,log(pgapB$upr),log(0.015))
								pgapA=pgap[pgap$exper=='a',]
								#pgapA$pred=ifelse(pgapA$pred>0,log(pgapA$pred),log(0.015))
								#pgapA$lwr=log(0.015)#log(pgapA$lwr)
								#pgapA$upr=log(pgapA$upr)
					}			
				}
			}
			{# run first - prepare aggreated datasets
				{# incubation
					x=ddply(inc,.(nest,sex,exper, treated_bird), summarise, med_bout=median(bout_length),mean_bout=mean(bout_length), med_eff=median(inc_eff),mean_eff=mean(inc_eff))
								ggplot(x,aes(x=med_eff,y=mean_eff))+geom_point()
								plot(med_eff~mean_eff,x, xlim=c(0,1), ylim=c(0,1))
								abline(a=0,b=1)
								
								x[x$med_eff<0.2,]
								densityplot(~x$med_eff)
					
					xb=x[x$exper=='b',]
					xa=x[x$exper=='a',]
					xb$med_eff_after=xa$med_eff[match(paste(xb$nest,xb$sex),paste(xa$nest,xa$sex))]
					xb$med_bout_after=xa$med_bout[match(paste(xb$nest,xb$sex),paste(xa$nest,xa$sex))]
					xb$col_=ifelse(xb$treated_bird=="y", "#FCB42C","#535F7C")
				}
				{# gap
					y=ddply(gapp,.(nest,sex,exper, treated_bird), summarise, med_gap=median(bout_length), num_gaps=sum(present))
				
					yb=y[y$exper=='b',] s622 s527 
					ya=y[y$exper=='a',] s527 s519 
					yb$med_gap_after=ya$med_gap[match(paste(yb$nest,yb$sex),paste(ya$nest,ya$sex))]
					yb$col_=ifelse(yb$treated_bird=="y", "#FCB42C","#535F7C")
					yb = yb[!is.na(yb$med_gap_after),]
				
			}
			}
			{# run first - prepare for boxplots
				inc$treated_sex_exper=as.factor(ifelse(inc$treated_bird=='y', ifelse(inc$sex=='f', ifelse(inc$exper=='b', 1, 2), ifelse(inc$exper=='b', 3, 4)), ifelse(inc$sex=='f', ifelse(inc$exper=='b', 5, 6), ifelse(inc$exper=='b', 7, 8))))
				
				gapp$treated_sex_exper=as.factor(ifelse(gapp$treated_bird=='y', ifelse(gapp$sex=='f', ifelse(gapp$exper=='b', 1, 2), ifelse(gapp$exper=='b', 3, 4)), ifelse(gapp$sex=='f', ifelse(gapp$exper=='b', 5, 6), ifelse(gapp$exper=='b', 7, 8))))
			}
		
			{# plot a,b,c, even less labels
				  dev.new(width=3.5,height=1.85*2.15)
				  #png(paste(out2,"Figure_5_.png", sep=""), width=3.5,height=1.85*2.15,units="in",res=600)
				   par(mfrow=c(3,3),mar=c(1,0,0,0.5),oma = c(2, 1.8, 0.7, 0.1),ps=12, mgp=c(1.2,0.35,0), las=1, cex=1, col.axis="grey30",font.main = 1, col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt
				{# constancy
					{# medians
					par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
					plot(xb$med_eff_after~xb$med_eff ,pch=19,xlim=c(0.6,1), ylim=c(0.6,1), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
						
						axis(1, at=c(0.7,0.9),labels=FALSE,col.ticks="grey90")						
						axis(1, at=seq(0.6,1,by=0.2),labels=c('0.6','0.8','1.0'),cex.axis=0.5,mgp=c(0,-0.20,0))
							mtext('Before experiment',side=1,line=-1, cex=0.5, las=1, col='grey45')
							#mtext('Nest attendance',side=1,line=0.3, cex=0.6, las=1, col='grey30')
							
						axis(2, at=c(0.7,0.9),labels=FALSE,col.ticks="grey90")	
						axis(2, at=seq(0.6,1,by=0.2),labels=c('0.6','0.8','1.0'))
							mtext('Post experiment',side=2,line=-0.6, cex=0.5, las=3, col='grey45')
							mtext('Nest attendance',side=2,line=1.1, cex=0.6, las=3, col='grey30')
						
					lines(c(0.6,1),c(0.6,1), lty=3, col="lightgrey")
						#mtext(expression(bold("a")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
						
					# data
						points(xb$med_eff_after~xb$med_eff, col=xb$col_,bg=adjustcolor(xb$col_, alpha.f = 0.4), pch=21, cex=0.5)	
						
						text(x=0.65,y=1, labels='Focal', col='#FCB42C', cex=0.5,adj=0)
						text(x=0.65,y=0.96, labels='Removed', col='#535F7C', cex=0.5,adj=0)
						
						mtext(expression(bold("a")),side=3,line=0, cex=0.7, las=1, col="grey30")
					# predictions	
						# 95%CI for before
						arrows(y0=pcsB$lwr, x0=pcsA$pred,y1=pcsB$upr, x1=pcsA$pred,
							code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)
						# 95%CI for after
						arrows(y0=pcsB$pred, x0=pcsA$lwr,y1=pcsB$pred, x1=pcsA$upr,
							code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)	
						
						points(x=pcsA$pred,y=pcsB$pred, pch=20, cex=0.5,col="red")
					}
					{# boxplot - no points, color box outline - before after
						par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
					
						boxplot(inc_eff ~ treated_sex_exper, data = inc,
											ylab = NULL,xaxt='n', yaxt='n',
											ylim=c(0.58,1),								
											#at=c(1,2,3.5,4.5),
											at=c(1,2,3.5,4.5,7,8,9.5,10.5),
											type='n',
											outcex=0.3,outpch=20,boxwex=0.6,whisklty=1,staplelty=0,#medlwd=1, 
											lwd = 1,
											border=c('#FCB42C','#535F7C','#FCB42C','#535F7C','#FCB42C','#535F7C','#FCB42C','#535F7C'),
											col = adjustcolor("white", alpha.f = 0), # trick for PNGs, to show what is underneath the boxplot else can be taken out
											#outcol="darkgrey",boxcol='darkgrey',whiskcol='darkgrey',staplecol='darkgrey',medcol='darkgrey', 
											par(bty='l')	
											#add=TRUE
											)					
						axis(1, at=c(2.75,8.75),labels=c('Focal','Removed'),cex.axis=0.5,mgp=c(0,-0.20,0))
						text(c(1.5,4,7.5,10), (par("usr")[3])*1.05326705, labels = c('\u2640','\u2642','\u2640','\u2642'), font=4, xpd = TRUE, cex=0.5, col='grey45')
						#mtext('Parent',side=1,line=0.3, cex=0.6, las=1, col='grey30')
						
						axis(2, at=c(0.7,0.9),labels=FALSE,col.ticks="grey90")	
						axis(2, at=seq(0.6,1,by=0.2), labels=FALSE)
						
						text(x=5.75,y=0.71, labels='Before', col='#FCB42C', cex=0.5)
						text(x=5.75,y=0.67, labels='After', col='#535F7C', cex=0.5)
						
						#text(x=10.5,y=0.58, labels='*', col=col_2011, cex=0.7) #labels='0, 0.3, 0.4, 0.6'
						#text(x=10.5,y=0.58, labels='*', col='grey30', cex=0.6) #labels='0, 0.3, 0.4, 0.6'
						points(x=10.5,y=0.58, pch=20 ,col=col_2011_light, cex=0.5) #labels='0, 0.3, 0.4, 0.6'
						#outliers 0.0129349,0.2786890,0.3592990, 0.5777780
						
						mtext(expression(bold("b")),side=3,line=0, cex=0.7, las=1, col="grey30")
						}
					{# day
						par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
						plot(pcon$pred~pcon$bout_start_j_c ,pch=19,xlim=c(-4.25,4), ylim=c(0.6,1), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
											
							axis(1, at=seq(-4,4,by=2),labels=seq(-4,4,by=2),cex.axis=0.5,mgp=c(0,-0.20,0))
								#mtext('Day',side=1,line=0.3, cex=0.6, las=1, col='grey30')
							
							axis(2, at=c(0.7,0.9),labels=FALSE,col.ticks="grey90")		
							axis(2, at=seq(0.6,1,by=0.2),labels=FALSE)
							lines(c(0,0),c(0.58,1), lty=3, col="red")							
						# data
							points(inc$inc_eff[inc$inc_eff>0.58]~inc$bout_start_j_c[inc$inc_eff>0.58], col=adjustcolor(inc$col_[inc$inc_eff>0.58], alpha.f = 0.5), pch=20, cex=0.2)	
							#points(inc$inc_eff~inc$bout_start_j_c, col=inc$col_,bg=adjustcolor(inc$col_, alpha.f = 0.4), pch=21, cex=0.5)	
							
						# predictions
							# before
							polygon(c(pconB$bout_start_j_c, rev(pconB$bout_start_j_c)), c(pconB$lwr, 
								rev(pconB$upr)), border=NA, col=adjustcolor(col_t ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pconB$bout_start_j_c, pconB$pred, col=col_t,lwd=1)
							
							# before
							polygon(c(pconA$bout_start_j_c, rev(pconA$bout_start_j_c)), c(pconA$lwr, 
								rev(pconA$upr)), border=NA, col=adjustcolor(col_c ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pconA$bout_start_j_c, pconA$pred, col=col_c,lwd=1)
												
						# text	
							mtext(expression(bold("c")),side=3,line=0, cex=0.7, las=1, col="grey30")
							#text(x=inc$bout_start_j_c[inc$inc_eff<0.6],y=0.6, labels='*', col=col_2011_light, cex=0.6) #labels='0, 0.3, 0.4, 0.6' #round(inc$inc_eff[inc$inc_eff<0.6],2 )
							points(x=inc$bout_start_j_c[inc$inc_eff<0.6],y=rep(0.6,length(inc$bout_start_j_c[inc$inc_eff<0.6])), pch=20, col=col_2011_light, cex=0.5) 
							text(x=1.20,y=0.66, labels='0.28', col='grey30', cex=0.4)
							text(x=1.20,y=0.63, labels='0.01', col='grey30', cex=0.4)
							
							text(x=2.75,y=0.63, labels='0.36', col='grey30', cex=0.4)
							text(x=3.5,y=0.66, labels='0.58', col='grey30', cex=0.4)
							
							text(x=-2,y=0.725, labels='Before', col='#FCB42C', cex=0.5)
							text(x=2,y=0.725, labels='After', col='#535F7C', cex=0.5)
					}
				}
				{# bout
					{# medians
					par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
					plot(xb$med_bout_after~xb$med_bout ,pch=19,xlim=c(0,16.5), ylim=c(0,16.5), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
						
						axis(1, at=c(4,12),labels=FALSE,col.ticks="grey90")							
						axis(1, at=seq(0,16,by=8),labels=c('0','8','16'),cex.axis=0.5,mgp=c(0,-0.20,0))
							mtext('Before experiment',side=1,line=-1, cex=0.5, las=1, col='grey45')
							#mtext('Bout [h]',side=1,line=0.3, cex=0.6, las=1, col='grey30')
						
						axis(2, at=c(4,12),labels=FALSE,col.ticks="grey90")		
						axis(2, at=seq(0,16,by=8),labels=c('0','8','16'))
							mtext('Post experiment',side=2,line=-0.6, cex=0.5, las=3, col='grey45')
							mtext('Bout [h]',side=2,line=1.1, cex=0.6, las=3, col='grey30')
						
					lines(c(0,16),c(0,16), lty=3, col="lightgrey")
						#mtext(expression(bold("a")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
						
					# data
						points(xb$med_bout_after~xb$med_bout, col=xb$col_,bg=adjustcolor(xb$col_, alpha.f = 0.4), pch=21, cex=0.5)	
						
						#text(x=3,y=1, labels='Focal', col='#FCB42C', cex=0.5,adj=0)
						#text(x=0.3,y=0.96, labels='Removed', col='#535F7C', cex=0.5,adj=0)
					
					# predictions	
						# 95%CI for before
						arrows(y0=pbsB$lwr, x0=pbsA$pred,y1=pbsB$upr, x1=pbsA$pred,
							code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)
						# 95%CI for after
						arrows(y0=pbsB$pred, x0=pbsA$lwr,y1=pbsB$pred, x1=pbsA$upr,
							code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)	
						
						points(x=pbsA$pred,y=pbsB$pred, pch=20, cex=0.5,col="red")
					}
					{# boxplot - no points, color box outline - before after
						par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
					
						boxplot(bout_length ~ treated_sex_exper, data = inc,
											ylab = NULL,xaxt='n', yaxt='n',
											ylim=c(0,16.5),								
											#at=c(1,2,3.5,4.5),
											at=c(1,2,3.5,4.5,7,8,9.5,10.5),
											type='n',
											outcex=0.3,outpch=20,boxwex=0.6,whisklty=1,staplelty=0,#medlwd=1, 
											lwd = 1,
											border=c('#FCB42C','#535F7C','#FCB42C','#535F7C','#FCB42C','#535F7C','#FCB42C','#535F7C'),
											col = adjustcolor("white", alpha.f = 0), # trick for PNGs, to show what is underneath the boxplot else can be taken out
											#outcol="darkgrey",boxcol='darkgrey',whiskcol='darkgrey',staplecol='darkgrey',medcol='darkgrey', 
											par(bty='l')	
											#add=TRUE
											)					
						axis(1, at=c(2.75,8.75),labels=c('Focal','Removed'),cex.axis=0.5,mgp=c(0,-0.20,0))
						text(c(1.5,4,7.5,10), par("usr")[3]+1.3, labels = c('\u2640','\u2642','\u2640','\u2642'), font=4, xpd = TRUE, cex=0.5, col='grey45')
						#mtext('Parent',side=1,line=0.3, cex=0.6, las=1, col='grey30')
						
						axis(2, at=c(4,12),labels=FALSE,col.ticks="grey90")		
						axis(2, at=seq(0,16,by=8), labels=FALSE)
						
						#text(x=5.75,y=0.71, labels='Before', col='#FCB42C', cex=0.5)
						#text(x=5.75,y=0.67, labels='After', col='#535F7C', cex=0.5)
						
						#text(x=10.9,y=0.58, labels='*', col='grey30', cex=0.6) #labels='0, 0.3, 0.4, 0.6'
						
						#outliers 0.0129349,0.2786890,0.3592990, 0.5777780
						
						
						}
					{# day
						par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
						plot(pbout$pred~pbout$bout_start_j_c ,pch=19,xlim=c(-4.25,4), ylim=c(0,16.5), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
											
							axis(1, at=seq(-4,4,by=2),labels=seq(-4,4,by=2),cex.axis=0.5,mgp=c(0,-0.20,0))
								#mtext('Day',side=1,line=0.3, cex=0.6, las=1, col='grey30')
								
							axis(2, at=c(4,12),labels=FALSE,col.ticks="grey90")		
							axis(2, at=seq(0,16,by=8), labels=FALSE)
							lines(c(0,0),c(0,16.5), lty=3, col="red")							
						# data
							points(inc$bout_length~inc$bout_start_j_c, col=adjustcolor(inc$col_, alpha.f = 0.5), pch=20, cex=0.2)	
							#points(inc$inc_eff~inc$bout_start_j_c, col=inc$col_,bg=adjustcolor(inc$col_, alpha.f = 0.4), pch=21, cex=0.5)	
						
						
						# predictions
							# before
							polygon(c(pboutB$bout_start_j_c, rev(pboutB$bout_start_j_c)), c(pboutB$lwr, 
								rev(pboutB$upr)), border=NA, col=adjustcolor(col_t ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pboutB$bout_start_j_c, pboutB$pred, col=col_t,lwd=1)
							
							# before
							polygon(c(pboutA$bout_start_j_c, rev(pboutA$bout_start_j_c)), c(pboutA$lwr, 
								rev(pboutA$upr)), border=NA, col=adjustcolor(col_c ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pboutA$bout_start_j_c, pboutA$pred, col=col_c,lwd=1)
													
					
							#text(x=-2,y=0.725, labels='Before', col='#FCB42C', cex=0.5)
							#text(x=2,y=0.725, labels='After', col='#535F7C', cex=0.5)
					}
				}
				{# gap
					{# medians
					par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
					plot(yb$med_gap_after~yb$med_gap ,pch=19,xlim=log(c(0.015,700)), ylim=log(c(0.015,700)), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
										
						axis(1, at=log(c(seq(0.1,0.9,by=0.1),seq(1,10,by=1), seq(20,100, by=10), seq(200,700, by=100))),labels=FALSE,col.ticks="grey90")
							axis(1,  at=log(c(0.1,1,10,100,700)),labels=c('0.1','1','10','100','700'),mgp=c(0,-0.20,0))
							mtext('Before experiment',side=1,line=-1, cex=0.5, las=1, col='grey45')
							#mtext('Exchange gap [min]',side=1,line=0.3, cex=0.6, las=1, col='grey30')
							
						axis(2,  at=log(c(seq(0.1,0.9,by=0.1),seq(1,10,by=1), seq(20,100, by=10), seq(200,700, by=100))),labels=FALSE,col.ticks="grey90")
						axis(2,  at=log(c(0.1,1,10,100,700)),labels=c('0.1','1','10','100','700'))
							mtext('Post experiment',side=2,line=-0.6, cex=0.5, las=3, col='grey45')
							mtext('Exchange gap [min]',side=2,line=1.1, cex=0.6, las=3, col='grey30')
						
					lines(log(c(0.001,700)),log(c(0.001,700)), lty=3, col="lightgrey")
						#mtext(expression(bold("a")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
						
					# data
						points(log(yb$med_gap_after)~log(yb$med_gap), col=yb$col_,bg=adjustcolor(yb$col_, alpha.f = 0.4), pch=21, cex=0.5)	
						
						#text(x=3,y=1, labels='Focal', col='#FCB42C', cex=0.5,adj=0)
						#text(x=0.3,y=0.96, labels='Removed', col='#535F7C', cex=0.5,adj=0)
					
					# predictions	
						# 95%CI for before
						arrows(x0=pgsA$pred,y0=pgsA$pred,x1=pgsB$upr, y1=pgsA$pred,
							code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)
						# 95%CI for after
						arrows(x0=pgsB$pred, y0=pgsA$lwr,x1=pgsB$pred, y1=pgsA$upr,
							code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)	
						
						points(y=pgsA$pred,x=pgsB$pred, pch=20, cex=0.5,col="red")
					}
					{# boxplot - no points, color box outline - before after
						par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
					
						boxplot(log(bout_length) ~ treated_sex_exper, data = gapp,
											ylab = NULL,xaxt='n', yaxt='n',
											ylim=log(c(0.015,700)),								
											#at=c(1,2,3.5,4.5),
											at=c(1,2,3.5,4.5,7,8,9.5,10.5),
											type='n',
											outcex=0.3,outpch=20,boxwex=0.6,whisklty=1,staplelty=0,#medlwd=1, 
											lwd = 1,
											border=c('#FCB42C','#535F7C','#FCB42C','#535F7C','#FCB42C','#535F7C','#FCB42C','#535F7C'),
											col = adjustcolor("white", alpha.f = 0), # trick for PNGs, to show what is underneath the boxplot else can be taken out
											#outcol="darkgrey",boxcol='darkgrey',whiskcol='darkgrey',staplecol='darkgrey',medcol='darkgrey', 
											par(bty='l')	
											#add=TRUE
											)					
						axis(1, at=c(2.75,8.75),labels=c('Focal','Removed'),cex.axis=0.5,mgp=c(0,-0.20,0))
						text(c(1.5,4,7.5,10), par("usr")[3]+log(2.3), labels = c('\u2640','\u2642','\u2640','\u2642'), font=4, xpd = TRUE, cex=0.5, col='grey45')
						mtext('Parent',side=1,line=0.5, cex=0.6, las=1, col='grey30')
						
						axis(2,  at=log(c(seq(0.1,0.9,by=0.1),seq(1,10,by=1), seq(20,100, by=10), seq(200,700, by=100))),labels=FALSE,col.ticks="grey90")
						axis(2,  at=log(c(0.1,1,10,100,700)),labels=FALSE)
						
						#text(x=5.75,y=0.71, labels='Before', col='#FCB42C', cex=0.5)
						#text(x=5.75,y=0.67, labels='After', col='#535F7C', cex=0.5)
						
						#text(x=10.9,y=0.58, labels='*', col='grey30', cex=0.6) #labels='0, 0.3, 0.4, 0.6'
						
						#outliers 0.0129349,0.2786890,0.3592990, 0.5777780
						
						}
					{# day
						par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
						plot(pgap$pred~pgap$bout_start_j_c ,pch=19,xlim=c(-4.25,4), ylim=log(c(0.015,700)),	 xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
											
							axis(1, at=seq(-4,4,by=2),labels=seq(-4,4,by=2),cex.axis=0.5,mgp=c(0,-0.20,0))
								mtext('Day',side=1,line=0.5, cex=0.6, las=1, col='grey30')
								
							axis(2,  at=log(c(seq(0.1,0.9,by=0.1),seq(1,10,by=1), seq(20,100, by=10), seq(200,700, by=100))),labels=FALSE,col.ticks="grey90")
							axis(2,  at=log(c(0.1,1,10,100,700)),labels=FALSE)			

							lines(c(0,0),log(c(0.015,700)), lty=3, col="red")								
						# data
							points(log(gapp$bout_length)~gapp$bout_start_j_c, col=adjustcolor(gapp$col_, alpha.f = 0.5), pch=20, cex=0.2)	
							#points(inc$inc_eff~inc$bout_start_j_c, col=inc$col_,bg=adjustcolor(inc$col_, alpha.f = 0.4), pch=21, cex=0.5)	
						
						
						# predictions
							# before
							polygon(c(pgapB$bout_start_j_c, rev(pgapB$bout_start_j_c)), c(pgapB$lwr, 
								rev(pgapB$upr)), border=NA, col=adjustcolor(col_t ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pgapB$bout_start_j_c, pgapB$pred, col=col_t,lwd=1)
							
							# before
							polygon(c(pgapA$bout_start_j_c, rev(pgapA$bout_start_j_c)), c(pgapA$lwr, 
								rev(pgapA$upr)), border=NA, col=adjustcolor(col_c ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pgapA$bout_start_j_c, pgapA$pred, col=col_c,lwd=1)
													
							#text(x=-2,y=0.725, labels='Before', col='#FCB42C', cex=0.5)
							#text(x=2,y=0.725, labels='After', col='#535F7C', cex=0.5)
					}
				}
												
				dev.off()
			}	
				
		
									
			{# not used	
				{# plot
				  dev.new(width=3.5,height=1.85*2.15)
				  #png(paste(out2,"Figure_5_all_.png", sep=""), width=3.5,height=1.85*2.15,units="in",res=600)
				  par(mfrow=c(3,3),mar=c(1.6,0,0,0.5),oma = c(1, 1.8, 0.2, 0.5),ps=12, mgp=c(1.2,0.35,0), las=1, cex=1, col.axis="grey30",font.main = 1, col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt
				{# constancy
					{# medians
					par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
					plot(xb$med_eff_after~xb$med_eff ,pch=19,xlim=c(0.6,1), ylim=c(0.6,1), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
										
						axis(1, at=seq(0.6,1,by=0.2),labels=c('0.6','0.8','1.0'),cex.axis=0.5,mgp=c(0,-0.20,0))
							mtext('Before experiment',side=1,line=-1, cex=0.5, las=1, col='grey45')
							mtext('Nest attendance',side=1,line=0.3, cex=0.6, las=1, col='grey30')
							
						axis(2, at=seq(0.6,1,by=0.2),labels=c('0.6','0.8','1.0'))
							mtext('After experiment',side=2,line=-0.6, cex=0.5, las=3, col='grey45')
							mtext('Nest attendance',side=2,line=1, cex=0.6, las=3, col='grey30')
						
					lines(c(0.6,1),c(0.6,1), lty=3, col="lightgrey")
						#mtext(expression(bold("a")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
						
					# data
						points(xb$med_eff_after~xb$med_eff, col=xb$col_,bg=adjustcolor(xb$col_, alpha.f = 0.4), pch=21, cex=0.5)	
						
						text(x=0.65,y=1, labels='Focal', col='#FCB42C', cex=0.5,adj=0)
						text(x=0.65,y=0.96, labels='Removed', col='#535F7C', cex=0.5,adj=0)
						
						mtext(expression(bold("a")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
					# predictions	
						# 95%CI for before
						arrows(y0=pcsB$lwr, x0=pcsA$pred,y1=pcsB$upr, x1=pcsA$pred,
							code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)
						# 95%CI for after
						arrows(y0=pcsB$pred, x0=pcsA$lwr,y1=pcsB$pred, x1=pcsA$upr,
							code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)	
						
						points(x=pcsA$pred,y=pcsB$pred, pch=20, cex=0.9,col="red")
					}
					{# boxplot - no points, color box outline - before after
						par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
					
						boxplot(inc_eff ~ treated_sex_exper, data = inc,
											ylab = NULL,xaxt='n', yaxt='n',
											ylim=c(0.58,1),								
											#at=c(1,2,3.5,4.5),
											at=c(1,2,3.5,4.5,7,8,9.5,10.5),
											type='n',
											outcex=0.3,outpch=20,boxwex=0.6,whisklty=1,staplelty=0,#medlwd=1, 
											lwd = 1,
											border=c('#FCB42C','#535F7C','#FCB42C','#535F7C','#FCB42C','#535F7C','#FCB42C','#535F7C'),
											col = adjustcolor("white", alpha.f = 0), # trick for PNGs, to show what is underneath the boxplot else can be taken out
											#outcol="darkgrey",boxcol='darkgrey',whiskcol='darkgrey',staplecol='darkgrey',medcol='darkgrey', 
											par(bty='l')	
											#add=TRUE
											)					
						axis(1, at=c(2.75,8.75),labels=c('Focal','Removed'),cex.axis=0.5,mgp=c(0,-0.20,0))
						text(c(1.5,4,7.5,10), (par("usr")[3])*1.05326705, labels = c('\u2640','\u2642','\u2640','\u2642'), font=4, xpd = TRUE, cex=0.5, col='grey45')
						mtext('Parent',side=1,line=0.3, cex=0.6, las=1, col='grey30')
						
						axis(2, at=seq(0.6,1,by=0.2), labels=FALSE)
						
						text(x=5.75,y=0.71, labels='Before', col='#FCB42C', cex=0.5)
						text(x=5.75,y=0.67, labels='After', col='#535F7C', cex=0.5)
						
						text(x=10.9,y=0.58, labels='*', col='grey30', cex=0.6) #labels='0, 0.3, 0.4, 0.6'
						
						#outliers 0.0129349,0.2786890,0.3592990, 0.5777780
						
						mtext(expression(bold("b")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
						}
					{# day
						par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
						plot(pcon$pred~pcon$bout_start_j_c ,pch=19,xlim=c(-4.25,4), ylim=c(0.6,1), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
											
							axis(1, at=seq(-4,4,by=2),labels=seq(-4,4,by=2),cex.axis=0.5,mgp=c(0,-0.20,0))
								mtext('Day',side=1,line=0.3, cex=0.6, las=1, col='grey30')
								
							axis(2, at=seq(0.6,1,by=0.2),labels=FALSE)
							lines(c(0,0),c(0.58,1), lty=3, col="red")							
						# data
							points(inc$inc_eff[inc$inc_eff>0.58]~inc$bout_start_j_c[inc$inc_eff>0.58], col=adjustcolor(inc$col_[inc$inc_eff>0.58], alpha.f = 0.5), pch=20, cex=0.2)	
							#points(inc$inc_eff~inc$bout_start_j_c, col=inc$col_,bg=adjustcolor(inc$col_, alpha.f = 0.4), pch=21, cex=0.5)	
							
						# predictions
							# before
							polygon(c(pconB$bout_start_j_c, rev(pconB$bout_start_j_c)), c(pconB$lwr, 
								rev(pconB$upr)), border=NA, col=adjustcolor(col_t ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pconB$bout_start_j_c, pconB$pred, col=col_t,lwd=1)
							
							# before
							polygon(c(pconA$bout_start_j_c, rev(pconA$bout_start_j_c)), c(pconA$lwr, 
								rev(pconA$upr)), border=NA, col=adjustcolor(col_c ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pconA$bout_start_j_c, pconA$pred, col=col_c,lwd=1)
												
						# text	
							mtext(expression(bold("c")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
							text(x=inc$bout_start_j_c[inc$inc_eff<0.6],y=0.59, labels='*', col='grey30', cex=0.6) #labels='0, 0.3, 0.4, 0.6' #round(inc$inc_eff[inc$inc_eff<0.6],2 )
							text(x=1.20,y=0.635, labels='0.28', col='grey30', cex=0.4)
							text(x=1.20,y=0.61, labels='0.01', col='grey30', cex=0.4)
							
							text(x=2.75,y=0.61, labels='0.36', col='grey30', cex=0.4)
							text(x=3.5,y=0.635, labels='0.58', col='grey30', cex=0.4)
							
							text(x=-2,y=0.725, labels='Before', col='#FCB42C', cex=0.5)
							text(x=2,y=0.725, labels='After', col='#535F7C', cex=0.5)
					}
				}
				{# bout
					{# medians
					par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
					plot(xb$med_bout_after~xb$med_bout ,pch=19,xlim=c(0,16.5), ylim=c(0,16.5), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
										
						axis(1, at=seq(0,16,by=4),labels=c('0','','8','','16'),cex.axis=0.5,mgp=c(0,-0.20,0))
							mtext('Before experiment',side=1,line=-1, cex=0.5, las=1, col='grey45')
							mtext('Bout [h]',side=1,line=0.3, cex=0.6, las=1, col='grey30')
							
						axis(2, at=seq(0,16,by=4),labels=c('0','','8','','16'))
							mtext('After experiment',side=2,line=-0.6, cex=0.5, las=3, col='grey45')
							mtext('Bout [h]',side=2,line=1, cex=0.6, las=3, col='grey30')
						
					lines(c(0,16),c(0,16), lty=3, col="lightgrey")
						#mtext(expression(bold("a")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
						
					# data
						points(xb$med_bout_after~xb$med_bout, col=xb$col_,bg=adjustcolor(xb$col_, alpha.f = 0.4), pch=21, cex=0.5)	
						
						#text(x=3,y=1, labels='Focal', col='#FCB42C', cex=0.5,adj=0)
						#text(x=0.3,y=0.96, labels='Removed', col='#535F7C', cex=0.5,adj=0)
					
					# text	
						mtext(expression(bold("d")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
					# predictions	
						# 95%CI for before
						arrows(y0=pbsB$lwr, x0=pbsA$pred,y1=pbsB$upr, x1=pbsA$pred,
							code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)
						# 95%CI for after
						arrows(y0=pbsB$pred, x0=pbsA$lwr,y1=pbsB$pred, x1=pbsA$upr,
							code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)	
						
						points(x=pbsA$pred,y=pbsB$pred, pch=20, cex=0.9,col="red")
					}
					{# boxplot - no points, color box outline - before after
						par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
					
						boxplot(bout_length ~ treated_sex_exper, data = inc,
											ylab = NULL,xaxt='n', yaxt='n',
											ylim=c(0,16.5),								
											#at=c(1,2,3.5,4.5),
											at=c(1,2,3.5,4.5,7,8,9.5,10.5),
											type='n',
											outcex=0.3,outpch=20,boxwex=0.6,whisklty=1,staplelty=0,#medlwd=1, 
											lwd = 1,
											border=c('#FCB42C','#535F7C','#FCB42C','#535F7C','#FCB42C','#535F7C','#FCB42C','#535F7C'),
											col = adjustcolor("white", alpha.f = 0), # trick for PNGs, to show what is underneath the boxplot else can be taken out
											#outcol="darkgrey",boxcol='darkgrey',whiskcol='darkgrey',staplecol='darkgrey',medcol='darkgrey', 
											par(bty='l')	
											#add=TRUE
											)					
						axis(1, at=c(2.75,8.75),labels=c('Focal','Removed'),cex.axis=0.5,mgp=c(0,-0.20,0))
						text(c(1.5,4,7.5,10), par("usr")[3]+1.3, labels = c('\u2640','\u2642','\u2640','\u2642'), font=4, xpd = TRUE, cex=0.5, col='grey45')
						mtext('Parent',side=1,line=0.3, cex=0.6, las=1, col='grey30')
						
						axis(2, at=seq(0,16,by=4), labels=FALSE)
						
						#text(x=5.75,y=0.71, labels='Before', col='#FCB42C', cex=0.5)
						#text(x=5.75,y=0.67, labels='After', col='#535F7C', cex=0.5)
						
						#text(x=10.9,y=0.58, labels='*', col='grey30', cex=0.6) #labels='0, 0.3, 0.4, 0.6'
						
						#outliers 0.0129349,0.2786890,0.3592990, 0.5777780
						
						mtext(expression(bold("e")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
						}
					{# day
						par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
						plot(pbout$pred~pbout$bout_start_j_c ,pch=19,xlim=c(-4.25,4), ylim=c(0,16.5), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
											
							axis(1, at=seq(-4,4,by=2),labels=seq(-4,4,by=2),cex.axis=0.5,mgp=c(0,-0.20,0))
								mtext('Day',side=1,line=0.3, cex=0.6, las=1, col='grey30')
								
							axis(2, at=seq(0,16,by=4),labels=FALSE)
							lines(c(0,0),c(0,16.5), lty=3, col="red")							
						# data
							points(inc$bout_length~inc$bout_start_j_c, col=adjustcolor(inc$col_, alpha.f = 0.5), pch=20, cex=0.2)	
							#points(inc$inc_eff~inc$bout_start_j_c, col=inc$col_,bg=adjustcolor(inc$col_, alpha.f = 0.4), pch=21, cex=0.5)	
						
						
						# predictions
							# before
							polygon(c(pboutB$bout_start_j_c, rev(pboutB$bout_start_j_c)), c(pboutB$lwr, 
								rev(pboutB$upr)), border=NA, col=adjustcolor(col_t ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pboutB$bout_start_j_c, pboutB$pred, col=col_t,lwd=1)
							
							# before
							polygon(c(pboutA$bout_start_j_c, rev(pboutA$bout_start_j_c)), c(pboutA$lwr, 
								rev(pboutA$upr)), border=NA, col=adjustcolor(col_c ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pboutA$bout_start_j_c, pboutA$pred, col=col_c,lwd=1)
													
						# text	
							mtext(expression(bold("f")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
							#text(x=-2,y=0.725, labels='Before', col='#FCB42C', cex=0.5)
							#text(x=2,y=0.725, labels='After', col='#535F7C', cex=0.5)
					}
				}
				{# gap
					{# medians
					par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
					plot(yb$med_gap_after~yb$med_gap ,pch=19,xlim=log(c(0.015,700)), ylim=log(c(0.015,700)), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
										
						axis(1, at=log(c(seq(0.1,0.9,by=0.1),seq(1,10,by=1), seq(20,100, by=10), seq(200,700, by=100))),labels=FALSE,col.ticks="grey90")
							axis(1,  at=log(c(0.1,1,10,100,700)),labels=c('0.1','1','10','100','700'),mgp=c(0,-0.20,0))
							mtext('Before experiment',side=1,line=-1, cex=0.5, las=1, col='grey45')
							mtext('Exchange gap [min]',side=1,line=0.3, cex=0.6, las=1, col='grey30')
							
						axis(2,  at=log(c(seq(0.1,0.9,by=0.1),seq(1,10,by=1), seq(20,100, by=10), seq(200,700, by=100))),labels=FALSE,col.ticks="grey90")
						axis(2,  at=log(c(0.1,1,10,100,700)),labels=c('0.1','1','10','100','700'))
							mtext('After experiment',side=2,line=-0.6, cex=0.5, las=3, col='grey45')
							mtext('Exchange gap [min]',side=2,line=1, cex=0.6, las=3, col='grey30')
						
					lines(log(c(0.001,700)),log(c(0.001,700)), lty=3, col="lightgrey")
						#mtext(expression(bold("a")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
						
					# data
						points(log(yb$med_gap_after)~log(yb$med_gap), col=yb$col_,bg=adjustcolor(yb$col_, alpha.f = 0.4), pch=21, cex=0.5)	
						
						#text(x=3,y=1, labels='Focal', col='#FCB42C', cex=0.5,adj=0)
						#text(x=0.3,y=0.96, labels='Removed', col='#535F7C', cex=0.5,adj=0)
					
					# text	
						mtext(expression(bold("g")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
					# predictions	
						# 95%CI for before
						arrows(x0=pgsA$pred,y0=pgsA$pred,x1=pgsB$upr, y1=pgsA$pred,
							code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)
						# 95%CI for after
						arrows(x0=pgsB$pred, y0=pgsA$lwr,x1=pgsB$pred, y1=pgsA$upr,
							code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)	
						
						points(y=pgsA$pred,x=pgsB$pred, pch=20, cex=0.9,col="red")
					}
					{# boxplot - no points, color box outline - before after
						par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
					
						boxplot(log(bout_length) ~ treated_sex_exper, data = gapp,
											ylab = NULL,xaxt='n', yaxt='n',
											ylim=log(c(0.015,700)),								
											#at=c(1,2,3.5,4.5),
											at=c(1,2,3.5,4.5,7,8,9.5,10.5),
											type='n',
											outcex=0.3,outpch=20,boxwex=0.6,whisklty=1,staplelty=0,#medlwd=1, 
											lwd = 1,
											border=c('#FCB42C','#535F7C','#FCB42C','#535F7C','#FCB42C','#535F7C','#FCB42C','#535F7C'),
											col = adjustcolor("white", alpha.f = 0), # trick for PNGs, to show what is underneath the boxplot else can be taken out
											#outcol="darkgrey",boxcol='darkgrey',whiskcol='darkgrey',staplecol='darkgrey',medcol='darkgrey', 
											par(bty='l')	
											#add=TRUE
											)					
						axis(1, at=c(2.75,8.75),labels=c('Focal','Removed'),cex.axis=0.5,mgp=c(0,-0.20,0))
						text(c(1.5,4,7.5,10), par("usr")[3]+log(2.3), labels = c('\u2640','\u2642','\u2640','\u2642'), font=4, xpd = TRUE, cex=0.5, col='grey45')
						mtext('Parent',side=1,line=0.3, cex=0.6, las=1, col='grey30')
						
						axis(2,  at=log(c(seq(0.1,0.9,by=0.1),seq(1,10,by=1), seq(20,100, by=10), seq(200,700, by=100))),labels=FALSE,col.ticks="grey90")
						axis(2,  at=log(c(0.1,1,10,100,700)),labels=FALSE)
						
						#text(x=5.75,y=0.71, labels='Before', col='#FCB42C', cex=0.5)
						#text(x=5.75,y=0.67, labels='After', col='#535F7C', cex=0.5)
						
						#text(x=10.9,y=0.58, labels='*', col='grey30', cex=0.6) #labels='0, 0.3, 0.4, 0.6'
						
						#outliers 0.0129349,0.2786890,0.3592990, 0.5777780
						
						mtext(expression(bold("h")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
						}
					{# day
						par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
						plot(pgap$pred~pgap$bout_start_j_c ,pch=19,xlim=c(-4.25,4), ylim=log(c(0.015,700)),	 xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
											
							axis(1, at=seq(-4,4,by=2),labels=seq(-4,4,by=2),cex.axis=0.5,mgp=c(0,-0.20,0))
								mtext('Day',side=1,line=0.3, cex=0.6, las=1, col='grey30')
								
							axis(2,  at=log(c(seq(0.1,0.9,by=0.1),seq(1,10,by=1), seq(20,100, by=10), seq(200,700, by=100))),labels=FALSE,col.ticks="grey90")
							axis(2,  at=log(c(0.1,1,10,100,700)),labels=FALSE)			

							lines(c(0,log(0.015),c(0,log(700)), lty=3, col="red")								
						# data
							points(log(gapp$bout_length)~gapp$bout_start_j_c, col=adjustcolor(gapp$col_, alpha.f = 0.5), pch=20, cex=0.2)	
							#points(inc$inc_eff~inc$bout_start_j_c, col=inc$col_,bg=adjustcolor(inc$col_, alpha.f = 0.4), pch=21, cex=0.5)	
						
						
						# predictions
							# before
							polygon(c(pgapB$bout_start_j_c, rev(pgapB$bout_start_j_c)), c(pgapB$lwr, 
								rev(pgapB$upr)), border=NA, col=adjustcolor(col_t ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pgapB$bout_start_j_c, pgapB$pred, col=col_t,lwd=1)
							
							# before
							polygon(c(pgapA$bout_start_j_c, rev(pgapA$bout_start_j_c)), c(pgapA$lwr, 
								rev(pgapA$upr)), border=NA, col=adjustcolor(col_c ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pgapA$bout_start_j_c, pgapA$pred, col=col_c,lwd=1)
													
						# text	
							mtext(expression(bold("i")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
							#text(x=-2,y=0.725, labels='Before', col='#FCB42C', cex=0.5)
							#text(x=2,y=0.725, labels='After', col='#535F7C', cex=0.5)
					}
				}
												
				dev.off()
			}	
			{# plot a,b,c
				  dev.new(width=3.5,height=1.85*2.15)
				  #png(paste(out2,"Figure_5_abc.png", sep=""), width=3.5,height=1.85*2.15,units="in",res=600)
				  par(mfrow=c(3,3),mar=c(1.6,0,0,0.5),oma = c(0.5, 1.8, 0.7, 0.5),ps=12, mgp=c(1.2,0.35,0), las=1, cex=1, col.axis="grey30",font.main = 1, col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt
				{# constancy
					{# medians
					par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
					plot(xb$med_eff_after~xb$med_eff ,pch=19,xlim=c(0.6,1), ylim=c(0.6,1), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
										
						axis(1, at=seq(0.6,1,by=0.2),labels=c('0.6','0.8','1.0'),cex.axis=0.5,mgp=c(0,-0.20,0))
							mtext('Before experiment',side=1,line=-1, cex=0.5, las=1, col='grey45')
							mtext('Nest attendance',side=1,line=0.3, cex=0.6, las=1, col='grey30')
							
						axis(2, at=seq(0.6,1,by=0.2),labels=c('0.6','0.8','1.0'))
							mtext('Post experiment',side=2,line=-0.6, cex=0.5, las=3, col='grey45')
							mtext('Nest attendance',side=2,line=1, cex=0.6, las=3, col='grey30')
						
					lines(c(0.6,1),c(0.6,1), lty=3, col="lightgrey")
						#mtext(expression(bold("a")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
						
					# data
						points(xb$med_eff_after~xb$med_eff, col=xb$col_,bg=adjustcolor(xb$col_, alpha.f = 0.4), pch=21, cex=0.5)	
						
						text(x=0.65,y=1, labels='Focal', col='#FCB42C', cex=0.5,adj=0)
						text(x=0.65,y=0.96, labels='Removed', col='#535F7C', cex=0.5,adj=0)
						
						mtext(expression(bold("a")),side=3,line=0, cex=0.7, las=1, col="grey30")
					# predictions	
						# 95%CI for before
						arrows(y0=pcsB$lwr, x0=pcsA$pred,y1=pcsB$upr, x1=pcsA$pred,
							code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)
						# 95%CI for after
						arrows(y0=pcsB$pred, x0=pcsA$lwr,y1=pcsB$pred, x1=pcsA$upr,
							code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)	
						
						points(x=pcsA$pred,y=pcsB$pred, pch=20, cex=0.5,col="red")
					}
					{# boxplot - no points, color box outline - before after
						par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
					
						boxplot(inc_eff ~ treated_sex_exper, data = inc,
											ylab = NULL,xaxt='n', yaxt='n',
											ylim=c(0.58,1),								
											#at=c(1,2,3.5,4.5),
											at=c(1,2,3.5,4.5,7,8,9.5,10.5),
											type='n',
											outcex=0.3,outpch=20,boxwex=0.6,whisklty=1,staplelty=0,#medlwd=1, 
											lwd = 1,
											border=c('#FCB42C','#535F7C','#FCB42C','#535F7C','#FCB42C','#535F7C','#FCB42C','#535F7C'),
											col = adjustcolor("white", alpha.f = 0), # trick for PNGs, to show what is underneath the boxplot else can be taken out
											#outcol="darkgrey",boxcol='darkgrey',whiskcol='darkgrey',staplecol='darkgrey',medcol='darkgrey', 
											par(bty='l')	
											#add=TRUE
											)					
						axis(1, at=c(2.75,8.75),labels=c('Focal','Removed'),cex.axis=0.5,mgp=c(0,-0.20,0))
						text(c(1.5,4,7.5,10), (par("usr")[3])*1.05326705, labels = c('\u2640','\u2642','\u2640','\u2642'), font=4, xpd = TRUE, cex=0.5, col='grey45')
						mtext('Parent',side=1,line=0.3, cex=0.6, las=1, col='grey30')
						
						axis(2, at=seq(0.6,1,by=0.2), labels=FALSE)
						
						text(x=5.75,y=0.71, labels='Before', col='#FCB42C', cex=0.5)
						text(x=5.75,y=0.67, labels='Post', col='#535F7C', cex=0.5)
						
						text(x=10.9,y=0.58, labels='*', col='grey30', cex=0.6) #labels='0, 0.3, 0.4, 0.6'
						
						#outliers 0.0129349,0.2786890,0.3592990, 0.5777780
						
						mtext(expression(bold("b")),side=3,line=0, cex=0.7, las=1, col="grey30")
						}
					{# day
						par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
						plot(pcon$pred~pcon$bout_start_j_c ,pch=19,xlim=c(-4.25,4), ylim=c(0.6,1), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
											
							axis(1, at=seq(-4,4,by=2),labels=seq(-4,4,by=2),cex.axis=0.5,mgp=c(0,-0.20,0))
								mtext('Day',side=1,line=0.3, cex=0.6, las=1, col='grey30')
								
							axis(2, at=seq(0.6,1,by=0.2),labels=FALSE)
							lines(c(0,0),c(0.58,1), lty=3, col="red")							
						# data
							points(inc$inc_eff[inc$inc_eff>0.58]~inc$bout_start_j_c[inc$inc_eff>0.58], col=adjustcolor(inc$col_[inc$inc_eff>0.58], alpha.f = 0.5), pch=20, cex=0.2)	
							#points(inc$inc_eff~inc$bout_start_j_c, col=inc$col_,bg=adjustcolor(inc$col_, alpha.f = 0.4), pch=21, cex=0.5)	
							
						# predictions
							# before
							polygon(c(pconB$bout_start_j_c, rev(pconB$bout_start_j_c)), c(pconB$lwr, 
								rev(pconB$upr)), border=NA, col=adjustcolor(col_t ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pconB$bout_start_j_c, pconB$pred, col=col_t,lwd=1)
							
							# before
							polygon(c(pconA$bout_start_j_c, rev(pconA$bout_start_j_c)), c(pconA$lwr, 
								rev(pconA$upr)), border=NA, col=adjustcolor(col_c ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pconA$bout_start_j_c, pconA$pred, col=col_c,lwd=1)
												
						# text	
							mtext(expression(bold("c")),side=3,line=0, cex=0.7, las=1, col="grey30")
							text(x=inc$bout_start_j_c[inc$inc_eff<0.6],y=0.59, labels='*', col='grey30', cex=0.6) #labels='0, 0.3, 0.4, 0.6' #round(inc$inc_eff[inc$inc_eff<0.6],2 )
							text(x=1.20,y=0.635, labels='0.28', col='grey30', cex=0.4)
							text(x=1.20,y=0.61, labels='0.01', col='grey30', cex=0.4)
							
							text(x=2.75,y=0.61, labels='0.36', col='grey30', cex=0.4)
							text(x=3.5,y=0.635, labels='0.58', col='grey30', cex=0.4)
							
							text(x=-2,y=0.725, labels='Before', col='#FCB42C', cex=0.5)
							text(x=2,y=0.725, labels='Post', col='#535F7C', cex=0.5)
					}
				}
				{# bout
					{# medians
					par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
					plot(xb$med_bout_after~xb$med_bout ,pch=19,xlim=c(0,16.5), ylim=c(0,16.5), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
										
						axis(1, at=seq(0,16,by=4),labels=c('0','','8','','16'),cex.axis=0.5,mgp=c(0,-0.20,0))
							mtext('Before experiment',side=1,line=-1, cex=0.5, las=1, col='grey45')
							mtext('Bout [h]',side=1,line=0.3, cex=0.6, las=1, col='grey30')
							
						axis(2, at=seq(0,16,by=4),labels=c('0','','8','','16'))
							mtext('Post experiment',side=2,line=-0.6, cex=0.5, las=3, col='grey45')
							mtext('Bout [h]',side=2,line=1, cex=0.6, las=3, col='grey30')
						
					lines(c(0,16),c(0,16), lty=3, col="lightgrey")
						#mtext(expression(bold("a")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
						
					# data
						points(xb$med_bout_after~xb$med_bout, col=xb$col_,bg=adjustcolor(xb$col_, alpha.f = 0.4), pch=21, cex=0.5)	
						
						#text(x=3,y=1, labels='Focal', col='#FCB42C', cex=0.5,adj=0)
						#text(x=0.3,y=0.96, labels='Removed', col='#535F7C', cex=0.5,adj=0)
					
					# predictions	
						# 95%CI for before
						arrows(y0=pbsB$lwr, x0=pbsA$pred,y1=pbsB$upr, x1=pbsA$pred,
							code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)
						# 95%CI for after
						arrows(y0=pbsB$pred, x0=pbsA$lwr,y1=pbsB$pred, x1=pbsA$upr,
							code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)	
						
						points(x=pbsA$pred,y=pbsB$pred, pch=20, cex=0.5,col="red")
					}
					{# boxplot - no points, color box outline - before after
						par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
					
						boxplot(bout_length ~ treated_sex_exper, data = inc,
											ylab = NULL,xaxt='n', yaxt='n',
											ylim=c(0,16.5),								
											#at=c(1,2,3.5,4.5),
											at=c(1,2,3.5,4.5,7,8,9.5,10.5),
											type='n',
											outcex=0.3,outpch=20,boxwex=0.6,whisklty=1,staplelty=0,#medlwd=1, 
											lwd = 1,
											border=c('#FCB42C','#535F7C','#FCB42C','#535F7C','#FCB42C','#535F7C','#FCB42C','#535F7C'),
											col = adjustcolor("white", alpha.f = 0), # trick for PNGs, to show what is underneath the boxplot else can be taken out
											#outcol="darkgrey",boxcol='darkgrey',whiskcol='darkgrey',staplecol='darkgrey',medcol='darkgrey', 
											par(bty='l')	
											#add=TRUE
											)					
						axis(1, at=c(2.75,8.75),labels=c('Focal','Removed'),cex.axis=0.5,mgp=c(0,-0.20,0))
						text(c(1.5,4,7.5,10), par("usr")[3]+1.3, labels = c('\u2640','\u2642','\u2640','\u2642'), font=4, xpd = TRUE, cex=0.5, col='grey45')
						mtext('Parent',side=1,line=0.3, cex=0.6, las=1, col='grey30')
						
						axis(2, at=seq(0,16,by=4), labels=FALSE)
						
						#text(x=5.75,y=0.71, labels='Before', col='#FCB42C', cex=0.5)
						#text(x=5.75,y=0.67, labels='After', col='#535F7C', cex=0.5)
						
						#text(x=10.9,y=0.58, labels='*', col='grey30', cex=0.6) #labels='0, 0.3, 0.4, 0.6'
						
						#outliers 0.0129349,0.2786890,0.3592990, 0.5777780
						
						
						}
					{# day
						par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
						plot(pbout$pred~pbout$bout_start_j_c ,pch=19,xlim=c(-4.25,4), ylim=c(0,16.5), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
											
							axis(1, at=seq(-4,4,by=2),labels=seq(-4,4,by=2),cex.axis=0.5,mgp=c(0,-0.20,0))
								mtext('Day',side=1,line=0.3, cex=0.6, las=1, col='grey30')
								
							axis(2, at=seq(0,16,by=4),labels=FALSE)
							lines(c(0,0),c(0,16.5), lty=3, col="red")							
						# data
							points(inc$bout_length~inc$bout_start_j_c, col=adjustcolor(inc$col_, alpha.f = 0.5), pch=20, cex=0.2)	
							#points(inc$inc_eff~inc$bout_start_j_c, col=inc$col_,bg=adjustcolor(inc$col_, alpha.f = 0.4), pch=21, cex=0.5)	
						
						
						# predictions
							# before
							polygon(c(pboutB$bout_start_j_c, rev(pboutB$bout_start_j_c)), c(pboutB$lwr, 
								rev(pboutB$upr)), border=NA, col=adjustcolor(col_t ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pboutB$bout_start_j_c, pboutB$pred, col=col_t,lwd=1)
							
							# before
							polygon(c(pboutA$bout_start_j_c, rev(pboutA$bout_start_j_c)), c(pboutA$lwr, 
								rev(pboutA$upr)), border=NA, col=adjustcolor(col_c ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pboutA$bout_start_j_c, pboutA$pred, col=col_c,lwd=1)
													
					
							#text(x=-2,y=0.725, labels='Before', col='#FCB42C', cex=0.5)
							#text(x=2,y=0.725, labels='After', col='#535F7C', cex=0.5)
					}
				}
				{# gap
					{# medians
					par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
					plot(yb$med_gap_after~yb$med_gap ,pch=19,xlim=log(c(0.015,700)), ylim=log(c(0.015,700)), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
										
						axis(1, at=log(c(seq(0.1,0.9,by=0.1),seq(1,10,by=1), seq(20,100, by=10), seq(200,700, by=100))),labels=FALSE,col.ticks="grey90")
							axis(1,  at=log(c(0.1,1,10,100,700)),labels=c('0.1','1','10','100','700'),mgp=c(0,-0.20,0))
							mtext('Before experiment',side=1,line=-1, cex=0.5, las=1, col='grey45')
							mtext('Exchange gap [min]',side=1,line=0.3, cex=0.6, las=1, col='grey30')
							
						axis(2,  at=log(c(seq(0.1,0.9,by=0.1),seq(1,10,by=1), seq(20,100, by=10), seq(200,700, by=100))),labels=FALSE,col.ticks="grey90")
						axis(2,  at=log(c(0.1,1,10,100,700)),labels=c('0.1','1','10','100','700'))
							mtext('Post experiment',side=2,line=-0.6, cex=0.5, las=3, col='grey45')
							mtext('Exchange gap [min]',side=2,line=1, cex=0.6, las=3, col='grey30')
						
					lines(log(c(0.001,700)),log(c(0.001,700)), lty=3, col="lightgrey")
						#mtext(expression(bold("a")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
						
					# data
						points(log(yb$med_gap_after)~log(yb$med_gap), col=yb$col_,bg=adjustcolor(yb$col_, alpha.f = 0.4), pch=21, cex=0.5)	
						
						#text(x=3,y=1, labels='Focal', col='#FCB42C', cex=0.5,adj=0)
						#text(x=0.3,y=0.96, labels='Removed', col='#535F7C', cex=0.5,adj=0)
					
					# predictions	
						# 95%CI for before
						arrows(x0=pgsA$pred,y0=pgsA$pred,x1=pgsB$upr, y1=pgsA$pred,
							code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)
						# 95%CI for after
						arrows(x0=pgsB$pred, y0=pgsA$lwr,x1=pgsB$pred, y1=pgsA$upr,
							code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)	
						
						points(y=pgsA$pred,x=pgsB$pred, pch=20, cex=0.5,col="red")
					}
					{# boxplot - no points, color box outline - before after
						par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
					
						boxplot(log(bout_length) ~ treated_sex_exper, data = gapp,
											ylab = NULL,xaxt='n', yaxt='n',
											ylim=log(c(0.015,700)),								
											#at=c(1,2,3.5,4.5),
											at=c(1,2,3.5,4.5,7,8,9.5,10.5),
											type='n',
											outcex=0.3,outpch=20,boxwex=0.6,whisklty=1,staplelty=0,#medlwd=1, 
											lwd = 1,
											border=c('#FCB42C','#535F7C','#FCB42C','#535F7C','#FCB42C','#535F7C','#FCB42C','#535F7C'),
											col = adjustcolor("white", alpha.f = 0), # trick for PNGs, to show what is underneath the boxplot else can be taken out
											#outcol="darkgrey",boxcol='darkgrey',whiskcol='darkgrey',staplecol='darkgrey',medcol='darkgrey', 
											par(bty='l')	
											#add=TRUE
											)					
						axis(1, at=c(2.75,8.75),labels=c('Focal','Removed'),cex.axis=0.5,mgp=c(0,-0.20,0))
						text(c(1.5,4,7.5,10), par("usr")[3]+log(2.3), labels = c('\u2640','\u2642','\u2640','\u2642'), font=4, xpd = TRUE, cex=0.5, col='grey45')
						mtext('Parent',side=1,line=0.3, cex=0.6, las=1, col='grey30')
						
						axis(2,  at=log(c(seq(0.1,0.9,by=0.1),seq(1,10,by=1), seq(20,100, by=10), seq(200,700, by=100))),labels=FALSE,col.ticks="grey90")
						axis(2,  at=log(c(0.1,1,10,100,700)),labels=FALSE)
						
						#text(x=5.75,y=0.71, labels='Before', col='#FCB42C', cex=0.5)
						#text(x=5.75,y=0.67, labels='After', col='#535F7C', cex=0.5)
						
						#text(x=10.9,y=0.58, labels='*', col='grey30', cex=0.6) #labels='0, 0.3, 0.4, 0.6'
						
						#outliers 0.0129349,0.2786890,0.3592990, 0.5777780
						
						}
					{# day
						par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
						plot(pgap$pred~pgap$bout_start_j_c ,pch=19,xlim=c(-4.25,4), ylim=log(c(0.015,700)),	 xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
											
							axis(1, at=seq(-4,4,by=2),labels=seq(-4,4,by=2),cex.axis=0.5,mgp=c(0,-0.20,0))
								mtext('Day',side=1,line=0.3, cex=0.6, las=1, col='grey30')
								
							axis(2,  at=log(c(seq(0.1,0.9,by=0.1),seq(1,10,by=1), seq(20,100, by=10), seq(200,700, by=100))),labels=FALSE,col.ticks="grey90")
							axis(2,  at=log(c(0.1,1,10,100,700)),labels=FALSE)			

							lines(c(0,0),log(c(0.015,700)), lty=3, col="red")								
						# data
							points(log(gapp$bout_length)~gapp$bout_start_j_c, col=adjustcolor(gapp$col_, alpha.f = 0.5), pch=20, cex=0.2)	
							#points(inc$inc_eff~inc$bout_start_j_c, col=inc$col_,bg=adjustcolor(inc$col_, alpha.f = 0.4), pch=21, cex=0.5)	
						
						
						# predictions
							# before
							polygon(c(pgapB$bout_start_j_c, rev(pgapB$bout_start_j_c)), c(pgapB$lwr, 
								rev(pgapB$upr)), border=NA, col=adjustcolor(col_t ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pgapB$bout_start_j_c, pgapB$pred, col=col_t,lwd=1)
							
							# before
							polygon(c(pgapA$bout_start_j_c, rev(pgapA$bout_start_j_c)), c(pgapA$lwr, 
								rev(pgapA$upr)), border=NA, col=adjustcolor(col_c ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pgapA$bout_start_j_c, pgapA$pred, col=col_c,lwd=1)
													
							#text(x=-2,y=0.725, labels='Before', col='#FCB42C', cex=0.5)
							#text(x=2,y=0.725, labels='After', col='#535F7C', cex=0.5)
					}
				}
												
				dev.off()
			}	
			{# plot a,b,c, less labels
				  dev.new(width=3.5,height=1.85*2.15)
				  #png(paste(out2,"Figure_5_abc_less_labels.png", sep=""), width=3.5,height=1.85*2.15,units="in",res=600)
				  par(mfrow=c(3,3),mar=c(1.6,0,0,0.5),oma = c(0.5, 1.8, 0.7, 0.5),ps=12, mgp=c(1.2,0.35,0), las=1, cex=1, col.axis="grey30",font.main = 1, col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt
				{# constancy
					{# medians
					par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
					plot(xb$med_eff_after~xb$med_eff ,pch=19,xlim=c(0.6,1), ylim=c(0.6,1), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
										
						axis(1, at=seq(0.6,1,by=0.2),labels=c('0.6','0.8','1.0'),cex.axis=0.5,mgp=c(0,-0.20,0))
							mtext('Before experiment',side=1,line=-1, cex=0.5, las=1, col='grey45')
							#mtext('Nest attendance',side=1,line=0.3, cex=0.6, las=1, col='grey30')
							
						axis(2, at=seq(0.6,1,by=0.2),labels=c('0.6','0.8','1.0'))
							mtext('Post experiment',side=2,line=-0.6, cex=0.5, las=3, col='grey45')
							mtext('Nest attendance',side=2,line=1, cex=0.6, las=3, col='grey30')
						
					lines(c(0.6,1),c(0.6,1), lty=3, col="lightgrey")
						#mtext(expression(bold("a")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
						
					# data
						points(xb$med_eff_after~xb$med_eff, col=xb$col_,bg=adjustcolor(xb$col_, alpha.f = 0.4), pch=21, cex=0.5)	
						
						text(x=0.65,y=1, labels='Focal', col='#FCB42C', cex=0.5,adj=0)
						text(x=0.65,y=0.96, labels='Removed', col='#535F7C', cex=0.5,adj=0)
						
						mtext(expression(bold("a")),side=3,line=0, cex=0.7, las=1, col="grey30")
					# predictions	
						# 95%CI for before
						arrows(y0=pcsB$lwr, x0=pcsA$pred,y1=pcsB$upr, x1=pcsA$pred,
							code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)
						# 95%CI for after
						arrows(y0=pcsB$pred, x0=pcsA$lwr,y1=pcsB$pred, x1=pcsA$upr,
							code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)	
						
						points(x=pcsA$pred,y=pcsB$pred, pch=20, cex=0.5,col="red")
					}
					{# boxplot - no points, color box outline - before after
						par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
					
						boxplot(inc_eff ~ treated_sex_exper, data = inc,
											ylab = NULL,xaxt='n', yaxt='n',
											ylim=c(0.58,1),								
											#at=c(1,2,3.5,4.5),
											at=c(1,2,3.5,4.5,7,8,9.5,10.5),
											type='n',
											outcex=0.3,outpch=20,boxwex=0.6,whisklty=1,staplelty=0,#medlwd=1, 
											lwd = 1,
											border=c('#FCB42C','#535F7C','#FCB42C','#535F7C','#FCB42C','#535F7C','#FCB42C','#535F7C'),
											col = adjustcolor("white", alpha.f = 0), # trick for PNGs, to show what is underneath the boxplot else can be taken out
											#outcol="darkgrey",boxcol='darkgrey',whiskcol='darkgrey',staplecol='darkgrey',medcol='darkgrey', 
											par(bty='l')	
											#add=TRUE
											)					
						axis(1, at=c(2.75,8.75),labels=c('Focal','Removed'),cex.axis=0.5,mgp=c(0,-0.20,0))
						text(c(1.5,4,7.5,10), (par("usr")[3])*1.05326705, labels = c('\u2640','\u2642','\u2640','\u2642'), font=4, xpd = TRUE, cex=0.5, col='grey45')
						mtext('Parent',side=1,line=0.3, cex=0.6, las=1, col='grey30')
						
						axis(2, at=seq(0.6,1,by=0.2), labels=FALSE)
						
						text(x=5.75,y=0.71, labels='Before', col='#FCB42C', cex=0.5)
						text(x=5.75,y=0.67, labels='Post', col='#535F7C', cex=0.5)
						
						text(x=10.9,y=0.58, labels='*', col='grey30', cex=0.6) #labels='0, 0.3, 0.4, 0.6'
						
						#outliers 0.0129349,0.2786890,0.3592990, 0.5777780
						
						mtext(expression(bold("b")),side=3,line=0, cex=0.7, las=1, col="grey30")
						}
					{# day
						par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
						plot(pcon$pred~pcon$bout_start_j_c ,pch=19,xlim=c(-4.25,4), ylim=c(0.6,1), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
											
							axis(1, at=seq(-4,4,by=2),labels=seq(-4,4,by=2),cex.axis=0.5,mgp=c(0,-0.20,0))
								mtext('Day',side=1,line=0.3, cex=0.6, las=1, col='grey30')
								
							axis(2, at=seq(0.6,1,by=0.2),labels=FALSE)
							lines(c(0,0),c(0.58,1), lty=3, col="red")							
						# data
							points(inc$inc_eff[inc$inc_eff>0.58]~inc$bout_start_j_c[inc$inc_eff>0.58], col=adjustcolor(inc$col_[inc$inc_eff>0.58], alpha.f = 0.5), pch=20, cex=0.2)	
							#points(inc$inc_eff~inc$bout_start_j_c, col=inc$col_,bg=adjustcolor(inc$col_, alpha.f = 0.4), pch=21, cex=0.5)	
							
						# predictions
							# before
							polygon(c(pconB$bout_start_j_c, rev(pconB$bout_start_j_c)), c(pconB$lwr, 
								rev(pconB$upr)), border=NA, col=adjustcolor(col_t ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pconB$bout_start_j_c, pconB$pred, col=col_t,lwd=1)
							
							# before
							polygon(c(pconA$bout_start_j_c, rev(pconA$bout_start_j_c)), c(pconA$lwr, 
								rev(pconA$upr)), border=NA, col=adjustcolor(col_c ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pconA$bout_start_j_c, pconA$pred, col=col_c,lwd=1)
												
						# text	
							mtext(expression(bold("c")),side=3,line=0, cex=0.7, las=1, col="grey30")
							text(x=inc$bout_start_j_c[inc$inc_eff<0.6],y=0.59, labels='*', col='grey30', cex=0.6) #labels='0, 0.3, 0.4, 0.6' #round(inc$inc_eff[inc$inc_eff<0.6],2 )
							text(x=1.20,y=0.635, labels='0.28', col='grey30', cex=0.4)
							text(x=1.20,y=0.61, labels='0.01', col='grey30', cex=0.4)
							
							text(x=2.75,y=0.61, labels='0.36', col='grey30', cex=0.4)
							text(x=3.5,y=0.635, labels='0.58', col='grey30', cex=0.4)
							
							text(x=-2,y=0.725, labels='Before', col='#FCB42C', cex=0.5)
							text(x=2,y=0.725, labels='Post', col='#535F7C', cex=0.5)
					}
				}
				{# bout
					{# medians
					par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
					plot(xb$med_bout_after~xb$med_bout ,pch=19,xlim=c(0,16.5), ylim=c(0,16.5), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
										
						axis(1, at=seq(0,16,by=4),labels=c('0','','8','','16'),cex.axis=0.5,mgp=c(0,-0.20,0))
							mtext('Before experiment',side=1,line=-1, cex=0.5, las=1, col='grey45')
							#mtext('Bout [h]',side=1,line=0.3, cex=0.6, las=1, col='grey30')
							
						axis(2, at=seq(0,16,by=4),labels=c('0','','8','','16'))
							mtext('Post experiment',side=2,line=-0.6, cex=0.5, las=3, col='grey45')
							mtext('Bout [h]',side=2,line=1, cex=0.6, las=3, col='grey30')
						
					lines(c(0,16),c(0,16), lty=3, col="lightgrey")
						#mtext(expression(bold("a")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
						
					# data
						points(xb$med_bout_after~xb$med_bout, col=xb$col_,bg=adjustcolor(xb$col_, alpha.f = 0.4), pch=21, cex=0.5)	
						
						#text(x=3,y=1, labels='Focal', col='#FCB42C', cex=0.5,adj=0)
						#text(x=0.3,y=0.96, labels='Removed', col='#535F7C', cex=0.5,adj=0)
					
					# predictions	
						# 95%CI for before
						arrows(y0=pbsB$lwr, x0=pbsA$pred,y1=pbsB$upr, x1=pbsA$pred,
							code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)
						# 95%CI for after
						arrows(y0=pbsB$pred, x0=pbsA$lwr,y1=pbsB$pred, x1=pbsA$upr,
							code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)	
						
						points(x=pbsA$pred,y=pbsB$pred, pch=20, cex=0.5,col="red")
					}
					{# boxplot - no points, color box outline - before after
						par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
					
						boxplot(bout_length ~ treated_sex_exper, data = inc,
											ylab = NULL,xaxt='n', yaxt='n',
											ylim=c(0,16.5),								
											#at=c(1,2,3.5,4.5),
											at=c(1,2,3.5,4.5,7,8,9.5,10.5),
											type='n',
											outcex=0.3,outpch=20,boxwex=0.6,whisklty=1,staplelty=0,#medlwd=1, 
											lwd = 1,
											border=c('#FCB42C','#535F7C','#FCB42C','#535F7C','#FCB42C','#535F7C','#FCB42C','#535F7C'),
											col = adjustcolor("white", alpha.f = 0), # trick for PNGs, to show what is underneath the boxplot else can be taken out
											#outcol="darkgrey",boxcol='darkgrey',whiskcol='darkgrey',staplecol='darkgrey',medcol='darkgrey', 
											par(bty='l')	
											#add=TRUE
											)					
						axis(1, at=c(2.75,8.75),labels=c('Focal','Removed'),cex.axis=0.5,mgp=c(0,-0.20,0))
						text(c(1.5,4,7.5,10), par("usr")[3]+1.3, labels = c('\u2640','\u2642','\u2640','\u2642'), font=4, xpd = TRUE, cex=0.5, col='grey45')
						mtext('Parent',side=1,line=0.3, cex=0.6, las=1, col='grey30')
						
						axis(2, at=seq(0,16,by=4), labels=FALSE)
						
						#text(x=5.75,y=0.71, labels='Before', col='#FCB42C', cex=0.5)
						#text(x=5.75,y=0.67, labels='After', col='#535F7C', cex=0.5)
						
						#text(x=10.9,y=0.58, labels='*', col='grey30', cex=0.6) #labels='0, 0.3, 0.4, 0.6'
						
						#outliers 0.0129349,0.2786890,0.3592990, 0.5777780
						
						
						}
					{# day
						par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
						plot(pbout$pred~pbout$bout_start_j_c ,pch=19,xlim=c(-4.25,4), ylim=c(0,16.5), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
											
							axis(1, at=seq(-4,4,by=2),labels=seq(-4,4,by=2),cex.axis=0.5,mgp=c(0,-0.20,0))
								mtext('Day',side=1,line=0.3, cex=0.6, las=1, col='grey30')
								
							axis(2, at=seq(0,16,by=4),labels=FALSE)
							lines(c(0,0),c(0,16.5), lty=3, col="red")							
						# data
							points(inc$bout_length~inc$bout_start_j_c, col=adjustcolor(inc$col_, alpha.f = 0.5), pch=20, cex=0.2)	
							#points(inc$inc_eff~inc$bout_start_j_c, col=inc$col_,bg=adjustcolor(inc$col_, alpha.f = 0.4), pch=21, cex=0.5)	
						
						
						# predictions
							# before
							polygon(c(pboutB$bout_start_j_c, rev(pboutB$bout_start_j_c)), c(pboutB$lwr, 
								rev(pboutB$upr)), border=NA, col=adjustcolor(col_t ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pboutB$bout_start_j_c, pboutB$pred, col=col_t,lwd=1)
							
							# before
							polygon(c(pboutA$bout_start_j_c, rev(pboutA$bout_start_j_c)), c(pboutA$lwr, 
								rev(pboutA$upr)), border=NA, col=adjustcolor(col_c ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pboutA$bout_start_j_c, pboutA$pred, col=col_c,lwd=1)
													
					
							#text(x=-2,y=0.725, labels='Before', col='#FCB42C', cex=0.5)
							#text(x=2,y=0.725, labels='After', col='#535F7C', cex=0.5)
					}
				}
				{# gap
					{# medians
					par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
					plot(yb$med_gap_after~yb$med_gap ,pch=19,xlim=log(c(0.015,700)), ylim=log(c(0.015,700)), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
										
						axis(1, at=log(c(seq(0.1,0.9,by=0.1),seq(1,10,by=1), seq(20,100, by=10), seq(200,700, by=100))),labels=FALSE,col.ticks="grey90")
							axis(1,  at=log(c(0.1,1,10,100,700)),labels=c('0.1','1','10','100','700'),mgp=c(0,-0.20,0))
							mtext('Before experiment',side=1,line=-1, cex=0.5, las=1, col='grey45')
							#mtext('Exchange gap [min]',side=1,line=0.3, cex=0.6, las=1, col='grey30')
							
						axis(2,  at=log(c(seq(0.1,0.9,by=0.1),seq(1,10,by=1), seq(20,100, by=10), seq(200,700, by=100))),labels=FALSE,col.ticks="grey90")
						axis(2,  at=log(c(0.1,1,10,100,700)),labels=c('0.1','1','10','100','700'))
							mtext('Post experiment',side=2,line=-0.6, cex=0.5, las=3, col='grey45')
							mtext('Exchange gap [min]',side=2,line=1, cex=0.6, las=3, col='grey30')
						
					lines(log(c(0.001,700)),log(c(0.001,700)), lty=3, col="lightgrey")
						#mtext(expression(bold("a")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
						
					# data
						points(log(yb$med_gap_after)~log(yb$med_gap), col=yb$col_,bg=adjustcolor(yb$col_, alpha.f = 0.4), pch=21, cex=0.5)	
						
						#text(x=3,y=1, labels='Focal', col='#FCB42C', cex=0.5,adj=0)
						#text(x=0.3,y=0.96, labels='Removed', col='#535F7C', cex=0.5,adj=0)
					
					# predictions	
						# 95%CI for before
						arrows(x0=pgsA$pred,y0=pgsA$pred,x1=pgsB$upr, y1=pgsA$pred,
							code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)
						# 95%CI for after
						arrows(x0=pgsB$pred, y0=pgsA$lwr,x1=pgsB$pred, y1=pgsA$upr,
							code = 0, col="red", angle = 90, length = .025, lwd=1, lty=1)	
						
						points(y=pgsA$pred,x=pgsB$pred, pch=20, cex=0.5,col="red")
					}
					{# boxplot - no points, color box outline - before after
						par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
					
						boxplot(log(bout_length) ~ treated_sex_exper, data = gapp,
											ylab = NULL,xaxt='n', yaxt='n',
											ylim=log(c(0.015,700)),								
											#at=c(1,2,3.5,4.5),
											at=c(1,2,3.5,4.5,7,8,9.5,10.5),
											type='n',
											outcex=0.3,outpch=20,boxwex=0.6,whisklty=1,staplelty=0,#medlwd=1, 
											lwd = 1,
											border=c('#FCB42C','#535F7C','#FCB42C','#535F7C','#FCB42C','#535F7C','#FCB42C','#535F7C'),
											col = adjustcolor("white", alpha.f = 0), # trick for PNGs, to show what is underneath the boxplot else can be taken out
											#outcol="darkgrey",boxcol='darkgrey',whiskcol='darkgrey',staplecol='darkgrey',medcol='darkgrey', 
											par(bty='l')	
											#add=TRUE
											)					
						axis(1, at=c(2.75,8.75),labels=c('Focal','Removed'),cex.axis=0.5,mgp=c(0,-0.20,0))
						text(c(1.5,4,7.5,10), par("usr")[3]+log(2.3), labels = c('\u2640','\u2642','\u2640','\u2642'), font=4, xpd = TRUE, cex=0.5, col='grey45')
						mtext('Parent',side=1,line=0.3, cex=0.6, las=1, col='grey30')
						
						axis(2,  at=log(c(seq(0.1,0.9,by=0.1),seq(1,10,by=1), seq(20,100, by=10), seq(200,700, by=100))),labels=FALSE,col.ticks="grey90")
						axis(2,  at=log(c(0.1,1,10,100,700)),labels=FALSE)
						
						#text(x=5.75,y=0.71, labels='Before', col='#FCB42C', cex=0.5)
						#text(x=5.75,y=0.67, labels='After', col='#535F7C', cex=0.5)
						
						#text(x=10.9,y=0.58, labels='*', col='grey30', cex=0.6) #labels='0, 0.3, 0.4, 0.6'
						
						#outliers 0.0129349,0.2786890,0.3592990, 0.5777780
						
						}
					{# day
						par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
						plot(pgap$pred~pgap$bout_start_j_c ,pch=19,xlim=c(-4.25,4), ylim=log(c(0.015,700)),	 xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
											
							axis(1, at=seq(-4,4,by=2),labels=seq(-4,4,by=2),cex.axis=0.5,mgp=c(0,-0.20,0))
								mtext('Day',side=1,line=0.3, cex=0.6, las=1, col='grey30')
								
							axis(2,  at=log(c(seq(0.1,0.9,by=0.1),seq(1,10,by=1), seq(20,100, by=10), seq(200,700, by=100))),labels=FALSE,col.ticks="grey90")
							axis(2,  at=log(c(0.1,1,10,100,700)),labels=FALSE)			

							lines(c(0,0),log(c(0.015,700)), lty=3, col="red")								
						# data
							points(log(gapp$bout_length)~gapp$bout_start_j_c, col=adjustcolor(gapp$col_, alpha.f = 0.5), pch=20, cex=0.2)	
							#points(inc$inc_eff~inc$bout_start_j_c, col=inc$col_,bg=adjustcolor(inc$col_, alpha.f = 0.4), pch=21, cex=0.5)	
						
						
						# predictions
							# before
							polygon(c(pgapB$bout_start_j_c, rev(pgapB$bout_start_j_c)), c(pgapB$lwr, 
								rev(pgapB$upr)), border=NA, col=adjustcolor(col_t ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pgapB$bout_start_j_c, pgapB$pred, col=col_t,lwd=1)
							
							# before
							polygon(c(pgapA$bout_start_j_c, rev(pgapA$bout_start_j_c)), c(pgapA$lwr, 
								rev(pgapA$upr)), border=NA, col=adjustcolor(col_c ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pgapA$bout_start_j_c, pgapA$pred, col=col_c,lwd=1)
													
							#text(x=-2,y=0.725, labels='Before', col='#FCB42C', cex=0.5)
							#text(x=2,y=0.725, labels='After', col='#535F7C', cex=0.5)
					}
				}
												
				dev.off()
			}	
			
				
				{# boxplot - no points, color box outline - sex
					par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
				
					boxplot(inc_eff ~ treated_sex_exper, data = inc, 
											 ylab =NULL, xaxt='n',yaxt='n',
											 ylim=c(0.60,1),								
											par(bty='n'),
										#at=c(1,2,3.5,4.5),at=c(1,2,3.7,4.7), 
										at=c(1,2,3.5,4.5,7,8,9.5,10.5), type='n',
										outcex=0.3, outpch=20,boxwex=0.6,whisklty=1,staplelty=0,medlwd=1, 
										lwd = 0.25, 
												
										outcol="darkgrey",boxcol='darkgrey',whiskcol='darkgrey',staplecol='darkgrey',medcol='darkgrey'
										) # col=z_g$cols, border=z_g$cols
										
						boxplot(inc_eff ~ treated_sex_exper, data = inc,
										ylab = NULL,xaxt='n', yaxt='n',
										ylim=c(0.60,1),								
										#at=c(1,2,3.5,4.5),
										at=c(1,2,3.5,4.5,7,8,9.5,10.5),
										type='n',
										outcex=0.3,outpch=20,boxwex=0.6,whisklty=1,staplelty=0,#medlwd=1, 
										lwd = 1,
										border=c('#FCB42C','#FCB42C','#535F7C','#535F7C','#FCB42C','#FCB42C','#535F7C','#535F7C'),
										col = c('grey','white','grey','white','grey','white','grey','white'), # trick for PNGs, to show what is underneath the boxplot else can be taken out
										#outcol="darkgrey",boxcol='darkgrey',whiskcol='darkgrey',staplecol='darkgrey',medcol='darkgrey', 
										par(bty='l'),		
										add=TRUE
										)					
					axis(1, at=c(2.75,8.75),labels=c('Focal','Removed'),cex.axis=0.5,mgp=c(0,-0.15,0))
					text(c(1.5,4,7.5,10), par("usr")[3]+0.03, labels = c('\u2640','\u2642','\u2640','\u2642'), font=4, xpd = TRUE, cex=0.5, col=c('#FCB42C','#535F7C'))
					mtext('Parent',side=1,line=0.3, cex=0.6, las=1, col='grey30')
					
					axis(2, at=seq(0.6,1,by=0.2), labels=FALSE)
					
					mtext(expression(bold("b")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
					}
				{# boxplot - no points, color box outline - sex
					par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
				
					boxplot(inc_eff ~ treated_sex_exper, data = inc, 
											 ylab =NULL, xaxt='n',yaxt='n',
											 ylim=c(0.60,1),								
											par(bty='n'),
										#at=c(1,2,3.5,4.5),at=c(1,2,3.7,4.7), 
										at=c(1,2,3.5,4.5,7,8,9.5,10.5), type='n',
										outcex=0.3, outpch=20,boxwex=0.6,whisklty=1,staplelty=0,medlwd=1, 
										lwd = 0.25, 
												
										outcol="darkgrey",boxcol='darkgrey',whiskcol='darkgrey',staplecol='darkgrey',medcol='darkgrey'
										) # col=z_g$cols, border=z_g$cols
										
					stripchart(inc_eff ~ factor(treated_sex_exper), vertical = TRUE, data = inc, method = "jitter", add = TRUE, 
										at=c(1,2,3.5,4.5,7,8,9.5,10.5),
										pch = 20,cex=0.3, 
										col="gray63"
										#bg=adjustcolor("gray63", alpha.f = 0.4)
										)
										
					boxplot(inc_eff ~ treated_sex_exper, data = inc,
										ylab = NULL,xaxt='n', yaxt='n',
										ylim=c(0.60,1),								
										#at=c(1,2,3.5,4.5),
										at=c(1,2,3.5,4.5,7,8,9.5,10.5),
										type='n',
										outcex=0.3,outpch=20,boxwex=0.6,whisklty=1,staplelty=0,#medlwd=1, 
										lwd = 1,
										border=c('#FCB42C','#FCB42C','#535F7C','#535F7C','#FCB42C','#FCB42C','#535F7C','#535F7C'),
										col = c('white','grey','white','grey','white','grey','white','grey'), # trick for PNGs, to show what is underneath the boxplot else can be taken out
										#outcol="darkgrey",boxcol='darkgrey',whiskcol='darkgrey',staplecol='darkgrey',medcol='darkgrey', 
										par(bty='l'),		
										add=TRUE
										)					
					axis(1, at=c(2.75,8.75),labels=c('Focal','Removed'),cex.axis=0.5,mgp=c(0,-0.15,0))
					text(c(1.5,4,7.5,10), par("usr")[3]+0.03, labels = c('\u2640','\u2642','\u2640','\u2642'), font=4, xpd = TRUE, cex=0.5, col=c('#FCB42C','#535F7C'))
					mtext('Parent',side=1,line=0.4, cex=0.6, las=1, col='grey30')
					
					axis(2, at=seq(0.6,1,by=0.2), labels=FALSE)
					
					mtext(expression(bold("b")),side=3,line=-0.4, cex=0.7, las=1,adj=1, col="grey30")
					}
				{# old		
			{# old constancy
				ms=lmer(inc_eff~exper+(1|bird_ID),inc, REML=FALSE)
				ms_=lmer(inc_eff~scale(bout_length)+exper+(1|bird_ID),inc, REML=FALSE)
				msw=lmer(inc_eff~exper+(1|bird_ID),inc,weights = bout_length, REML=FALSE)
				mst=lmer(inc_eff~exper*treated_bird+(1|bird_ID),inc, REML=FALSE)
				m=lmer(inc_eff~exper*bout_start_j_c+(bout_start_j_c|bird_ID),inc, REML=FALSE)
				m_int=lmer(inc_eff~exper*treated_bird*bout_start_j_c+(bout_start_j_c|bird_ID),inc, REML=FALSE)
				o=data.frame(model=c('m','ms_','ms','msw','mst','m_int'), AIC=c(AICc(m, nobs=25),AICc(ms_, nobs=25),AICc(ms,nobs=25),AICc(msw,nobs=25),AICc(mst,nobs=25),AICc(m_int,nobs=25)))
						o$delta=o$AIC-min(o$AIC)
						o$prob=round(exp(-0.5*o$delta)/sum(exp(-0.5*o$delta)),4)
						o$ER=max(o$prob)/o$prob
						o[order(o$delta),]
				
				}
			{#old  bout length
				ms=lmer(bout_length~exper+(1|bird_ID),inc, REML=FALSE)
				m=lmer(bout_length~exper*bout_start_j_c+(bout_start_j_c|bird_ID),inc, REML=FALSE)
				m_int=lmer(bout_length~exper*treated_bird*bout_start_j_c+(bout_start_j_c|bird_ID),inc, REML=FALSE)
				o=data.frame(model=c('m','ms','m_int'), AIC=c(AICc(m, nobs=25),AICc(ms,nobs=25),AICc(m_int,nobs=25)))
						o$delta=o$AIC-min(o$AIC)
						o$prob=round(exp(-0.5*o$delta)/sum(exp(-0.5*o$delta)),4)
						o$ER=max(o$prob)/o$prob
						o[order(o$delta),]
			}
			{#old  gap presence
				gap$present=ifelse(gap$bout_length==0,0,1)
				
				ms=glmer(present~exper+(1|bird_ID),gap, family='binomial', REML=FALSE)
				mst=glmer(present~exper*treated_bird+(1|bird_ID),gap, family='binomial',REML=FALSE)
				m=glmer(present~exper*bout_start_j_c+(bout_start_j_c|bird_ID),gap, family='binomial',control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) # uses different optimizer (results same as without it, but without it model does not converge)
				m_int=glmer(present~exper*treated_bird*bout_start_j_c+(bout_start_j_c|bird_ID),gap, family='binomial',control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) # uses different optimizer (results same as without it, but without it model does not converge)
					o=data.frame(model=c('m','mst','ms','m_int'), AIC=c(AICc(m, nobs=25),AICc(mst, nobs=25),AICc(ms,nobs=25),AICc(m_int,nobs=25)))
						o$delta=o$AIC-min(o$AIC)
						o$prob=round(exp(-0.5*o$delta)/sum(exp(-0.5*o$delta)),4)
						o$ER=max(o$prob)/o$prob
						o[order(o$delta),]
			}
			{# old gap length 
				gapp=gap[gap$bout_length>0,]
				densityplot(~gap$bout_length)
				ms=lmer(bout_length~exper+(1|bird_ID),gapp, REML=FALSE)
				m=lmer(bout_length~exper*bout_start_j_c+(bout_start_j_c|bird_ID),gapp, REML=FALSE)
				m_int=lmer(bout_length~exper*treated_bird*bout_start_j_c+(bout_start_j_c|bird_ID),gapp, REML=FALSE)
				o=data.frame(model=c('m','ms','m_int'), AIC=c(AICc(m, nobs=25),AICc(ms,nobs=25),AICc(m_int,nobs=25)))
						o$delta=o$AIC-min(o$AIC)
						o$prob=round(exp(-0.5*o$delta)/sum(exp(-0.5*o$delta)),4)
						o$ER=max(o$prob)/o$prob
						o[order(o$delta),]
			
			}
		}
				
				legend("bottomleft", legend=c('treated_bird','control_bird'), pch=19, col=c( '#FCB42C','#535F7C'), cex=0.8, bty='n')
					#dev.off()
			{# incubation
				{# boxplot per nest
				 # incubation
					ggplot(inc, aes(x=interaction(treated_bird, sex, nest), y=bout_length, col=exper))+geom_point()+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2))
					
					ggplot(inc, aes(x=interaction(sex, nest), y=bout_length, col=exper, fill=treated_bird))+geom_boxplot()+ stat_summary(data=inc, aes(x=interaction(sex, nest), y=bout_length, col=exper,fill=treated_bird),fun.y=mean, geom="point")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2))
					
					ggplot(inc, aes(x=nest, y=bout_length, col=exper, fill=treated_bird))+geom_boxplot()+ stat_summary(data=inc, aes(x=nest, y=bout_length, col=exper,fill=treated_bird),fun.y=mean, geom="point")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2))+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_colour_manual(values = c("black", "grey"))
					
					ggplot(inc, aes(x=interaction(sex, nest), y=bout_length, col=exper, fill=treated_bird))+geom_boxplot()+ stat_summary(data=inc, aes(x=interaction(sex, nest), y=bout_length, col=exper,fill=treated_bird),fun.y=mean, geom="point")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2))+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_colour_manual(values = c("black", "grey"))
					
					ggplot(inc, aes(x=interaction(treated_bird, sex, nest), y=bout_length, col=exper))+geom_boxplot()+ stat_summary(data=inc, aes(x=interaction(treated_bird, sex, nest), y=bout_length, col=exper),fun.y=mean, geom="point")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2))
					
				# inc_eff
					ggplot(inc, aes(x=interaction(treated_bird, sex, nest), y=inc_eff, col=exper))+geom_point()+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2))
					ggplot(inc, aes(x=interaction(treated_bird, sex, nest), y=inc_eff, col=exper))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2))
					ggplot(inc, aes(x=interaction(treated_bird, sex, nest), y=inc_eff, col=exper))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2))+stat_summary(data=inc, aes(x=interaction(treated_bird, sex, nest), y=inc_eff, col=exper),fun.y=mean, geom="point")
					
					
					ggplot(inc, aes(x=nest, y=inc_eff, col=exper, fill=treated_bird))+geom_boxplot()+ stat_summary(data=inc, aes(x=nest, y=inc_eff, col=exper,fill=treated_bird),fun.y=mean, geom="point")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2))+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_colour_manual(values = c("black", "grey"))
					
					ggplot(inc, aes(x=nest, y=inc_eff, col=treated_bird, fill=exper))+geom_boxplot()+ stat_summary(data=inc, aes(x=nest, y=inc_eff, col=treated_bird, fill=exper),fun.y=mean, geom="point")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2))+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_colour_manual(values = c("black", "grey"))
						
				}
				{# boxplot per sex
					# incubation
					ggplot(inc, aes(x=sex, y=bout_length, col=treated_bird, fill=exper))+geom_boxplot()+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_colour_manual(values = c("black", "grey"))
					ggplot(inc, aes(x=treated_bird, y=bout_length, col=exper, fill=sex))+geom_boxplot()+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_colour_manual(values = c("black", "grey"))
					
					# inc_eff
					ggplot(inc, aes(x=sex, y=inc_eff, col=treated_bird, fill=exper))+geom_boxplot()+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_colour_manual(values = c("black", "grey"))
					ggplot(inc, aes(x=sex, y=inc_eff, col=treated_bird, fill=exper))+geom_boxplot()+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_colour_manual(values = c("black", "grey"))+geom_boxplot()+coord_cartesian(ylim = c(0.8,1))
					ggplot(inc, aes(x=treated_bird, y=inc_eff, col=exper, fill=sex))+geom_boxplot()+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_colour_manual(values = c("black", "grey"))
					
					ggplot(inc, aes(x=interaction(treated_bird, sex), y=inc_eff, col=exper))+geom_boxplot()
					ggplot(inc, aes(x=interaction(treated_bird, sex), y=inc_eff, col=exper))+geom_boxplot()+coord_cartesian(ylim = c(0.8,1))
					#ggplot(inc[inc$treated_bird=='y',], aes(x=interaction(exper,sex), y=bout_length, col=interaction(exper,sex)))+geom_boxplot()
					#ggplot(inc[inc$treated_bird=='y',], aes(x=interaction(exper,sex), y=inc_eff, col=interaction(exper,sex)))+geom_boxplot()
					#ggplot(inc[inc$treated_bird=='y',], aes(x=interaction(exper,sex), y=inc_eff, col=interaction(exper,sex)))+geom_boxplot()+coord_cartesian(ylim = c(0.8,1))
				
				}
				{# all inc eff
						{# plot
						x=ddply(inc,.(nest,sex,exper, treated_bird), summarise, med_bout=median(bout_length),mean_bout=mean(bout_length), med_eff=median(inc_eff),mean_eff=mean(inc_eff))
									ggplot(x,aes(x=med_eff,y=mean_eff))+geom_point()
									plot(med_eff~mean_eff,x, xlim=c(0,1), ylim=c(0,1))
									abline(a=0,b=1)
									
									x[x$med_eff<0.2,]
									densityplot(~x$med_eff)
						
						xb=x[x$exper=='b',]
						xa=x[x$exper=='a',]
						xb$med_eff_after=xa$med_eff[match(paste(xb$nest,xb$sex),paste(xa$nest,xa$sex))]
						xb$col_=ifelse(xb$treated_bird=="y", "#FCB42C","#535F7C")
						dev.new(width=2.5,height=2.5)
						#png(paste(out2,"constancy_treated_control_relationship.png", sep=""), width=2.5,height=2.5,units="in",res=600)
						par(mar=c(2,2,0.7,0.7), ps=12, mgp=c(1.2,0.35,0), las=1, cex.lab=0.7, cex.axis=0.6, tcl=-0.15,bty="l",xpd=TRUE, col.axis="grey30", col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt

						plot(xb$med_eff_after[xb$nest!='s806']~xb$med_eff[xb$nest!='s806'], col=xb$col_[xb$nest!='s806'],bg=adjustcolor( xb$col_[xb$nest!='s806'], alpha.f = 0.4), pch=21,xlim=c(0,1), ylim=c(0,1), ylab='Nest attendance - after', xlab=NA, xaxt='n')
						lines(c(0,1),c(0,1), lty=3, col="lightgrey")
						axis(1, at=seq(0,1,by=0.2), labels=FALSE)
						text(seq(0,1,by=0.2), par("usr")[3]-0.05, labels = seq(0,1,by=0.2),  xpd = TRUE, cex=0.6, col="grey30") #labels = z_g$genus
						mtext('Nest attendance - before',side=1,line=0.6, cex=0.7, las=1, col='grey30')
						
						legend("bottomleft", legend=c('treated_bird','control_bird'), pch=19, col=c( '#FCB42C','#535F7C'), cex=0.8, bty='n')
						#dev.off()
						}
						{# model
							m=lmer(inc_eff~scale(bout_length)+exper+(1|nest),inc, REML=FALSE)
							m_=lmer(log(inc_eff)~scale(bout_length)+exper+(1|nest),inc, REML=FALSE)
							m1=lmer(inc_eff~scale(bout_length)+exper*treated_bird+(1|nest),inc, REML=FALSE)
							m2=lmer(inc_eff~scale(bout_length)+exper*sex+(1|nest),inc, REML=FALSE)
							m3=lmer(inc_eff~scale(bout_length)+exper*treated_bird*sex+(1|nest),inc, REML=FALSE)
								aic=AIC(m,m_,m1,m2,m3)
								aic[order(aic$AIC),]
								
							plot(allEffects(m))
							summary(m)
							# simulation		
							nsim <- 2000
							bsim <- sim(m, n.sim=nsim)  
							apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							fixef(m)
						}
				}
				{# limited inc eff
					x=ddply(inc,.(nest,sex,exper, treated_bird), summarise, med_bout=median(bout_length), med_eff=median(inc_eff))
					x=x[x$med_eff>0.5,]
					xb=x[x$exper=='b',]
					xa=x[x$exper=='a',]
					xb$med_eff_after=xa$med_eff[match(paste(xb$nest,xb$sex),paste(xa$nest,xa$sex))]
					xb$col_=ifelse(xb$treated_bird=="y", "#FCB42C","#535F7C")
				dev.new(width=2.5,height=2.5)
				#png(paste(out2,"constancy_treated_control_relationship.png", sep=""), width=2.5,height=2.5,units="in",res=600)
				par(mar=c(2,2,0.7,0.7), ps=12, mgp=c(1.2,0.35,0), las=1, cex.lab=0.7, cex.axis=0.6, tcl=-0.15,bty="l",xpd=TRUE, col.axis="grey30", col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt

				plot(xb$med_eff_after[xb$nest!='s806']~xb$med_eff[xb$nest!='s806'], col=xb$col_[xb$nest!='s806'],bg=adjustcolor( xb$col_[xb$nest!='s806'], alpha.f = 0.4), pch=21,xlim=c(0.8,1), ylim=c(0.7,1), ylab='Nest attendance - after', xlab=NA, xaxt='n')
				lines(c(0.8,1),c(0.8,1), lty=3, col="lightgrey")
				axis(1, at=seq(0.8,1,by=0.05), labels=FALSE)
				text(seq(0.8,1,by=0.05), par("usr")[3]-0.015, labels = seq(0.8,1,by=0.05),  xpd = TRUE, cex=0.6, col="grey30") #labels = z_g$genus
				mtext('Nest attendance - before',side=1,line=0.6, cex=0.7, las=1, col='grey30')
				
				legend("bottomleft", legend=c('treated_bird','control_bird'), pch=19, col=c( '#FCB42C','#535F7C'), cex=0.8, bty='n')
				#dev.off()
				}
				
				{# bout
					{# all bout b vs a
						x=ddply(inc,.(nest,sex,exper, treated_bird), summarise, med_bout=median(bout_length),mean_bout=mean(bout_length), med_eff=median(inc_eff),mean_eff=mean(inc_eff))
						
							plot(med_bout~mean_bout,x, xlim=c(0,15), ylim=c(0,15))
							abline(a=0, b=1)
								
						xb=x[x$exper=='b',]
						xa=x[x$exper=='a',]
						xb$med_bout_after=xa$med_bout[match(paste(xb$nest,xb$sex),paste(xa$nest,xa$sex))]
						xb$col_=ifelse(xb$treated_bird=="y", "#FCB42C","#535F7C")
						dev.new(width=2.5,height=2.5)
						#png(paste(out2,"constancy_treated_control_relationship.png", sep=""), width=2.5,height=2.5,units="in",res=600)
						par(mar=c(2,2,0.7,0.7), ps=12, mgp=c(1.2,0.35,0), las=1, cex.lab=0.7, cex.axis=0.6, tcl=-0.15,bty="l",xpd=TRUE, col.axis="grey30", col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt

						plot(xb$med_bout_after~xb$med_bout, col=xb$col_,bg=adjustcolor( xb$col_, alpha.f = 0.4), pch=21,xlim=c(0,15), ylim=c(0,15), ylab='Bout - after [h]', xlab=NA, xaxt='n')
						lines(c(0,15),c(0,15), lty=3, col="lightgrey")
						axis(1, at=seq(0,15,by=5), labels=FALSE)
						text(seq(0,15,by=5), par("usr")[3]-0.7, labels = seq(0,15,by=5),  xpd = TRUE, cex=0.6, col="grey30") #labels = z_g$genus
						mtext('bout - before [h]',side=1,line=0.6, cex=0.7, las=1, col='grey30')
						
						legend("bottomleft", legend=c('treated_bird','control_bird'), pch=19, col=c( '#FCB42C','#535F7C'), cex=0.6, bty='n')
						
						ggplot(x[x$treated_bird=='y',], aes(x=interaction(exper,sex), y=med_bout, col=interaction(exper,sex)))+geom_boxplot()
						ggplot(inc[inc$treated_bird=='y',], aes(x=interaction(exper,sex), y=bout_length, col=interaction(exper,sex)))+geom_boxplot()
						#dev.off()
					}
					{# all bout f vs m
						x=ddply(inc,.(nest,sex,exper, treated_bird), summarise, med_bout=median(bout_length), med_eff=median(inc_eff))
						xf=x[x$sex=='f',]
						xm=x[x$sex=='m',]
						xf$med_bout_male=xm$med_bout[match(paste(xf$nest,xf$exper),paste(xm$nest,xm$exper))]
						xf$col_=ifelse(xf$exper=="b", "#FCB42C","#535F7C")
						
						dev.new(width=2.5,height=2.5)
						#png(paste(out2,"constancy_treated_control_relationship.png", sep=""), width=2.5,height=2.5,units="in",res=600)
						par(mar=c(2,2,0.7,0.7), ps=12, mgp=c(1.2,0.35,0), las=1, cex.lab=0.7, cex.axis=0.6, tcl=-0.15,bty="l",xpd=TRUE, col.axis="grey30", col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt

						plot(xf$med_bout_male~xf$med_bout, col=xf$col_,bg=adjustcolor( xf$col_, alpha.f = 0.4), pch=21,xlim=c(0,15), ylim=c(0,15), ylab='Bout - female [h]', xlab=NA, xaxt='n')
						lines(c(0,15),c(0,15), lty=3, col="lightgrey")
						axis(1, at=seq(0,15,by=5), labels=FALSE)
						text(seq(0,15,by=5), par("usr")[3]-0.7, labels = seq(0,15,by=5),  xpd = TRUE, cex=0.6, col="grey30") #labels = z_g$genus
						mtext('bout - male [h]',side=1,line=0.6, cex=0.7, las=1, col='grey30')
						
						legend("bottomleft", legend=c('before','after'), pch=19, col=c( '#FCB42C','#535F7C'), cex=0.6, bty='n')
						
						ggplot(x[x$treated_bird=='y',], aes(x=interaction(exper,sex), y=med_bout, col=interaction(exper,sex)))+geom_boxplot()
						ggplot(x[x$treated_bird=='n',], aes(x=interaction(exper,sex), y=med_bout, col=interaction(exper,sex)))+geom_boxplot()
						ggplot(inc[inc$treated_bird=='y',], aes(x=interaction(exper,sex), y=bout_length, col=interaction(exper,sex)))+geom_boxplot()
						ggplot(inc[inc$treated_bird=='n',], aes(x=interaction(exper,sex), y=bout_length, col=interaction(exper,sex)))+geom_boxplot()
						#dev.off()
					
						ggplot(inc, aes(x=interaction(exper,sex, treated_bird), y=bout_length, col=interaction(exper,sex,treated_bird)))+geom_boxplot()
					}
					{# change in bout length in captive and treated over time 0.5,1.5,2.5 reflect captive
						inc_a=inc[inc$exper=='a',]
						inc_a=inc_a[order(inc_a$nest, inc_a$treated_bird, inc_a$bout_start),]
						inc_a=ddply(inc_a,.(nest, treated_bird), transform, bout_order=c(1:length(bout_start)))
						head(inc_a)
						inc_a$bout_order[inc_a$treated_bird=='n']=inc_a$bout_order[inc_a$treated_bird=='n']-0.5
						inc_a_n=inc_a[inc_a$treated_bird=='n',]
						inc_a_y=inc_a[inc_a$treated_bird=='y',]
						inc_a_y$bout_captive=inc_a_n$bout_length[match(paste(inc_a_y$nest,inc_a_y$sex),paste(inc_a_n$nest,inc_a_n$sex))]
						
						ggplot(inc_a, aes(x=bout_order,y=bout_length, col=nest))+geom_point()+geom_line(aes(group = nest))
					
						gap_a=gap[gap$exper=='a',]
						gap_a=gap_a[order(gap_a$nest, gap_a$treated_bird, gap_a$bout_start),]
						gap_a=ddply(gap_a,.(nest, treated_bird), transform, bout_order=c(1:length(bout_start)))
						head(gap_a)
						gap_a$bout_order[gap_a$treated_bird=='n']=gap_a$bout_order[gap_a$treated_bird=='n']-0.5
						gap_a_n=gap_a[gap_a$treated_bird=='n',]
						gap_a_y=gap_a[gap_a$treated_bird=='y',]
						gap_a_y$bout_captive=gap_a_n$bout_length[match(paste(gap_a_y$nest,gap_a_y$sex),paste(gap_a_n$nest,gap_a_n$sex))]
						
						ggplot(gap_a, aes(x=bout_order,y=bout_length, col=nest))+geom_point()+geom_line(aes(group = nest))
						ggplot(gap_a, aes(x=bout_order,y=log(bout_length+0.01), col=nest))+geom_point()+geom_line(aes(group = nest))
					}
					
					
					{# model
								{# mean center bout - not necessary
									#within ID and with before
									incsplit=split(inc,paste(inc$nest, inc$sex))
										foo=lapply(incsplit,function(x) {
																		#x=incsplit$"s404 f"
																		#x$t_a=c(x$treat[-1], NA) 	
																		x$bout_length_b=x$bout_length - mean(x$bout_length[x$exper=='b'])
																		return(x)
																		}
																		)
									
									inc=do.call(rbind, foo)	
									m=lmer(bout_length_b~exper+(1|nest),inc)
								}
								
							
							inc$bird_ID=paste(inc$nest,inc$sex)
							inc$int=as.factor(interaction(inc$exper, inc$treated,inc$sex))
							unique(inc$int)
							densityplot(~inc$bout_length)
							
							m=lmer(bout_length~exper+(1|bird_ID),inc, REML=FALSE)
							m_=lmer(bout_length~exper*treated_bird+(1|bird_ID),inc, REML=FALSE)
							mp=lmer(bout_length~exper*bout_start_j_c+(bout_start_j_c|bird_ID),inc, REML=FALSE)
							mp_=lmer(bout_length~exper*treated_bird*bout_start_j_c+(bout_start_j_c|bird_ID),inc, REML=FALSE)

							o=data.frame(model=c('m','m_','mp','mp_'), AIC=c(AICc(m, nobs=25),AICc(m_,nobs=25),AICc(mp,nobs=25),AICc(mp_,nobs=25)))
							o$delta=o$AIC-min(o$AIC)
							o$prob=round(exp(-0.5*o$delta)/sum(exp(-0.5*o$delta)),4)
							o$ER=max(o$prob)/o$prob
							o[order(o$delta),]
						
							
							m=lmer(bout_length~exper+(1|nest),inc)
							
							m_=lmer(log(bout_length)~exper+(1|nest),inc)
							m1=lmer(bout_length~exper*treated_bird+(1|nest),inc)
							m2=lmer(bout_length~exper*sex+(1|nest),inc)
							m3=lmer(bout_length~exper*treated_bird*sex+(1|nest),inc)
							m4=lmer(bout_length~int+(1|nest),inc)
								aic=AIC(m,m_,m1,m2,m3)
								aic[order(aic$AIC),]
								aic=AIC(m,m1,m2,m3,m4)
								aic[order(aic$AIC),]
								
							plot(allEffects(mps))
							summary(m)
							# simulation		
							nsim <- 2000
							bsim <- sim(m, n.sim=nsim)  
							apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							fixef(m)
							
							plot(allEffects(m3))
							plot(allEffects(m4))
							# simulation		
							nsim <- 2000
							bsim <- sim(m3, n.sim=nsim)  
							apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							fixef(m3)
							
							summary(glht(m4, linfct = mcp(int = c(	"b.y.m - a.y.m = 0",
																	"b.y.f - a.y.f = 0", #
																	"b.n.m - a.n.m = 0", #
																	"b.n.f - a.n.f = 0",
																	
																	"a.n.m - a.n.f = 0",
																	"a.y.m - a.y.f = 0",#
																	
																	"a.y.m - a.n.m = 0",#
																	"a.y.f - a.n.f = 0"
																))))
							
					}
					{# inc$bout_start_j 
					x=ddply(inc,.(nest,sex,exper, treated_bird), summarise, med_bout=median(bout_length), med_eff=median(inc_eff),med_bout_start_j=median(bout_start_j))
					
					ggplot(inc, aes(x=bout_start_j, y=bout_length, col=interaction(exper,treated_bird)))+geom_point()+stat_smooth()
					ggplot(inc, aes(x=bout_start_j, y=bout_length, col=interaction(exper,treated_bird)))+geom_point()+stat_smooth(method='lm')
					ggplot(x, aes(x=med_bout_start_j, y=med_bout, col=interaction(exper,treated_bird)))+geom_point()+stat_smooth(method='lm')
					
					ggplot(inc, aes(x=bout_start_j, y=inc_eff, col=interaction(exper,treated_bird)))+geom_point()+stat_smooth(method='lm')
					ggplot(inc, aes(x=bout_start_j, y=inc_eff, col=interaction(exper,treated_bird)))+geom_point()+stat_smooth(method='lm',se = FALSE)
					ggplot(inc, aes(x=bout_start_j, y=inc_eff, col=interaction(exper,treated_bird)))+geom_point()+stat_smooth(method='lm')+coord_cartesian(ylim = c(0.6,1))
					
					
					ggplot(gap, aes(x=bout_start_j, y=bout_length, col=interaction(exper,treated_bird)))+geom_point()+stat_smooth()
					ggplot(gap, aes(x=bout_start_j, y=bout_length, col=interaction(exper,treated_bird)))+geom_point()+stat_smooth(method='lm')
					ggplot(gapp, aes(x=bout_start_j, y=bout_length, col=interaction(exper,treated_bird)))+geom_point()+stat_smooth(method='lm')
					ggplot(gapp, aes(x=bout_start_j, y=log(bout_length), col=interaction(exper,treated_bird)))+geom_point()+stat_smooth(method='lm')
					ggplot(gap, aes(x=bout_start_j, y=log(bout_length+0.01), col=interaction(exper,treated_bird)))+geom_point()+stat_smooth(method='lm')
					ggplot(x, aes(x=med_bout_start_j, y=med_bout, col=interaction(exper,treated_bird)))+geom_point()+stat_smooth(method='lm')
						
						m=lmer(bout_length~scale(bout_start_j)*exper+(1|nest), gap)
						m_=lmer(log(bout_length+0.01)~scale(bout_start_j)*exper+(1|nest), gap)
						m2=lmer(bout_length~scale(bout_start_j)*exper*treated_bird+(1|nest), gap)
						m2_=lmer(log(bout_length+0.01)~scale(bout_start_j)*exper*treated_bird+(1|nest), gap)
						aic=AIC(m,m_,m2,m2_)
								aic[order(aic$AIC),]
								
							plot(allEffects(m2_))
							summary(m)
							# simulation		
							nsim <- 2000
							bsim <- sim(m2_, n.sim=nsim)  
							apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							fixef(m)
					
				}
			}	
			}
			{# gap
				gap$present=ifelse(gap$bout_length==0,0,1)
				gapy=gap[gap$present==1,]
				
				densityplot(~gap$bout_length)
				densityplot(~log(gap$bout_length+0.1))
				densityplot(~log(gap$bout_length+0.1), groups=gap$sex, auto.key=TRUE)
				densityplot(~log(gap$bout_length+0.1), groups=gap$exper, auto.key=TRUE)
				
				densityplot(~gapy$bout_length)
				densityplot(~log(gapy$bout_length+0.1))
				densityplot(~log(gapy$bout_length+0.1), groups=gapy$sex, auto.key=TRUE)
				densityplot(~log(gapy$bout_length+0.1), groups=gapy$exper, auto.key=TRUE)
				
				ggplot(gap, aes(x=interaction(sex, nest), y=bout_length, col=exper, fill=treated_bird))+geom_boxplot()+ stat_summary(data=inc, aes(x=interaction(sex, nest), y=bout_length, col=exper,fill=treated_bird),fun.y=mean, geom="point")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2))
					
				ggplot(gap, aes(x=nest, y=bout_length, col=exper, fill=treated_bird))+geom_boxplot()+ stat_summary(data=inc, aes(x=nest, y=bout_length, col=exper,fill=treated_bird),fun.y=mean, geom="point")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2))+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_colour_manual(values = c("black", "grey"))
				
				ggplot(gap, aes(x=interaction(sex, nest), y=bout_length, col=exper, fill=treated_bird))+geom_boxplot()+ stat_summary(data=gap, aes(x=interaction(sex, nest), y=bout_length, col=exper,fill=treated_bird),fun.y=mean, geom="point")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2))+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_colour_manual(values = c("black", "grey"))			
				ggplot(gap, aes(x=interaction(sex, nest), y=log(bout_length+0.1), col=exper, fill=treated_bird))+geom_boxplot()+ stat_summary(data=gap, aes(x=interaction(sex, nest), y=log(bout_length+0.1), col=exper,fill=treated_bird),fun.y=mean, geom="point")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2))+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_colour_manual(values = c("black", "grey"))			
				
				ggplot(gap, aes(x=interaction(treated_bird, sex), y=bout_length, col=exper))+geom_boxplot()
				ggplot(gap, aes(col=interaction(exper,treated_bird, sex), y=bout_length, x=bout_start_j_c))+geom_point()+stat_smooth()
				ggplot(gap, aes(col=interaction(exper,treated_bird), y=bout_length, x=bout_start_j_c))+geom_point()+stat_smooth(method='lm')
				ggplot(gap, aes(col=exper, y=bout_length, x=bout_start_j_c))+geom_point()+stat_smooth(method='lm')
				ggplot(gap, aes(x=interaction(treated_bird, sex), y=log(bout_length+0.1), col=exper))+geom_boxplot()
				ggplot(gapp, aes(x=interaction(treated_bird, sex), y=log(bout_length+0.1), col=exper))+geom_boxplot()
				ggplot(gapy, aes(x=interaction(treated_bird, sex), y=log(bout_length+0.1), col=exper))+geom_boxplot()
				ggplot(gap[gap$treated_bird=='y',], aes(x=sex, y=log(bout_length+0.1), col=exper))+geom_boxplot()
				
				x=ddply(gap,.(nest,sex,exper, treated_bird), summarise, med_gap=median(bout_length), num_gaps=sum(present))
					{# gap plot	
						xb=x[x$exper=='b',]
						xa=x[x$exper=='a',]
						xb$med_gap_after=xa$med_gap[match(paste(xb$nest,xb$sex),paste(xa$nest,xa$sex))]
						xb$col_=ifelse(xb$treated_bird=="y", "#FCB42C","#535F7C")
						dev.new(width=2.5,height=2.5)
						#png(paste(out2,"constancy_treated_control_relationship.png", sep=""), width=2.5,height=2.5,units="in",res=600)
						par(mar=c(2,2,0.7,0.7), ps=12, mgp=c(1.2,0.35,0), las=1, cex.lab=0.7, cex.axis=0.6, tcl=-0.15,bty="l",xpd=TRUE, col.axis="grey30", col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt

						plot(log(xb$med_gap_after+0.1)~log(xb$med_gap+0.1), col=xb$col_,bg=adjustcolor( xb$col_, alpha.f = 0.4), pch=21, ylab='Exchange gap - after [min]', xlab=NA, yaxt='n',xaxt='n',xlim=log(c(0,260)+0.1), ylim=log(c(0,260)+0.1))
						lines(log(c(0,260)+0.1),log(c(0,260)+0.1), lty=3, col="lightgrey")
						axis(2, at=log(c(0,1,10,100,260)+0.1), labels=c(0,1,10,100,260))
						axis(1, at=log(c(0,1,10,100,260)+0.1), labels=FALSE)
						text(log(c(0,1,10,100,260)+0.1), par("usr")[3]-0.4, labels = c(0,1,10,100,260),  xpd = TRUE, cex=0.6, col="grey30") #labels = z_g$genus
						mtext('Exchange gap - before [min]',side=1,line=0.6, cex=0.7, las=1, col='grey30')
						
						legend("topright", legend=c('treated bird','control bird'), pch=19, col=c( '#FCB42C','#535F7C'), cex=0.6, bty='n')
					}
					{# gap presence plot - makes little sense
						densityplot(~x$num_gaps)
						xb=x[x$exper=='b',]
						xa=x[x$exper=='a',]
						xb$num_gaps_after=xa$num_gaps[match(paste(xb$nest,xb$sex),paste(xa$nest,xa$sex))]
						xb$col_=ifelse(xb$treated_bird=="y", "#FCB42C","#535F7C")
						dev.new(width=2.5,height=2.5)
						#png(paste(out2,"constancy_treated_control_relationship.png", sep=""), width=2.5,height=2.5,units="in",res=600)
						par(mar=c(2,2,0.7,0.7), ps=12, mgp=c(1.2,0.35,0), las=1, cex.lab=0.7, cex.axis=0.6, tcl=-0.15,bty="l",xpd=TRUE, col.axis="grey30", col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt

						plot(xb$num_gaps_after~xb$num_gaps, col=xb$col_,bg=adjustcolor( xb$col_, alpha.f = 0.4), pch=21, ylab='Exchange gap - after [min]', xlab=NA, xaxt='n',xlim=c(0,3), ylim=c(0,3))
					
						axis(1, at=seq(0,3,by=1), labels=FALSE)
						text(seq(0,3,by=1), par("usr")[3]-0.05, labels = seq(0,3,by=1),  xpd = TRUE, cex=0.6, col="grey30") #labels = z_g$genus
						mtext('bout - before [h]',side=1,line=0.6, cex=0.7, las=1, col='grey30')
						
						legend("bottomleft", legend=c('treated_bird','control_bird'), pch=19, col=c( '#FCB42C','#535F7C'), cex=0.6, bty='n')
					}
				{# model
					gap$int=as.factor(interaction(gap$exper, gap$treated,gap$sex))
							unique(gap$int)
			
							densityplot(~gap$bout_length)
							densityplot(~log(gap$bout_length))
							
							m=glmer(present~exper+(1|nest),gap, family='binomial')
							m1=glmer(present~exper*treated_bird+(1|nest),gap, family='binomial')
							m2=glmer(present~exper*sex+(1|nest),gap, family='binomial')
							m3=glmer(present~int+(1|nest),gap, family='binomial')
								aic=AIC(m,m1,m2,m3)
								aic[order(aic$AIC),]
									plot(allEffects(m))
									plot(allEffects(m3))
							plot(allEffects(m))
							summary(m)
							summary(m3)
							# simulation		
							nsim <- 2000
							bsim <- sim(msw, n.sim=nsim)  
							apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							fixef(m)
							
							m=lmer(bout_length~exper+(1|nest),gap)
							m_=lmer(log(bout_length+0.01)~exper+(1|nest),gap)
							m1=lmer(bout_length~exper*treated_bird+(1|nest),gap)
							m2=lmer(bout_length~exper*sex+(1|nest),gap)
							m3=lmer(bout_length~exper*treated_bird*sex+(1|nest),gap)
							m4=lmer(bout_length~int+(1|nest),gap)
							m5=lmer(log(bout_length+0.01)~int+(1|nest),gap)
								aic=AIC(m,m_,m1,m2,m3,m4,m5)
								aic[order(aic$AIC),]
								aic=AIC(m,m1,m2,m3,m4)
								aic[order(aic$AIC),]
								
							plot(allEffects(m))
							plot(allEffects(m_))
							summary(m)
							# simulation		
							nsim <- 2000
							bsim <- sim(m, n.sim=nsim)  
							apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							fixef(m)
							
							plot(allEffects(m3))
							plot(allEffects(m4))
							plot(allEffects(m5))
							# simulation		
							nsim <- 2000
							bsim <- sim(m3, n.sim=nsim)  
							apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							fixef(m3)
							
							summary(glht(m5, linfct = mcp(int = c(	"b.y.m - a.y.m = 0",
																	"b.y.f - a.y.f = 0", #
																	"b.n.m - a.n.m = 0", 
																	"b.n.f - a.n.f = 0",
																	
																	"a.n.m - a.n.f = 0",
																	"a.y.m - a.y.f = 0",
																	
																	"a.y.m - a.n.m = 0",
																	"a.y.f - a.n.f = 0"#
																))))
							
				
				}
			}
			{# after bout of captive bird in relation to mass loss
					
					g=dbGetQuery(conMy, "select*from barrow_2013.captivity")	
					gg=ddply(g[!is.na(g$mass),],.(ring_num),summarise,mass_loss=(mass[phase=='start']-mass[phase=='end']), mass_end=mass[phase=='end'])
							ggplot(gg,aes(x=mass_loss))+geom_density()
							ggplot(gg,aes(x=mass_loss, y=mass_end))+geom_point()+stat_smooth()
					a=inc[inc$exper=='a' & inc$treated_bird=='n',]
					a$mass_loss=gg$mass_loss[match(a$bird_ID_filled, gg$ring_num)]
					x=ddply(a,.(nest,sex), summarise, med_bout=median(bout_length), med_eff=median(inc_eff), mass_loss=median(mass_loss))
					y=ddply(a,.(nest,sex), summarise, first_bout=bout_length[bout_start==min(bout_start)], first_eff=inc_eff[bout_start==min(bout_start)], mass_loss=median(mass_loss))
				  {# bout				
					# three bouts
						ggplot(a, aes(x=nest, y=mass_loss, fill=sex))+geom_point()
						ggplot(a, aes(x=mass_loss, y=bout_length))+geom_point()
						ggplot(a, aes(x=mass_loss, y=bout_length))+geom_point()+stat_smooth(method='lm')
							m=lmer(bout_length~mass_loss+(1|nest),a)
							nsim <- 2000
							bsim <- sim(m, n.sim=nsim)  
							apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							summary(glht(m))
							plot(allEffects(m))
						
						ggplot(a, aes(x=mass_loss, y=bout_length, col=sex))+geom_point()
						ggplot(a, aes(x=mass_loss, y=bout_length, col=sex))+geom_point()+stat_smooth(method='lm')
							ms=lmer(bout_length~mass_loss*sex+(1|nest),a)
							nsim <- 2000
							bsim <- sim(ms, n.sim=nsim)  
							apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							summary(glht(ms))
							plot(allEffects(ms))
							AIC(m,ms)
						
					# median bout
						ggplot(x, aes(x=mass_loss, y=med_bout))+geom_point()+stat_smooth(method='lm')
								m=lm(med_bout~mass_loss,x)
								nsim <- 2000
								bsim <- sim(m, n.sim=nsim)  
								apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
								summary(glht(m))
								plot(allEffects(m))
						
						ggplot(x, aes(x=mass_loss, y=med_bout, col=sex))+geom_point()+stat_smooth(method='lm')
							ms=lm(med_bout~mass_loss*sex,x)
							nsim <- 2000
							bsim <- sim(ms, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
							summary(glht(ms))
							plot(allEffects(ms))
							AIC(m,ms)
							
					# first bout
						ggplot(y, aes(x=mass_loss, y=first_bout))+geom_point()+stat_smooth(method='lm')
								m=lm(first_bout~mass_loss,y)
								nsim <- 2000
								bsim <- sim(m, n.sim=nsim)  
								apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
								summary(glht(m))
								plot(allEffects(m))
						
						ggplot(y, aes(x=mass_loss, y=first_bout, col=sex))+geom_point()+stat_smooth(method='lm')
							ms=lm(first_bout~mass_loss*sex,y)
							nsim <- 2000
							bsim <- sim(ms, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
							summary(glht(ms))
							plot(allEffects(ms))
							AIC(m,ms)
				}				
				  {# efficiency	
					######CHECK what is going on here
						
								a[a$inc_eff<0.1,]

								
					densityplot(~a$inc_eff)
					densityplot(~a$inc_eff)
					# three bouts
						ggplot(a, aes(x=mass_loss, y=inc_eff))+geom_point()
						ggplot(a, aes(x=mass_loss, y=inc_eff))+geom_point()+stat_smooth(method='lm')
							m=lmer(inc_eff~scale(bout_length)+scale(mass_loss)+(1|nest),a)
							nsim <- 2000
							bsim <- sim(m, n.sim=nsim)  
							apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							summary(glht(m))
							plot(allEffects(m))
						
						ggplot(a, aes(x=mass_loss, y=inc_eff, col=sex))+geom_point()
						ggplot(a, aes(x=mass_loss, y=inc_eff, col=sex))+geom_point()+stat_smooth(method='lm')
							ms=lmer(inc_eff~scale(bout_length)+scale(mass_loss)*sex+(1|nest),a)
							ms=lmer(bout_length~scale(mass_loss)*sex+(1|nest),a)
							nsim <- 2000
							bsim <- sim(ms, n.sim=nsim)  
							apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							summary(glht(ms))
							plot(allEffects(ms))
							AIC(m,ms)
						
					# median bout
						ggplot(x, aes(x=mass_loss, y=med_bout))+geom_point()+stat_smooth(method='lm')
								m=lm(med_bout~mass_loss,x)
								nsim <- 2000
								bsim <- sim(m, n.sim=nsim)  
								apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
								summary(glht(m))
								plot(allEffects(m))
						
						ggplot(x, aes(x=mass_loss, y=med_bout, col=sex))+geom_point()+stat_smooth(method='lm')
							ms=lm(med_bout~mass_loss*sex,x)
							nsim <- 2000
							bsim <- sim(ms, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
							summary(glht(ms))
							plot(allEffects(ms))
							AIC(m,ms)
							
					# first bout
						ggplot(y, aes(x=mass_loss, y=first_bout))+geom_point()+stat_smooth(method='lm')
								m=lm(first_bout~mass_loss,y)
								nsim <- 2000
								bsim <- sim(m, n.sim=nsim)  
								apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
								summary(glht(m))
								plot(allEffects(m))
						
						ggplot(y, aes(x=mass_loss, y=first_bout, col=sex))+geom_point()+stat_smooth(method='lm')
							ms=lm(first_bout~mass_loss*sex,y)
							nsim <- 2000
							bsim <- sim(ms, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
							summary(glht(ms))
							plot(allEffects(ms))
							AIC(m,ms)
				}				
				}	
			{# after bout constancy/length of treated bird in relation to its compensation
					load(file=paste(wd,'experimental.Rdata', sep=""))
					bb=b[b$exper=='c',]
					bt=b[b$exper=='t',]
					bb$inc_eff_treated=bt$inc_eff[match(bb$nest, bt$nest)]
					bb$compensation=bb$inc_eff_treated/bb$inc_eff
									
					a=inc[inc$exper=='a' & inc$treated_bird=='y',]
					a$compensation=bb$compensation[match(a$nest, bb$nest)]
					x=ddply(a,.(nest,sex), summarise, med_bout=median(bout_length), mean_bout=mean(bout_length),med_eff=median(inc_eff), mean_eff=mean(inc_eff), compensation=median(compensation))
					y=ddply(a,.(nest,sex), summarise, first_bout=bout_length[bout_start==min(bout_start)], first_eff=inc_eff[bout_start==min(bout_start)], compensation=median(compensation))
					
					densityplot(~a$compensation)
					densityplot(~a$bout_length)
				  {# bout				
					# three bouts
						ggplot(a, aes(x=compensation, y=bout_length))+geom_point()
						ggplot(a, aes(x=compensation, y=bout_length))+geom_point()+stat_smooth(method='lm')
						
						
							m=lmer(bout_length~scale(compensation)+(1|nest),a)
							nsim <- 2000
							bsim <- sim(m, n.sim=nsim)  
							apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							summary(glht(m))
							plot(allEffects(m))
										
						ggplot(a, aes(x=compensation, y=bout_length, col=sex))+geom_point()
						ggplot(a, aes(x=compensation, y=bout_length, col=sex))+geom_point()+stat_smooth(method='lm')
							ms=lmer(bout_length~scale(compensation)*sex+(1|nest),a)
							nsim <- 2000
							bsim <- sim(ms, n.sim=nsim)  
							apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							summary(ms)
							summary(glht(ms))
							plot(allEffects(ms))
							AIC(m,ms)
						
					# median bout
						ggplot(x, aes(x=compensation, y=med_bout))+geom_point()+stat_smooth(method='lm')
								m=lm(med_bout~scale(compensation),x)
								nsim <- 2000
								bsim <- sim(m, n.sim=nsim)  
								apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
								summary(glht(m))
								plot(allEffects(m))
						
						ggplot(x, aes(x=compensation, y=med_bout, col=sex))+geom_point()+stat_smooth(method='lm')
							ms=lm(med_bout~compensation*sex,x)
							nsim <- 2000
							bsim <- sim(ms, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
							summary(glht(ms))
							plot(allEffects(ms))
							AIC(m,ms)
							
					# mean bout
						ggplot(x, aes(x=compensation, y=mean_bout))+geom_point()+stat_smooth(method='lm')
								m=lm(mean_bout~scale(compensation),x)
								nsim <- 2000
								bsim <- sim(m, n.sim=nsim)  
								apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
								summary(glht(m))
								plot(allEffects(m))
						
						ggplot(x, aes(x=compensation, y=mean_bout, col=sex))+geom_point()+stat_smooth(method='lm')
							ms=lm(mean_bout~compensation*sex,x)
							nsim <- 2000
							bsim <- sim(ms, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
							summary(glht(ms))
							plot(allEffects(ms))
							AIC(m,ms)
							
					# first bout
						ggplot(y, aes(x=compensation, y=first_bout))+geom_point()+stat_smooth(method='lm')
								m=lm(first_bout~compensation,y)
								nsim <- 2000
								bsim <- sim(m, n.sim=nsim)  
								apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
								summary(glht(m))
								plot(allEffects(m))
						
						ggplot(y, aes(x=compensation, y=first_bout, col=sex))+geom_point()+stat_smooth(method='lm')
							ms=lm(first_bout~compensation*sex,y)
							nsim <- 2000
							bsim <- sim(ms, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
							summary(glht(ms))
							plot(allEffects(ms))
							AIC(m,ms)
				}				
				  {# efficiency	
						a[a$inc_eff<0.1,]
						densityplot(~a$inc_eff)
					
					# three bouts
						ggplot(a, aes(x=compensation, y=inc_eff))+geom_point()
						ggplot(a, aes(x=compensation, y=inc_eff))+geom_point()+stat_smooth(method='lm')
							m=lmer(inc_eff~scale(bout_length)+compensation+(1|nest),a)
							nsim <- 2000
							bsim <- sim(m, n.sim=nsim)  
							apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							summary(glht(m))
							plot(allEffects(m))
						
						ggplot(a, aes(x=compensation, y=inc_eff, col=sex))+geom_point()
						ggplot(a, aes(x=compensation, y=inc_eff, col=sex))+geom_point()+stat_smooth(method='lm')
							ms=lmer(inc_eff~scale(bout_length)+compensation*sex+(1|nest),a)
							nsim <- 2000
							bsim <- sim(ms, n.sim=nsim)  
							apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							summary(glht(ms))
							plot(allEffects(ms))
							AIC(m,ms)
						
					# med_eff
						ggplot(x, aes(x=compensation, y=med_eff))+geom_point()+stat_smooth(method='lm')
								m=lm(med_eff~scale(med_bout)+compensation,x)
								nsim <- 2000
								bsim <- sim(m, n.sim=nsim)  
								apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
								summary(glht(m))
								plot(allEffects(m))
						
						ggplot(x, aes(x=compensation, y=med_eff, col=sex))+geom_point()+stat_smooth(method='lm')
							ms=lm(med_eff~scale(med_bout)+compensation*sex,x)
							nsim <- 2000
							bsim <- sim(ms, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
							summary(glht(ms))
							plot(allEffects(ms))
							AIC(m,ms)
					# mean bout
						ggplot(x, aes(x=compensation, y=mean_eff))+geom_point()+stat_smooth(method='lm')
								m=lm(mean_eff~scale(med_bout)+compensation,x)
								nsim <- 2000
								bsim <- sim(m, n.sim=nsim)  
								apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
								summary(glht(m))
								plot(allEffects(m))
						
						ggplot(x, aes(x=compensation, y=mean_eff, col=sex))+geom_point()+stat_smooth(method='lm')
							ms=lm(mean_eff~scale(med_bout)+compensation*sex,x)
							nsim <- 2000
							bsim <- sim(ms, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
							summary(glht(ms))
							plot(allEffects(ms))
							AIC(m,ms)
									
					# first bout
						ggplot(y, aes(x=compensation, y=first_eff))+geom_point()+stat_smooth(method='lm')
								m=lm(first_eff~scale(first_bout)+compensation,y)
								nsim <- 2000
								bsim <- sim(m, n.sim=nsim)  
								apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
								summary(glht(m))
								plot(allEffects(m))
						
						ggplot(y, aes(x=compensation, y=first_eff, col=sex))+geom_point()+stat_smooth(method='lm')
							ms=lm(first_eff~scale(first_bout)+compensation*sex,y)
							nsim <- 2000
							bsim <- sim(ms, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
							summary(glht(ms))
							plot(allEffects(ms))
							AIC(m,ms)
				}				
				}	
			}
		}
	
		{# Supplementary Table 6 - after bout of captive bird in relation to mass loss
			{# run first - add mass loss to data
				u=read.csv(paste(wd,'experiment_metadata.csv', sep=""), stringsAsFactors=FALSE)	
					#uu = u[9:nrow(u),]
					#summary(-uu$mass_loss)
				a=inc[inc$exper=='a' & inc$treated_bird=='n',]
				a$mass_loss=-u$mass_loss[match(a$bird_ID_filled, u$ID_taken)]
				a$bird_ID=as.factor(a$bird_ID)
				x=ddply(a,.(nest,sex), summarise, med_bout=median(bout_length), med_eff=median(inc_eff), mass_loss=median(mass_loss))
				y=ddply(a,.(nest,sex), summarise, first_bout=bout_length[bout_start==min(bout_start)], first_eff=inc_eff[bout_start==min(bout_start)], mass_loss=median(mass_loss))
			 }
		
			{# constancy
				{# 1 - simple model - convergence problem, because of singularity, we should remove the randome slope, but since the results do not change and leaving the slope in makes comparioson with other models easier, we keep it in
				 m=lmer(inc_eff~scale(mass_loss)+(scale(mass_loss)|bird_ID),a, REML=FALSE)
				 pred=c('Intercept','Mass_loss')
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
				# Fixed effects
					v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))
					ci=apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
					oi=data.frame(model='1',response='constancy',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
					rownames(oi) = NULL
						oi$estimate_r=round(oi$estimate,3)
						oi$lwr_r=round(oi$lwr,3)
						oi$upr_r=round(oi$upr,3)
						#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
					oii=oi[c('model','response','type',"effect", "estimate_r","lwr_r",'upr_r')]	
				# Random effects
					l=data.frame(summary(m)$varcor)
						l=l[is.na(l$var2),]
						l$var2=NULL
						ri=data.frame(model='1',response='constancy',type='random (var)',effect=l$var, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA, stringsAsFactors=FALSE)
							ri$effect[nrow(ri)]='residual'
					o1=rbind(oii,ri)
			
			}
				{# 2 - with sex
				m2=lmer(inc_eff~scale(mass_loss)*sex+(mass_loss|bird_ID),a, REML=FALSE)
				pred=c('Intercept (female)','Mass_loss','Sex (male)','Mass:Sex')
						nsim <- 2000
						bsim <- sim(m2, n.sim=nsim)  
				# Fixed effects
					v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))
					ci=apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
					oi=data.frame(model='1',response='constancy',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
					rownames(oi) = NULL
						oi$estimate_r=round(oi$estimate,3)
						oi$lwr_r=round(oi$lwr,3)
						oi$upr_r=round(oi$upr,3)
						#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
					oii=oi[c('model','response','type',"effect", "estimate_r","lwr_r",'upr_r')]	
				# Random effects
					l=data.frame(summary(m2)$varcor)
						l=l[is.na(l$var2),]
						l$var2=NULL
						ri=data.frame(model='1',response='constancy',type='random (var)',effect=l$var, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA, stringsAsFactors=FALSE)
						ri$effect[nrow(ri)]='residual'
					o2=rbind(oii,ri)
			
						
			}	
					o_con=rbind(o2,o1)
				{# AICc # number of observations = number of birds*2 (2 = before and after period) length(unique(paste(inc$bird_ID,inc$exper)))
						o=data.frame(response='constancy',model=c('simple','with sex'), AIC=c(AICc(m, nobs=18),AICc(m2,nobs=18)))
						o$delta=o$AIC-min(o$AIC)
						o$prob=exp(-0.5*o$delta)/sum(exp(-0.5*o$delta))
						o$ER=max(o$prob)/o$prob
						o$AIC=round(o$AIC,1)
						o$delta=round(o$delta,2)
						o$prob=round(o$prob,3)
						o$ER=round(o$ER,2)
						o_1=o
						#o[order(o$delta),]
					}
			}
			{# bout
				{# 1 - simple model 
				m=lmer(bout_length~scale(mass_loss)+(mass_loss|nest),a, REML=FALSE)
				pred=c('Intercept','Mass_loss')
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
				# Fixed effects
					v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))
					ci=apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
					oi=data.frame(model='1',response='constancy',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
					rownames(oi) = NULL
						oi$estimate_r=round(oi$estimate,3)
						oi$lwr_r=round(oi$lwr,3)
						oi$upr_r=round(oi$upr,3)
						#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
					oii=oi[c('model','response','type',"effect", "estimate_r","lwr_r",'upr_r')]	
				# Random effects
					l=data.frame(summary(m)$varcor)
						l=l[is.na(l$var2),]
						l$var2=NULL
						ri=data.frame(model='1',response='constancy',type='random (var)',effect=l$var, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA, stringsAsFactors=FALSE)
							ri$effect[nrow(ri)]='residual'
					o1=rbind(oii,ri)
			
			}
				{# 2 - with sex
			 m2=lmer(bout_length~scale(mass_loss)*sex+(mass_loss|nest),a, REML=FALSE)
			 pred=c('Intercept (female)','Mass_loss','Sex (male)','Mass:Sex')
						nsim <- 2000
						bsim <- sim(m2, n.sim=nsim)  
				# Fixed effects
					v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))
					ci=apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
										oi=data.frame(model='1',response='constancy',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
					rownames(oi) = NULL
						oi$estimate_r=round(oi$estimate,3)
						oi$lwr_r=round(oi$lwr,3)
						oi$upr_r=round(oi$upr,3)
						#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
					oii=oi[c('model','response','type',"effect", "estimate_r","lwr_r",'upr_r')]	
				# Random effects
					l=data.frame(summary(m2)$varcor)
						l=l[is.na(l$var2),]
						l$var2=NULL
						ri=data.frame(model='1',response='constancy',type='random (var)',effect=l$var, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA, stringsAsFactors=FALSE)
							ri$effect[nrow(ri)]='residual'
					o2=rbind(oii,ri)
			
						
			}	
					o_bout=rbind(o2,o1)
				{# AICc # number of observations = number of birds*2 (2 = before and after period) length(unique(paste(inc$bird_ID,inc$exper)))
						o=data.frame(response='constancy',model=c('simple','with sex'), AIC=c(AICc(m, nobs=18),AICc(m2,nobs=18)))
						o$delta=o$AIC-min(o$AIC)
						o$prob=exp(-0.5*o$delta)/sum(exp(-0.5*o$delta))
						o$ER=max(o$prob)/o$prob
						o$AIC=round(o$AIC,1)
						o$delta=round(o$delta,2)
						o$prob=round(o$prob,3)
						o$ER=round(o$ER,2)
						o_2=o
						#o[order(o$delta),]
					}
			}			
			{# combine and export to excel table
							sname = tempfile(fileext='.xls')
							wb = loadWorkbook(sname,create = TRUE)	
							createSheet(wb, name = "output")
							writeWorksheet(wb, rbind(o_con,o_bout), sheet = "output")
							createSheet(wb, name = "output_AIC")
							writeWorksheet(wb, rbind(o_1,o_2), sheet = "output_AIC")
							saveWorkbook(wb)
							shell(sname)
			}
		
			{# model assumptions
				{# constancy
					{# 1
						dev.new(width=6,height=9)
									
						m=lmer(inc_eff~scale(mass_loss)+(scale(mass_loss)|bird_ID),a, REML=FALSE) #convergence problem, because of singularity, we should remove the randome slope, but since the results do not change and leaving the slope in makes comparioson with other models easier, we keep it in
						par(mfrow=c(4,3))
						
						scatter.smooth(fitted(m),resid(m),col='red');abline(h=0, lty=2)
						scatter.smooth(fitted(m),sqrt(abs(resid(m))), col='red')
													
						qqnorm(resid(m), main=list("Normal Q-Q Plot: residuals", cex=0.8),col='red') 
						qqline(resid(m))
			 
						qqnorm(unlist(ranef(m)$bird_ID[1]), main = " Bird")
						qqline(unlist(ranef(m)$bird_ID[1]))
						
						qqnorm(unlist(ranef(m)$bird_ID[2]), main = " Slope")
						qqline(unlist(ranef(m)$bird_ID[2]))
						
						scatter.smooth(resid(m)~scale(a$mass_loss));abline(h=0, lty=2, col='red')
											
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=a$lon, y=a$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
					}
					{# 2
						dev.new(width=6,height=9)
									
						m2=lmer(inc_eff~scale(mass_loss)*sex+(scale(mass_loss)|bird_ID),a, REML=FALSE)
						par(mfrow=c(4,3))
						
						scatter.smooth(fitted(m),resid(m),col='red');abline(h=0, lty=2)
						scatter.smooth(fitted(m),sqrt(abs(resid(m))), col='red')
													
						qqnorm(resid(m), main=list("Normal Q-Q Plot: residuals", cex=0.8),col='red') 
						qqline(resid(m))
			 
						qqnorm(unlist(ranef(m)$bird_ID[1]), main = " Bird")
						qqline(unlist(ranef(m)$bird_ID[1]))
						
						scatter.smooth(resid(m)~scale(a$mass_loss));abline(h=0, lty=2, col='red')
						plot(resid(m)~factor(a$sex)); abline(h=0, lty=2, col='red')
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=a$lon, y=a$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
					}
				}	
				{# bout
					{# 1
						dev.new(width=6,height=9)
									
						m=lmer(bout_length~scale(mass_loss)+(mass_loss|bird_ID),a, REML=FALSE)
						par(mfrow=c(4,3))
						
						scatter.smooth(fitted(m),resid(m),col='red');abline(h=0, lty=2)
						scatter.smooth(fitted(m),sqrt(abs(resid(m))), col='red')
													
						qqnorm(resid(m), main=list("Normal Q-Q Plot: residuals", cex=0.8),col='red') 
						qqline(resid(m))
			 
						qqnorm(unlist(ranef(m)$bird_ID[1]), main = " Bird")
						qqline(unlist(ranef(m)$bird_ID[1]))
						
						qqnorm(unlist(ranef(m)$bird_ID[2]), main = " Slope")
						qqline(unlist(ranef(m)$bird_ID[2]))
						
						scatter.smooth(resid(m)~scale(a$mass_loss));abline(h=0, lty=2, col='red')
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=a$lon, y=a$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
					}
					{# 2
						dev.new(width=6,height=9)
									
						m2=lmer(bout_length~scale(mass_loss)*sex+(mass_loss|bird_ID),a, REML=FALSE)
						par(mfrow=c(4,3))
						
						scatter.smooth(fitted(m),resid(m),col='red');abline(h=0, lty=2)
						scatter.smooth(fitted(m),sqrt(abs(resid(m))), col='red')
													
						qqnorm(resid(m), main=list("Normal Q-Q Plot: residuals", cex=0.8),col='red') 
						qqline(resid(m))
			 
						qqnorm(unlist(ranef(m)$bird_ID[1]), main = " Bird")
						qqline(unlist(ranef(m)$bird_ID[1]))
						
						scatter.smooth(resid(m)~scale(a$mass_loss));abline(h=0, lty=2, col='red')
						plot(resid(m)~factor(a$sex)); abline(h=0, lty=2, col='red')
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=a$lon, y=a$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
					}
				}
			}					
				
		}	
		{# Supplementary Table 7 - after bouts of focal bird in relation to amount of compensation
			{# run first - add compensatin to data
				bb$compensation=bb$inc_eff_treated/bb$inc_eff
								
				a=inc[inc$exper=='a' & inc$treated_bird=='y',]
				a$compensation=bb$compensation[match(a$nest, bb$nest)]
				x=ddply(a,.(nest,sex), summarise, med_bout=median(bout_length), mean_bout=mean(bout_length),med_eff=median(inc_eff), mean_eff=mean(inc_eff), compensation=median(compensation))
				y=ddply(a,.(nest,sex), summarise, first_bout=bout_length[bout_start==min(bout_start)], first_eff=inc_eff[bout_start==min(bout_start)], compensation=median(compensation))
			}					
			
			{# constancy
				{# 1 - simple model 
				 m=lmer(inc_eff~scale(compensation)+(scale(compensation)|bird_ID),a, REML=FALSE)
				 pred=c('Intercept','compensation')
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
				# Fixed effects
					v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))
					ci=apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
					oi=data.frame(model='1',response='constancy',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
					rownames(oi) = NULL
						oi$estimate_r=round(oi$estimate,3)
						oi$lwr_r=round(oi$lwr,3)
						oi$upr_r=round(oi$upr,3)
						#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
					oii=oi[c('model','response','type',"effect", "estimate_r","lwr_r",'upr_r')]	
				# Random effects
					l=data.frame(summary(m)$varcor)
						l=l[is.na(l$var2),]
						l$var2=NULL
						ri=data.frame(model='1',response='constancy',type='random (var)',effect=l$var, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA, stringsAsFactors=FALSE)
							ri$effect[nrow(ri)]='residual'
					o1=rbind(oii,ri)
			
			}
				{# 2 - with sex
				m2=lmer(inc_eff~scale(compensation)*sex+(scale(compensation)|bird_ID),a, REML=FALSE)
				pred=c('Intercept (female)','compensation','Sex (male)','Mass:Sex')
						nsim <- 2000
						bsim <- sim(m2, n.sim=nsim)  
				# Fixed effects
					v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))
					ci=apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
					oi=data.frame(model='sex',response='constancy',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
					rownames(oi) = NULL
						oi$estimate_r=round(oi$estimate,3)
						oi$lwr_r=round(oi$lwr,3)
						oi$upr_r=round(oi$upr,3)
						#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
					oii=oi[c('model','response','type',"effect", "estimate_r","lwr_r",'upr_r')]	
				# Random effects
					l=data.frame(summary(m2)$varcor)
						l=l[is.na(l$var2),]
						l$var2=NULL
						ri=data.frame(model='sex',response='constancy',type='random (var)',effect=l$var, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA, stringsAsFactors=FALSE)
						ri$effect[nrow(ri)]='residual'
					o2=rbind(oii,ri)
			
						
			}	
					o_con=rbind(o2,o1)
				{# AICc # number of observations = number of birds*2 (2 = before and after period) length(unique(paste(inc$bird_ID,inc$exper)))
						o=data.frame(response='constancy',model=c('simple','with sex'), AIC=c(AICc(m, nobs=18),AICc(m2,nobs=18)))
						o$delta=o$AIC-min(o$AIC)
						o$prob=exp(-0.5*o$delta)/sum(exp(-0.5*o$delta))
						o$ER=max(o$prob)/o$prob
						o$AIC=round(o$AIC,1)
						o$delta=round(o$delta,2)
						o$prob=round(o$prob,3)
						o$ER=round(o$ER,2)
						o_1=o
						#o[order(o$delta),]
					}
			}
			{# bout
				{# 1 - simple model 
				m=lmer(bout_length~scale(compensation)+(scale(compensation)|nest),a, REML=FALSE)
				pred=c('Intercept','compensation')
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
				# Fixed effects
					v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))
					ci=apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
					oi=data.frame(model='1',response='constancy',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
					rownames(oi) = NULL
						oi$estimate_r=round(oi$estimate,3)
						oi$lwr_r=round(oi$lwr,3)
						oi$upr_r=round(oi$upr,3)
						#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
					oii=oi[c('model','response','type',"effect", "estimate_r","lwr_r",'upr_r')]	
				# Random effects
					l=data.frame(summary(m)$varcor)
						l=l[is.na(l$var2),]
						l$var2=NULL
						ri=data.frame(model='1',response='constancy',type='random (var)',effect=l$var, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA, stringsAsFactors=FALSE)
							ri$effect[nrow(ri)]='residual'
					o1=rbind(oii,ri)
			
			}
				{# 2 - with sex
			 m2=lmer(bout_length~scale(compensation)*sex+(scale(compensation)|nest),a, REML=FALSE)
			 pred=c('Intercept (female)','compensation','Sex (male)','Mass:Sex')
						nsim <- 2000
						bsim <- sim(m2, n.sim=nsim)  
				# Fixed effects
					v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))
					ci=apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
										oi=data.frame(model='sex',response='constancy',type='fixed',effect=pred,estimate=v, lwr=ci[1,], upr=ci[2,])
					rownames(oi) = NULL
						oi$estimate_r=round(oi$estimate,3)
						oi$lwr_r=round(oi$lwr,3)
						oi$upr_r=round(oi$upr,3)
						#oi$CI=paste("(", oi$lwr_r, "-", oi$upr_r, ")", sep = "", collapse = NULL)
					oii=oi[c('model','response','type',"effect", "estimate_r","lwr_r",'upr_r')]	
				# Random effects
					l=data.frame(summary(m2)$varcor)
						l=l[is.na(l$var2),]
						l$var2=NULL
						ri=data.frame(model='sex',response='constancy',type='random (var)',effect=l$var, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA, stringsAsFactors=FALSE)
							ri$effect[nrow(ri)]='residual'
					o2=rbind(oii,ri)
			
						
			}	
					o_bout=rbind(o2,o1)
				{# AICc # number of observations = number of birds*2 (2 = before and after period) length(unique(paste(inc$bird_ID,inc$exper)))
						o=data.frame(response='constancy',model=c('simple','with sex'), AIC=c(AICc(m, nobs=18),AICc(m2,nobs=18)))
						o$delta=o$AIC-min(o$AIC)
						o$prob=exp(-0.5*o$delta)/sum(exp(-0.5*o$delta))
						o$ER=max(o$prob)/o$prob
						o$AIC=round(o$AIC,1)
						o$delta=round(o$delta,2)
						o$prob=round(o$prob,3)
						o$ER=round(o$ER,2)
						o_2=o
						#o[order(o$delta),]
					}
			}			
			{# combine and export to excel table
							sname = tempfile(fileext='.xls')
							wb = loadWorkbook(sname,create = TRUE)	
							createSheet(wb, name = "output")
							writeWorksheet(wb, rbind(o_con,o_bout), sheet = "output")
							createSheet(wb, name = "output_AIC")
							writeWorksheet(wb, rbind(o_1,o_2), sheet = "output_AIC")
							saveWorkbook(wb)
							shell(sname)
			}
		
			{# model assumptions
				{# constancy
					{# 1
						dev.new(width=6,height=9)
									
						 m=lmer(inc_eff~scale(compensation)+(scale(compensation)|bird_ID),a, REML=FALSE)
						par(mfrow=c(4,3))
						
						scatter.smooth(fitted(m),resid(m),col='red');abline(h=0, lty=2)
						scatter.smooth(fitted(m),sqrt(abs(resid(m))), col='red')
													
						qqnorm(resid(m), main=list("Normal Q-Q Plot: residuals", cex=0.8),col='red') 
						qqline(resid(m))
			 
						qqnorm(unlist(ranef(m)$bird_ID[1]), main = " Bird")
						qqline(unlist(ranef(m)$bird_ID[1]))
						
						qqnorm(unlist(ranef(m)$bird_ID[2]), main = " Slope")
						qqline(unlist(ranef(m)$bird_ID[2]))
						
						scatter.smooth(resid(m)~a$compensation);abline(h=0, lty=2, col='red')
												
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=a$lon, y=a$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
					}
					{# 2
						dev.new(width=6,height=9)
									
						m=lmer(inc_eff~scale(compensation)*sex+(scale(compensation)|bird_ID),a, REML=FALSE)
						par(mfrow=c(4,3))
						
						scatter.smooth(fitted(m),resid(m),col='red');abline(h=0, lty=2)
						scatter.smooth(fitted(m),sqrt(abs(resid(m))), col='red')
													
						qqnorm(resid(m), main=list("Normal Q-Q Plot: residuals", cex=0.8),col='red') 
						qqline(resid(m))
			 
						qqnorm(unlist(ranef(m)$bird_ID[1]), main = " Bird")
						qqline(unlist(ranef(m)$bird_ID[1]))
						
						scatter.smooth(resid(m)~a$compensation);abline(h=0, lty=2, col='red')
						plot(resid(m)~factor(a$sex)); abline(h=0, lty=2, col='red')
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=a$lon, y=a$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
					}
				}	
				{# bout
					{# 1
						dev.new(width=6,height=9)
									
						m=lmer(bout_length~scale(compensation)+(compensation|bird_ID),a, REML=FALSE)
						par(mfrow=c(4,3))
						
						scatter.smooth(fitted(m),resid(m),col='red');abline(h=0, lty=2)
						scatter.smooth(fitted(m),sqrt(abs(resid(m))), col='red')
													
						qqnorm(resid(m), main=list("Normal Q-Q Plot: residuals", cex=0.8),col='red') 
						qqline(resid(m))
			 
						qqnorm(unlist(ranef(m)$bird_ID[1]), main = " Bird")
						qqline(unlist(ranef(m)$bird_ID[1]))
						
						qqnorm(unlist(ranef(m)$bird_ID[2]), main = " Slope")
						qqline(unlist(ranef(m)$bird_ID[2]))
						
						scatter.smooth(resid(m)~a$compensation);abline(h=0, lty=2, col='red')
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=a$lon, y=a$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
					}
					{# 2
						dev.new(width=6,height=9)
									
						m=lmer(bout_length~scale(compensation)*sex+(compensation|nest),a, REML=FALSE)
						par(mfrow=c(4,3))
						
						scatter.smooth(fitted(m),resid(m),col='red');abline(h=0, lty=2)
						scatter.smooth(fitted(m),sqrt(abs(resid(m))), col='red')
													
						qqnorm(resid(m), main=list("Normal Q-Q Plot: residuals", cex=0.8),col='red') 
						qqline(resid(m))
			 
						qqnorm(unlist(ranef(m)$bird_ID[1]), main = " Bird")
						qqline(unlist(ranef(m)$bird_ID[1]))
						
						scatter.smooth(resid(m)~a$compensation);abline(h=0, lty=2, col='red')
						plot(resid(m)~factor(a$sex)); abline(h=0, lty=2, col='red')
						
						acf(resid(m), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
						
						# spatial autocorrelations - nest location
							spdata=data.frame(resid=resid(m), x=a$lon, y=a$lat)
								spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
								#cex_=c(1,2,3,3.5,4)
								cex_=c(1,1.5,2,2.5,3)
								spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
								plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
								
								plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
								plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
					}
				}
			}
		}	
			{# not used
			{# bout				
				# three bouts
					ggplot(a, aes(x=nest, y=mass_loss, fill=sex))+geom_point()
					ggplot(a, aes(x=mass_loss, y=bout_length))+geom_point()
					ggplot(a, aes(x=mass_loss, y=bout_length))+geom_point()+stat_smooth(method='lm')
						m=lmer(bout_length~mass_loss+(1|nest),a)
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
						apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
						summary(glht(m))
						plot(allEffects(m))
					
					ggplot(a, aes(x=mass_loss, y=bout_length, col=sex))+geom_point()
					ggplot(a, aes(x=mass_loss, y=bout_length, col=sex))+geom_point()+stat_smooth(method='lm')
						ms=lmer(bout_length~mass_loss*sex+(1|nest),a)
						nsim <- 2000
						bsim <- sim(ms, n.sim=nsim)  
						apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
						summary(glht(ms))
						plot(allEffects(ms))
						AIC(m,ms)
					
				# median bout
					ggplot(x, aes(x=mass_loss, y=med_bout))+geom_point()+stat_smooth(method='lm')
							m=lm(med_bout~mass_loss,x)
							nsim <- 2000
							bsim <- sim(m, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
							summary(glht(m))
							plot(allEffects(m))
					
					ggplot(x, aes(x=mass_loss, y=med_bout, col=sex))+geom_point()+stat_smooth(method='lm')
						ms=lm(med_bout~mass_loss*sex,x)
						nsim <- 2000
						bsim <- sim(ms, n.sim=nsim)  
						apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
						summary(glht(ms))
						plot(allEffects(ms))
						AIC(m,ms)
						
				# first bout
					ggplot(y, aes(x=mass_loss, y=first_bout))+geom_point()+stat_smooth(method='lm')
							m=lm(first_bout~mass_loss,y)
							nsim <- 2000
							bsim <- sim(m, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
							summary(glht(m))
							plot(allEffects(m))
					
					ggplot(y, aes(x=mass_loss, y=first_bout, col=sex))+geom_point()+stat_smooth(method='lm')
						ms=lm(first_bout~mass_loss*sex,y)
						nsim <- 2000
						bsim <- sim(ms, n.sim=nsim)  
						apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
						summary(glht(ms))
						plot(allEffects(ms))
						AIC(m,ms)
			}				
			{# efficiency	
				######CHECK what is going on here
					
							a[a$inc_eff<0.1,]

							
				densityplot(~a$inc_eff)
				densityplot(~a$inc_eff)
				# three bouts
					ggplot(a, aes(x=mass_loss, y=inc_eff))+geom_point()
					ggplot(a, aes(x=mass_loss, y=inc_eff))+geom_point()+stat_smooth(method='lm')
						m=lmer(inc_eff~scale(bout_length)+scale(mass_loss)+(1|nest),a)
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
						apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
						summary(glht(m))
						plot(allEffects(m))
					
					ggplot(a, aes(x=mass_loss, y=inc_eff, col=sex))+geom_point()
					ggplot(a, aes(x=mass_loss, y=inc_eff, col=sex))+geom_point()+stat_smooth(method='lm')
						ms=lmer(inc_eff~scale(bout_length)+scale(mass_loss)*sex+(1|nest),a)
						ms=lmer(bout_length~scale(mass_loss)*sex+(1|nest),a)
						nsim <- 2000
						bsim <- sim(ms, n.sim=nsim)  
						apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
						summary(glht(ms))
						plot(allEffects(ms))
						AIC(m,ms)
					
				# median bout
					ggplot(x, aes(x=mass_loss, y=med_bout))+geom_point()+stat_smooth(method='lm')
							m=lm(med_bout~mass_loss,x)
							nsim <- 2000
							bsim <- sim(m, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
							summary(glht(m))
							plot(allEffects(m))
					
					ggplot(x, aes(x=mass_loss, y=med_bout, col=sex))+geom_point()+stat_smooth(method='lm')
						ms=lm(med_bout~mass_loss*sex,x)
						nsim <- 2000
						bsim <- sim(ms, n.sim=nsim)  
						apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
						summary(glht(ms))
						plot(allEffects(ms))
						AIC(m,ms)
						
				# first bout
					ggplot(y, aes(x=mass_loss, y=first_bout))+geom_point()+stat_smooth(method='lm')
							m=lm(first_bout~mass_loss,y)
							nsim <- 2000
							bsim <- sim(m, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
							summary(glht(m))
							plot(allEffects(m))
					
					ggplot(y, aes(x=mass_loss, y=first_bout, col=sex))+geom_point()+stat_smooth(method='lm')
						ms=lm(first_bout~mass_loss*sex,y)
						nsim <- 2000
						bsim <- sim(ms, n.sim=nsim)  
						apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
						summary(glht(ms))
						plot(allEffects(ms))
						AIC(m,ms)
			}				
			{# after bout constancy/length of treated bird in relation to its compensation
				load(file=paste(wd,'experimental.Rdata', sep=""))
				bb=b[b$exper=='c',]
				bt=b[b$exper=='t',]
				bb$inc_eff_treated=bt$inc_eff[match(bb$nest, bt$nest)]
				bb$compensation=bb$inc_eff_treated/bb$inc_eff
								
				a=inc[inc$exper=='a' & inc$treated_bird=='y',]
				a$compensation=bb$compensation[match(a$nest, bb$nest)]
				x=ddply(a,.(nest,sex), summarise, med_bout=median(bout_length), mean_bout=mean(bout_length),med_eff=median(inc_eff), mean_eff=mean(inc_eff), compensation=median(compensation))
				y=ddply(a,.(nest,sex), summarise, first_bout=bout_length[bout_start==min(bout_start)], first_eff=inc_eff[bout_start==min(bout_start)], compensation=median(compensation))
				
				densityplot(~a$compensation)
				densityplot(~a$bout_length)
			  {# bout				
				# three bouts
					ggplot(a, aes(x=compensation, y=bout_length))+geom_point()
					ggplot(a, aes(x=compensation, y=bout_length))+geom_point()+stat_smooth(method='lm')
					
					
						m=lmer(bout_length~scale(compensation)+(1|nest),a)
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
						apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
						summary(glht(m))
						plot(allEffects(m))
									
					ggplot(a, aes(x=compensation, y=bout_length, col=sex))+geom_point()
					ggplot(a, aes(x=compensation, y=bout_length, col=sex))+geom_point()+stat_smooth(method='lm')
						ms=lmer(bout_length~scale(compensation)*sex+(1|nest),a)
						nsim <- 2000
						bsim <- sim(ms, n.sim=nsim)  
						apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
						summary(ms)
						summary(glht(ms))
						plot(allEffects(ms))
						AIC(m,ms)
					
				# median bout
					ggplot(x, aes(x=compensation, y=med_bout))+geom_point()+stat_smooth(method='lm')
							m=lm(med_bout~scale(compensation),x)
							nsim <- 2000
							bsim <- sim(m, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
							summary(glht(m))
							plot(allEffects(m))
					
					ggplot(x, aes(x=compensation, y=med_bout, col=sex))+geom_point()+stat_smooth(method='lm')
						ms=lm(med_bout~compensation*sex,x)
						nsim <- 2000
						bsim <- sim(ms, n.sim=nsim)  
						apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
						summary(glht(ms))
						plot(allEffects(ms))
						AIC(m,ms)
						
				# mean bout
					ggplot(x, aes(x=compensation, y=mean_bout))+geom_point()+stat_smooth(method='lm')
							m=lm(mean_bout~scale(compensation),x)
							nsim <- 2000
							bsim <- sim(m, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
							summary(glht(m))
							plot(allEffects(m))
					
					ggplot(x, aes(x=compensation, y=mean_bout, col=sex))+geom_point()+stat_smooth(method='lm')
						ms=lm(mean_bout~compensation*sex,x)
						nsim <- 2000
						bsim <- sim(ms, n.sim=nsim)  
						apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
						summary(glht(ms))
						plot(allEffects(ms))
						AIC(m,ms)
						
				# first bout
					ggplot(y, aes(x=compensation, y=first_bout))+geom_point()+stat_smooth(method='lm')
							m=lm(first_bout~compensation,y)
							nsim <- 2000
							bsim <- sim(m, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
							summary(glht(m))
							plot(allEffects(m))
					
					ggplot(y, aes(x=compensation, y=first_bout, col=sex))+geom_point()+stat_smooth(method='lm')
						ms=lm(first_bout~compensation*sex,y)
						nsim <- 2000
						bsim <- sim(ms, n.sim=nsim)  
						apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
						summary(glht(ms))
						plot(allEffects(ms))
						AIC(m,ms)
			}				
			  {# efficiency	
					a[a$inc_eff<0.1,]
					densityplot(~a$inc_eff)
				
				# three bouts
					ggplot(a, aes(x=compensation, y=inc_eff))+geom_point()
					ggplot(a, aes(x=compensation, y=inc_eff))+geom_point()+stat_smooth(method='lm')
						m=lmer(inc_eff~scale(bout_length)+compensation+(1|nest),a)
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
						apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
						summary(glht(m))
						plot(allEffects(m))
					
					ggplot(a, aes(x=compensation, y=inc_eff, col=sex))+geom_point()
					ggplot(a, aes(x=compensation, y=inc_eff, col=sex))+geom_point()+stat_smooth(method='lm')
						ms=lmer(inc_eff~scale(bout_length)+compensation*sex+(1|nest),a)
						nsim <- 2000
						bsim <- sim(ms, n.sim=nsim)  
						apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
						summary(glht(ms))
						plot(allEffects(ms))
						AIC(m,ms)
					
				# med_eff
					ggplot(x, aes(x=compensation, y=med_eff))+geom_point()+stat_smooth(method='lm')
							m=lm(med_eff~scale(med_bout)+compensation,x)
							nsim <- 2000
							bsim <- sim(m, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
							summary(glht(m))
							plot(allEffects(m))
					
					ggplot(x, aes(x=compensation, y=med_eff, col=sex))+geom_point()+stat_smooth(method='lm')
						ms=lm(med_eff~scale(med_bout)+compensation*sex,x)
						nsim <- 2000
						bsim <- sim(ms, n.sim=nsim)  
						apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
						summary(glht(ms))
						plot(allEffects(ms))
						AIC(m,ms)
				# mean bout
					ggplot(x, aes(x=compensation, y=mean_eff))+geom_point()+stat_smooth(method='lm')
							m=lm(mean_eff~scale(med_bout)+compensation,x)
							nsim <- 2000
							bsim <- sim(m, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
							summary(glht(m))
							plot(allEffects(m))
					
					ggplot(x, aes(x=compensation, y=mean_eff, col=sex))+geom_point()+stat_smooth(method='lm')
						ms=lm(mean_eff~scale(med_bout)+compensation*sex,x)
						nsim <- 2000
						bsim <- sim(ms, n.sim=nsim)  
						apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
						summary(glht(ms))
						plot(allEffects(ms))
						AIC(m,ms)
								
				# first bout
					ggplot(y, aes(x=compensation, y=first_eff))+geom_point()+stat_smooth(method='lm')
							m=lm(first_eff~scale(first_bout)+compensation,y)
							nsim <- 2000
							bsim <- sim(m, n.sim=nsim)  
							apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
							summary(glht(m))
							plot(allEffects(m))
					
					ggplot(y, aes(x=compensation, y=first_eff, col=sex))+geom_point()+stat_smooth(method='lm')
						ms=lm(first_eff~scale(first_bout)+compensation*sex,y)
						nsim <- 2000
						bsim <- sim(ms, n.sim=nsim)  
						apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
						summary(glht(ms))
						plot(allEffects(ms))
						AIC(m,ms)
			}				
			}	
		}
		
		{# Figure 6
			{# run first - add mass loss to data
				u=read.csv(paste(wd,'experiment_metadata.csv', sep=""), stringsAsFactors=FALSE)	
				a=inc[inc$exper=='a' & inc$treated_bird=='n',]
				a$mass_loss=u$mass_loss[match(a$bird_ID_filled, u$ID_taken)]
				a$col_=ifelse(a$sex=="f", "#FCB42C","#535F7C")
				
			}
			{# run first - add compensatin to data
				bb$compensation=bb$inc_eff_treated/bb$inc_eff
				f=inc[inc$exper=='a' & inc$treated_bird=='y',]
				f$compensation=bb$compensation[match(f$nest, bb$nest)]
				f$col_=ifelse(f$sex=="f", "#FCB42C","#535F7C")
			}					
			{# run first - prepare predictions
				{# mass loss
					{# constancy	
						m=lmer(inc_eff~mass_loss*sex+(mass_loss|bird_ID),a, REML=FALSE)
								# simulation		
								nsim <- 2000
								bsim <- sim(m, n.sim=nsim)  
								apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							
							# coefficients
								v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))	
							# predicted values		
								newD=data.frame(mass_loss=seq(min(a$mass_loss),max(a$mass_loss),length.out=100),
												sex=c('f','m')
												)
								
							# exactly the model which was used has to be specified here 
								X <- model.matrix(~ mass_loss*sex,data=newD)	
											
							# calculate predicted values and creditability intervals
								newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
										predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
										for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@fixef[i,]
										newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
										newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
										#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
										#newD=newD[order(newD$t_tundra),]
								pcon=newD
								pconF=pcon[pcon$sex=='f',]
								pconM=pcon[pcon$sex=='m',]
					}			
					{# bout	
						m=lmer(bout_length~mass_loss*sex+(mass_loss|bird_ID),a,REML=FALSE)
								# simulation		
								nsim <- 2000
								bsim <- sim(m, n.sim=nsim)  
								apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							
							# coefficients
								v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))	
							# predicted values		
								newD=data.frame(mass_loss=seq(min(a$mass_loss),max(a$mass_loss),length.out=100),
												sex=c('f','m')
												)
								
							# exactly the model which was used has to be specified here 
								X <- model.matrix(~ mass_loss*sex,data=newD)	
											
							# calculate predicted values and creditability intervals
								newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
										predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
										for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@fixef[i,]
										newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
										newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
										#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
										#newD=newD[order(newD$t_tundra),]
								pbout=newD
								pboutF=pbout[pbout$sex=='f',]
								pboutM=pbout[pbout$sex=='m',]
					}			
				
				}
				{# compensation
					{# constancy	
						m=lmer(inc_eff~compensation*sex+(scale(compensation)|bird_ID),f,REML=FALSE)
								# simulation		
								nsim <- 2000
								bsim <- sim(m, n.sim=nsim)  
								apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							
							# coefficients
								v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))	
							# predicted values		
								newD=data.frame(compensation=seq(min(f$compensation),max(f$compensation),length.out=100),
												sex=c('f','m')
												)
								
							# exactly the model which was used has to be specified here 
								X <- model.matrix(~ compensation*sex,data=newD)	
											
							# calculate predicted values and creditability intervals
								newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
										predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
										for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@fixef[i,]
										newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
										newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
										#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
										#newD=newD[order(newD$t_tundra),]
								pccon=newD
								pcconF=pccon[pccon$sex=='f',]
								pcconM=pccon[pccon$sex=='m',]
					}			
					{# bout	
						m=lmer(bout_length~compensation*sex+(compensation|nest),f,REML=FALSE)
								# simulation		
								nsim <- 2000
								bsim <- sim(m, n.sim=nsim)  
								apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
							
							# coefficients
								v <- apply(bsim@fixef, 2, quantile, prob=c(0.5))	
							# predicted values		
								newD=data.frame(compensation=seq(min(f$compensation),max(f$compensation),length.out=100),
												sex=c('f','m')
												)
								
							# exactly the model which was used has to be specified here 
								X <- model.matrix(~ compensation*sex,data=newD)	
											
							# calculate predicted values and creditability intervals
								newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
										predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
										for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@fixef[i,]
										newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
										newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
										#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
										#newD=newD[order(newD$t_tundra),]
								pcbout=newD
								pcboutF=pcbout[pcbout$sex=='f',]
								pcboutM=pcbout[pcbout$sex=='m',]
					}			
				
				}
			}
				{# plot a,b,c, even less labels
				  dev.new(width=3.5-1.06299,height=(1.85*2.15)-1.14173)
				  #png(paste(out2,"Figure_6_abcd.png", sep=""), width=3.5-1.06299,height=(1.85*2.15)-1.14173,units="in",res=600)
				   par(mfrow=c(2,2),mar=c(0.5,0,0,0.5),oma = c(2.2, 1.8, 1.2, 0.1),ps=12, mgp=c(1.2,0.35,0), las=1, cex=1, col.axis="grey30",font.main = 1, col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt
				{# constancy
					{# mass loss
						par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
						plot(pcon$pred~pcon$mass_loss ,pch=19,ylim=c(0.58,1.05), xlim=c(-0.2,3.2), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
										
						axis(1, at=-c(-2.5,-1.5,-0.5),labels=FALSE,col.ticks="grey90")				
						axis(1, at=-seq(-3,0,by=1),labels=FALSE,cex.axis=0.5,mgp=c(0,-0.20,0))
							#mtext('Mass loss [g]',side=1,line=0.3, cex=0.6, las=1, col='grey30')
						
						axis(2, at=c(0.7,0.9),labels=FALSE,col.ticks="grey90")						
						axis(2, at=seq(0.6,1,by=0.2),labels=c('0.6','0.8','1.0'))	
							mtext('Nest attendance',side=2,line=1.1, cex=0.6, las=3, col='grey30')
																	
						mtext(expression(bold("a")),side=3,line=0, cex=0.7, las=1, col="grey30")
						
						text(x=1.3,y=0.70, labels='\u2640', col='#FCB42C', cex=0.6)
						text(x=1.6,y=0.72, labels='\u2642', col='#535F7C', cex=0.6)
											
					
					# predictions
							# female
							polygon(c(pconF$mass_loss, rev(pconF$mass_loss)), c(pconF$lwr, 
								rev(pconF$upr)), border=NA, col=adjustcolor(col_t ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pconF$mass_loss, pconF$pred, col=col_t,lwd=1)
							
							# male
							polygon(c(pconM$mass_loss, rev(pconM$mass_loss)), c(pconM$lwr, 
								rev(pconM$upr)), border=NA, col=adjustcolor(col_c ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
								lines(pconM$mass_loss, pconM$pred, col=col_c,lwd=1)
					
					# data
						#points(a$inc_eff~a$mass_loss, col=a$col_,bg=adjustcolor(a$col_, alpha.f = 0.4), pch=21, cex=0.5)	
						points(a$inc_eff~a$mass_loss, col=adjustcolor(a$col_, alpha.f = 0.7), pch=20, cex=0.6)	
					
					}
					{# compensation
						par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
						plot(pccon$pred~pccon$compensation ,pch=19,ylim=c(0.58,1.05), xlim=c(0,1.05), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
										
						axis(1, at=c(0.25,0.75),labels=FALSE,col.ticks="grey90")
						axis(1, at=seq(0,1,by=0.5),labels=FALSE)
							#mtext('Mass loss [g]',side=1,line=0.3, cex=0.6, las=1, col='grey30')
							
						axis(2, at=c(0.7,0.9),labels=FALSE,col.ticks="grey90")						
						axis(2, at=seq(0.6,1,by=0.2),labels=FALSE)	
							#mtext('Nest attendance',side=2,line=1.1, cex=0.6, las=3, col='grey30')
						
						mtext(expression(bold("b")),side=3,line=0, cex=0.7, las=1, col="grey30")
						
						#text(x=-3,y=0.60, labels='\u2640', col='#FCB42C', cex=0.6)
						#text(x=-2.7,y=0.62, labels='\u2642', col='#535F7C', cex=0.6)
											
					
					# predictions
							# female
							polygon(c(pcconF$compensation, rev(pcconF$compensation)), c(pcconF$lwr, 
								rev(pcconF$upr)), border=NA, col=adjustcolor(col_t ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pcconF$compensation, pcconF$pred, col=col_t,lwd=1)
							
							# male
							polygon(c(pcconM$compensation, rev(pcconM$compensation)), c(pcconM$lwr, 
								rev(pcconM$upr)), border=NA, col=adjustcolor(col_c ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
								lines(pcconM$compensation, pcconM$pred, col=col_c,lwd=1)
					
					# data
						#points(a$inc_eff~a$compensation, col=a$col_,bg=adjustcolor(a$col_, alpha.f = 0.4), pch=21, cex=0.5)	
						points(f$inc_eff~f$compensation, col=adjustcolor(f$col_, alpha.f = 0.7), pch=20, cex=0.6)	
					
					}
				}
				{# bout
					{# mass loss
					par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
					plot(pbout$pred~pbout$mass_loss,pch=19,xlim=c(-0.2,3.2), ylim=c(0,16.5), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
										
						axis(1, at=-c(-2.5,-1.5,-0.5),labels=FALSE,col.ticks="grey90")				
						axis(1, at=seq(0,3,by=1),labels=c('0','1',"2","3"),cex.axis=0.5,mgp=c(0,-0.20,0))
							mtext('Mass loss [g]',side=1,line=0.5, cex=0.6, las=1, col='grey30')			
						
						axis(2, at=c(4,8),labels=FALSE,col.ticks="grey90")
						axis(2, at=seq(0,16,by=8),labels=seq(0,16,by=8))
							mtext('Bout [h]',side=2,line=1.1, cex=0.6, las=3, col='grey30')
							
					# predictions
							# female
							polygon(c(pboutF$mass_loss, rev(pboutF$mass_loss)), c(pboutF$lwr, 
								rev(pboutF$upr)), border=NA, col=adjustcolor(col_t ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pboutF$mass_loss, pboutF$pred, col=col_t,lwd=1)
							
							# male
							polygon(c(pboutM$mass_loss, rev(pboutM$mass_loss)), c(pboutM$lwr, 
								rev(pboutM$upr)), border=NA, col=adjustcolor(col_c ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
								lines(pboutM$mass_loss, pboutM$pred, col=col_c,lwd=1)
					
					# data
						#points(a$inc_eff~a$mass_loss, col=a$col_,bg=adjustcolor(a$col_, alpha.f = 0.4), pch=21, cex=0.5)	
						points(a$bout_length~a$mass_loss, col=adjustcolor(a$col_, alpha.f = 0.7), pch=20, cex=0.6)	
					}
					{# compensation
						par(ps=12,cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.15,bty="l",xpd=TRUE)
						plot(pcbout$pred~pcbout$compensation,pch=19,xlim=c(0,1.05), ylim=c(0,16.5), xlab=NA, ylab=NA, yaxt='n',xaxt='n', type='n')
										
						axis(1, at=c(0.25,0.75),labels=FALSE,col.ticks="grey90")
						axis(1, at=seq(0,1,by=0.5),labels=c('0','50',"100"),cex.axis=0.5,mgp=c(0,-0.20,0))
							mtext('Compensation [%]',side=1,line=0.5, cex=0.6, las=1, col='grey30')			
						
						axis(2, at=c(4,8),labels=FALSE,col.ticks="grey90")
						axis(2, at=seq(0,16,by=8),labels=FALSE)
							#mtext('Bout [h]',side=2,line=1.1, cex=0.6, las=3, col='grey30')
							
					# predictions
							# female
							polygon(c(pcboutF$compensation, rev(pcboutF$compensation)), c(pcboutF$lwr, 
								rev(pcboutF$upr)), border=NA, col=adjustcolor(col_t ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
							lines(pcboutF$compensation, pcboutF$pred, col=col_t,lwd=1)
							
							# male
							polygon(c(pcboutM$compensation, rev(pcboutM$compensation)), c(pcboutM$lwr, 
								rev(pcboutM$upr)), border=NA, col=adjustcolor(col_c ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
								lines(pcboutM$compensation, pcboutM$pred, col=col_c,lwd=1)
					
					# data
						#points(a$inc_eff~a$mass_loss, col=a$col_,bg=adjustcolor(a$col_, alpha.f = 0.4), pch=21, cex=0.5)	
						points(f$bout_length~f$compensation, col=adjustcolor(f$col_, alpha.f = 0.7), pch=20, cex=0.6)	
					}
				}
															
				dev.off()
			}	
		
		}
			
	}	
	
	{# not in the paper
		{# control vs treated period length
			densityplot(~bout_length,groups=exper, data=b, auto.key=TRUE)
		
			dev.new(width=2.5,height=2.5)
					#png(paste(out,"figs/tree_period_high_.png", sep=""), width=2.3,height=3.5,units="in",res=600)
			par(mar=c(2,2,0.7,0.7), ps=12, mgp=c(1.1,0.35,0), las=1, cex.lab=0.7, cex.axis=0.6, tcl=-0.15,bty="l",xpd=TRUE, col.axis="grey30", col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt
		
			plot(bt$bout_length~bb$bout_length, xlab=NA, ylab='Treatment bout length [h]', col='grey20', xaxt='n')
				axis(1, at=seq(10,14,by=1), labels=FALSE)
				text(seq(10,14,by=1), par("usr")[3]-0.35, labels = seq(10,14,by=1),  xpd = TRUE, cex=0.6, col="grey30") #labels = z_g$genus
				mtext('Control bout length [h]',side=1,line=0.6, cex=0.7, las=1, col='grey30')
		}
		{# fat a change in removed parents
		d =read.csv(paste(wd,'Data/captivity.csv', sep=""),stringsAsFactors=FALSE)
		dd = d[!is.na(d$mass) & d$phase %in%c('start','end'),]
		d = d[!is.na(d$fat) & d$phase %in%c('start','end'),]
		d$n = 1
		v = ddply(d,.(ring_num),summarise, n = sum(n))
		v = ddply(d,.(ring_num, phase),summarise, fat = median(fat))
		pd <- position_dodge(0.4)
		v$phase = factor(v$phase, levels = c('start','end'))
		ggplot(v,aes(x=phase, y = fat, col = factor(ring_num)))+geom_point( position = pd) +geom_line(position = pd, aes(group = factor(ring_num), col = factor(ring_num))) +
		ggplot(v,aes(x=phase, y = fat, col = factor(ring_num)))+geom_point( position = pd) +geom_line(position = pd, aes(group = factor(ring_num), col = factor(ring_num))) +
			theme(legend.position = "none") + ylab('fat score') + xlab('captivity phase')
			xyplot(fat~numeric(phase), group = factor(ring_num), data = v)
	}
	}



{# script to create metadata with what is happening on the nest
	{# metadata - preparation
		#con = dbcon(user='root',host='127.0.0.1', password='',database='barrow_2013')						   
		#id=dbq(con, "select distinct(ring_num) as ring_num from barrow_2013.captivity")							   

	id=dbGetQuery(conMy, "select distinct(ring_num) as ID, s.sex from barrow_2013.captivity c
							LEFT JOIN ( select ID, (case when (sex = 1) then 'm' else 'f' end) as sex from avesatbarrow.sex) s on c.ring_num=s.ID")							
	length(unique(id$ID))
	n=dbGetQuery(conMy,	"select ring_num as ID, nest  from barrow_2013.captures where  bird_taken='y' and ring_num is not null and species='sesa'")
    length(unique(n$nest))
	n[!n$ID%in%id$ID,]
		
	n$sex=id$sex[match(n$ID,id$ID)]
	summary(factor(id$sex))	# 17 female, 12 male
	
	nn=dbGetQuery(conMy, "select lower(nest) nest, IDfemale, IDmale from avesatbarrow.nests where year_=2013 and species='SESA'")
	
	nn=nn[nn$nest%in%n$nest,]	
	nn$taken=n$sex[match(nn$nest,n$nest)]
	nn$ID_taken=ifelse(nn$taken=="f", nn$IDfemale, nn$IDmale)
	nn$ID_treated=ifelse(nn$taken=="f", nn$IDmale, nn$IDfemale)
	
	nn$back=c('n','y',
	
	write.table(nn,file=paste(wd,'nests_birds.csv', sep=""), sep=",",col.names = TRUE, row.names=FALSE)
		#nn$taken_ID=n$ID[match(nn$nest,n$nest)]
	d=read.csv(paste(wd,'nests_birds.csv', sep=""), stringsAsFactors=FALSE)	
	p=dbGetQuery(conMy, "select distinct(ring_num) as ID, partner_present from barrow_2013.captivity 
								where partner_present is not null")
	d$present=p$partner_present[match(d$ID_taken, p$ID)]
	
	write.table(d,file=paste(wd,'experiment_metadata_new.csv', sep=""), sep=",",col.names = TRUE, row.names=FALSE)		
	

	# add release time		
	h=wet(conMy, "SELECT*FROM extract_2013.mb_temp1_2013 where taken='free'")	
	o=read.csv(paste(wd,'experiment_metadata.csv', sep=""), stringsAsFactors=FALSE)	
	o$release=h$release_time[match(o$nest,h$nest)]
	write.table(o,file=paste(wd,'experiment_metadata.csv', sep=""), sep=",",col.names = TRUE, row.names=FALSE)	
	
	# add capture time
		h=wet(conMy, "SELECT*FROM extract_2013.mb_temp1_2013 where taken='y'")	
		o=read.csv(paste(wd,'experiment_metadata.csv', sep=""), stringsAsFactors=FALSE)	
		o$taken_time=h$caught_date_time[match(o$nest,h$nest)]
		write.table(o,file=paste(wd,'experiment_metadata.csv', sep=""), sep=",",col.names = TRUE, row.names=FALSE)	
	# add nest location
		u=read.csv(paste(wd,'experiment_metadata_.csv', sep=""), stringsAsFactors=FALSE)	
		h=wet(conMy, "SELECT lower(nest) as nest, latit as lat, longit as lon FROM avesatbarrow.nests where species='sesa' and year_=2013")
		u$lat=h$lat[match(u$nest,h$nest)]
		u$lon=h$lon[match(u$nest,h$nest)]
		write.table(u,file=paste(wd,'experiment_metadata__.csv', sep=""), sep=",",col.names = TRUE, row.names=FALSE)	
		
	# add mass loss
		u=read.csv(paste(wd,'experiment_metadata__.csv', sep=""), stringsAsFactors=FALSE)	
		g=dbGetQuery(conMy, "select*from barrow_2013.captivity")
		gg=ddply(g[!is.na(g$mass),],.(ring_num),summarise,mass_loss=(mass[phase=='start']-mass[phase=='end']), mass_end=mass[phase=='end'])
						#ggplot(gg,aes(x=mass_loss))+geom_density()
						#ggplot(gg,aes(x=mass_loss, y=mass_end))+geom_point()+stat_smooth()
		u$mass_loss=gg$mass_loss[match(u$ID_taken, gg$ring_num)]
		write.table(u,file=paste(wd,'experiment_metadata.csv', sep=""), sep=",",col.names = TRUE, row.names=FALSE)	
		
	}

			{# prepare bout lengths - for definition of control treatment
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
			nn$bt=round(bt$inc_eff[match(nn$nest,bt$nest)],2)
			nn$bb=round(bb$inc_eff[match(nn$nest,bb$nest)],2)
			nn$bt=ifelse(nchar(nn$bt)==3, paste(nn$bt,"0",sep=""), nn$bt)
			nn$bb=ifelse(nchar(nn$bb)==3, paste(nn$bb,"0",sep=""), nn$bb)
			nn$label=paste(nn$bb,nn$bt,"              ",sep="       ")
			
			nn$sex=bt$sex[match(nn$nest,bt$nest)]
			nn$col_=ifelse(nn$sex=="f", "#FCB42C","#535F7C")
			nnn=nn$label
			names(nnn)=nn$n
			nnnest=nn$nest
			names(nnnest)=nn$n
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
												if(tst=='f'){x=data.frame(nest=x$nest[1],start_=min(x$bout_start[!is.na(x$exper) & x$exper=='t']), end_=min(x$bout_start[x$bout_start>tst2 & x$sex=='m']),stringsAsFactors=FALSE)}
												if(tst=='m'){x=data.frame(nest=x$nest[1], start_=min(x$bout_start[!is.na(x$exper) & x$exper=='t']),end_=min(x$bout_start[x$bout_start>tst2 & x$sex=='f']),stringsAsFactors=FALSE)}
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
                       n = e$n, sex=e$sex, nest=e$nest)
				
				o=read.csv(paste(wd,'experiment_metadata.csv', sep=""), stringsAsFactors=FALSE)	
				e$taken=as.POSIXct(o$taken_time[match(e$nest,o$nest)])
				write.table(e,file=paste(wd,'experiment_times.csv', sep=""), sep=",",col.names = TRUE, row.names=FALSE)	
		}
}
				
		
		
		
		
		
{# results experimental

	{# constancy plot
			dev.new(width=2.5,height=2.5)
			#png(paste(out2,"constancy_treated_control_relationship.png", sep=""), width=2.5,height=2.5,units="in",res=600)
			par(mar=c(1.5,2,0.7,0.7), ps=12, mgp=c(1.2,0.35,0), las=1, cex.lab=0.7, cex.axis=0.6, tcl=-0.15,bty="l",xpd=TRUE, col.axis="grey30", col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt

	  plot(bb$inc_eff_treated~bb$inc_eff, col=bb$col_,bg=adjustcolor( bb$col_, alpha.f = 0.4), pch=21,xlim=c(0,1), ylim=c(0,1), ylab='Nest attendance - treated bout', xlab=NA, xaxt='n')
		lines(c(0,1),c(0,1), lty=3, col="lightgrey")
		axis(1, at=seq(0,1,by=0.2), labels=FALSE)
		text(seq(0,1,by=0.2), par("usr")[3]-0.05, labels = seq(0,1,by=0.2),  xpd = TRUE, cex=0.6, col="grey30") #labels = z_g$genus
		mtext('Nest attendance - before bout',side=1,line=0.6, cex=0.7, las=1, col='grey30')
		
		legend("bottomleft", legend=c('\u2640','\u2642'), pch=19, col=c( '#FCB42C','#535F7C'), cex=0.8, bty='n')
		#dev.off()
		estimate_m=data.frame(y = c(0.91273,  (0.91273-0.36650)),x=c('c','t'), stringsAsFactors=FALSE)
		estimate_t=data.frame(y = c( 0.9270140, 0.5314018),x=c('c','t'), stringsAsFactors=FALSE)
	
		ggplot(b,aes(y=inc_eff, x=exper, fill=exper))+geom_point(aes(color=exper),alpha=0.5,position = position_jitterdodge(dodge.width=0.9,jitter.width=2))+geom_boxplot(fill=NA,alpha=0.5,position = position_dodge(width=0.9))+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_color_manual(values=c("#E69F00", "#56B4E9"))+ geom_point(data=estimate_t,aes(y = y,x=x),fill="red", colour="red",size=2)+ geom_point(data=estimate_m,aes(y = y,x=x),fill="pink", colour="pink",size=2)+theme(legend.position="none")
		
		ggplot(b,aes(y=inc_eff, x=sex, col=exper))+geom_boxplot(fill=NA,alpha=0.5,position = position_dodge(width=0.9))
		ggplot(b,aes(y=inc_eff, x=exper, col=sex))+
							geom_boxplot(aes(fill=NA,col=sex),alpha=0.5,position = position_dodge(width=0.9))+
							geom_point(aes(col=sex),alpha=0.5,position = position_jitterdodge(dodge.width=0.9,jitter.width=1))+
							scale_fill_manual(values = c("#E69F00", "#56B4E9"))+
							scale_color_manual(values=c("#E69F00", "#56B4E9"))
		
		+geom_boxplot(fill=NA,alpha=0.5,position = position_dodge(width=0.9))+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_color_manual(values=c("#E69F00", "#56B4E9"))+ geom_point(data=estimate_t,aes(y = y,x=x),fill="red", colour="red",size=2)+ geom_point(data=estimate_m,aes(y = y,x=x),fill="pink", colour="pink",size=2)+theme(legend.position="none")
		ggplot(b,aes(y=inc_eff, x=sex, fill=exper))+geom_point(aes(color=exper),alpha=0.5,position = position_jitterdodge(dodge.width=0.9,jitter.width=2))+geom_boxplot(fill=NA,alpha=0.5,position = position_dodge(width=0.9))+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_color_manual(values=c("#E69F00", "#56B4E9"))+ geom_point(data=estimate_t,aes(y = y,x=x),fill="red", colour="red",size=2)+ geom_point(data=estimate_m,aes(y = y,x=x),fill="pink", colour="pink",size=2)+theme(legend.position="none")
		
	}
	{# compensation distribution plot
		summary(100*bt$inc_eff/bb$inc_eff)
		pr=data.frame(difference_in_constancy_percent=(bb$inc_eff-bt$inc_eff)/(bb$inc_eff/100), boxp='boxplot')
		ggplot(pr,aes(x=boxp, y=difference_in_constancy_percent))+geom_boxplot()+
		ylab("Difference in constancy\nbetween control and treated bout\n[%]") +theme(axis.title.x=element_blank())
		+coord_cartesian(ylim = c(-5,60)) 
		pr=data.frame(difference_in_constancy=(bb$inc_eff-bt$inc_eff), boxp='boxplot')
		ggplot(pr,aes(x=boxp, y=difference_in_constancy))+geom_boxplot()+coord_cartesian(ylim = c(-0.05,0.6)) 
	}
	{# inc eff only for treated
			ms=lm(inc_eff~sex,bt)
			summary(ms)
			nsim <- 2000
				bsim <- sim(ms, n.sim=nsim)  
				apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
	}		

	{# constancy over experimental period graph
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
			{# all in plot with facets
		l=list()
			{# prepare dataset	
				for (j in 1:length(nests_)){
						print(nests_[j])
						pp=paste("extract_2013.mb_inc_2013_",nests_[j],sep="")	
						b=wet(conMy, paste("SELECT pk, datetime_, inc_t, incubation FROM  ", pp," where bout_ID_tr in ('50')"))
						b$datetime_= as.POSIXct(b$datetime_)
						b$datetime_c=as.numeric(difftime(b$datetime_,b$datetime_[1]+60*dd$bout[dd$nest==nests_[j]], units='hours'))
						
						if(nests_[j]=='s806'){b$inc_t=b$incubation
											  b$inc_t[is.na(b$inc_t)]=0}
						
						b$roll_con=rollmean(b$inc_t, 12*60, align="center", fill="extend") # 12=1min
						b$nest=nests_[j]
						l[[j]]=b
						}
				u=do.call(rbind,l)
				# add order for plotting
				u$n=nn$n[match(u$nest,nn$nest)]
				u$sex=nn$sex[match(u$nest,nn$nest)]
				u$col_=nn$col_[match(u$nest,nn$nest)]
				
				}
			{# plot
				summary(u)
				ggplot(b,aes(y=inc_eff, x=exper, fill=exper))+geom_point(aes(color=exper),alpha=0.5,position = position_jitterdodge(dodge.width=0.9,jitter.width=2))+geom_boxplot(fill=NA,alpha=0.5,position = position_dodge(width=0.9))+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_color_manual(values=c("#E69F00", "#56B4E9"))+ geom_point(data=estimate_t,aes(y = y,x=x),fill="red", colour="red",size=2)+ geom_point(data=estimate_m,aes(y = y,x=x),fill="pink", colour="pink",size=2)+theme(legend.position="none")
				
				nn$label=paste(nn$bb,nn$bt,sep="       ")
				nnn=nn$label
				names(nnn)=nn$n
				
				{# 7x4 with labels -with constancy lables and sex color
				dev.new(width=7.4, height=4.35) # for 4*7
				ggplot(u,aes(x=datetime_c,y=roll_con, fill=sex))+
						#stat_smooth(method="lm", fill="grey95",alpha=1, aes(colour=genus_))+
						#scale_colour_manual(values=go$cols, name="Genus")+
						facet_wrap(~n, nrow=5, labeller=as_labeller(nnn))+
						geom_point(aes(color=sex), shape=20,size=0.25)+#colour="grey38",aes(fill=genus_))+
						scale_color_manual(values=c("#FCB42C","#535F7C"))+
						#scale_fill_manual(values=go$cols, name="Genus")+ 
						labs(x="Time in experimental period [h]", y="Nest attendance")+
						scale_x_continuous(limits=c(-16,16), breaks=seq(-15,15,by=5), labels=c(seq(-15,15,by=5)))+
						geom_vline(xintercept=0, lty=3,lwd=0.5, col='red')+
						#scale_x_continuous(breaks=15, labels=15)+
						scale_y_continuous(limits=c(-0.05,1.05),breaks=seq(0,1,by=0.25), labels=seq(0,1,by=0.25))+
						#coord_cartesian(ylim=c(0,1), xlim=c(-15,15))+ #coord_cartesian(ylim=c(-5,61))+ # needs to be 10 if it is not to touch x-axis
						theme_light()+
						theme(	axis.line=element_line(colour="grey70", size=0.25),
								#panel.border=element_rect(colour="white"),
								panel.border=element_rect(colour="grey70", size=0.25),
								panel.grid = element_blank(),
								
								axis.title=element_text(size=7, colour="grey50"),
								axis.title.y = element_text(vjust=1),
								axis.title.x = element_text(vjust=0.2),
								axis.text=element_text(size=6),# margin=units(0.5,"mm")),
								axis.ticks.length=unit(0.5,"mm"),
								#axis.ticks.margin,
								
								strip.text.x =element_text(size = 6, color="grey30",  margin=margin(1,1,1,1,"mm")), #grey50
								strip.background=element_rect(fill="grey99",colour="grey70", size=0.25),
								#strip.background = element_blank(), 
								axis.text=element_text(size=6),
								panel.margin = unit(0, "mm"),
								
								#legend.background=element_rect(colour="white"),
								legend.key=element_rect(fill="grey99", colour="white"),
								legend.text=element_text(size=7, colour="grey30"),
								legend.title=element_text(size=8, colour="grey30")
								)
						
						ggsave(paste(out2,"constancy_per_nest_with_constancy_labels_tick_sex_color.png", sep=""))
				}
					{# 7x4 with labels -with constancy lables
					dev.new(width=7, height=4.35) # for 4*7
					ggplot(u,aes(x=datetime_c,y=roll_con))+
						#stat_smooth(method="lm", fill="grey95",alpha=1, aes(colour=genus_))+
						#scale_colour_manual(values=go$cols, name="Genus")+
						facet_wrap(~n, nrow=5, labeller=as_labeller(nnn))+
						geom_point(shape=20,size=0.25)+#colour="grey38",aes(fill=genus_))+
						#scale_color_manual(values=c("#FCB42C","#535F7C"))+
						#scale_fill_manual(values=go$cols, name="Genus")+ 
						labs(x="Time in experimental period [h]", y="Incubation constancy")+
						scale_x_continuous(limits=c(-16,16), breaks=seq(-15,15,by=5), labels=c(seq(-15,15,by=5)))+
						geom_vline(xintercept=0, lty=3,lwd=0.5, col='red')+
						#scale_x_continuous(breaks=15, labels=15)+
						scale_y_continuous(limits=c(-0.05,1.05),breaks=seq(0,1,by=0.25), labels=seq(0,1,by=0.25))+
						#coord_cartesian(ylim=c(0,1), xlim=c(-15,15))+ #coord_cartesian(ylim=c(-5,61))+ # needs to be 10 if it is not to touch x-axis
						theme_light()+
						theme(	axis.line=element_line(colour="grey70", size=0.25),
								#panel.border=element_rect(colour="white"),
								panel.border=element_rect(colour="grey70", size=0.25),
								panel.grid = element_blank(),
								
								axis.title=element_text(size=7, colour="grey50"),
								axis.title.y = element_text(vjust=1),
								axis.title.x = element_text(vjust=0.2),
								axis.text=element_text(size=6),# margin=units(0.5,"mm")),
								axis.ticks.length=unit(0.5,"mm"),
								#axis.ticks.margin,
								
								strip.text.x =element_text(size = 6, color="grey30",  margin=margin(1,1,1,1,"mm")), #grey50
								strip.background=element_rect(fill="grey99",colour="grey70", size=0.25),
								#strip.background = element_blank(), 
								axis.text=element_text(size=6),
								panel.margin = unit(0, "mm"),
								
								#legend.background=element_rect(colour="white"),
								legend.key=element_rect(fill="grey99", colour="white"),
								legend.text=element_text(size=7, colour="grey30"),
								legend.title=element_text(size=8, colour="grey30")
								)
						
						ggsave(paste(out2,"constancy_per_nest_with_constancy_labels__.png", sep=""))
				}
							
				{# 7x4 with nest labels
				dev.new(width=7.4, height=4.35) # for 4*7
				ggplot(u,aes(x=datetime_c,y=roll_con))+
						#stat_smooth(method="lm", fill="grey95",alpha=1, aes(colour=genus_))+
						#scale_colour_manual(values=go$cols, name="Genus")+
						facet_wrap(~nest, nrow=4)+
						geom_point(shape=20,size=0.25)+#colour="grey38",aes(fill=genus_))+
						#scale_fill_manual(values=go$cols, name="Genus")+ 
						labs(x="Time in experimental period [h]", y="Incubation constancy")+
						scale_x_continuous(limits=c(-16,16), breaks=seq(-15,15,by=5), labels=c(seq(-15,15,by=5)))+
						geom_vline(xintercept=0, lty=3,lwd=0.5, col='red')+
						#scale_x_continuous(breaks=15, labels=15)+
						scale_y_continuous(limits=c(-0.05,1.05),breaks=seq(0,1,by=0.25), labels=seq(0,1,by=0.25))+
						#coord_cartesian(ylim=c(0,1), xlim=c(-15,15))+ #coord_cartesian(ylim=c(-5,61))+ # needs to be 10 if it is not to touch x-axis
						theme_light()+
						theme(	axis.line=element_line(colour="grey70", size=0.25),
								#panel.border=element_rect(colour="white"),
								panel.border=element_rect(colour="grey70", size=0.25),
								panel.grid = element_blank(),
								
								axis.title=element_text(size=7, colour="grey50"),
								axis.title.y = element_text(vjust=1),
								axis.title.x = element_text(vjust=0.2),
								axis.text=element_text(size=6),# margin=units(0.5,"mm")),
								axis.ticks.length=unit(0.5,"mm"),
								#axis.ticks.margin,
								
								strip.text.x =element_text(size = 6, color="grey30",  margin=margin(1,1,1,1,"mm")), #grey50
								strip.background=element_rect(fill="grey99",colour="grey70", size=0.25),
								#strip.background = element_blank(), 
								#strip.text = element_blank(),
								panel.margin = unit(0, "mm"),
								
								#legend.background=element_rect(colour="white"),
								legend.key=element_rect(fill="grey99", colour="white"),
								legend.text=element_text(size=7, colour="grey30"),
								legend.title=element_text(size=8, colour="grey30")
								)
						
						ggsave(paste(out2,"constancy_per_nest_with_nest_labels_tick.png", sep=""))
				}
				
				{# 7x4 without labels
				dev.new(width=7.4, height=4.35) # for 4*7
				ggplot(u,aes(x=datetime_c,y=roll_con))+
						#stat_smooth(method="lm", fill="grey95",alpha=1, aes(colour=genus_))+
						#scale_colour_manual(values=go$cols, name="Genus")+
						facet_wrap(~n, nrow=4)+
						geom_point(shape=20,size=0.25)+#colour="grey38",aes(fill=genus_))+
						#scale_fill_manual(values=go$cols, name="Genus")+ 
						labs(x="Time in experimental period [h]", y="Incubation constancy")+
						scale_x_continuous(limits=c(-16,16), breaks=seq(-15,15,by=5), labels=c(seq(-15,15,by=5)))+
						geom_vline(xintercept=0, lty=3,lwd=0.5, col='red')+
						#scale_x_continuous(breaks=15, labels=15)+
						scale_y_continuous(limits=c(-0.05,1.05),breaks=seq(0,1,by=0.25), labels=seq(0,1,by=0.25))+
						#coord_cartesian(ylim=c(0,1), xlim=c(-15,15))+ #coord_cartesian(ylim=c(-5,61))+ # needs to be 10 if it is not to touch x-axis
						theme_light()+
						theme(	axis.line=element_line(colour="grey70", size=0.25),
								#panel.border=element_rect(colour="white"),
								panel.border=element_rect(colour="grey70", size=0.25),
								panel.grid = element_blank(),
								
								axis.title=element_text(size=7, colour="grey50"),
								axis.title.y = element_text(vjust=1),
								axis.title.x = element_text(vjust=0.2),
								axis.text=element_text(size=6),# margin=units(0.5,"mm")),
								axis.ticks.length=unit(0.5,"mm"),
								#axis.ticks.margin,
								
								#strip.text.x =element_text(size = 6, color="grey30",  margin=margin(1,1,1,1,"mm")), #grey50
								#strip.background=element_rect(fill="grey99",colour="grey70", size=0.25),
								strip.background = element_blank(), 
								strip.text = element_blank(),
								panel.margin = unit(0, "mm"),
								
								#legend.background=element_rect(colour="white"),
								legend.key=element_rect(fill="grey99", colour="white"),
								legend.text=element_text(size=7, colour="grey30"),
								legend.title=element_text(size=8, colour="grey30")
								)
						
						ggsave(paste(out2,"constancy_per_nest_ordered.png", sep=""))
				}
				{# 3x9 without labels
				dev.new(width=3.3, height=9.5) # for 4*7
				ggplot(u,aes(x=datetime_c,y=roll_con))+
						#stat_smooth(method="lm", fill="grey95",alpha=1, aes(colour=genus_))+
						#scale_colour_manual(values=go$cols, name="Genus")+
						facet_wrap(~n, nrow=9)+
						geom_point(shape=20,size=0.25)+#colour="grey38",aes(fill=genus_))+
						#scale_fill_manual(values=go$cols, name="Genus")+ 
						labs(x="Time in experimental period [h]", y="Incubation constancy")+
						scale_x_continuous(limits=c(-16,16), breaks=seq(-15,15,by=5), labels=c(seq(-15,15,by=5)))+
						geom_vline(xintercept=0, lty=3,lwd=0.5, col='red')+
						#scale_x_continuous(breaks=15, labels=15)+
						scale_y_continuous(limits=c(-0.05,1.05),breaks=seq(0,1,by=0.25), labels=seq(0,1,by=0.25))+
						#coord_cartesian(ylim=c(0,1), xlim=c(-15,15))+ #coord_cartesian(ylim=c(-5,61))+ # needs to be 10 if it is not to touch x-axis
						theme_light()+
						theme(	axis.line=element_line(colour="grey70", size=0.25),
								#panel.border=element_rect(colour="white"),
								panel.border=element_rect(colour="grey70", size=0.25),
								panel.grid = element_blank(),
								
								axis.title=element_text(size=7, colour="grey50"),
								axis.title.y = element_text(vjust=1),
								axis.title.x = element_text(vjust=0.2),
								axis.text=element_text(size=6),# margin=units(0.5,"mm")),
								axis.ticks.length=unit(0.5,"mm"),
								#axis.ticks.margin,
								
								#strip.text.x =element_text(size = 6, color="grey30",  margin=margin(1,1,1,1,"mm")), #grey50
								#strip.background=element_rect(fill="grey99",colour="grey70", size=0.25),
								strip.background = element_blank(), 
								strip.text = element_blank(),
								panel.margin = unit(0, "mm"),
								
								#legend.background=element_rect(colour="white"),
								legend.key=element_rect(fill="grey99", colour="white"),
								legend.text=element_text(size=7, colour="grey30"),
								legend.title=element_text(size=8, colour="grey30")
								)
				 
				 	ggsave(paste(out2,"constancy_per_nest_9x3_ascending.png", sep=""))
				}
			}
		}
		
		# old
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
		{# prepare color palet
			n <- 30
			
			dev.new(width=12, height=4)
			par(mfrow=c(1,3))
				qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
				col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
			col_1=sample(col_vector, n)
			pie(rep(1,n), col=col_1)
			
			col_2 <- distinctColorPalette(n)
			pie(rep(1, n), col=col_2)
			
			col_3=iwanthue(n=30,hmin=0, hmax=360, cmin=30, cmax=80, lmin=35, lmax=80)
			pie(rep(1,n),col=iwanthue(n=30,hmin=0, hmax=360, cmin=30, cmax=80, lmin=35, lmax=80))
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
			nn$bt=round(bt$inc_eff[match(nn$nest,bt$nest)],2)
			nn$bb=round(bb$inc_eff[match(nn$nest,bb$nest)],2)
			nn$bt=ifelse(nchar(nn$bt)==3, paste(nn$bt,"0",sep=""), nn$bt)
			nn$bb=ifelse(nchar(nn$bb)==3, paste(nn$bb,"0",sep=""), nn$bb)
			nn$label=paste(nn$bb,nn$bt,sep="       ")
			nn$sex=bt$sex[match(nn$nest,bt$nest)]
			nn$col_=ifelse(nn$sex=="f", "#FCB42C","#535F7C")
			nnn=nn$label
			names(nnn)=nn$n
		}	
		{# all in one or in separate plots
		type='png1'		#type='dev', type='png'		
		for (j in 1:length(nests_)){
						print(nests_[j])
						pp=paste("extract_2013.mb_inc_2013_",nests_[j],sep="")	
						b=wet(conMy, paste("SELECT pk, datetime_, inc_t, incubation FROM  ", pp," where bout_ID_tr in ('50')"))
						b$datetime_= as.POSIXct(b$datetime_)
						b$datetime_c=as.numeric(difftime(b$datetime_,b$datetime_[1]+60*dd$bout[dd$nest==nests_[j]], units='hours'))
						
						if(nests_[j]=='s806'){b$inc_t=b$incubation
											  b$inc_t[is.na(b$inc_t)]=0}
						
						b$roll_con=rollmean(b$inc_t, 12*60, align="center", fill="extend") # 12=1min
						
						if(type=='png'){
							png(paste(out,nests_[j],".png", sep=""), width=3.5,height=2.5,units="in",res=600)
							par(mar=c(2,2,0.7,0.7), ps=12, mgp=c(1.2,0.35,0), las=1, cex.lab=0.7, cex.axis=0.6, tcl=-0.15,bty="l",xpd=TRUE, col.axis="grey30", col.lab="grey30", col.main="grey30", fg="grey70")
							plot(b$roll_con~b$datetime_c, type='n', pch=19,xlim=c(-15,15), ylim=c(0,1), ylab='Incubation constancy', xlab=NA, xaxt='n')
							axis(1, at=seq(-15,15,by=5), labels=FALSE)
							text(seq(-15,15,by=5), par("usr")[3]-0.05, labels = seq(-15,15,by=5),  xpd = TRUE, cex=0.6, col="grey30") #labels = z_g$genus
							mtext('Time [h]',side=1,line=0.6, cex=0.7, las=1, col='grey30')
							points(b$roll_con~b$datetime_c, cex=0.1, pch=20, col=col_1[j])
							abline(v=0, lty=3, col='grey',xpd=FALSE)
							dev.off()
							print(nests_[j])
						}else{
							if(j==1){
								if(type=='png1'){
									png(paste(out2,"constancy.png", sep=""), width=3.5,height=2.5,units="in",res=600)}else{dev.new(width=3.5,height=2.5)}
								par(mar=c(2,2,0.7,0.7), ps=12, mgp=c(1.2,0.35,0), las=1, cex.lab=0.7, cex.axis=0.6, tcl=-0.15,bty="l",xpd=TRUE, col.axis="grey30", col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt
								plot(b$roll_con~b$datetime_c, type='n', pch=19,xlim=c(-15,15), ylim=c(0,1), ylab='Incubation constancy', xlab=NA, xaxt='n')
								axis(1, at=seq(-15,15,by=5), labels=FALSE)
								text(seq(-15,15,by=5), par("usr")[3]-0.05, labels = seq(-15,15,by=5),  xpd = TRUE, cex=0.6, col="grey30") 
								mtext('Time [h]',side=1,line=0.6, cex=0.7, las=1, col='grey30')
								}
							points(b$roll_con~b$datetime_c, cex=0.1, pch=20, col=col_1[j])
							print(nests_[j])
							}	
						}						
			if(type=='png1'){abline(v=0, lty=3, col='grey',xpd=FALSE)
							 dev.off()}else{abline(v=0, lty=3, col='grey',xpd=FALSE)}
		} 
}
    {# constancy over experimental period until bird's returngraph
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

		{# all in plot with facets
		
			{# prepare nest datasets	
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
				
								
				}
			{# plot
				summary(u)
				{# 7x4 with labels -with constancy lables and sex label, black and grey line
					nn$label=paste(nn$bb,nn$bt,"   ",nn$sex_symbol,sep="        ")
					nn$label=paste("     ",nn$bb,"        ",nn$bt,"                    ",nn$sex_symbol,sep="")
					
					# label constancy
					nnn=nn$label
					names(nnn)=nn$n
					
					# label nest
					nnnest=nn$nest
					names(nnnest)=nn$n
				
				
				#dev.new(width=7.4, height=4.35) # for 4*7
				dev.new(width=7, height=4.35) # for 4*7
				  ggplot(u,aes(x=datetime_c,y=roll_con, fill=col_line))+
						#stat_smooth(method="lm", fill="grey95",alpha=1, aes(colour=genus_))+
						#scale_colour_manual(values=go$cols, name="Genus")+
						#facet_wrap(~n, nrow=5, labeller=as_labeller(nnnest))+ #label is nests
						facet_wrap(~n, nrow=5, labeller=as_labeller(nnn))+ #label is constancy
						geom_point(aes(color=col_line), shape=20,size=0.25)+#colour="grey38",aes(fill=genus_))+
						scale_color_manual(values=c("black","gray90"), guide=FALSE)+
						#scale_fill_manual(values=go$cols, name="Genus")+ 
						labs(x="Time in experimental period [h]", y="Incubation constancy")+
						scale_x_continuous(limits=c(-16,33), breaks=seq(-15,30,by=15), labels=c(seq(-15,30,by=15)))+
						geom_vline(xintercept=0, lty=3,lwd=0.5, col='red')+
						#geom_vline(data=e,aes(xintercept=r_time), lty=3,lwd=0.5, col='red')+
						geom_text(data = ann_text,aes(label =lab),size=1.9, angle = 90, colour='grey50',hjust='middle',vjust='middle')+
						#scale_x_continuous(breaks=15, labels=15)+
						scale_y_continuous(limits=c(-0.05,1.05),breaks=seq(0,1,by=0.25), labels=seq(0,1,by=0.25))+
						#coord_cartesian(ylim=c(0,1), xlim=c(-15,15))+ #coord_cartesian(ylim=c(-5,61))+ # needs to be 10 if it is not to touch x-axis
						theme_light()+
						theme(	axis.line=element_line(colour="grey70", size=0.25),
								#panel.border=element_rect(colour="white"),
								panel.border=element_rect(colour="grey70", size=0.25),
								panel.grid = element_blank(),
								
								axis.title=element_text(size=7, colour="grey50"),
								axis.title.y = element_text(vjust=1),
								axis.title.x = element_text(vjust=0.2),
								axis.text=element_text(size=6),# margin=units(0.5,"mm")),
								axis.ticks.length=unit(0.5,"mm"),
								#axis.ticks.margin,
								
								strip.text.x =element_text(size = 6, color="grey30",  margin=margin(1,1,1,1,"mm")), #grey50
								strip.background=element_rect(fill="grey99",colour="grey70", size=0.25),
								#strip.background = element_blank(), 
								#strip.text = element_blank(),
								panel.margin = unit(0, "mm"),
								legend.position="none"
								#legend.background=element_rect(colour="white"),
								#legend.key=element_rect(fill="grey99", colour="white"),
								#legend.text=element_text(size=7, colour="grey30"),
								#legend.title=element_text(size=8, colour="grey30")
								)
						
						ggsave(paste(out2,"constancy_per_nest_with_constancy_sex_labels_after_period_in_grey_smaller.png", sep=""))
						#ggsave(paste(out2,"constancy_per_nest_with_nest_labels_sex_after_period_3.png", sep=""))
				}
											
				{# 7x4 with labels -with constancy lables and color sex
				nn$label=paste(nn$bb,nn$bt,"              ",sep="       ")
				# label constancy
				nnn=nn$label
				names(nnn)=nn$n
				
				# label nest
				nnnest=nn$nest
				names(nnnest)=nn$n
				
				dev.new(width=7.4, height=4.35) # for 4*7
				  ggplot(u,aes(x=datetime_c,y=roll_con, fill=sex))+
						#stat_smooth(method="lm", fill="grey95",alpha=1, aes(colour=genus_))+
						#scale_colour_manual(values=go$cols, name="Genus")+
						facet_wrap(~n, nrow=5, labeller=as_labeller(nnnest))+ #label is nests
						#facet_wrap(~n, nrow=5, labeller=as_labeller(nnn))+ #label is constancy
						geom_point(aes(color=sex), shape=20,size=0.25)+#colour="grey38",aes(fill=genus_))+
						scale_color_manual(values=c("#FCB42C","#535F7C"), guide=FALSE)+
						#scale_fill_manual(values=go$cols, name="Genus")+ 
						labs(x="Time in experimental period [h]", y="Incubation constancy")+
						scale_x_continuous(limits=c(-16,33), breaks=seq(-15,30,by=15), labels=c(seq(-15,30,by=15)))+
						geom_vline(xintercept=0, lty=3,lwd=0.5, col='red')+
						geom_vline(data=e,aes(xintercept=r_time), lty=3,lwd=0.5, col='red')+
						geom_text(data = ann_text,aes(label =lab),size=1.9, angle = 90, colour='grey50',hjust='middle',vjust='middle')+
						#scale_x_continuous(breaks=15, labels=15)+
						scale_y_continuous(limits=c(-0.05,1.05),breaks=seq(0,1,by=0.25), labels=seq(0,1,by=0.25))+
						#coord_cartesian(ylim=c(0,1), xlim=c(-15,15))+ #coord_cartesian(ylim=c(-5,61))+ # needs to be 10 if it is not to touch x-axis
						theme_light()+
						theme(	axis.line=element_line(colour="grey70", size=0.25),
								#panel.border=element_rect(colour="white"),
								panel.border=element_rect(colour="grey70", size=0.25),
								panel.grid = element_blank(),
								
								axis.title=element_text(size=7, colour="grey50"),
								axis.title.y = element_text(vjust=1),
								axis.title.x = element_text(vjust=0.2),
								axis.text=element_text(size=6),# margin=units(0.5,"mm")),
								axis.ticks.length=unit(0.5,"mm"),
								#axis.ticks.margin,
								
								strip.text.x =element_text(size = 6, color="grey30",  margin=margin(1,1,1,1,"mm")), #grey50
								strip.background=element_rect(fill="grey99",colour="grey70", size=0.25),
								#strip.background = element_blank(), 
								#strip.text = element_blank(),
								panel.margin = unit(0, "mm"),
								legend.position="none"
								#legend.background=element_rect(colour="white"),
								#legend.key=element_rect(fill="grey99", colour="white"),
								#legend.text=element_text(size=7, colour="grey30"),
								#legend.title=element_text(size=8, colour="grey30")
								)
						
						ggsave(paste(out2,"constancy_per_nest_with_constancy_labels_sex_after_period_3.png", sep=""))
						ggsave(paste(out2,"constancy_per_nest_with_nest_labels_sex_after_period_3.png", sep=""))
				}
			
				{# 7x4 with nest labels
				dev.new(width=7.4, height=4.35) # for 4*7
				ggplot(u,aes(x=datetime_c,y=roll_con))+
						#stat_smooth(method="lm", fill="grey95",alpha=1, aes(colour=genus_))+
						#scale_colour_manual(values=go$cols, name="Genus")+
						facet_wrap(~nest, nrow=4)+
						geom_point(shape=20,size=0.25)+#colour="grey38",aes(fill=genus_))+
						#scale_fill_manual(values=go$cols, name="Genus")+ 
						labs(x="Time in experimental period [h]", y="Incubation constancy")+
						scale_x_continuous(limits=c(-16,16), breaks=seq(-15,15,by=5), labels=c(seq(-15,15,by=5)))+
						geom_vline(xintercept=0, lty=3,lwd=0.5, col='red')+
						#scale_x_continuous(breaks=15, labels=15)+
						scale_y_continuous(limits=c(-0.05,1.05),breaks=seq(0,1,by=0.25), labels=seq(0,1,by=0.25))+
						#coord_cartesian(ylim=c(0,1), xlim=c(-15,15))+ #coord_cartesian(ylim=c(-5,61))+ # needs to be 10 if it is not to touch x-axis
						theme_light()+
						theme(	axis.line=element_line(colour="grey70", size=0.25),
								#panel.border=element_rect(colour="white"),
								panel.border=element_rect(colour="grey70", size=0.25),
								panel.grid = element_blank(),
								
								axis.title=element_text(size=7, colour="grey50"),
								axis.title.y = element_text(vjust=1),
								axis.title.x = element_text(vjust=0.2),
								axis.text=element_text(size=6),# margin=units(0.5,"mm")),
								axis.ticks.length=unit(0.5,"mm"),
								#axis.ticks.margin,
								
								strip.text.x =element_text(size = 6, color="grey30",  margin=margin(1,1,1,1,"mm")), #grey50
								strip.background=element_rect(fill="grey99",colour="grey70", size=0.25),
								#strip.background = element_blank(), 
								#strip.text = element_blank(),
								panel.margin = unit(0, "mm"),
								
								#legend.background=element_rect(colour="white"),
								legend.key=element_rect(fill="grey99", colour="white"),
								legend.text=element_text(size=7, colour="grey30"),
								legend.title=element_text(size=8, colour="grey30")
								)
						
						ggsave(paste(out2,"constancy_per_nest_with_nest_labels_tick.png", sep=""))
				}
				
			}
		}
		}

	{# effect of time
	  {# run first
		b$bout_start=as.POSIXct(b$bout_start)
		b$exper=as.factor(b$exper)
		# time as start of bout
			#b$time = as.numeric(difftime(b$bout_start,trunc(b$bout_start,"day"),   #DO NOT USE rfid_T$day here!! (the time will be wrong)
			#					units = "hours"))
		# time as midpoint of bout
			b$midbout=b$bout_start+60*60*b$bout_length/2
			b$time = as.numeric(difftime(b$midbout,trunc(b$midbout,"day"),   #DO NOT USE rfid_T$day here!! (the time will be wrong)
								units = "hours"))
		b$rad=b$time*pi/12
		b$sin_=sin(b$rad)
		b$cos_=cos(b$rad)
		
		bt_=b[b$exper=='t',]
		bc_=b[b$exper=='c',]
		col_t='#FCB42C'
		col_t_light='#fef0d4'
		col_c='#535F7C'
		col_c_light='#cbcfd7'
		
		ggplot(b,aes(x=rad, y=inc_eff, col=sex,shape=exper))+geom_point()
		}
	  
	  {# c+t in one
	    {# c+t in one model
			#m=lm(inc_eff~sin(rad)+cos(rad), b[b$exper=='t',])
			m=lmer(inc_eff~sin_*exper+cos_*exper+(1|nest), b)
			summary(m)
			plot(allEffects(m))
		
			# simulation		
				nsim <- 2000
				bsim <- sim(m, n.sim=nsim)  
				apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
			
			# coefficients
				v <- fixef(m)
			# predicted values		
				newD=data.frame(time_=seq(0,24,0.25))
						newD$rad=2*pi*newD$time_ / 24
						newD$sin_=sin(newD$rad)
						newD$cos_=cos(newD$rad)
						newD$exper='t'
				n1=newD		
				newD=data.frame(time_=seq(0,24,0.25))
						newD$rad=2*pi*newD$time_ / 24
						newD$sin_=sin(newD$rad)
						newD$cos_=cos(newD$rad)
						newD$exper='c'	
				newD=rbind(n1,newD)		
			# exactly the model which was used has to be specified here 
				X <- model.matrix(~ sin_*exper + cos_*exper,data=newD)	
							
			# calculate predicted values and creditability intervals
				newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
						predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
						for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@fixef[i,]
						newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
						newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
						#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
						#newD=newD[order(newD$t_tundra),]
					pc=newD[newD$exper=='c',]
					pp=newD[newD$exper=='t',]
		}			
	    {# plot			
			dev.new(width=2.5,height=2.5)
			png(paste(out2,"efficiency_time_start_single_model.png", sep=""), width=2.3,height=3.5,units="in",res=600)
			par(mar=c(2,2,0.7,0.7), ps=12, mgp=c(1.2,0.35,0), las=1, cex.lab=0.7, cex.axis=0.6, tcl=-0.15,bty="l",xpd=TRUE, col.axis="grey30", col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt
			
	  plot(pp$pred~pp$time_, pch=19,xlim=c(0,24), ylim=c(0,1), ylab='Incubation constancy', xlab=NA, xaxt='n', type='n')
		axis(1, at=seq(0,24,by=4), labels=FALSE)
		text(seq(0,24,by=4), par("usr")[3]-0.05, labels = seq(0,24,by=4),  xpd = TRUE, cex=0.6, col="grey30") #labels = z_g$genus
		mtext('Time [h]',side=1,line=0.6, cex=0.7, las=1, col='grey30')
		
		polygon(c(pc$time_, rev(pc$time_)), c(pc$lwr, 
							rev(pc$upr)), border=NA, col=col_c_light) #0,0,0 black 0.5 is transparents RED
		lines(pc$time_, pc$pred, col=col_c,lwd=2)
		


		polygon(c(pp$time_, rev(pp$time_)), c(pp$lwr, 
							rev(pp$upr)), border=NA, col=col_t_light) #0,0,0 black 0.5 is transparents RED
		lines(pp$time_, pp$pred, col=col_t,lwd=2)
		
		#points(bt_$time, bt_$inc_eff, col=col_t, pch=21, bg=col_t_light, cex=0.8)
		#points(bc_$time, bc_$inc_eff, col=col_c,  pch=21, bg=col_c_light, cex=0.8)
		#legend("bottom", legend=c('treated','control'), pch=21, cex=0.6,pt.cex=0.8, bty='n', text.col="grey50", bg=c(col_t_light, col_c_light), col=c(col_t,col_c))	
		
		points(bt_$time, bt_$inc_eff, col=col_t, pch=20)
		points(bc_$time, bc_$inc_eff, col=col_c, pch=20)
		legend("bottom", legend=c('treated','control'), pch=20, col=c( '#FCB42C','#535F7C'), cex=0.6,pt.cex=1, bty='n', text.col="grey50")	
		
		dev.off()
			}
	  }
	  {# separate models
		{# treated
			#m=lm(inc_eff~sin(rad)+cos(rad), b[b$exper=='t',])
			m=lm(inc_eff~sin_+cos_, b[b$exper=='t',])
			ms=lm(inc_eff~sin_*sex+cos_*sex, b[b$exper=='t',])
			AIC(m,ms)
			summary(m)
			plot(allEffects(m))
		
			# simulation		
				nsim <- 2000
				bsim <- sim(m, n.sim=nsim)  
				apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
			# coefficients
				v <- coef(m)
			# predicted values		
				newD=data.frame(time_=seq(0,24,0.25))
						newD$rad=2*pi*newD$time_ / 24
						newD$sin_=sin(newD$rad)
						newD$cos_=cos(newD$rad)
			# exactly the model which was used has to be specified here 
				X <- model.matrix(~ sin_ + cos_,data=newD)	
							
			# calculate predicted values and creditability intervals
				newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
						predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
						for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@coef[i,]
						newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
						newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
						#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
						#newD=newD[order(newD$t_tundra),]
					pp=newD
		}			
		{# control
			#m=lm(inc_eff~sin(rad)+cos(rad), b[b$exper=='t',])
			m=lm(inc_eff~sin_+cos_, b[b$exper=='c',])
			summary(m)
			plot(allEffects(m))
		
			# simulation		
				nsim <- 2000
				bsim <- sim(m, n.sim=nsim)  
				apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
			# coefficients
				v <- coef(m)
			# predicted values		
				newD=data.frame(time_=seq(0,24,0.25))
						newD$rad=2*pi*newD$time_ / 24
						newD$sin_=sin(newD$rad)
						newD$cos_=cos(newD$rad)
			# exactly the model which was used has to be specified here 
				X <- model.matrix(~ sin_ + cos_,data=newD)	
							
			# calculate predicted values and creditability intervals
				newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
						predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
						for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@coef[i,]
						newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
						newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
						#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
						#newD=newD[order(newD$t_tundra),]
					pc=newD
		}			
		{# plot			
			dev.new(width=2.5,height=2.5)
			png(paste(out2,"efficiency_time_start.png", sep=""), width=2.3,height=3.5,units="in",res=600)
			par(mar=c(2,2,0.7,0.7), ps=12, mgp=c(1.2,0.35,0), las=1, cex.lab=0.7, cex.axis=0.6, tcl=-0.15,bty="l",xpd=TRUE, col.axis="grey30", col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt
			
	  plot(pp$pred~pp$time_, pch=19,xlim=c(0,24), ylim=c(0,1), ylab='Incubation constancy', xlab=NA, xaxt='n', type='n')
		axis(1, at=seq(0,24,by=4), labels=FALSE)
		text(seq(0,24,by=4), par("usr")[3]-0.05, labels = seq(0,24,by=4),  xpd = TRUE, cex=0.6, col="grey30") #labels = z_g$genus
		mtext('Time [h]',side=1,line=0.6, cex=0.7, las=1, col='grey30')
		
		polygon(c(pp$time_, rev(pp$time_)), c(pp$lwr, 
							rev(pp$upr)), border=NA, col=col_t_light) #0,0,0 black 0.5 is transparents RED
		lines(pp$time_, pp$pred, col=col_t,lwd=2)
	
		polygon(c(pc$time_, rev(pc$time_)), c(pc$lwr, 
							rev(pc$upr)), border=NA, col=col_c_light) #0,0,0 black 0.5 is transparents RED
		lines(pc$time_, pc$pred, col=col_c,lwd=2)
	
		#points(bt_$time, bt_$inc_eff, col=col_t, pch=21, bg=col_t_light, cex=0.8)
		#points(bc_$time, bc_$inc_eff, col=col_c,  pch=21, bg=col_c_light, cex=0.8)
		#legend("bottom", legend=c('treated','control'), pch=21, cex=0.6,pt.cex=0.8, bty='n', text.col="grey50", bg=c(col_t_light, col_c_light), col=c(col_t,col_c))	
	
		points(bt_$time, bt_$inc_eff, col=col_t, pch=20)
		points(bc_$time, bc_$inc_eff, col=col_c, pch=20)
		legend("bottom", legend=c('treated','control'), pch=20, col=c( '#FCB42C','#535F7C'), cex=0.6,pt.cex=1, bty='n', text.col="grey50")	
		
		dev.off()
			}
		}			
	  {# constancy from 2011
		g=read.csv(paste(wd,'B-bout_length_constancy.csv', sep=""), stringsAsFactors=FALSE)
		g=g[-which(is.na(g$inc_eff)),]
		g$bout_start=as.POSIXct(g$bout_start)
		g$midbout=g$bout_start+60*g$bout_length/2
			g$time = as.numeric(difftime(g$midbout,trunc(g$midbout,"day"),   #DO NOT USE rfid_T$day here!! (the time will be wrong)
								units = "hours"))
		g$rad=g$time*pi/12
		g$sin_=sin(g$rad)
		g$cos_=cos(g$rad)
		
		{# model
			m=lmer(inc_eff~sin_+cos_+(1|nest), g)
			summary(m)
			plot(allEffects(m))
		
			# simulation		
				nsim <- 2000
				bsim <- sim(m, n.sim=nsim)  
				apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
			
			# coefficients
				v <- fixef(m)
			# predicted values		
				newD=data.frame(time_=seq(0,24,0.25))
						newD$rad=2*pi*newD$time_ / 24
						newD$sin_=sin(newD$rad)
						newD$cos_=cos(newD$rad)
				
			# exactly the model which was used has to be specified here 
				X <- model.matrix(~ sin_ + cos_,data=newD)	
							
			# calculate predicted values and creditability intervals
				newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
						predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
						for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@fixef[i,]
						newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
						newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
						#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
						#newD=newD[order(newD$t_tundra),]
				pp=newD
		}
		{# plot
			dev.new(width=2.5,height=2.5)
			png(paste(out2,"efficiency_time_2011.png", sep=""), width=2.3,height=3.5,units="in",res=600)
			par(mar=c(2,2,0.7,0.7), ps=12, mgp=c(1.2,0.35,0), las=1, cex.lab=0.7, cex.axis=0.6, tcl=-0.15,bty="l",xpd=TRUE, col.axis="grey30", col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt
			
		plot(pp$pred~pp$time_, pch=19,xlim=c(0,24), ylim=c(0,1), ylab='Incubation constancy', xlab=NA, xaxt='n', type='n')
		axis(1, at=seq(0,24,by=4), labels=FALSE)
		text(seq(0,24,by=4), par("usr")[3]-0.05, labels = seq(0,24,by=4),  xpd = TRUE, cex=0.6, col="grey30") #labels = z_g$genus
		mtext('Time [h]',side=1,line=0.6, cex=0.7, las=1, col='grey30')
		
		#points(g$time, g$inc_eff, col=adjustcolor(col_c ,alpha.f = 0.4),bg=adjustcolor(col_c ,alpha.f = 0.1), pch=21)# pch=20)
		points(g$time, g$inc_eff, col=adjustcolor(col_c ,alpha.f = 0.4), pch=20,cex=0.1)
		
		polygon(c(pp$time_, rev(pp$time_)), c(pp$lwr, 
							rev(pp$upr)), border=NA, col=adjustcolor(col_c ,alpha.f = 0.5)) #0,0,0 black 0.5 is transparents RED
		lines(pp$time_, pp$pred, col=col_c,lwd=2)
	
		
		#points(bt_$time, bt_$inc_eff, col=col_t, pch=21, bg=col_t_light, cex=0.8)
		#points(bc_$time, bc_$inc_eff, col=col_c,  pch=21, bg=col_c_light, cex=0.8)
		#legend("bottom", legend=c('treated','control'), pch=21, cex=0.6,pt.cex=0.8, bty='n', text.col="grey50", bg=c(col_t_light, col_c_light), col=c(col_t,col_c))	
	
		
				
		dev.off()
		}
	}
	  {# treatment and constancy 2011
		{# prepare 2011
			g=read.csv(paste(wd,'B-bout_length_constancy.csv', sep=""), stringsAsFactors=FALSE)
			g=g[-which(is.na(g$inc_eff)),]
			g$bout_start=as.POSIXct(g$bout_start)
			g$midbout=g$bout_start+60*g$bout_length/2
				g$time = as.numeric(difftime(g$midbout,trunc(g$midbout,"day"),   #DO NOT USE rfid_T$day here!! (the time will be wrong)
									units = "hours"))
			g$rad=g$time*pi/12
			g$sin_=sin(g$rad)
			g$cos_=cos(g$rad)
			
		}
		{# all 2011
		  {# model
			m=lmer(inc_eff~sin_+cos_+(1|nest), g)
			#summary(m)
			#plot(allEffects(m))
		
			# simulation		
				nsim <- 2000
				bsim <- sim(m, n.sim=nsim)  
				apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
			
			# coefficients
				v <- fixef(m)
			# predicted values		
				newD=data.frame(time_=seq(0,24,0.25))
						newD$rad=2*pi*newD$time_ / 24
						newD$sin_=sin(newD$rad)
						newD$cos_=cos(newD$rad)
				
			# exactly the model which was used has to be specified here 
				X <- model.matrix(~ sin_ + cos_,data=newD)	
							
			# calculate predicted values and creditability intervals
				newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
						predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
						for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@fixef[i,]
						newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
						newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
						#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
						#newD=newD[order(newD$t_tundra),]
				pb=newD
		}
	      {# treated
			#m=lm(inc_eff~sin(rad)+cos(rad), b[b$exper=='t',])
			m=lm(inc_eff~sin_+cos_, b[b$exper=='t',])
			summary(m)
			plot(allEffects(m))
		
			# simulation		
				nsim <- 2000
				bsim <- sim(m, n.sim=nsim)  
				apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
			# coefficients
				v <- coef(m)
			# predicted values		
				newD=data.frame(time_=seq(0,24,0.25))
						newD$rad=2*pi*newD$time_ / 24
						newD$sin_=sin(newD$rad)
						newD$cos_=cos(newD$rad)
			# exactly the model which was used has to be specified here 
				X <- model.matrix(~ sin_ + cos_,data=newD)	
							
			# calculate predicted values and creditability intervals
				newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
						predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
						for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@coef[i,]
						newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
						newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
						#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
						#newD=newD[order(newD$t_tundra),]
					pp=newD
		}			
	      {# plot			
			dev.new(width=2.5,height=2.5)
			#png(paste(out2,"efficiency_time_treated_2011.png", sep=""), width=2.3,height=3.5,units="in",res=600)
			par(mar=c(2,2,0.7,0.7), ps=12, mgp=c(1.2,0.35,0), las=1, cex.lab=0.7, cex.axis=0.6, tcl=-0.15,bty="l",xpd=TRUE, col.axis="grey30", col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt
			
			plot(pp$pred~pp$time_, pch=19,xlim=c(0,24), ylim=c(0,1), ylab='Incubation constancy', xlab=NA, xaxt='n', type='n')
			axis(1, at=seq(0,24,by=4), labels=FALSE)
			text(seq(0,24,by=4), par("usr")[3]-0.05, labels = seq(0,24,by=4),  xpd = TRUE, cex=0.6, col="grey30") #labels = z_g$genus
			mtext('Time [h]',side=1,line=0.6, cex=0.7, las=1, col='grey30')
			
			points(g$time, g$inc_eff, col=adjustcolor(col_c ,alpha.f = 0.4), pch=20,cex=0.1)
			polygon(c(pb$time_, rev(pb$time_)), c(pb$lwr, 
			rev(pb$upr)), border=NA, col=adjustcolor(col_c ,alpha.f = 0.4)) #0,0,0 black 0.5 is transparents RED
			lines(pb$time_, pb$pred, col=col_c,lwd=2)
					
			polygon(c(pp$time_, rev(pp$time_)), c(pp$lwr, 
								rev(pp$upr)), border=NA, col=adjustcolor(col_t ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
			lines(pp$time_, pp$pred, col=col_t,lwd=2)
		
			points(bt_$time, bt_$inc_eff, col=col_t, pch=20)
			
			legend("bottom", legend=c('treated','control'), pch=20, col=c( '#FCB42C','#535F7C'), cex=0.6,pt.cex=1, bty='n', text.col="grey50")	
			
			dev.off()
				}
		}
		{# 2011 limited to mid incubation
			g=g[which(g$bout_start_j>9 & g$bout_start_j<13),] # limit to mid incubation
		  {# model
			m=lmer(inc_eff~sin_+cos_+(1|nest), g)
			#summary(m)
			#plot(allEffects(m))
		
			# simulation		
				nsim <- 2000
				bsim <- sim(m, n.sim=nsim)  
				apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
			
			# coefficients
				v <- fixef(m)
			# predicted values		
				newD=data.frame(time_=seq(0,24,0.25))
						newD$rad=2*pi*newD$time_ / 24
						newD$sin_=sin(newD$rad)
						newD$cos_=cos(newD$rad)
				
			# exactly the model which was used has to be specified here 
				X <- model.matrix(~ sin_ + cos_,data=newD)	
							
			# calculate predicted values and creditability intervals
				newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
						predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
						for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@fixef[i,]
						newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
						newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
						#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
						#newD=newD[order(newD$t_tundra),]
				pb=newD
		  }
	      {# treated
			#m=lm(inc_eff~sin(rad)+cos(rad), b[b$exper=='t',])
			m=lm(inc_eff~sin_+cos_, b[b$exper=='t',])
			summary(m)
			plot(allEffects(m))
		
			# simulation		
				nsim <- 2000
				bsim <- sim(m, n.sim=nsim)  
				apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))	
			# coefficients
				v <- coef(m)
			# predicted values		
				newD=data.frame(time_=seq(0,24,0.25))
						newD$rad=2*pi*newD$time_ / 24
						newD$sin_=sin(newD$rad)
						newD$cos_=cos(newD$rad)
			# exactly the model which was used has to be specified here 
				X <- model.matrix(~ sin_ + cos_,data=newD)	
							
			# calculate predicted values and creditability intervals
				newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
						predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
						for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@coef[i,]
						newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
						newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
						#newD$other <- apply(predmatrix, 1, quantile, prob=0.5)
						#newD=newD[order(newD$t_tundra),]
					pp=newD
		}			
	      {# plot			
			dev.new(width=2.5,height=2.5)
			#png(paste(out2,"efficiency_time_treated_2011.png", sep=""), width=2.3,height=3.5,units="in",res=600)
			#png(paste(out2,"efficiency_time_treated_2011_limited_mid_inc_2.png", sep=""), width=2.3,height=3.5,units="in",res=600)
			par(mar=c(2,2,0.7,0.7), ps=12, mgp=c(1.2,0.35,0), las=1, cex.lab=0.7, cex.axis=0.6, tcl=-0.15,bty="l",xpd=TRUE, col.axis="grey30", col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt
			
			plot(pp$pred~pp$time_, pch=19,xlim=c(0,24), ylim=c(0,1), ylab='Incubation constancy', xlab=NA, xaxt='n', type='n')
			axis(1, at=seq(0,24,by=4), labels=FALSE)
			text(seq(0,24,by=4), par("usr")[3]-0.05, labels = seq(0,24,by=4),  xpd = TRUE, cex=0.6, col="grey30") #labels = z_g$genus
			mtext('Time [h]',side=1,line=0.6, cex=0.7, las=1, col='grey30')
			
			points(g$time, g$inc_eff, col=adjustcolor(col_c ,alpha.f = 0.4), pch=20)#points(g$time, g$inc_eff, col=adjustcolor(col_c ,alpha.f = 0.4), pch=20,cex=0.1)
			polygon(c(pb$time_, rev(pb$time_)), c(pb$lwr, 
			rev(pb$upr)), border=NA, col=adjustcolor(col_c ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
			lines(pb$time_, pb$pred, col=col_c,lwd=2)
					
			polygon(c(pp$time_, rev(pp$time_)), c(pp$lwr, 
								rev(pp$upr)), border=NA, col=adjustcolor(col_t ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
			lines(pp$time_, pp$pred, col=col_t,lwd=2)
		
			points(bt_$time, bt_$inc_eff, col=col_t, pch=20)
			
			legend("bottom", legend=c('treated','control'), pch=20, col=c( '#FCB42C','#535F7C'), cex=0.6,pt.cex=1, bty='n', text.col="grey50")	
			
			dev.off()
			}
		}
		}
	}
	{# share of incubation in last three days prior to treatment
			{# prepare data
			load(file=paste(wd,'experimental.Rdata', sep=""))
			b=b[-which(b$nest%in%c('s410','s514','s524','s624')),]
			b$inc_eff[b$nest=='s806']=b$inc_eff_2[b$nest=='s806']
			b$bout_start=as.POSIXct(b$bout_start)
			b$exper=as.factor(b$exper)
			bc=b[b$exper=='c',]
			bt=b[b$exper=='t',]
			l=list()
			for(i in 1:nrow(bc)){
								x=bc[bc$nest==bc$nest[i],]
								day=as.Date(trunc(x$bout_start, "day"))
								pp=paste("extract_2013.mb_inc_2013_",bc$nest[i],sep="")	
								if('s402'==b$nest[i]){
								y=wet(conMy, paste("SELECT pk, datetime_, inc_t, incubation, sex_inc_filled FROM  ", pp," where exchange_gap is NULL and datetime_<'",day,"' and datetime_>='",day-4,"'", sep=""))
								y$datetime_= as.POSIXct(y$datetime_)
								y$day=as.Date(trunc(y$datetime_, "day"))
								y=y[y$day!=as.Date("2013-06-20"),]
								}else{
									y=wet(conMy, paste("SELECT pk, datetime_, inc_t, incubation, sex_inc_filled FROM  ", pp," where exchange_gap is NULL and datetime_<'",day,"' and datetime_>='",day-3,"'", sep=""))
									y$datetime_= as.POSIXct(y$datetime_)
									y$day=as.Date(trunc(y$datetime_, "day"))
									}
								l[[i]]=ddply(y,.(day), summarise, nest=x$nest, p_f=length(sex_inc_filled[sex_inc_filled=='f'])/length(sex_inc_filled))
									#length(y$sex_inc_filled[y$sex_inc_filled=='f' & y$day=='2013-06-19'])/length(y$sex_inc_filled[y$day=='2013-06-19'])
								
								print(bc$nest[i])
								}
			p_f=do.call(rbind,l)
			p_ff=ddply(p_f,.(nest),summarise, mean_=mean(p_f), med_=median(p_f))
			bt$p_f=p_ff$med_[match(bt$nest,p_ff$nest)]
			bt$prop=ifelse(bt$sex=='f', bt$p_f, 1-bt$p_f)
			}
			{# analyze  proportion of inc of the treated bird
				ggplot(p_f,aes(x=nest, y=p_f, col=nest))+geom_boxplot()+geom_point(col='darkgrey')+stat_smooth()+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2))+stat_summary(fun.y=mean, colour="darkred", geom="point", shape=18, size=2,show_guide = FALSE)
				densityplot(~bt$prop)
				m=lm(inc_eff~prop,bt)
				  summary(m)
				  plot(allEffects(m))
					nsim <- 2000
					bsim <- sim(m, n.sim=nsim)  
					apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))
				
				ggplot(bt,aes(x=prop, y=inc_eff, col=sex))+geom_point()+stat_smooth()
				
				
				
			}
	}
	{# escape distance
		{# run first
		load(file=paste(wd,'escape.Rdata',sep="")) #vv - all estimates, es - aggregated per nest and sex
		load(file=paste(wd,'experimental.Rdata', sep=""))
			b=b[-which(b$nest%in%c('s410','s514','s524','s624')),]
			b$inc_eff[b$nest=='s806']=b$inc_eff_2[b$nest=='s806']
			b$bout_start=as.POSIXct(b$bout_start)
			b$exper=as.factor(b$exper)
			b$sex=as.factor(b$sex)
			bc=b[b$exper=='c',]
			bt=b[b$exper=='t',]
			bt$esc=es$esc[match(paste(bt$nest, bt$sex),paste(es$nest,es$sex))]
			bt=bt[-which(is.na(bt$esc)),]
		}
		{# analyze
			# distribution of escape
			densityplot(~bt$esc)
			densityplot(~bt$esc, groups=bt$sex, auto.key=TRUE)
			summary(bt$esc)
			densityplot(~bt$esc)
			densityplot(~log(bt$esc))
			densityplot(~log(bt$esc),groups=bt$sex, auto.key=TRUE)
			bt$esc_log=log(bt$esc)
			m=lm(inc_eff~esc,bt)#m=lm(inc_eff~esc,bt[bt$esc<40,])
			m2=lm(inc_eff~esc*sex,bt)
			m3=lm(inc_eff~log(esc),bt)
			m4=lm(inc_eff~log(esc)*sex,bt)
			AIC(m,m2,m3,m4)
				summary(m)
				  plot(allEffects(m))
				  plot(allEffects(m2))
				  plot(allEffects(m3))
				  plot(allEffects(m4))
					nsim <- 2000
					bsim <- sim(m, n.sim=nsim)  
					apply(bsim@coef, 2, quantile, prob=c(0.025,0.975))
				
				ggplot(bt,aes(x=esc, y=inc_eff, col=sex))+geom_point()
				ggplot(bt,aes(x=esc, y=inc_eff, col=sex))+geom_point()+stat_smooth(method='lm')
				ggplot(bt,aes(x=log(esc), y=inc_eff, col=sex))+geom_point()
				ggplot(bt,aes(x=log(esc), y=inc_eff, col=sex))+geom_point()+stat_smooth()
				ggplot(bt,aes(x=log(esc), y=inc_eff, col=sex))+geom_point()+stat_smooth(method='lm')
				
				+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2))+stat_summary(fun.y=mean, colour="darkred", geom="point", shape=18, size=3,show_guide = FALSE)
		
		}
	}
	{# season
		b$bout_start=as.POSIXct(b$bout_start)
		b$b_st_season_j = as.numeric(format(b$bout_start,"%j")) - as.numeric(format(as.POSIXct('2013-06-01'),"%j"))+1 
		b$exper=as.factor(b$exper)
		b$sex=as.factor(b$sex)		
		bt=b[b$exper=='t',]

		
		densityplot(~bt$b_st_season_j)
		densityplot(~log(bt$b_st_season_j-20))
		densityplot(~bt$b_st_season_j, groups=bt$sex, auto.key=TRUE)
		densityplot(~log(bt$b_st_season_j-20), groups=bt$sex, auto.key=TRUE)
		ggplot(bt,aes(x=b_st_season_j, y=inc_eff, col=sex))+geom_point()
		ggplot(bt,aes(x=b_st_season_j, y=inc_eff, col=sex))+geom_point()+stat_smooth(method='lm')
		ggplot(bt,aes(x=log(b_st_season_j-20), y=inc_eff, col=sex))+geom_point()
		ggplot(bt,aes(x=log(b_st_season_j-20), y=inc_eff, col=sex))+geom_point()+stat_smooth(method='lm')
		
			m=lm(inc_eff~b_st_season_j,bt)#m=lm(inc_eff~esc,bt[bt$esc<40,])
			m2=lm(inc_eff~b_st_season_j*sex,bt)
			m3=lm(inc_eff~log(b_st_season_j-20),bt)
			m4=lm(inc_eff~log(b_st_season_j-20)*sex,bt)
			AIC(m,m2,m3,m4)
				summary(m)
				plot(allEffects(m))
	}
	{# incubation period
		load(file=paste(wd,'experimental.Rdata', sep=""))
		b$inc_eff[b$nest=='s806']=b$inc_eff_2[b$nest=='s806']
		b$bout_start=as.POSIXct(b$bout_start)
		b$exper=as.factor(b$exper)
		
		load(file = paste(wd,"inc_start_estimation_2013_new.Rdata",sep=""))	
		e$nest=tolower(e$nest)
		b$inc_start=e$inc_start[match(b$nest,e$nest)]
		b$bout_start=as.POSIXct(b$bout_start)
		b$inc_start=as.POSIXct(b$inc_start)
		
			# start of bout within incubation
				b$bout_start_j  = as.numeric(format(b$bout_start ,"%j")) - as.numeric(format(b$inc_start,"%j")) +1
			# start of bout within season	
				b$b_st_season_j = as.numeric(format(b$bout_start,"%j")) - min(as.numeric(format(b$inc_start,"%j")))+1 
			# start of incubation within season
				b$inc_start_j = as.numeric(format(b$inc_start,"%j")) - min(as.numeric(format(b$inc_start,"%j")))+1
				
		b$exper=as.factor(b$exper)
		b$sex=as.factor(b$sex)		
		bt=b[b$exper=='t',]
		# start of bout within incubation
		densityplot(~bt$bout_start_j)
		densityplot(~bt$bout_start_j, groups=bt$sex, auto.key=TRUE)
		
		ggplot(bt,aes(x=bout_start_j, y=inc_eff, col=sex))+geom_point()
		ggplot(bt,aes(x=bout_start_j, y=inc_eff, col=sex))+geom_point()+stat_smooth(method='lm')
		
			m=lm(inc_eff~bout_start_j,bt)#m=lm(inc_eff~esc,bt[bt$esc<40,])
			m2=lm(inc_eff~poly(bout_start_j,2),bt)#m=lm(inc_eff~esc,bt[bt$esc<40,])
			m3=lm(inc_eff~bout_start_j*sex,bt)
			AIC(m,m2,m3)
				summary(m)
				plot(allEffects(m))
		# start of incubation within season		
		densityplot(~bt$inc_start_j)
		densityplot(~bt$inc_start_j, groups=bt$sex, auto.key=TRUE)
		
		ggplot(bt,aes(x=inc_start_j, y=inc_eff, col=sex))+geom_point()
		ggplot(bt,aes(x=inc_start_j, y=inc_eff, col=sex))+geom_point()+stat_smooth(method='lm')
		
			m=lm(inc_eff~inc_start_j,bt)#m=lm(inc_eff~esc,bt[bt$esc<40,])
			m2=lm(inc_eff~poly(inc_start_j,2),bt)#m=lm(inc_eff~esc,bt[bt$esc<40,])
			m3=lm(inc_eff~inc_start_j*sex,bt)
			AIC(m,m2,m3)
				summary(m)
				plot(allEffects(m))		
	}
	
	{# post treatment - only for nests where somebody returned
  #### #S515 NEED EXCEPTON AS BOUT 16,18 ARE TAKEN AS SEPARATE BOUTS, BUT THEY ARE ONE
  
		{# prepare metadata
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
			nn$bt=round(bt$inc_eff[match(nn$nest,bt$nest)],2)
			nn$bb=round(bb$inc_eff[match(nn$nest,bb$nest)],2)
			nn$bt=ifelse(nchar(nn$bt)==3, paste(nn$bt,"0",sep=""), nn$bt)
			nn$bb=ifelse(nchar(nn$bb)==3, paste(nn$bb,"0",sep=""), nn$bb)
			nn$label=paste(nn$bb,nn$bt,"              ",sep="       ")
			
			nn$sex=bt$sex[match(nn$nest,bt$nest)]
			nn$col_=ifelse(nn$sex=="f", "#FCB42C","#535F7C")
			nnn=nn$label
			names(nnn)=nn$n
			nnnest=nn$nest
			names(nnnest)=nn$n
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
												if(tst=='f'){x=data.frame(nest=x$nest[1],start_=min(x$bout_start[!is.na(x$exper) & x$exper=='t']), end_=min(x$bout_start[x$bout_start>tst2 & x$sex=='m']),stringsAsFactors=FALSE)}
												if(tst=='m'){x=data.frame(nest=x$nest[1], start_=min(x$bout_start[!is.na(x$exper) & x$exper=='t']),end_=min(x$bout_start[x$bout_start>tst2 & x$sex=='f']),stringsAsFactors=FALSE)}
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
                       n = e$n, sex=e$sex, nest=e$nest)
				
				o=read.csv(paste(wd,'experiment_metadata.csv', sep=""), stringsAsFactors=FALSE)	
				e$taken=as.POSIXct(o$taken_time[match(e$nest,o$nest)])
				write.table(e,file=paste(wd,'experiment_times.csv', sep=""), sep=",",col.names = TRUE, row.names=FALSE)	
		}
		}
		{# prepare data	
			load(file=paste(wd,'bout.Rdata', sep=""))
				#ggplot(b[!b$nest=='s806' & b$bout_type=='inc',], aes(x=inc_eff_3, y=inc_eff))+geom_point()+stat_smooth()
				#ggplot(b[!b$nest=='s806' & b$bout_type=='inc',], aes(x=inc_eff_3, y=inc_eff))+geom_point()+stat_smooth()+coord_cartesian(ylim = c(0.8,1), xlim = c(0.8,1))
			b$inc_eff[b$nest=='s806']=b$inc_eff_2[b$nest=='s806']
			b$inc_eff[b$nest=='s515' & b$bout_ID==16 | b$nest=='s702' & b$bout_ID==12]=b$inc_eff_3[b$nest=='s515' & b$bout_ID==16 | b$nest=='s702' & b$bout_ID==12] # bird came just for a second - so to have bout with some constancy, we assign constancy based on tag_ID
			d=b[-which(b$nest%in%c('s410','s514','s524','s624')),]
			d$bout_start=as.POSIXct(d$bout_start)
			d$bout_end=as.POSIXct(d$bout_end)
			# limit to nests where a bird returned
			d=d[-which(d$nest%in%c('s402','s409','s510','s516','s526','s711','s806')),] 
			e=e[-which(e$nest%in%c('s402','s409','s510','s516',"s526",'s711','s806')),] 
			# create zero exchange gaps where no exchange gap detected
			d=d[order(d$nest, d$bout_start),]
			d=ddply(d,.(nest), transform, before=as.character((c(NA,bout_type[-length(bout_type)]))))
			d$pk=c(1:nrow(d))
			l=list()
			for(i in 1:nrow(d)){
								di=d[i,]
								if(di$bout_type=='inc' & !is.na(di$before) & di$before=='inc'){
										di$pk=di$pk-0.5
										di$bout_ID=di$bout_ID-0.5
										di$bout_type='gap'
										di$bout_end=di$bout_start
										di$bout_length=di$inc_eff=di$inc_eff_2=di$inc_eff_3=di$disturb=di$disturb_log=di$dist_visit=di$dist_v_log=di$capture=0
										l[[i]]=di
										}
								}		
				d2=do.call(rbind, l)
			d=rbind(d,d2)
			d$before=NULL
			d=d[order(d$nest, d$bout_ID, d$bout_start),]
				# test correct assignment - YES
					#d_=ddply(d,.(nest), transform, before=as.character((c(NA,bout_type[-length(bout_type)]))))
					#d_=d_[which(!is.na(d_$before)),]
					#d_[d_$before==d_$bout_type,]
			
			li=list()
			lg=list()
			for(i in 1:nrow(e)){
							print(e$nest[i])
							ei=e[e$nest==e$nest[i],]
							di=d[d$nest==e$nest[i],]
							di=di[1:(nrow(di)-1),]
							{# treated bird 
								# before 
									db=di[di$sex==ei$sex & di$bout_end<(ei$taken-30*60),]
									db$pk=c(1:nrow(db))
									
								  # incubation bout
									dbt=db[db$bout_type=='inc' ,] # limit to bouts ending minimum 30 minutes before treatment capture
									dbt=dbt[(nrow(dbt)-2):nrow(dbt),] # take 3 bouts prior to treatment
									dbt$exper='b'
									dbt$treated_bird='y'
								  # gaps
											btg=db[db$bout_type=='gap' ,]
											btg=btg[(nrow(btg)-2):nrow(btg),] 
											btg$exper='b'
											btg$treated_bird='y'
								# after # if long gap and then still treated incubates, it is considered as his new bout
									da=di[di$sex==ei$sex & di$bout_start>(ei$release),]									
									da$pk=c(1:nrow(da))
								 # incubation bout
									dat=da[da$bout_type=='inc' ,] 
									dat=dat[1:3,] # take 3 bouts after to treatment
									dat$exper='a'
									dat$treated_bird='y'
								  # gaps
										atg=da[da$bout_type=='gap' ,] 
										atg=atg[1:3,]
										atg$exper='a'
										atg$treated_bird='y'
								inc_t=rbind(dbt, dat)
								gap_t=rbind(btg, atg)
								}
							{# control bird 
								# before 
									db=di[di$sex!=ei$sex & di$bout_end<(ei$taken-30*60),]
									db$pk=c(1:nrow(db))
								  # incubation bout
									dbt=db[db$bout_type=='inc' ,] # limit to bouts ending minimum 30 minutes before treatment capture
									dbt=dbt[(nrow(dbt)-2):nrow(dbt),] # take 3 bouts prior to treatment
									dbt$exper='b'
									dbt$treated_bird='n'
								 # gaps
											btg=db[db$bout_type=='gap' ,]
											btg=btg[(nrow(btg)-2):nrow(btg),] 
											btg$exper='b'
											btg$treated_bird='n'
								# after
									da=di[di$sex!=ei$sex & di$bout_start>(ei$release),]
									da$pk=c(1:nrow(da))
									
								  # incubation bout
									dat=da[da$bout_type=='inc' ,] # limit to bouts ending minimum 30 minutes before treatment capture
									dat=dat[1:3,] # take 3 bouts after to treatment
									dat$exper='a'
									dat$treated_bird='n'
								  # gaps
										atg=da[da$bout_type=='gap' ,] 
										atg=atg[1:3,]
										atg$exper='a'
										atg$treated_bird='n'
										
								inc_c=rbind(dbt, dat)
								gap_c=rbind(btg, atg)
								}
							li[[i]]=rbind(inc_t,inc_c)
							lg[[i]]=rbind(gap_t,gap_c)
							}
			inc=do.call(rbind,li)
			inc=inc[which(!is.na(inc$bout_length)),]

			gap=do.call(rbind,lg)
			gap=gap[which(!is.na(gap$bout_length)),]
			
			load(file = paste(wd,"inc_start_estimation_2013_new.Rdata",sep=""))	
			e$nest=tolower(e$nest)
			
			inc$inc_start=e$inc_start[match(inc$nest,e$nest)]
			inc$bout_start=as.POSIXct(inc$bout_start)
			inc$inc_start=as.POSIXct(inc$inc_start)
			inc$bout_length=inc$bout_length/60
			inc$sex=as.factor(inc$sex)
			inc$exper=as.factor(inc$exper)
			inc$treated_bird=as.factor(inc$treated_bird)
			
			# start of bout within incubation
				inc$bout_start_j  = as.numeric(format(inc$bout_start ,"%j")) - as.numeric(format(inc$inc_start,"%j")) +1
			# start of bout within season	
				inc$b_st_season_j = as.numeric(format(inc$bout_start,"%j")) - min(as.numeric(format(inc$inc_start,"%j")))+1 
			# start of incubation within season
				inc$inc_start_j = as.numeric(format(inc$inc_start,"%j")) - min(as.numeric(format(inc$inc_start,"%j")))+1
				
			
			gap$inc_start=e$inc_start[match(gap$nest,e$nest)]
			gap$bout_start=as.POSIXct(gap$bout_start)
			gap$inc_start=as.POSIXct(gap$inc_start)
			gap$sex=as.factor(gap$sex)
			gap$sex=as.factor(gap$sex)
			gap$exper=as.factor(gap$exper)
			gap$treated_bird=as.factor(gap$treated_bird)
			
			# start of bout within oncpubation
				gap$bout_start_j  = as.numeric(format(gap$bout_start ,"%j")) - as.numeric(format(gap$inc_start,"%j")) +1
			# start of bout within season	
				gap$b_st_season_j = as.numeric(format(gap$bout_start,"%j")) - min(as.numeric(format(gap$inc_start,"%j")))+1 
			# start of gapubation within season
				gap$inc_start_j = as.numeric(format(gap$inc_start,"%j")) - min(as.numeric(format(gap$inc_start,"%j")))+1
		}
		
		{# incubation
			{# boxplot per nest
			 # incubation
				ggplot(inc, aes(x=interaction(treated_bird, sex, nest), y=bout_length, col=exper))+geom_point()+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2))
				
				ggplot(inc, aes(x=interaction(sex, nest), y=bout_length, col=exper, fill=treated_bird))+geom_boxplot()+ stat_summary(data=inc, aes(x=interaction(sex, nest), y=bout_length, col=exper,fill=treated_bird),fun.y=mean, geom="point")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2))
				
				ggplot(inc, aes(x=nest, y=bout_length, col=exper, fill=treated_bird))+geom_boxplot()+ stat_summary(data=inc, aes(x=nest, y=bout_length, col=exper,fill=treated_bird),fun.y=mean, geom="point")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2))+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_colour_manual(values = c("black", "grey"))
				
				ggplot(inc, aes(x=interaction(sex, nest), y=bout_length, col=exper, fill=treated_bird))+geom_boxplot()+ stat_summary(data=inc, aes(x=interaction(sex, nest), y=bout_length, col=exper,fill=treated_bird),fun.y=mean, geom="point")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2))+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_colour_manual(values = c("black", "grey"))
				
				ggplot(inc, aes(x=interaction(treated_bird, sex, nest), y=bout_length, col=exper))+geom_boxplot()+ stat_summary(data=inc, aes(x=interaction(treated_bird, sex, nest), y=bout_length, col=exper),fun.y=mean, geom="point")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2))
				
			# inc_eff
				ggplot(inc, aes(x=interaction(treated_bird, sex, nest), y=inc_eff, col=exper))+geom_point()+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2))
				ggplot(inc, aes(x=interaction(treated_bird, sex, nest), y=inc_eff, col=exper))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2))
				ggplot(inc, aes(x=interaction(treated_bird, sex, nest), y=inc_eff, col=exper))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2))+stat_summary(data=inc, aes(x=interaction(treated_bird, sex, nest), y=inc_eff, col=exper),fun.y=mean, geom="point")
				
				
				ggplot(inc, aes(x=nest, y=inc_eff, col=exper, fill=treated_bird))+geom_boxplot()+ stat_summary(data=inc, aes(x=nest, y=inc_eff, col=exper,fill=treated_bird),fun.y=mean, geom="point")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2))+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_colour_manual(values = c("black", "grey"))
				
				ggplot(inc, aes(x=nest, y=inc_eff, col=treated_bird, fill=exper))+geom_boxplot()+ stat_summary(data=inc, aes(x=nest, y=inc_eff, col=treated_bird, fill=exper),fun.y=mean, geom="point")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2))+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_colour_manual(values = c("black", "grey"))
					
			}
			{# boxplot per sex
				# incubation
				ggplot(inc, aes(x=sex, y=bout_length, col=treated_bird, fill=exper))+geom_boxplot()+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_colour_manual(values = c("black", "grey"))
				ggplot(inc, aes(x=treated_bird, y=bout_length, col=exper, fill=sex))+geom_boxplot()+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_colour_manual(values = c("black", "grey"))
				
				# inc_eff
				ggplot(inc, aes(x=sex, y=inc_eff, col=treated_bird, fill=exper))+geom_boxplot()+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_colour_manual(values = c("black", "grey"))
				ggplot(inc, aes(x=sex, y=inc_eff, col=treated_bird, fill=exper))+geom_boxplot()+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_colour_manual(values = c("black", "grey"))+geom_boxplot()+coord_cartesian(ylim = c(0.8,1))
				ggplot(inc, aes(x=treated_bird, y=inc_eff, col=exper, fill=sex))+geom_boxplot()+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_colour_manual(values = c("black", "grey"))
				
				ggplot(inc, aes(x=interaction(treated_bird, sex), y=inc_eff, col=exper))+geom_boxplot()
				ggplot(inc, aes(x=interaction(treated_bird, sex), y=inc_eff, col=exper))+geom_boxplot()+coord_cartesian(ylim = c(0.8,1))
				#ggplot(inc[inc$treated_bird=='y',], aes(x=interaction(exper,sex), y=bout_length, col=interaction(exper,sex)))+geom_boxplot()
				#ggplot(inc[inc$treated_bird=='y',], aes(x=interaction(exper,sex), y=inc_eff, col=interaction(exper,sex)))+geom_boxplot()
				#ggplot(inc[inc$treated_bird=='y',], aes(x=interaction(exper,sex), y=inc_eff, col=interaction(exper,sex)))+geom_boxplot()+coord_cartesian(ylim = c(0.8,1))
			
			}
			{# all inc eff
					{# plot
					x=ddply(inc,.(nest,sex,exper, treated_bird), summarise, med_bout=median(bout_length),mean_bout=mean(bout_length), med_eff=median(inc_eff),mean_eff=mean(inc_eff))
								ggplot(x,aes(x=med_eff,y=mean_eff))+geom_point()
								plot(med_eff~mean_eff,x, xlim=c(0,1), ylim=c(0,1))
								abline(a=0,b=1)
								
								x[x$med_eff<0.2,]
								densityplot(~x$med_eff)
					
					xb=x[x$exper=='b',]
					xa=x[x$exper=='a',]
					xb$med_eff_after=xa$med_eff[match(paste(xb$nest,xb$sex),paste(xa$nest,xa$sex))]
					xb$col_=ifelse(xb$treated_bird=="y", "#FCB42C","#535F7C")
					dev.new(width=2.5,height=2.5)
					#png(paste(out2,"constancy_treated_control_relationship.png", sep=""), width=2.5,height=2.5,units="in",res=600)
					par(mar=c(2,2,0.7,0.7), ps=12, mgp=c(1.2,0.35,0), las=1, cex.lab=0.7, cex.axis=0.6, tcl=-0.15,bty="l",xpd=TRUE, col.axis="grey30", col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt

					plot(xb$med_eff_after[xb$nest!='s806']~xb$med_eff[xb$nest!='s806'], col=xb$col_[xb$nest!='s806'],bg=adjustcolor( xb$col_[xb$nest!='s806'], alpha.f = 0.4), pch=21,xlim=c(0,1), ylim=c(0,1), ylab='Incubation constancy - after', xlab=NA, xaxt='n')
					lines(c(0,1),c(0,1), lty=3, col="lightgrey")
					axis(1, at=seq(0,1,by=0.2), labels=FALSE)
					text(seq(0,1,by=0.2), par("usr")[3]-0.05, labels = seq(0,1,by=0.2),  xpd = TRUE, cex=0.6, col="grey30") #labels = z_g$genus
					mtext('Incubation constancy - before',side=1,line=0.6, cex=0.7, las=1, col='grey30')
					
					legend("bottomleft", legend=c('treated_bird','control_bird'), pch=19, col=c( '#FCB42C','#535F7C'), cex=0.8, bty='n')
					#dev.off()
					}
					{# model
						m=lmer(inc_eff~scale(bout_length)+exper+(1|nest),inc)
						m_=lmer(log(inc_eff)~scale(bout_length)+exper+(1|nest),inc)
						m1=lmer(inc_eff~scale(bout_length)+exper*treated_bird+(1|nest),inc)
						m2=lmer(inc_eff~scale(bout_length)+exper*sex+(1|nest),inc)
						m3=lmer(inc_eff~scale(bout_length)+exper*treated_bird*sex+(1|nest),inc)
							aic=AIC(m,m_,m1,m2,m3)
							aic[order(aic$AIC),]
							
						plot(allEffects(m))
						summary(m)
						# simulation		
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
						apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
						fixef(m)
					}
			}
			{# limited inc eff
				x=ddply(inc,.(nest,sex,exper, treated_bird), summarise, med_bout=median(bout_length), med_eff=median(inc_eff))
				x=x[x$med_eff>0.5,]
				xb=x[x$exper=='b',]
				xa=x[x$exper=='a',]
				xb$med_eff_after=xa$med_eff[match(paste(xb$nest,xb$sex),paste(xa$nest,xa$sex))]
				xb$col_=ifelse(xb$treated_bird=="y", "#FCB42C","#535F7C")
			dev.new(width=2.5,height=2.5)
			#png(paste(out2,"constancy_treated_control_relationship.png", sep=""), width=2.5,height=2.5,units="in",res=600)
			par(mar=c(2,2,0.7,0.7), ps=12, mgp=c(1.2,0.35,0), las=1, cex.lab=0.7, cex.axis=0.6, tcl=-0.15,bty="l",xpd=TRUE, col.axis="grey30", col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt

			plot(xb$med_eff_after[xb$nest!='s806']~xb$med_eff[xb$nest!='s806'], col=xb$col_[xb$nest!='s806'],bg=adjustcolor( xb$col_[xb$nest!='s806'], alpha.f = 0.4), pch=21,xlim=c(0.8,1), ylim=c(0.7,1), ylab='Incubation constancy - after', xlab=NA, xaxt='n')
			lines(c(0.8,1),c(0.8,1), lty=3, col="lightgrey")
			axis(1, at=seq(0.8,1,by=0.05), labels=FALSE)
			text(seq(0.8,1,by=0.05), par("usr")[3]-0.015, labels = seq(0.8,1,by=0.05),  xpd = TRUE, cex=0.6, col="grey30") #labels = z_g$genus
			mtext('Incubation constancy - before',side=1,line=0.6, cex=0.7, las=1, col='grey30')
			
			legend("bottomleft", legend=c('treated_bird','control_bird'), pch=19, col=c( '#FCB42C','#535F7C'), cex=0.8, bty='n')
			#dev.off()
			}
			
			{# bout
				{# all bout b vs a
					x=ddply(inc,.(nest,sex,exper, treated_bird), summarise, med_bout=median(bout_length),mean_bout=mean(bout_length), med_eff=median(inc_eff),mean_eff=mean(inc_eff))
					
						plot(med_bout~mean_bout,x, xlim=c(0,15), ylim=c(0,15))
						abline(a=0, b=1)
							
					xb=x[x$exper=='b',]
					xa=x[x$exper=='a',]
					xb$med_bout_after=xa$med_bout[match(paste(xb$nest,xb$sex),paste(xa$nest,xa$sex))]
					xb$col_=ifelse(xb$treated_bird=="y", "#FCB42C","#535F7C")
					dev.new(width=2.5,height=2.5)
					#png(paste(out2,"constancy_treated_control_relationship.png", sep=""), width=2.5,height=2.5,units="in",res=600)
					par(mar=c(2,2,0.7,0.7), ps=12, mgp=c(1.2,0.35,0), las=1, cex.lab=0.7, cex.axis=0.6, tcl=-0.15,bty="l",xpd=TRUE, col.axis="grey30", col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt

					plot(xb$med_bout_after~xb$med_bout, col=xb$col_,bg=adjustcolor( xb$col_, alpha.f = 0.4), pch=21,xlim=c(0,15), ylim=c(0,15), ylab='Bout - after [h]', xlab=NA, xaxt='n')
					lines(c(0,15),c(0,15), lty=3, col="lightgrey")
					axis(1, at=seq(0,15,by=5), labels=FALSE)
					text(seq(0,15,by=5), par("usr")[3]-0.7, labels = seq(0,15,by=5),  xpd = TRUE, cex=0.6, col="grey30") #labels = z_g$genus
					mtext('bout - before [h]',side=1,line=0.6, cex=0.7, las=1, col='grey30')
					
					legend("bottomleft", legend=c('treated_bird','control_bird'), pch=19, col=c( '#FCB42C','#535F7C'), cex=0.6, bty='n')
					
					ggplot(x[x$treated_bird=='y',], aes(x=interaction(exper,sex), y=med_bout, col=interaction(exper,sex)))+geom_boxplot()
					ggplot(inc[inc$treated_bird=='y',], aes(x=interaction(exper,sex), y=bout_length, col=interaction(exper,sex)))+geom_boxplot()
					#dev.off()
				}
				{# all bout f vs m
					x=ddply(inc,.(nest,sex,exper, treated_bird), summarise, med_bout=median(bout_length), med_eff=median(inc_eff))
					xf=x[x$sex=='f',]
					xm=x[x$sex=='m',]
					xf$med_bout_male=xm$med_bout[match(paste(xf$nest,xf$exper),paste(xm$nest,xm$exper))]
					xf$col_=ifelse(xf$exper=="b", "#FCB42C","#535F7C")
					
					dev.new(width=2.5,height=2.5)
					#png(paste(out2,"constancy_treated_control_relationship.png", sep=""), width=2.5,height=2.5,units="in",res=600)
					par(mar=c(2,2,0.7,0.7), ps=12, mgp=c(1.2,0.35,0), las=1, cex.lab=0.7, cex.axis=0.6, tcl=-0.15,bty="l",xpd=TRUE, col.axis="grey30", col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt

					plot(xf$med_bout_male~xf$med_bout, col=xf$col_,bg=adjustcolor( xf$col_, alpha.f = 0.4), pch=21,xlim=c(0,15), ylim=c(0,15), ylab='Bout - female [h]', xlab=NA, xaxt='n')
					lines(c(0,15),c(0,15), lty=3, col="lightgrey")
					axis(1, at=seq(0,15,by=5), labels=FALSE)
					text(seq(0,15,by=5), par("usr")[3]-0.7, labels = seq(0,15,by=5),  xpd = TRUE, cex=0.6, col="grey30") #labels = z_g$genus
					mtext('bout - male [h]',side=1,line=0.6, cex=0.7, las=1, col='grey30')
					
					legend("bottomleft", legend=c('before','after'), pch=19, col=c( '#FCB42C','#535F7C'), cex=0.6, bty='n')
					
					ggplot(x[x$treated_bird=='y',], aes(x=interaction(exper,sex), y=med_bout, col=interaction(exper,sex)))+geom_boxplot()
					ggplot(x[x$treated_bird=='n',], aes(x=interaction(exper,sex), y=med_bout, col=interaction(exper,sex)))+geom_boxplot()
					ggplot(inc[inc$treated_bird=='y',], aes(x=interaction(exper,sex), y=bout_length, col=interaction(exper,sex)))+geom_boxplot()
					ggplot(inc[inc$treated_bird=='n',], aes(x=interaction(exper,sex), y=bout_length, col=interaction(exper,sex)))+geom_boxplot()
					#dev.off()
				
					ggplot(inc, aes(x=interaction(exper,sex, treated_bird), y=bout_length, col=interaction(exper,sex,treated_bird)))+geom_boxplot()
				}
				{# change in bout length in captive and treated over time 0.5,1.5,2.5 reflect captive
					inc_a=inc[inc$exper=='a',]
					inc_a=inc_a[order(inc_a$nest, inc_a$treated_bird, inc_a$bout_start),]
					inc_a=ddply(inc_a,.(nest, treated_bird), transform, bout_order=c(1:length(bout_start)))
					head(inc_a)
					inc_a$bout_order[inc_a$treated_bird=='n']=inc_a$bout_order[inc_a$treated_bird=='n']-0.5
					inc_a_n=inc_a[inc_a$treated_bird=='n',]
					inc_a_y=inc_a[inc_a$treated_bird=='y',]
					inc_a_y$bout_captive=inc_a_n$bout_length[match(paste(inc_a_y$nest,inc_a_y$sex),paste(inc_a_n$nest,inc_a_n$sex))]
					
					ggplot(inc_a, aes(x=bout_order,y=bout_length, col=nest))+geom_point()+geom_line(aes(group = nest))
				
					gap_a=gap[gap$exper=='a',]
					gap_a=gap_a[order(gap_a$nest, gap_a$treated_bird, gap_a$bout_start),]
					gap_a=ddply(gap_a,.(nest, treated_bird), transform, bout_order=c(1:length(bout_start)))
					head(gap_a)
					gap_a$bout_order[gap_a$treated_bird=='n']=gap_a$bout_order[gap_a$treated_bird=='n']-0.5
					gap_a_n=gap_a[gap_a$treated_bird=='n',]
					gap_a_y=gap_a[gap_a$treated_bird=='y',]
					gap_a_y$bout_captive=gap_a_n$bout_length[match(paste(gap_a_y$nest,gap_a_y$sex),paste(gap_a_n$nest,gap_a_n$sex))]
					
					ggplot(gap_a, aes(x=bout_order,y=bout_length, col=nest))+geom_point()+geom_line(aes(group = nest))
					ggplot(gap_a, aes(x=bout_order,y=log(bout_length+0.01), col=nest))+geom_point()+geom_line(aes(group = nest))
				}
				
				
				{# model
						
						
						
						inc$int=as.factor(interaction(inc$exper, inc$treated,inc$sex))
						unique(inc$int)
		
						densityplot(~inc$bout_length)
						densityplot(~log(inc$bout_length))
						m=lmer(bout_length~exper+(1|nest),inc)
						m_=lmer(log(bout_length)~exper+(1|nest),inc)
						m1=lmer(bout_length~exper*treated_bird+(1|nest),inc)
						m2=lmer(bout_length~exper*sex+(1|nest),inc)
						m3=lmer(bout_length~exper*treated_bird*sex+(1|nest),inc)
						m4=lmer(bout_length~int+(1|nest),inc)
							aic=AIC(m,m_,m1,m2,m3)
							aic[order(aic$AIC),]
							aic=AIC(m,m1,m2,m3,m4)
							aic[order(aic$AIC),]
							
						plot(allEffects(m))
						summary(m)
						# simulation		
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
						apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
						fixef(m)
						
						plot(allEffects(m3))
						plot(allEffects(m4))
						# simulation		
						nsim <- 2000
						bsim <- sim(m3, n.sim=nsim)  
						apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
						fixef(m3)
						
						summary(glht(m4, linfct = mcp(int = c(	"b.y.m - a.y.m = 0",
																"b.y.f - a.y.f = 0", #
																"b.n.m - a.n.m = 0", #
																"b.n.f - a.n.f = 0",
																
																"a.n.m - a.n.f = 0",
																"a.y.m - a.y.f = 0",#
																
																"a.y.m - a.n.m = 0",#
																"a.y.f - a.n.f = 0"
															))))
						
				}
				{# inc$bout_start_j 
				x=ddply(inc,.(nest,sex,exper, treated_bird), summarise, med_bout=median(bout_length), med_eff=median(inc_eff),med_bout_start_j=median(bout_start_j))
				
				ggplot(inc, aes(x=bout_start_j, y=bout_length, col=interaction(exper,treated_bird)))+geom_point()+stat_smooth()
				ggplot(inc, aes(x=bout_start_j, y=bout_length, col=interaction(exper,treated_bird)))+geom_point()+stat_smooth(method='lm')
				ggplot(x, aes(x=med_bout_start_j, y=med_bout, col=interaction(exper,treated_bird)))+geom_point()+stat_smooth(method='lm')
				
				ggplot(inc, aes(x=bout_start_j, y=inc_eff, col=interaction(exper,treated_bird)))+geom_point()+stat_smooth(method='lm')
				ggplot(inc, aes(x=bout_start_j, y=inc_eff, col=interaction(exper,treated_bird)))+geom_point()+stat_smooth(method='lm',se = FALSE)
				ggplot(inc, aes(x=bout_start_j, y=inc_eff, col=interaction(exper,treated_bird)))+geom_point()+stat_smooth(method='lm')+coord_cartesian(ylim = c(0.6,1))
				
				
				ggplot(gap, aes(x=bout_start_j, y=bout_length, col=interaction(exper,treated_bird)))+geom_point()+stat_smooth()
				ggplot(gap, aes(x=bout_start_j, y=bout_length, col=interaction(exper,treated_bird)))+geom_point()+stat_smooth(method='lm')
				ggplot(gap, aes(x=bout_start_j, y=log(bout_length+0.01), col=interaction(exper,treated_bird)))+geom_point()+stat_smooth(method='lm')
				ggplot(x, aes(x=med_bout_start_j, y=med_bout, col=interaction(exper,treated_bird)))+geom_point()+stat_smooth(method='lm')
					
					m=lmer(bout_length~scale(bout_start_j)*exper+(1|nest), gap)
					m_=lmer(log(bout_length+0.01)~scale(bout_start_j)*exper+(1|nest), gap)
					m2=lmer(bout_length~scale(bout_start_j)*exper*treated_bird+(1|nest), gap)
					m2_=lmer(log(bout_length+0.01)~scale(bout_start_j)*exper*treated_bird+(1|nest), gap)
					aic=AIC(m,m_,m2,m2_)
							aic[order(aic$AIC),]
							
						plot(allEffects(m2_))
						summary(m)
						# simulation		
						nsim <- 2000
						bsim <- sim(m2_, n.sim=nsim)  
						apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
						fixef(m)
				
			}
		}	
		}
		{# gap
			gap$present=ifelse(gap$bout_length==0,0,1)
			gapy=gap[gap$present==1,]
			
			densityplot(~gap$bout_length)
			densityplot(~log(gap$bout_length+0.1))
			densityplot(~log(gap$bout_length+0.1), groups=gap$sex, auto.key=TRUE)
			densityplot(~log(gap$bout_length+0.1), groups=gap$exper, auto.key=TRUE)
			
			densityplot(~gapy$bout_length)
			densityplot(~log(gapy$bout_length+0.1))
			densityplot(~log(gapy$bout_length+0.1), groups=gapy$sex, auto.key=TRUE)
			densityplot(~log(gapy$bout_length+0.1), groups=gapy$exper, auto.key=TRUE)
			
			ggplot(gap, aes(x=interaction(sex, nest), y=bout_length, col=exper, fill=treated_bird))+geom_boxplot()+ stat_summary(data=inc, aes(x=interaction(sex, nest), y=bout_length, col=exper,fill=treated_bird),fun.y=mean, geom="point")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2))
				
			ggplot(gap, aes(x=nest, y=bout_length, col=exper, fill=treated_bird))+geom_boxplot()+ stat_summary(data=inc, aes(x=nest, y=bout_length, col=exper,fill=treated_bird),fun.y=mean, geom="point")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2))+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_colour_manual(values = c("black", "grey"))
			
			ggplot(gap, aes(x=interaction(sex, nest), y=bout_length, col=exper, fill=treated_bird))+geom_boxplot()+ stat_summary(data=gap, aes(x=interaction(sex, nest), y=bout_length, col=exper,fill=treated_bird),fun.y=mean, geom="point")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2))+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_colour_manual(values = c("black", "grey"))			
			ggplot(gap, aes(x=interaction(sex, nest), y=log(bout_length+0.1), col=exper, fill=treated_bird))+geom_boxplot()+ stat_summary(data=gap, aes(x=interaction(sex, nest), y=log(bout_length+0.1), col=exper,fill=treated_bird),fun.y=mean, geom="point")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.2))+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_colour_manual(values = c("black", "grey"))			
			
			ggplot(gap, aes(x=interaction(treated_bird, sex), y=bout_length, col=exper))+geom_boxplot()
			ggplot(gap, aes(x=interaction(treated_bird, sex), y=log(bout_length+0.1), col=exper))+geom_boxplot()
			ggplot(gapy, aes(x=interaction(treated_bird, sex), y=log(bout_length+0.1), col=exper))+geom_boxplot()
			ggplot(gap[gap$treated_bird=='y',], aes(x=sex, y=log(bout_length+0.1), col=exper))+geom_boxplot()
			
			x=ddply(gap,.(nest,sex,exper, treated_bird), summarise, med_gap=median(bout_length), num_gaps=sum(present))
				{# gap plot	
					xb=x[x$exper=='b',]
					xa=x[x$exper=='a',]
					xb$med_gap_after=xa$med_gap[match(paste(xb$nest,xb$sex),paste(xa$nest,xa$sex))]
					xb$col_=ifelse(xb$treated_bird=="y", "#FCB42C","#535F7C")
					dev.new(width=2.5,height=2.5)
					#png(paste(out2,"constancy_treated_control_relationship.png", sep=""), width=2.5,height=2.5,units="in",res=600)
					par(mar=c(2,2,0.7,0.7), ps=12, mgp=c(1.2,0.35,0), las=1, cex.lab=0.7, cex.axis=0.6, tcl=-0.15,bty="l",xpd=TRUE, col.axis="grey30", col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt

					plot(log(xb$med_gap_after+0.1)~log(xb$med_gap+0.1), col=xb$col_,bg=adjustcolor( xb$col_, alpha.f = 0.4), pch=21, ylab='Exchange gap - after [min]', xlab=NA, yaxt='n',xaxt='n',xlim=log(c(0,260)+0.1), ylim=log(c(0,260)+0.1))
					lines(log(c(0,260)+0.1),log(c(0,260)+0.1), lty=3, col="lightgrey")
					axis(2, at=log(c(0,1,10,100,260)+0.1), labels=c(0,1,10,100,260))
					axis(1, at=log(c(0,1,10,100,260)+0.1), labels=FALSE)
					text(log(c(0,1,10,100,260)+0.1), par("usr")[3]-0.4, labels = c(0,1,10,100,260),  xpd = TRUE, cex=0.6, col="grey30") #labels = z_g$genus
					mtext('Exchange gap - before [min]',side=1,line=0.6, cex=0.7, las=1, col='grey30')
					
					legend("topright", legend=c('treated bird','control bird'), pch=19, col=c( '#FCB42C','#535F7C'), cex=0.6, bty='n')
				}
				{# gap presence plot - makes little sense
					densityplot(~x$num_gaps)
					xb=x[x$exper=='b',]
					xa=x[x$exper=='a',]
					xb$num_gaps_after=xa$num_gaps[match(paste(xb$nest,xb$sex),paste(xa$nest,xa$sex))]
					xb$col_=ifelse(xb$treated_bird=="y", "#FCB42C","#535F7C")
					dev.new(width=2.5,height=2.5)
					#png(paste(out2,"constancy_treated_control_relationship.png", sep=""), width=2.5,height=2.5,units="in",res=600)
					par(mar=c(2,2,0.7,0.7), ps=12, mgp=c(1.2,0.35,0), las=1, cex.lab=0.7, cex.axis=0.6, tcl=-0.15,bty="l",xpd=TRUE, col.axis="grey30", col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt

					plot(xb$num_gaps_after~xb$num_gaps, col=xb$col_,bg=adjustcolor( xb$col_, alpha.f = 0.4), pch=21, ylab='Exchange gap - after [min]', xlab=NA, xaxt='n',xlim=c(0,3), ylim=c(0,3))
				
					axis(1, at=seq(0,3,by=1), labels=FALSE)
					text(seq(0,3,by=1), par("usr")[3]-0.05, labels = seq(0,3,by=1),  xpd = TRUE, cex=0.6, col="grey30") #labels = z_g$genus
					mtext('bout - before [h]',side=1,line=0.6, cex=0.7, las=1, col='grey30')
					
					legend("bottomleft", legend=c('treated_bird','control_bird'), pch=19, col=c( '#FCB42C','#535F7C'), cex=0.6, bty='n')
				}
			{# model
				gap$int=as.factor(interaction(gap$exper, gap$treated,gap$sex))
						unique(gap$int)
		
						densityplot(~gap$bout_length)
						densityplot(~log(gap$bout_length))
						
						m=glmer(present~exper+(1|nest),gap, family='binomial')
						m1=glmer(present~exper*treated_bird+(1|nest),gap, family='binomial')
						m2=glmer(present~exper*sex+(1|nest),gap, family='binomial')
						m3=glmer(present~int+(1|nest),gap, family='binomial')
							aic=AIC(m,m1,m2,m3)
							aic[order(aic$AIC),]
								plot(allEffects(m))
								plot(allEffects(m3))
						plot(allEffects(m))
						summary(m)
						summary(m3)
						# simulation		
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
						apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
						fixef(m)
						
						m=lmer(bout_length~exper+(1|nest),gap)
						m_=lmer(log(bout_length+0.01)~exper+(1|nest),gap)
						m1=lmer(bout_length~exper*treated_bird+(1|nest),gap)
						m2=lmer(bout_length~exper*sex+(1|nest),gap)
						m3=lmer(bout_length~exper*treated_bird*sex+(1|nest),gap)
						m4=lmer(bout_length~int+(1|nest),gap)
						m5=lmer(log(bout_length+0.01)~int+(1|nest),gap)
							aic=AIC(m,m_,m1,m2,m3,m4,m5)
							aic[order(aic$AIC),]
							aic=AIC(m,m1,m2,m3,m4)
							aic[order(aic$AIC),]
							
						plot(allEffects(m))
						plot(allEffects(m_))
						summary(m)
						# simulation		
						nsim <- 2000
						bsim <- sim(m, n.sim=nsim)  
						apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
						fixef(m)
						
						plot(allEffects(m3))
						plot(allEffects(m4))
						plot(allEffects(m5))
						# simulation		
						nsim <- 2000
						bsim <- sim(m3, n.sim=nsim)  
						apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
						fixef(m3)
						
						summary(glht(m5, linfct = mcp(int = c(	"b.y.m - a.y.m = 0",
																"b.y.f - a.y.f = 0", #
																"b.n.m - a.n.m = 0", 
																"b.n.f - a.n.f = 0",
																
																"a.n.m - a.n.f = 0",
																"a.y.m - a.y.f = 0",
																
																"a.y.m - a.n.m = 0",
																"a.y.f - a.n.f = 0"#
															))))
						
			
			}
		}


	}
}	


}
{# sample size
	# sample size (17 females, 12 males taken)
	d=read.csv(paste(wd,'experiment_metadata.csv', sep=""), stringsAsFactors=FALSE)
	d$comments=NULL
	d=d[d$type=="exp",]
	head(d)
	summary(factor(d$taken[d$use_t=="y"])) # number of nests with f & m taken into captivity and useful for treatment
	summary(factor(d$taken[d$use_b=="y"]))# number of nests with f & m taken into captivity and useful for return
	summary(factor(d$back[d$use_b=="y"]))# number of nests where captive bird does not return upon release
			summary(factor(d$taken[d$use_b=="y" & d$back=="n"])) # sex of birds that never returned - 5 females
	
	d[d$use_b=="y" & d$back=="n",] # return of the captive bird does not seem to dependent on treated birds presence during catching
	summary(factor(d$present))
}	
{# number of mealworms
	
	# median mass of our mealworms
		m=wet(con, "select mass_worms mass, worms_num as num from barrow_2013.captivity where worms_num is not null and mass_worms is not null")
		m$mass[m$num==106]=8.9
		m$g=100*m$mass/m$num # for 100 worms
		ggplot(m[m$num>99,],aes(x=g, fill=factor(num)))+geom_density()
		ggplot(m[m$num>99,],aes(y=mass, x=factor(num), fill=factor(num)))+geom_boxplot()
		summary(m) #7.5g/100
	
	# energy metabolizable from 100 worms # 188kj 
		7.5*24.2 
	# SESA energetic requirements
	# based on Norton 1973 page 64 for  in 
				# given in cubic centimeter per gram and hour
				# cc = 0.001 liters (that is why we devide by 1000 to get it in leter
				# O2 in liters multiply by 20.10 to get KJ
		#27g SESA and 6.2C - median tundra temperature in Barrow
		27*((10.69+(-0.203)*6.2)/1000)*20.1 # kj/27g/h # 4.699782	
		24*27*((10.69+(-0.203)*6.2)/1000)*20.1 # kj/27g/d # 4.699782	
		
	# based on Ashkenazei and Safriel 1979
		# 19-59 duing incubation
		
	m=wet(con, "select*from barrow_2013.captivity")
	m$datetime_=as.POSIXct(m$datetime_)
	m=m[order(m$captivity_pk,m$datetime_),]	
	write.csv(m,paste(wd,'captivity.csv', sep=""),row.names=FALSE)
}	
{# nests with cages
    # nest without cage 
		length(c('s514','s515','s516','s519','s520', 's527','s624','s625','s628','s629','s712','s806','s808'))
	# nests with cage on
	d=dbGetQuery(conMy, "select nest, datetime_, cage from barrow_2013.visits where cage is not null and nest in 		('s402','s404','s406','s409','s410','s502','s509','s510','s512','s514','s515','s516','s519','s520','s524','s526','s527','s622','s624','s625','s628','s629','s702','s704','s711','s712','s805','s806','s808')")
	length(unique(d$nest))		
	
	d[order(d$nest,d$datetime_),]						
	summary(factor(d$cage))
	on=d$nest[d$cage=="on"]
	off=d$nest[d$cage=="off"]
	on[!on%in%off]
	off[!off%in%on]
	d$datetime_=as.POSIXct(d$datetime_)
	dd = ddply(d,.(nest),summarise, on = min(datetime_[cage=='on']), off = max(datetime_[cage=='off']))
	dd$off[dd$nest=='s622'] = as.POSIXct('2013-07-09 12:00:00')
	write.csv(dd,paste(wd,'cages.csv', sep=""),row.names=FALSE)
	
	}
{# mass change
	g=dbGetQuery(conMy, "select*from barrow_2013.captivity")	
    gg=ddply(g[!is.na(g$mass),],.(ring_num),summarise,mass_delta=(mass[phase=='start']-mass[phase=='end']), mass_end=mass[phase=='end'])
	ggplot(gg,aes(x=mass_delta))+geom_density()
	ggplot(gg,aes(x=mass_delta, y=mass_end))+geom_point()+stat_smooth()
}
{# nest with video on
	length(c('s402','s406','s409','s410','s510','s520','s524','s526','s527','s624','s702','s805'))
}	

{# cage effects 2012
	load("M:\\Science\\Projects\\PHD\\STATS\\DATASETS\\C_SESA_inc_start.Rdata")
	d1 = d
	load("M:\\Science\\Projects\\PHD\\STATS\\DATASETS\\C_SESA_short_inc_start.Rdata")
	d = rbind(d1,d)
	d$nest= tolower(d$nest)	
	d$bout_length = d$bout_length/60 
		{# DONE identify when cage was on
		require('RSQLite')
		db=paste("M:/Science/Projects/PHD/FIELD_METHODS/Barrow/2012/DATA/Database/Barrow_2012.sqlite",sep="") # database name
		conLite = dbConnect(dbDriver("SQLite"),dbname = db)
		cc=dbGetQuery(conLite,paste("SELECT nest, datetime_, cage FROM visits WHERE cage is not null"))
		dbDisconnect(conLite)
		length(unique(cc$nest))		
		cc$datetime_ = as.POSIXct(cc$datetime_)
		cc = cc[order(cc$nest,cc$datetime_),]	
		cc = cc[substring(cc$nest,1,1)== 's',]
		on=cc$nest[cc$cage==1]
		off=cc$nest[cc$cage==0]
	
		on[!on%in%off]
		off[!off%in%on]
		
		dd = ddply(cc,.(nest),summarise, on = min(datetime_[cage==1]), off = max(datetime_[cage==0]))
		dd$off= as.character(dd$off)
		dd$off[is.na(dd$off)] = as.character(dd$on[is.na(dd$off)] + 30*24*60)
		dd$off = as.POSIXct(dd$off)
		#dd$on[dd$nest == 's101'] = as.POSIXct('2012-06-11 21:17:00')
		dd = dd[!is.na(dd$on),]
		write.csv(dd,paste(wd,'cages_2012.csv', sep=""),row.names=FALSE)
		}
		cc =  read.csv(paste(wd,'cages_2012.csv', sep=""), stringsAsFactors = FALSE)
			cc$on = as.POSIXct(cc$on, tz='America/Anchorage')
			cc$off = as.POSIXct(cc$off, tz='America/Anchorage')
			d$cage  = 'n'
			for(i in 1:nrow(d)){
				n = d$nest[i]
				ci = cc[cc$nest == n,]
				if(nrow(ci)!=0){
				d$cage[i] = ifelse(d$bout_start[i]>ci$on & d$bout_start[i]<ci$off, 'y','n')
				}
				print(i)
			}
			summary(factor(d$cage))
			table(d$nest, d$cage)
			ggplot(d,aes(x = cage, y = bout_length)) + geom_boxplot()
			ggplot(d,aes(x = nest, y = bout_length, col = cage)) + geom_boxplot() #+ geom_dotplot(aes(fill = cage),binaxis = 'y', stackdir = 'center',  position = position_dodge())
			ggplot(d,aes(x = cage, y = inc_eff)) + geom_boxplot() 
			ggplot(d,aes(x = nest, y = inc_eff, col = cage)) + geom_boxplot()
			
		     d$cage = as.factor(d$cage)
			 m=lmer(bout_length~cage+(1|nest)+(1|bird_ID_filled),d, REML=FALSE)
			 summary(m)
			 summary(glht(m))
			 plot(allEffects(m))
			 m1=lmer(inc_eff~cage+(1|nest)+(1|bird_ID_filled),d, weights=sqrt(bout_length),REML=FALSE)
			 plot(allEffects(m1))
			 summary(glht(m1))
			
			{# median bout per bird, nest and cage 
			dd = ddply(d,.(nest, bird_ID_filled, sex), summarise, cage_y = median(bout_length[cage=='y']), cage_n = median(bout_length[cage=='n']), n_y = length(bout_length[cage=='y']), n_n = length(bout_length[cage=='n']))
			t.test(dd$cage_y,dd$cage_n,paired=TRUE)
			dd = dd[!is.na(dd$cage_y) & !is.na(dd$cage_n),]
			dd1 = dd
			dd1$cage_n = dd1$n_n = NULL
			names(dd1)[names(dd1) == 'cage_y'] = 'median_bout'
			names(dd1)[names(dd1) == 'n_y'] = 'n_bouts'
			dd1$cage = 'y'
			
			dd2 = dd
			dd2$cage_y = dd2$n_y = NULL
			names(dd2)[names(dd2) == 'cage_n'] = 'median_bout'
			names(dd2)[names(dd2) == 'n_n'] = 'n_bouts'
			dd2$cage = 'n'
			dd = rbind(dd1,dd2)
			dd$cage = as.factor(dd$cage)
			#dd = ddply(d,.(nest, bird_ID_filled, cage), summarise, median_bout = median(bout_length), cn_y = length(bout_length[cage=='y']), n_n = length(bout_length))
			
			pd <- position_dodge(0.4)
			ggplot(dd,aes(x=cage, y = median_bout, col = factor(bird_ID_filled)))+geom_point(aes(size = n_bouts), position = pd) +geom_line(position = pd, aes(group = factor(bird_ID_filled), col = factor(bird_ID_filled))) +
			theme(legend.position = "none") + ylab('median bird bout length [h]') + xlab('cage on the nest')
			
			m=lmer(median_bout~cage+(1|nest),dd, REML=FALSE)
			 summary(glht(m))
			 plot(allEffects(m))
			 }
			{# median bout per  nest and cage 
			dd = ddply(d,.(nest), summarise, cage_y = median(bout_length[cage=='y']), cage_n = median(bout_length[cage=='n']), n_y = length(bout_length[cage=='y']), n_n = length(bout_length[cage=='n']))
			t.test(dd$cage_y,dd$cage_n,paired=TRUE)
			dd = dd[!is.na(dd$cage_y) & !is.na(dd$cage_n),]
			dd1 = dd
			dd1$cage_n = dd1$n_n = NULL
			names(dd1)[2] = 'median_bout'
			names(dd1)[3] = 'n_bouts'
			dd1$cage = 'y'
			
			dd2 = dd
			dd2$cage_y = dd2$n_y = NULL
			names(dd2)[2] = 'median_bout'
			names(dd2)[3] = 'n_bouts'
			dd2$cage = 'n'
			dd = rbind(dd1,dd2)
			dd$cage = as.factor(dd$cage)
			#dd = ddply(d,.(nest, bird_ID_filled, cage), summarise, median_bout = median(bout_length), cn_y = length(bout_length[cage=='y']), n_n = length(bout_length))
			
			pd <- position_dodge(0.4)
			ggplot(dd,aes(x=cage, y = median_bout, col = factor(nest)))+geom_point(aes(size = n_bouts), position = pd) +geom_line(position = pd, aes(group = factor(nest), col = factor(nest))) +
			theme(legend.position = "none") + ylab('median nest bout length [h]') + xlab('cage on the nest')
			
			m=lmer(median_bout~cage+(1|nest),dd, REML=FALSE)
			 summary(glht(m))
			 plot(allEffects(m))
			 }
			
			{# median nest attentiveness per bird, nest and cage 
			dd = ddply(d,.(nest, bird_ID_filled), summarise, cage_y = median(inc_eff[cage=='y']), cage_n = median(inc_eff[cage=='n']), n_y = length(inc_eff[cage=='y']), n_n = length(inc_eff[cage=='n']))
			t.test(dd$cage_y,dd$cage_n,paired=TRUE)
			dd = dd[!is.na(dd$cage_y) & !is.na(dd$cage_n),]
			dd1 = dd
			dd1$cage_n = dd1$n_n = NULL
			names(dd1)[3] = 'median_attentiveness'
			names(dd1)[4] = 'n_bouts'
			dd1$cage = 'y'
			
			dd2 = dd
			dd2$cage_y = dd2$n_y = NULL
			names(dd2)[3] = 'median_attentiveness'
			names(dd2)[4] = 'n_bouts'
			dd2$cage = 'n'
			dd = rbind(dd1,dd2)
			dd$cage = as.factor(dd$cage)
			#dd = ddply(d,.(nest, bird_ID_filled, cage), summarise, median_bout = median(bout_length), cn_y = length(bout_length[cage=='y']), n_n = length(bout_length))
			
			pd <- position_dodge(0.4)
			ggplot(dd,aes(x=cage, y = median_attentiveness, col = factor(bird_ID_filled)))+geom_point(aes(size = n_bouts), position = pd) +geom_line(position = pd, aes(group = factor(bird_ID_filled), col = factor(bird_ID_filled))) +
			theme(legend.position = "none") + ylab('median bird attentiveness') + xlab('cage on the nest')
			
			m=lmer(median_attentiveness~cage+(1|nest),dd, REML=FALSE)
			 summary(glht(m))
			 plot(allEffects(m))
			 }
			{# median nest attentiveness per nest and cage 
			dd = ddply(d,.(nest), summarise, cage_y = median(inc_eff[cage=='y']), cage_n = median(inc_eff[cage=='n']), n_y = length(inc_eff[cage=='y']), n_n = length(inc_eff[cage=='n']))
			t.test(dd$cage_y,dd$cage_n,paired=TRUE)
			dd = dd[!is.na(dd$cage_y) & !is.na(dd$cage_n),]
			dd1 = dd
			dd1$cage_n = dd1$n_n = NULL
			names(dd1)[2] = 'median_attentiveness'
			names(dd1)[3] = 'n_bouts'
			dd1$cage = 'y'
			
			dd2 = dd
			dd2$cage_y = dd2$n_y = NULL
			names(dd2)[2] = 'median_attentiveness'
			names(dd2)[3] = 'n_bouts'
			dd2$cage = 'n'
			dd = rbind(dd1,dd2)
			dd$cage = as.factor(dd$cage)
			#dd = ddply(d,.(nest, bird_ID_filled, cage), summarise, median_bout = median(bout_length), cn_y = length(bout_length[cage=='y']), n_n = length(bout_length))
			
			pd <- position_dodge(0.4)
			ggplot(dd,aes(x=cage, y = median_attentiveness, col = factor(nest)))+geom_point(aes(size = n_bouts), position = pd) +geom_line(position = pd, aes(group = factor(nest), col = factor(nest))) +
			theme(legend.position = "none") + ylab('median nest attentiveness') + xlab('cage on the nest')
			
			m=lmer(median_attentiveness~cage+(1|nest),dd, REML=FALSE)
			 summary(glht(m))
			 plot(allEffects(m))
			 }

}
{# camera effects 2012
	load("M:\\Science\\Projects\\PHD\\STATS\\DATASETS\\C_SESA_inc_start.Rdata")
	#d1 = d
	#load("M:\\Science\\Projects\\PHD\\STATS\\DATASETS\\C_SESA_short_inc_start.Rdata")
	#d = rbind(d1,d)
	d$nest= tolower(d$nest)	
	d$bout_length = d$bout_length/60 
		{# DONE identify when recorder was on
		require('RSQLite')
		db=paste("M:/Science/Projects/PHD/FIELD_METHODS/Barrow/2012/DATA/Database/Barrow_2012.sqlite",sep="") # database name
		conLite = dbConnect(dbDriver("SQLite"),dbname = db)
		cc=dbGetQuery(conLite,paste("SELECT nest, datetime_, recorder FROM visits WHERE recorder is not null"))
		dbDisconnect(conLite)
		length(unique(cc$nest))		
		cc$datetime_ = as.POSIXct(cc$datetime_)
		cc = cc[order(cc$nest,cc$datetime_),]	
		cc = cc[substring(cc$nest,1,1)== 's',]
		on=cc$nest[cc$recorder==1]
		off=cc$nest[cc$recorder==0]
	
		on[!on%in%off]
		off[!off%in%on]
		
		dd = ddply(cc,.(nest),summarise, on = min(datetime_[recorder%in%c(1,2)]), off = max(datetime_[recorder%in%c(2,0)]))
		dd$off= as.character(dd$off)
		dd = dd[!is.na(dd$off),]
		dd$off = as.POSIXct(dd$off)
		#dd$on[dd$nest == 's101'] = as.POSIXct('2012-06-11 21:17:00')

		write.csv(dd,paste(wd,'recorders_2012.csv', sep=""),row.names=FALSE)
		}
		cc =  read.csv(paste(wd,'recorders_2012.csv', sep=""), stringsAsFactors = FALSE)
			cc$on = as.POSIXct(cc$on, tz='America/Anchorage')
			cc$off = as.POSIXct(cc$off, tz='America/Anchorage')
			d$recorder  = 'n'
			for(i in 1:nrow(d)){
				n = d$nest[i]
				ci = cc[cc$nest == n,]
				if(nrow(ci)!=0){
				d$recorder[i] = ifelse(d$bout_start[i]>ci$on & d$bout_start[i]<ci$off, 'y','n')
				}
				print(i)
			}
			summary(factor(d$recorder))
			unique(d$nest[d$recorder == 'y'])
			length(unique(d$nest[d$recorder == 'n']))
			table(d$nest, d$recorder)
			ggplot(d,aes(x = recorder, y = bout_length)) + geom_boxplot()
			ggplot(d,aes(x = nest, y = bout_length, col = recorder)) + geom_boxplot() #+ geom_dotplot(aes(fill = recorder),binaxis = 'y', stackdir = 'center',  position = position_dodge())
			ggplot(d,aes(x = recorder, y = inc_eff)) + geom_boxplot() 
			ggplot(d,aes(x = nest, y = inc_eff, col = recorder)) + geom_boxplot()
			
		     d$recorder = as.factor(d$recorder)
			 m=lmer(bout_length~recorder+(1|nest)+(1|bird_ID_filled),d, REML=FALSE)
			 m=lmer(bout_length~recorder+(1|nest)+(1|bird_ID_filled),d[d$nest%in%unique(d$nest[d$recorder == 'y']),], REML=FALSE)
			 nrow(d[d$nest%in%unique(d$nest[d$recorder == 'y']),])
			 summary(d$recorder[d$nest%in%unique(d$nest[d$recorder == 'y'])])
			 summary(m)
			 summary(glht(m))
			 plot(allEffects(m))
			 m1=lmer(inc_eff~recorder+(1|nest)+(1|bird_ID_filled),d, weights=sqrt(bout_length),REML=FALSE)
			  m1=lmer(inc_eff~recorder+(1|nest)+(1|bird_ID_filled),d[d$nest%in%unique(d$nest[d$recorder == 'y']),],weights=sqrt(bout_length), REML=FALSE)
			 plot(allEffects(m1))
			 summary(glht(m1))
			
			{# median bout per bird, nest and recorder 
			dd = ddply(d,.(nest, bird_ID_filled, sex), summarise, recorder_y = median(bout_length[recorder=='y']), recorder_n = median(bout_length[recorder=='n']), n_y = length(bout_length[recorder=='y']), n_n = length(bout_length[recorder=='n']))
			dd = dd[!is.na(dd$recorder_y) & !is.na(dd$recorder_n),]
			t.test(dd$recorder_y,dd$recorder_n,paired=TRUE)
			dd1 = dd
			dd1$recorder_n = dd1$n_n = NULL
			names(dd1)[names(dd1) == 'recorder_y'] = 'median_bout'
			names(dd1)[names(dd1) == 'n_y'] = 'n_bouts'
			dd1$recorder = 'y'
			
			dd2 = dd
			dd2$recorder_y = dd2$n_y = NULL
			names(dd2)[names(dd2) == 'recorder_n'] = 'median_bout'
			names(dd2)[names(dd2) == 'n_n'] = 'n_bouts'
			dd2$recorder = 'n'
			dd = rbind(dd1,dd2)
			dd$recorder = as.factor(dd$recorder)
			#dd = ddply(d,.(nest, bird_ID_filled, recorder), summarise, median_bout = median(bout_length), cn_y = length(bout_length[recorder=='y']), n_n = length(bout_length))
			
			pd <- position_dodge(0.4)
			ggplot(dd,aes(x=recorder, y = median_bout, col = factor(bird_ID_filled)))+geom_point(aes(size = n_bouts), position = pd) +geom_line(position = pd, aes(group = factor(bird_ID_filled), col = factor(bird_ID_filled))) +
			theme(legend.position = "none") + ylab('median bird bout length [h]') + xlab('recorder on the nest')
			
			m=lmer(median_bout~recorder+(1|nest),dd, REML=FALSE)
			 summary(glht(m))
			 plot(allEffects(m))
			 }
			{# median bout per  nest and recorder 
			dd = ddply(d,.(nest), summarise, recorder_y = median(bout_length[recorder=='y']), recorder_n = median(bout_length[recorder=='n']), n_y = length(bout_length[recorder=='y']), n_n = length(bout_length[recorder=='n']))
			dd = dd[!is.na(dd$recorder_y) & !is.na(dd$recorder_n),]
			t.test(dd$recorder_y,dd$recorder_n,paired=TRUE)
			dd1 = dd
			dd1$recorder_n = dd1$n_n = NULL
			names(dd1)[2] = 'median_bout'
			names(dd1)[3] = 'n_bouts'
			dd1$recorder = 'y'
			
			dd2 = dd
			dd2$recorder_y = dd2$n_y = NULL
			names(dd2)[2] = 'median_bout'
			names(dd2)[3] = 'n_bouts'
			dd2$recorder = 'n'
			dd = rbind(dd1,dd2)
			dd$recorder = as.factor(dd$recorder)
			#dd = ddply(d,.(nest, bird_ID_filled, recorder), summarise, median_bout = median(bout_length), cn_y = length(bout_length[recorder=='y']), n_n = length(bout_length))
			
			pd <- position_dodge(0.4)
			ggplot(dd,aes(x=recorder, y = median_bout, col = factor(nest)))+geom_point(aes(size = n_bouts), position = pd) +geom_line(position = pd, aes(group = factor(nest), col = factor(nest))) +
			theme(legend.position = "none") + ylab('median nest bout length [h]') + xlab('recorder on the nest')
			
			m=lmer(median_bout~recorder+(1|nest),dd, REML=FALSE)
			 summary(glht(m))
			 plot(allEffects(m))
			 }
			
			{# median nest attentiveness per bird, nest and recorder 
			dd = ddply(d,.(nest, bird_ID_filled), summarise, recorder_y = median(inc_eff[recorder=='y']), recorder_n = median(inc_eff[recorder=='n']), n_y = length(inc_eff[recorder=='y']), n_n = length(inc_eff[recorder=='n']))
			dd = dd[!is.na(dd$recorder_y) & !is.na(dd$recorder_n),]
			t.test(dd$recorder_y,dd$recorder_n,paired=TRUE)
			dd1 = dd
			dd1$recorder_n = dd1$n_n = NULL
			names(dd1)[3] = 'median_attentiveness'
			names(dd1)[4] = 'n_bouts'
			dd1$recorder = 'y'
			
			dd2 = dd
			dd2$recorder_y = dd2$n_y = NULL
			names(dd2)[3] = 'median_attentiveness'
			names(dd2)[4] = 'n_bouts'
			dd2$recorder = 'n'
			dd = rbind(dd1,dd2)
			dd$recorder = as.factor(dd$recorder)
			#dd = ddply(d,.(nest, bird_ID_filled, recorder), summarise, median_bout = median(bout_length), cn_y = length(bout_length[recorder=='y']), n_n = length(bout_length))
			
			pd <- position_dodge(0.4)
			ggplot(dd,aes(x=recorder, y = median_attentiveness, col = factor(bird_ID_filled)))+geom_point(aes(size = n_bouts), position = pd) +geom_line(position = pd, aes(group = factor(bird_ID_filled), col = factor(bird_ID_filled))) +
			theme(legend.position = "none") + ylab('median bird attentiveness') + xlab('recorder on the nest')
			
			m=lmer(median_attentiveness~recorder+(1|nest),dd, REML=FALSE)
			 summary(glht(m))
			 plot(allEffects(m))
			 }
			{# median nest attentiveness per nest and recorder 
			dd = ddply(d,.(nest), summarise, recorder_y = median(inc_eff[recorder=='y']), recorder_n = median(inc_eff[recorder=='n']), n_y = length(inc_eff[recorder=='y']), n_n = length(inc_eff[recorder=='n']))
			dd = dd[!is.na(dd$recorder_y) & !is.na(dd$recorder_n),]
			t.test(dd$recorder_y,dd$recorder_n,paired=TRUE)
			dd1 = dd
			dd1$recorder_n = dd1$n_n = NULL
			names(dd1)[2] = 'median_attentiveness'
			names(dd1)[3] = 'n_bouts'
			dd1$recorder = 'y'
			
			dd2 = dd
			dd2$recorder_y = dd2$n_y = NULL
			names(dd2)[2] = 'median_attentiveness'
			names(dd2)[3] = 'n_bouts'
			dd2$recorder = 'n'
			dd = rbind(dd1,dd2)
			dd$recorder = as.factor(dd$recorder)
			#dd = ddply(d,.(nest, bird_ID_filled, recorder), summarise, median_bout = median(bout_length), cn_y = length(bout_length[recorder=='y']), n_n = length(bout_length))
			
			pd <- position_dodge(0.4)
			ggplot(dd,aes(x=recorder, y = median_attentiveness, col = factor(nest)))+geom_point(aes(size = n_bouts), position = pd) +geom_line(position = pd, aes(group = factor(nest), col = factor(nest))) +
			theme(legend.position = "none") + ylab('median nest attentiveness') + xlab('recorder on the nest')
			
			m=lmer(median_attentiveness~recorder+(1|nest),dd, REML=FALSE)
			 summary(glht(m))
			 plot(allEffects(m))
			 }

}

{# escape distance - careful we have sometimes two estimates (our and gps) for a single event
	load("C:\\Users\\mbulla\\Documents\\Dropbox\\Science\\Projects\\MS\\Removal\\Analyses\\Data\\escape.Rdata")
	str(vv)
	table(vv$nest_sex)
	ggplot(vv,aes(x = nest_sex, y = distm)) + geom_boxplot() + geom_dotplot(aes(fill = nest_sex),binaxis = 'y', stackdir = 'center',  position = position_dodge())
	densityplot(~vv$distm)
	densityplot(~log(vv$distm+0.5))
	
	vg = vv[vv$type == 'gps',]	
	ve = vv[!paste(vv$nest_sex, vv$datetime_)%in%paste(vg$nest_sex, vg$datetime_),]
	v = rbind(vg,ve)
	x = data.frame(table(v$nest_sex))
	summary(factor(x$Freq))
	summary(x$Freq)
	nrow(x)+1 #1 bird had no escape distance estimate)
	
	m = lmer(distm ~ 1 + (1|nest_sex), v)
	summary(m)
	
	xx = x$Var1[x$Freq>1]
	m = lmer(distm ~ 1 + (1|nest_sex), v[v$nest_sex%in%xx,])
		summary(m)
	table(v$nest_sex[v$nest_sex%in%xx])
	
	vv = v
	vv$time = as.numeric(difftime(vv$datetime_,trunc(vv$datetime_,"day"), units = "hours"))
	vv$day_j=as.numeric(format(as.Date(trunc(vv$datetime_, "day")),"%j")) - as.numeric(format(as.Date(trunc(min(vv$datetime_), "day")),"%j"))+1
	ggplot(vv,aes(x = day_j, y = distm)) + geom_point() +stat_smooth(method = 'lm')
		 load("C:\\Users\\mbulla\\Documents\\Dropbox\\Science\\Projects\\MS\\Removal\\Analyses\\Data\\inc_start_estimation_2013_new.Rdata")
		 e$nest = tolower(e$nest)
		vv$inc_start = e$inc_start[match(vv$nest,e$nest)] 	
		vv$day_inc=as.numeric(difftime(vv$datetime_,vv$inc_start, units = "days"))	
	ggplot(vv,aes(x = day_inc, y = distm)) + geom_point() +stat_smooth(method = 'lm')
	ggplot(vv,aes(x = day_inc, y = distm)) + geom_point() +stat_smooth()
	
	m = lmer(distm ~ 1 + (1|nest_sex), vv)
	m = lmer(distm ~ day_j + (day_j|nest_sex), vv)
	m = lmer(distm ~ day_inc + (day_inc|nest_sex), vv)
	m = lmer(distm ~ 1 + (day_inc|nest_sex), vv)
	75.56/(75.56+0.08+269)
	pp = profile(m)
	confint(pp)
	m = lmer(distm ~ day_inc + (day_inc|nest_sex), vv[vv$nest_sex%in%xx,])
	plot(allEffects(m))
	summary(glht(m))
	summary(m)
	
	v = vv[paste(vv$datetime_,vv$nest_sex)%in%paste(vv$datetime_,vv$nest_sex)[duplicated(paste(vv$datetime_,vv$nest_sex))],]
	v = v[order(v$nest_sex, v$datetime_),]
	table(paste(v$nest_sex, v$datetime))
	pd <- position_dodge(0.4)
	ggplot(v,aes(x=type, y = distm, col = factor(paste(nest_sex, datetime_))))+geom_point(position = pd) +geom_line(position = pd, aes(group = factor(paste(nest_sex, datetime_)), col = factor(paste(nest_sex, datetime_)))) +
			theme(legend.position = "none") + ylab('escape distance [m]') + xlab('type of estimation')
	

		
	
}
{# escape distance 2011-2012
 {# tools	
	wet_=dbGetQuery
	
	
	
	db_barr12="M:\\Science\\Projects\\PHD\\FIELD_METHODS\\Barrow\\2012\\DATA\\Database\\Barrow_2012.sqlite"
	con_barr12 = dbConnect(dbDriver("SQLite"),dbname = db_barr12)
	
	db_barr11="M:\\Science\\Projects\\PHD\\FIELD_METHODS\\Barrow\\2011\\DATA\\Database\\Barrow_2011.sqlite"
	con_barr11 = dbConnect(dbDriver("SQLite"),dbname = db_barr11)
}	
 {# Barrow 2011
	v=dbGetQuery(con_barr11,"SELECT date_time as datetime_, nest, dist_left_nest dist_est FROM NESTS where dist_left_nest is not null and species='sesa'")
		v[nchar(v$datetime_)!=nchar('2012-06-12 20:36:00'),]
		v = v[nchar(v$datetime_)>nchar('2012-06-12 20'),]
		v$datetime_ = as.POSIXct(v$datetime_,format="%Y-%m-%d %H:%M:%S", tz = 'America/Anchorage')		
	
	load("M:\\Science\\Projects\\PHD\\STATS\\DATASETS\\weather_SESA_2011_MS_BeEc.Rdata")
		d$nest = tolower(d$nest)
	v = v[v$nest %in% unique(d$nest),]	
			l = list()
			for(i in 1:nrow(v)){
							fi=v[i,]
							di=d[which(d$nest==fi$nest & fi$datetime_>=d$bout_start & fi$datetime_<=d$bout_end+60),] # adding 1 min to the bout end to make the selection more robust
							if(nrow(di)==0){#fi$sex=NA
											 print(fi) #all these visits where before we had rfid on and hence sex cannot be determined
											}else{	fi$sex=di$sex
													fi$bird_ID_filled = di$bird_ID_filled
													l[[i]]=fi
													print(i)
												 }
							
							
							}
	v=do.call(rbind,l)		
	v$year = 2011	
	v11 =v
	}
 {# Barrow 2012
	require(sp)
	v=dbGetQuery(con_barr12,"SELECT visits_pk pk,datetime_, nest, gps_left, wp_left, dist_left_nest dist_est FROM VISITS where gps_left is not null")
	g=dbGetQuery(con_barr12,"SELECT*FROM GPS")
	n=dbGetQuery(conMy,"SELECT nest, species, latit, longit FROM avesatbarrow.nests where year_=2012 and species = 'sesa'")
	v$lat_g=g$latit[match(paste(v$gps_left,v$wp_left),paste(g$gps,g$wp))]
	v$lon_g=g$longit[match(paste(v$gps_left,v$wp_left),paste(g$gps,g$wp))]
	v=v[!is.na(v$lat_g),]
	#v=v[-which(v$pk%in%c(45,259,607,1385)),]
	v$lat=n$latit[match(v$nest,tolower(n$nest))]
	v$lon=n$longit[match(v$nest,tolower(n$nest))]
	v$sp=n$species[match(v$nest,tolower(n$nest))]
		for(i in 1:nrow(v) ) { 
			if(is.na(v$lat[i])){v$distm[i]==NA} else {
					v[i,'distm'] = spDists( as.matrix(v[i, c("lon_g", "lat_g")]), v[i, c("lon", "lat")], longlat=TRUE)*1000
						}
						}
	v=v[v$distm<1250,]					
	v$issues=0				
	v1=v
	# add nests where dist estimated
	v=dbGetQuery(con_barr12,"SELECT visits_pk pk,datetime_, nest, gps_left, wp_left, dist_left_nest dist_est FROM VISITS where dist_left_nest is not null and wp_left is NULL")
	g=dbGetQuery(con_barr12,"SELECT*FROM GPS")
	n=dbGetQuery(conMy,"SELECT nest, species, latit, longit FROM avesatbarrow.nests where year_=2012")
	v$lat_g=g$latit[match(paste(v$gps_left,v$wp_left),paste(g$gps,g$wp))]
	v$lon_g=g$longit[match(paste(v$gps_left,v$wp_left),paste(g$gps,g$wp))]
	
	v$lat=n$latit[match(v$nest,tolower(n$nest))]
	v$lon=n$longit[match(v$nest,tolower(n$nest))]
	v$sp=n$species[match(v$nest,tolower(n$nest))]
	
	v$distm=NA
	v$issues=0
	v2=v
	
		v=rbind(v1,v2)
		v$issues[!is.na(v$dist_est) & !is.na(v$distm)]=1
		v_=v[v$issues==1,]
		ggplot(v_,aes(y=v_$dist_est,x=v_$distm))+geom_point()+stat_smooth(method="lm")
	v2$distm=v2$dist_est
	v1$type="gps"
	v2$type="est"
	v=rbind(v1,v2)
	v = v[!is.na(v$datetime_),]
	v[nchar(v$datetime_)!=nchar('2012-06-12 20:36:00'),]
	v$datetime_ = as.POSIXct(v$datetime_, tz = 'America/Anchorage')		
	
	load("M:\\Science\\Projects\\PHD\\STATS\\DATASETS\\C_SESA_inc_start.Rdata")
		d1 = d
	load("M:\\Science\\Projects\\PHD\\STATS\\DATASETS\\C_SESA_short_inc_start.Rdata")
	d = rbind(d1,d)
	d$nest= tolower(d$nest)	
	v = v[v$nest %in% unique(d$nest),]	
			l = list()
			for(i in 1:nrow(v)){
							fi=v[i,]
							di=d[which(d$nest==fi$nest & fi$datetime_>=d$bout_start & fi$datetime_<=d$bout_end+60),] # adding 1 min to the bout end to make the selection more robust
							if(nrow(di)==0){#fi$sex=NA
											 print(fi) #all these visits where before we had rfid on and hence sex cannot be determined
											}else{	fi$sex=di$sex
													fi$bird_ID_filled = di$bird_ID_filled
													l[[i]]=fi
													print(i)
												 }
							
							
							}
	v=do.call(rbind,l)		
	v$year = 2012	
	v12 =v
}
 {# combine
	names(v11)[names(v11) == 'dist_est'] = 'distm'	
	v11$type = 'est'
	v12 = v12[v12$sp == 'SESA',]
	v = rbind(v12[c('year','datetime_','nest','sex', 'bird_ID_filled','distm','type')], v11)
	nrow(v)
	save(v, file = paste(wd,'escpate_2011-2012.Rdata', sep=''))
	}

 {# analyses
	load(file = paste(wd,'escpate_2011-2012.Rdata', sep=''))
	nrow(v)
	#nrow(v[is.na(v$bird_ID_filled),])
	#v[v$nest == 's510',]
	v$bird_ID_filled[is.na(v$bird_ID_filled)] = 'f'
	v$nest_sex = paste(v$nest, v$sex)
	nrow(v[duplicated(paste(v$nest_sex,v$datetime_)),])
	#v[paste(v$nest_sex,v$datetime_) %in% paste(v$nest_sex,v$datetime_)[duplicated(paste(v$nest_sex,v$datetime_))],]
	#v = v[!duplicated(paste(v$nest_sex,v$datetime_)),]
	#vg = v[v$type == 'gps',]	
	#ve = v[!paste(v$nest_sex, v$datetime_)%in%paste(vg$nest_sex, vg$datetime_),]
	#v = rbind(vg,ve)
	#nrow(v)
	
	x = data.frame(table(v$bird_ID_filled))
	summary(factor(x$Freq))
	summary(x$Freq)
	
	densityplot(~v$distm)
	
	v = v[v$distm<200,]
	
	m = lmer(distm ~ 1 + (1|bird_ID_filled), v)
	m = lmer(distm ~ 1 + (1|bird_ID_filled), v[v$type=='est',])
	summary(m)
	
	xx = x$Var1[x$Freq>1]
	m = lmer(distm ~ 1 + (1|bird_ID_filled), v[v$bird_ID_filled%in%xx,])
		summary(m)
	table(v$nest_sex[v$nest_sex%in%xx])
	
	vv = v
	vv$time = as.numeric(difftime(vv$datetime_,trunc(vv$datetime_,"day"), units = "hours"))
	vv$day_j=as.numeric(format(as.Date(trunc(vv$datetime_, "day")),"%j")) - as.numeric(format(as.Date(trunc(min(vv$datetime_), "day")),"%j"))+1
	ggplot(vv,aes(x = day_j, y = distm)) + geom_point() +stat_smooth(method = 'lm')
		 load("M:\\Science\\Projects\\PHD\\STATS\\DATASETS\\inc_start_est_2011_SESA_NON-SESA.Rdata")
		 e$year = 2011
		 e1 = e
		 
		 load("M:\\Science\\Projects\\PHD\\STATS\\DATASETS\\inc_start_est_2012_SESA_NON-SESA.Rdata")
		  e$year = 2012
		 e2 = e
		 e = rbind(e1,e2)
		 e = e[!is.na(e$inc_est),]
		 e$inc_est = as.POSIXct(e$inc_est,format="%Y-%m-%d %H:%M:%S", tz = 'America/Anchorage')		
		  e$nest = tolower(e$nest)
		vv$inc_start = e$inc_est[match(paste(vv$year, vv$nest),paste(e$year,e$nest))] 	
		vv$day_inc=as.numeric(difftime(vv$datetime_,vv$inc_start, units = "days"))	
	ggplot(vv,aes(x = day_inc, y = distm)) + geom_point() +stat_smooth(method = 'lm')
	ggplot(vv,aes(x = day_inc, y = distm)) + geom_point() +stat_smooth()
	
	m = lmer(distm ~ 1 + (1|bird_ID_filled), vv)
	m = lmer(distm ~ day_j + (day_j|bird_ID_filled), vv)
	m = lmer(distm ~ day_inc + (day_inc|bird_ID_filled), vv,na.action='na.omit')
	m = lmer(distm ~ 1 + (day_inc|bird_ID_filled), vv,na.action='na.omit')
	m = lmer(distm ~ day_inc + (day_inc|bird_ID_filled), vv[vv$bird_ID_filled%in%xx,])
	m = lmer(distm ~ 1 + (day_inc|bird_ID_filled), vv[vv$bird_ID_filled%in%xx,])
	plot(allEffects(m))
	summary(glht(m))
	summary(m)
	pp = profile(m)
	confint(pp)
	library(MCMCglmm)
		mm <- MCMCglmm(distm ~ 1,random = ~us(day_inc):bird_ID_filled, data=vv[vv$bird_ID_filled%in%xx & !is.na(vv$day_inc),])
summary(mm)
}
}	
	
{# hatching success
	d=read.csv(paste(wd,'experiment_metadata.csv', sep=""),stringsAsFactors=FALSE)
	d = d[d$type == 'exp',]
	u=readWorksheetFromFile(paste(wd,'endstates_actogram_based_2016-06-15.xls', sep=""), sheet="temp")
	u$nest = tolower(u$nest)
	u$state = ifelse(is.na(u$state), u$nest_endstate, u$state)
	u_ = ddply(u,.(nest), summarise, end_state = ifelse('hd'%in%state | 'hg'%in%state | '?hg'%in%state |'fl'%in%state | 'hd or fledged?'%in%state |  '?hd'%in%state |  'fledged?'%in%state |  'hs'%in%state |  'hg load MSR first'%in%state|  '1d, 3?hd'%in%state, 'success', ifelse( 'ud'%in%state |'ud?'%in%state |'w'%in%state | 'p? hg?'%in%state | '1b,3w'%in%state | 'un'%in%state | '1b, 1ho, 2'%in%state | 'ud - shall we include ud states?'%in%state, 'unknown', 'failed')))
	
	u_ = u_[substring(u_$nest,1,1) == 's',]
	summary(factor(u_$end_state))
	
	d$end_state = u$state[match(d$nest,u$nest)]
		summary(factor(d$end_state))
	d$end_s = ifelse(d$end_state%in%c('hd','hg','?hg'), 'success', ifelse(d$end_state%in%c('un','ud','un'), 'unknown', 'failed'))
		summary(factor(d$end_s))
	d_$end_s2 = u_$end_state[match(d$nest,u_$nest)]	
		summary(factor(d$end_s2))
		
	
	uc = u_[!u_$nest%in%d$nest,]	
		summary(factor(uc$end_state))
	uu = u
	u$
	
}
	
{# not used - makes little sense
{# bout length
	{# run first - prepare data	
		{# prepare before treatment bout lengths 
			load(file=paste(wd,'bout.Rdata', sep=""))
			d=b[-which(b$nest%in%c('s524')),]
			d$bout_start=as.POSIXct(d$bout_start)
			d$bout_length=d$bout_length/60
			d=d[-which(d$nest=='s402' & d$bout_start>as.POSIXct("2013-06-19 12:00:00")& d$bout_start<as.POSIXct("2013-06-21 02:00:00")),] 
			d[which(d$nest=='s404' & d$exper=='t'),3:9]
			d[which(d$nest=='s404'),c(4:9,12)]
			min(d$bout_ID[which(!is.na(d$exper) & d$exper=='t' & d$nest=='s404')])
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
													
				db=do.call(rbind, foo)
			ggplot(db,aes(x=nest,y=bout_length,fill=nest))+geom_boxplot()+geom_point()
			#dd[which(dd$nest=='s402'),]
			#dd[which(dd$nest=='s808'),]
			dm=ddply(db,.(nest,sex),summarise, bout=median(bout_length))
											# test if sex correctly assigned -YES
												g=read.csv(paste(wd,'experiment_metadata.csv', sep=""), stringsAsFactors=FALSE)
												g=g[g$type=='exp',]
												dm$taken=g$taken[match(dm$nest,g$nest)]
												dm[dm$sex==dm$taken,]
			}									

		{# preapare treated bout length
		load(file=paste(wd,'50_51.Rdata', sep=""))
		d=b[-which(b$nest%in%c('s524')),]
		d$bout_start=as.POSIXct(d$bout_start)
		dd=d
		dd$col_=ifelse(dd$sex=="f", "#FCB42C","#535F7C")
		dd$bout_length=dd$bout_length/60
		}

		{# bring them together
		dd_=dd[,c('nest', 'sex', 'exper', 'bout_start', 'bout_length')]
		db_=db[,c('nest', 'sex', 'exper', 'bout_start', 'bout_length')]
		dd_$exper='t'
		db_$exper='b'
		
		g=rbind(dd_,db_)
		
			w=read.csv(paste(wd,'experiment_metadata.csv', sep=""), stringsAsFactors=FALSE)
			ne=w$nest[w$type=='exp' & w$use_t=='y']
		g_=g[g$nest%in%ne,]	
		}
	}
	{# model
		g$exper=as.factor(g$exper)
		g$sex=as.factor(g$sex)
		g$exper_sex=as.factor(ifelse(g$exper=='b', 'b', ifelse(g$sex=='f', 't_female', 't_male')))
		g$exper_sex_2=relevel(g$exper_sex,ref="t_female")
		g$exper_sex_int=as.factor(ifelse(g$exper=='b' & g$sex=='f', 'b_female',ifelse(g$exper=='b' & g$sex=='m', 'b_male', ifelse(g$sex=='f', 't_female', 't_male'))))
		
		m=lmer(bout_length~exper+(1|nest),g)






		m2=lmer(bout_length~exper*sex+(1|nest),g)
		m3=lmer(bout_length~exper_sex_2+(1|nest),g)
		m4=lmer(bout_length~exper_sex_int+(1|nest),g)
		
		AIC(m,m2,m3,m4)
		AIC(m,m2)
		
		summary(m)
		summary(m2)
		summary(m3)
		summary(m4)
		plot(allEffects(m))
		plot(allEffects(m2))
		plot(allEffects(m3))
		plot(allEffects(m4))
		o=summary(glht(m3, linfct = mcp(exper_sex = c(	"t_female - t_male = 0",
															"t_female - c = 0",
															"t_male - c = 0"
															))))
		summary(glht(m4, linfct = mcp(exper_sex_int = c(	"t_female - t_male = 0",
															"b_female - b_male = 0"
															))))			
					oo=confint(o)
		t.test(bb$inc_eff, bt$inc_eff,paired=TRUE, alternative="greater") 
		t.test(bb$inc_eff, bt$inc_eff, alternative="greater") 
		}
	{# plot
		dev.new(width=2.5,height=2.5)
					#png(paste(out,"figs/tree_period_high_.png", sep=""), width=2.3,height=3.5,units="in",res=600)
		par(mar=c(2,2,0.7,0.7), ps=12, mgp=c(1.2,0.35,0), las=1, cex.lab=0.7, cex.axis=0.6, tcl=-0.15,bty="l",xpd=TRUE, col.axis="grey30", col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt
			
		plot(dd$bout_length~dm$bout, col=dd$col_,bg=adjustcolor(dd$col_ ,alpha.f = 0.4), pch=21, ylab='Bout length - treated bout [h]', xlab=NA, ylim=c(0,120),xlim=c(10,15), xaxt='n')
		lines(c(0,15),c(0,15), lty=3, col="lightgrey", xpd=FALSE)
		axis(1, at=seq(10,15,by=1), labels=FALSE)
		text(seq(10,15,by=1), par("usr")[3]-6, labels = seq(10,15,by=1),  xpd = TRUE, cex=0.6, col="grey30") #labels = z_g$genus
		mtext('Bout Length - before bout [h]',side=1,line=0.6, cex=0.7, las=1, col='grey30')
		
		legend("topleft", legend=c('\u2640','\u2642'), pch=21, col=c( '#FCB42C','#535F7C'), bg=adjustcolor(c( '#FCB42C','#535F7C'),alpha.f = 0.4),cex=0.8, bty='n')
		
		estimate_m=data.frame(y = c(0.91273,  (0.91273-0.36650)),x=c('c','t'), stringsAsFactors=FALSE)
		estimate_t=data.frame(y = c( 0.9270140, 0.5314018),x=c('c','t'), stringsAsFactors=FALSE)
		
		ggplot(g,aes(y=bout_length, x=exper, fill=exper))+geom_point(aes(color=exper),alpha=0.5,position = position_jitterdodge(dodge.width=0.9,jitter.width=2))+geom_boxplot(fill=NA,alpha=0.5,position = position_dodge(width=0.9))+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_color_manual(values=c("#E69F00", "#56B4E9"))+theme(legend.position="none")
		
		geom_point(data=estimate_t,aes(y = y,x=x),fill="red", colour="red",size=2)+ geom_point(data=estimate_m,aes(y = y,x=x),fill="pink", colour="pink",size=2)
		ggplot(g,aes(x=exper, y=bout_length))+geom_point()
	}
	{#limited
	{# model
		g_$exper=as.factor(g_$exper)
		g_$sex=as.factor(g_$sex)
		g_$exper_sex=as.factor(ifelse(g_$exper=='b', 'b', ifelse(g_$sex=='f', 't_female', 't_male')))
		g_$exper_sex_2=relevel(g_$exper_sex,ref="t_female")
		g_$exper_sex_int=as.factor(ifelse(g_$exper=='b' & g_$sex=='f', 'b_female',ifelse(g_$exper=='b' & g_$sex=='m', 'b_male', ifelse(g_$sex=='f', 't_female', 't_male'))))
		
		m=lmer(bout_length~exper+(1|nest),g_)
		m2=lmer(bout_length~exper*sex+(1|nest),g_)
		m3=lmer(bout_length~exper_sex_2+(1|nest),g_)
		m4=lmer(bout_length~exper_sex_int+(1|nest),g_)
		
		AIC(m,m2,m3,m4)
		AIC(m,m2)
		
		summary(m)
		summary(m2)
		summary(m3)
		summary(m4)
		plot(allEffects(m))
		plot(allEffects(m2))
		plot(allEffects(m3))
		plot(allEffects(m4))
		o=summary(glht(m3, linfct = mcp(exper_sex = c(	"t_female - t_male = 0",
															"t_female - c = 0",
															"t_male - c = 0"
															))))
		summary(glht(m4, linfct = mcp(exper_sex_int = c(	"t_female - t_male = 0",
															"b_female - b_male = 0"
															))))			
					oo=confint(o)
		t.test(bb$inc_eff, bt$inc_eff,paired=TRUE, alternative="greater") 
		t.test(bb$inc_eff, bt$inc_eff, alternative="greater") 
		}
	{# plot
		dev.new(width=2.5,height=2.5)
					#png(paste(out,"figs/tree_period_high_.png", sep=""), width=2.3,height=3.5,units="in",res=600)
		par(mar=c(2,2,0.7,0.7), ps=12, mgp=c(1.2,0.35,0), las=1, cex.lab=0.7, cex.axis=0.6, tcl=-0.15,bty="l",xpd=TRUE, col.axis="grey30", col.lab="grey30", col.main="grey30", fg="grey70") # 0.6 makes font 7pt, 0.7 8pt

		plot(dd$bout_length~dm$bout, col=dd$col_,bg=adjustcolor(dd$col_ ,alpha.f = 0.4), pch=21, ylab='Bout length - treated bout [h]', xlab=NA, ylim=c(0,120),xlim=c(10,15), xaxt='n')
		lines(c(0,15),c(0,15), lty=3, col="lightgrey", xpd=FALSE)
		axis(1, at=seq(10,15,by=1), labels=FALSE)
		text(seq(10,15,by=1), par("usr")[3]-6, labels = seq(10,15,by=1),  xpd = TRUE, cex=0.6, col="grey30") #labels = z_g_$genus
		mtext('Bout Length - before bout [h]',side=1,line=0.6, cex=0.7, las=1, col='grey30')
		
		legend("topleft", legend=c('\u2640','\u2642'), pch=21, col=c( '#FCB42C','#535F7C'), bg=adjustcolor(c( '#FCB42C','#535F7C'),alpha.f = 0.4),cex=0.8, bty='n')
		
		estimate_m=data.frame(y = c(0.91273,  (0.91273-0.36650)),x=c('c','t'), stringsAsFactors=FALSE)
		estimate_t=data.frame(y = c( 0.9270140, 0.5314018),x=c('c','t'), stringsAsFactors=FALSE)
		
		ggplot(g,aes(y=bout_length, x=exper, fill=exper))+geom_point(aes(color=exper),alpha=0.5,position = position_jitterdodge(dodge.width=0.9,jitter.width=2))+geom_boxplot(fill=NA,alpha=0.5,position = position_dodge(width=0.9))+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_color_manual(values=c("#E69F00", "#56B4E9"))+theme(legend.position="none")
		
		geom_point(data=estimate_t,aes(y = y,x=x),fill="red", colour="red",size=2)+ geom_point(data=estimate_m,aes(y = y,x=x),fill="pink", colour="pink",size=2)
		ggplot(g,aes(x=exper, y=bout_length))+geom_point()


	
}
}



{# results cons
load(paste(wd,'const.Rdata',sep=""))
	b=b[which(b$bout_ID_tr>5),]
	b=b[-which(b$nes%in%c('s806','s524')),]
	b=b[!(is.na(b$inc_eff)),]
	b$bout_length=b$bout_length/60
	b$exper=as.factor(b$exper)
	b$exper_sex=as.factor(ifelse(b$exper=='b', 'b', ifelse(b$sex=='f', 't_female', 't_male')))
	b$exper_sex_2=relevel(b$exper_sex,ref="t_female")
	b[b$bout_length<18,]
	b$col_=ifelse(b$exper=="t", 'orange','slateblue')
	bb=b[b$exper=='b',]
	# distribution of constancy in before treatment bouts
			ggplot(bb,aes(y=inc_eff, x=nest))+geom_boxplot()+ theme(axis.text.x=element_text(angle = 90, hjust =1, vjust=0.2))+coord_cartesian(ylim = c(0.68, 1)) 
	
	bt=b[b$exper=='t',]
	bb$inc_eff_t=bt$inc_eff[match(bb$nest, bt$nest)]
	ggplot(bb,aes(x=inc_eff, y=inc_eff_t))+geom_point()+stat_smooth()

		m=lmer(inc_eff~bout_length+exper+(1|nest),b)
		m2=lmer(inc_eff~bout_length+exper_sex+(1|nest),b)
		AIC(m,m2)
		summary(m)
		summary(m2)
		plot(allEffects(m))
		plot(allEffects(m2))
		o=summary(glht(m2, linfct = mcp(exper_sex = c(	"t_female - t_male = 0",
															"t_female - b = 0",
															"t_male - b = 0"
															))))
					oo=confint(o)
		m3=lmer(inc_eff~bout_length+exper_sex_2+(1|nest),b)
		summary(m3)
		ggplot(b,aes(y=inc_eff, x=exper, fill=exper))+geom_point(aes(color=exper),position = position_jitterdodge(dodge.width=0.9))+geom_boxplot(alpha=0.5,position = position_dodge(width=0.9))+scale_fill_manual(values = c("#FCB42C","#535F7C"))+scale_color_manual(values=c("#FCB42C","#535F7C"))  #"#E69F00", "#56B4E9"
		
		estimate=data.frame(y = c( 0.857663,  (0.857663-0.206913)),x=c('b','t'), stringsAsFactors=FALSE)
		
		ggplot(b,aes(y=inc_eff, x=exper, fill=exper))+geom_point(aes(color=exper),alpha=0.5,position = position_jitterdodge(dodge.width=0.9,jitter.width=2))+geom_boxplot(fill=NA,alpha=0.5,position = position_dodge(width=0.9))+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_color_manual(values=c("#E69F00", "#56B4E9"))+ geom_point(data=estimate,aes(y = y,x=x),fill="red", colour="red",size=2)
		size=2)
		#"#E69F00", "#56B4E9"
		
		ggplot(b,aes(y=inc_eff, x=exper, fill=exper))+geom_point(aes(color=exper),alpha=0.5,position = position_jitterdodge(dodge.width=0.9))+geom_boxplot(aes(colour=exper),fill=NA,position = position_dodge(width=0.9))+scale_fill_manual(values = c("#E69F00", "#56B4E9"))+scale_color_manual(values=c("#E69F00", "#56B4E9"))  #"#E69F00", "#56B4E9"
		+geom_point(position=geom_jitter(width=0.5))
		summary(m)
		
		# sex ggplot
		estimate_sex=data.frame(y = c( 0.872710,  (0.872710-0.181997),(0.872710-0.228490)),x=c('b','t_female','t_male'), stringsAsFactors=FALSE)
		ggplot(b,aes(y=inc_eff, x=exper_sex, fill=exper_sex))+geom_point(aes(color=exper_sex),alpha=0.5,position = position_jitterdodge(dodge.width=0.9,jitter.width=2))+geom_boxplot(fill=NA,alpha=0.5,position = position_dodge(width=0.9))+scale_fill_manual(values = c("#5eab2b","#E69F00", "#56B4E9"))+scale_color_manual(values=c("#5eab2b","#E69F00", "#56B4E9"))+ geom_point(data=estimate_sex,aes(y = y,x=x),fill="red", colour="red",size=2)
		bb=ddply(b,.(nest),summarise, inc_eff=median(inc_eff,na.rm=TRUE))
		bb$exper='b'
		bb=bb[!is.na(bb$inc_eff),]
		bt=b[which(b$exper=='t'),c('nest','inc_eff','exper','sex')]
		bt=bt[!is.na(bt$inc_eff),]
		bt$inc_eff_before=bb$inc_eff[match(bt$nest,bb$nest)]
		bt$col_=ifelse(bt$sex=="f", "#FCB42C","#535F7C")
		
		t.test(bb$inc_eff, bt$inc_eff,paired=TRUE, alternative="greater") 
		t.test(bb$inc_eff, bt$inc_eff, alternative="greater") 
				
		plot(bt$inc_eff~bt$inc_eff_before, xlim=c(0.4,1), ylim=c(0.4,1), col=bt$col_, pch=19, ylab='incubation constancy - experimental bout', xlab='incubation constancy - before bout')
		abline(a=0,b=1)
		0.1990582/(0.9300718/100)
		pr=data.frame(difference_in_constancy_percent=(bb$inc_eff-bt$inc_eff)/(bb$inc_eff/100), boxp='boxplot')
		ggplot(pr,aes(x=boxp, y=difference_in_constancy_percent))+geom_boxplot()+coord_cartesian(ylim = c(-5,60)) 
		
		pr=data.frame(difference_in_constancy=(bb$inc_eff-bt$inc_eff), boxp='boxplot')
		ggplot(pr,aes(x=boxp, y=difference_in_constancy))+geom_boxplot()+coord_cartesian(ylim = c(-0.05,0.6)) 
}
}