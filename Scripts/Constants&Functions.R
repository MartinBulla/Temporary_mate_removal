{# load packages
		require(AICcmodavg)
		require(Amelia)
		require(arm)
		devtools::source_gist('45b49da5e260a9fc1cd7') #loads IwantHue fuctions https://gist.github.com/johnbaums/45b49da5e260a9fc1cd7
		require(effects)
		require(ggplot2)
		require(lattice)
		require(MCMCglmm)
		require(multcomp)
		require(plyr)
		require(randomcoloR) #install.packages('randomcoloR')
		require(RColorBrewer)
		require(RMySQL)
		require(XLConnect)
		require(zoo)
}
{# connect to database
	conMy=dbConnect(MySQL(),user='root',host='127.0.0.1', password='',dbname='')
	wet=dbGetQuery	
}
# define time and load packages
		Sys.setenv(TZ="UTC")	
# define constants
		fo=72/25.4 # use to multiply size of font (in mm) in ggplot2 to get to point size from base package 
# functions
	Mode <- function(x) {
			ux <- unique(x)
			ux[which.max(tabulate(match(x, ux)))]
			}
		
 