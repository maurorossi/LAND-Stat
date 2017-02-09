#!/opt/anaconda/bin/Rscript --vanilla --slave --quiet


#########################################################################
#########################################################################
####                                                                 ####
####                                                                 ####
####                                                                 ####
####                          LAND-Stat                              ####
####            LANDslide size Statistics evaluation                 ####
####                           IRPI CNR                              ####
####                    MAURO ROSSI - IRPI CNR                       ####
####                  v1r0b15 - 16 November 2016                     ####
####                                                                 ####
#### Copyright (C) 2016 Mauro Rossi                                  ####
####                                                                 ####
#### This program is free software; you can redistribute it and/or   ####
#### modify it under the terms of the GNU General Public License     ####
#### as published by the Free Software Foundation; either version 2  ####
#### of the License, or (at your option) any later version. ###      ####
####                                                                 ####
#### This program is distributed in the hope that it will be useful, ####
#### but WITHOUT ANY WARRANTY; without even the implied warranty of  ####
#### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the ####
#### GNU General Public License for more details.                    ####
####                                                                 ####
####     Istituto di Ricerca per la Protezione Idrogeologica         ####
####              Consiglio Nazionale delle Ricerche                 ####
####                    Gruppo di Geomorfologia                      ####
####                  Via della Madonna Alta, 126                    ####
####                    06128 Perugia (Italia)                       ####
####                       +39 075 5014421                           ####
####                       +39 075 5014420                           ####
####                   mauro.rossi@irpi.cnr.it                       ####
####                  geomorfologia@irpi.cnr.it                      ####
####                                                                 ####
####           This script was prepared using R 3.0.1                ####
####         The script requires the following R packages:           ####
####                       1: rgdal                                  ####
####                       2: rgeos                                  ####
####                       3: bbmle                                  ####
####                       4: Matching                               ####
####                                                                 ####
####     INPUTS: 1) data.txt file (tab delimited)                    ####
####                one named column -> name will be appended        ####
####                to resultsidentification value                   ####
####                OR                                               ####
####                training.shp (polygons or points)                ####
####                with a column with the size values               ####
####                                                                 ####
####             2) configuration.txt file (tab delimited)           ####
####                1st column -> parameter name                     ####
####                2nd column -> min value of parameter             ####
####                3rd column -> max value of parameter             ####
####                4th column -> specified_by_user can be           ####
####                              YES or NO (default) to enable the  ####
####                              use of min and max parameter       ####
####                              values                             ####
####                                                                 ####
#########################################################################
#########################################################################

###################### WPS ESA modification
# rm(list=(ls()))
# graphics.off()
# workdir<-"X:/R/MLE_LandslideArea/Run_Volume_Sandra_20161109/Salzberg" # For windows
# #workdir<-"/media/disco_dati/R/MLE_LandslideArea/Run_Romania_20150525"  # For linux
# setwd(workdir)
# #setwd("/media/disco_dati/R/MLE_LandslideArea/TEST_KCL_20140211")
# #memory.limit(size=16000)

#workdir<-paste(TMPDIR,"output_landstat",sep="/")
#unlink(workdir,recursive=TRUE,force=TRUE)
#dir.create(workdir)
#setwd(workdir)

library("rciop")
param_configuration_file_name <- rciop.getparam("configuration_file_name")
res_configuration<-rciop.copy(param_configuration_file_name, TMPDIR, uncompress=TRUE)
if (res_configuration$exit.code==0) local.url.configuration <- res_configuration$output
configuration<-read.table(local.url.configuration,header=FALSE,skip=1,dec=".", sep="\t",as.is=TRUE)
file.remove(local.url.configuration)
# configuration<-read.table("configuration.txt",header = FALSE,skip=1,dec=".", sep="\t",as.is=TRUE)
######################
	
summary(configuration)
str(configuration)
names(configuration)


### Selcting the type of data to be analyzed
xlabel<-expression(A~~(m^{2})) # For area
#xlabel<-expression(V~~(m^{3})) # For volume
#xlabel<-expression(X~~(m)) # For a length  
#xlabel<-expression(frac(L, W)~~bgroup("(",frac(m,m),")")) # For a length/width ratio


ks_boot_samples<-100
use_shape<-FALSE
	

###################### WPS ESA modification
param_data_file_name <- rciop.getparam("data_file_name")
res_data<-rciop.copy(param_data_file_name, TMPDIR, uncompress=TRUE)
if (res_data$exit.code==0) local.url.data <- res_data$output
name_file_data<-local.url.data
# name_file_data<-"Data1_nocolname.txt" # .txt or .shp if use_shape=TRUE
######################
	
use_shape_field_stat<-TRUE
shape_field_stat<-"area_sqrm"

bin_method<-c("Sturges") # "Sturges" or "scott" or "FD"

cdf_percentiles<-c(0.01,0.05,0.25,0.5,0.75,0.95,0.99)
executing_CDF_sensitivity_analysis<-FALSE

if (use_shape==FALSE)
  {
  data.series<-read.table(name_file_data,header = TRUE,dec=".", sep="\t",as.is=TRUE)
  #data.series<-read.table("data.txt",header = TRUE,dec=".", sep="\t",as.is=TRUE)
  #data.series<-data.series[order(data.series),]
  rowlen<-dim(data.series)[1]
  data.series.results<-data.frame(data.series,hdedp_pdf=as.numeric(rep(NA,rowlen)),hdedp_cdf=as.numeric(rep(NA,rowlen)),hdedps_pdf=as.numeric(rep(NA,rowlen)),hdedps_cdf=as.numeric(rep(NA,rowlen)),hdeig_pdf=as.numeric(rep(NA,rowlen)),hdeig_cdf=as.numeric(rep(NA,rowlen)),kdedp_pdf=as.numeric(rep(NA,rowlen)),kdedp_cdf=as.numeric(rep(NA,rowlen)),kdedps_pdf=as.numeric(rep(NA,rowlen)),kdedps_cdf=as.numeric(rep(NA,rowlen)),kdeig_pdf=as.numeric(rep(NA,rowlen)),kdeig_cdf=as.numeric(rep(NA,rowlen)),mledp_pdf=as.numeric(rep(NA,rowlen)),mledp_cdf=as.numeric(rep(NA,rowlen)),mledps_pdf=as.numeric(rep(NA,rowlen)),mledps_cdf=as.numeric(rep(NA,rowlen)),mleig_pdf=as.numeric(rep(NA,rowlen)),mleig_cdf=as.numeric(rep(NA,rowlen)))
  } else
  {
  library(rgdal)
  data.shape<-readOGR(dsn=name_file_data,layer=gsub(".shp","",name_file_data))
  library(rgeos)
  data.shape@data<-data.frame(cbind(data.shape@data,areargeos=gArea(data.shape,byid=TRUE),stringsAsFactors = FALSE))
  if (use_shape_field_stat==TRUE)
    {
    data.series<-as.data.frame(data.shape@data[,which(names(data.shape@data)==shape_field_stat)])
    names(data.series)<-shape_field_stat
    } else
    {
    data.series<-data.frame(areargeos=data.shape@data$areargeos)
    }
  rowlen<-dim(data.series)[1]
  data.series.results<-data.frame(data.shape@data,hdedp_pdf=as.numeric(rep(NA,rowlen)),hdedp_cdf=as.numeric(rep(NA,rowlen)),hdedps_pdf=as.numeric(rep(NA,rowlen)),hdedps_cdf=as.numeric(rep(NA,rowlen)),hdeig_pdf=as.numeric(rep(NA,rowlen)),hdeig_cdf=as.numeric(rep(NA,rowlen)),kdedp_pdf=as.numeric(rep(NA,rowlen)),kdedp_cdf=as.numeric(rep(NA,rowlen)),kdedps_pdf=as.numeric(rep(NA,rowlen)),kdedps_cdf=as.numeric(rep(NA,rowlen)),kdeig_pdf=as.numeric(rep(NA,rowlen)),kdeig_cdf=as.numeric(rep(NA,rowlen)),mledp_pdf=as.numeric(rep(NA,rowlen)),mledp_cdf=as.numeric(rep(NA,rowlen)),mledps_pdf=as.numeric(rep(NA,rowlen)),mledps_cdf=as.numeric(rep(NA,rowlen)),mleig_pdf=as.numeric(rep(NA,rowlen)),mleig_cdf=as.numeric(rep(NA,rowlen)))
  }

summary(data.series)
str(data.series)
names(data.series)
range(data.series,na.rm=TRUE)

dir.create(paste(getwd(),"/Results",sep=""))
setwd(paste(getwd(),"/Results",sep=""))

#### Box Plot
centers<-1:dim(data.series)[2]
color.series<-heat.colors(dim(data.series)[2], alpha = 1)
if (Sys.info()[1] == "Linux") {x11()}
if (Sys.info()[1] == "Windows") {windows()}
if (Sys.info()[1] == "Darwin") {quartz()}
par(cex.axis=0.5)
for (center.count in 1:dim(data.series)[2])
   {
   if (center.count == 1)
    {
    boxplot(data.series[,center.count],main=paste("Box Plot comparison",sep=""),plot=TRUE,at=centers[center.count],log="y",xlim=c(0,dim(data.series)[2]+1),ylim=range(data.series,na.rm=TRUE),width=10,border="black",col=color.series[center.count])
    } else
    {
    boxplot(data.series[,center.count],add=TRUE, plot=TRUE,at=centers[center.count],border="black",col=color.series[center.count],width=10)
    }
    lines(c(centers[center.count]-0.4,centers[center.count]+0.4),rep(mean(na.omit(data.series[,center.count])),2),lty=2,lwd=1)
    points(centers[center.count],mean(na.omit(data.series[,center.count])),pch=4,lwd=1)
    points(centers[center.count],min(na.omit(data.series[,center.count])),pch=24,lwd=1,col="black",bg="darkgreen")
    points(centers[center.count],max(na.omit(data.series[,center.count])),pch=25,lwd=1,col="black",bg="darkgreen")	
    axis(side=1, at=center.count, labels = paste(names(data.series)[center.count],sep=""))
   }


pdf(file = paste("BoxPlot_Comparison.pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
par(cex.axis=0.4)
for (center.count in 1:dim(data.series)[2])
   {
   if (center.count == 1)
    {
    boxplot(data.series[,center.count],main=paste("Box Plot comparison",sep=""),plot=TRUE,at=centers[center.count],log="y",xlim=c(0,dim(data.series)[2]+1),ylim=range(data.series,na.rm=TRUE),width=10,border="black",col=color.series[center.count])
    } else
    {
    boxplot(data.series[,center.count],add=TRUE, plot=TRUE,at=centers[center.count],border="black",col=color.series[center.count],width=10)
    }
    lines(c(centers[center.count]-0.4,centers[center.count]+0.4),rep(mean(na.omit(data.series[,center.count])),2),lty=2,lwd=1)
    points(centers[center.count],mean(na.omit(data.series[,center.count])),pch=4,lwd=1)
    points(centers[center.count],min(na.omit(data.series[,center.count])),pch=24,lwd=1,col="black",bg="darkgreen")
    points(centers[center.count],max(na.omit(data.series[,center.count])),pch=25,lwd=1,col="black",bg="darkgreen")	
    axis(side=1, at=center.count, labels = paste(names(data.series)[center.count],sep=""))
   }
dev.off()

### Single Box Plot
for (center.count in 1:dim(data.series)[2])
   {
#   if (Sys.info()[1] == "Linux") {x11()}
#	if (Sys.info()[1] == "Windows") {windows()}
#	if (Sys.info()[1] == "Darwin") {quartz()}
#   boxplot(data.series[,center.count],main=paste("Box Plot: ",names(data.series)[center.count],sep=""),plot=TRUE,log="y",at=1,ylim=range(na.omit(data.series)),border="black",col=color.series[center.count])
#   lines(c(1-0.2,1+0.2),rep(mean(na.omit(data.series[,center.count])),2),lty=2,lwd=1)
#   points(1,mean(na.omit(data.series[,center.count])),pch=4,lwd=1)
#   axis(side=1, at=1, labels = paste(names(data.series)[center.count],sep=""))
#
   pdf(file = paste("BoxPlot_",names(data.series)[center.count],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
   par(cex.axis=0.6)
   boxplot(data.series[,center.count],main=paste("Box Plot: ",names(data.series)[center.count],sep=""),plot=TRUE,log="y",at=1,ylim=range(data.series,na.rm=TRUE),border="black",col=color.series[center.count])
   lines(c(1-0.2,1+0.2),rep(mean(na.omit(data.series[,center.count])),2),lty=2,lwd=1)
   points(1,mean(na.omit(data.series[,center.count])),pch=4,lwd=1)
   points(1,min(na.omit(data.series[,center.count])),pch=24,lwd=1,col="black",bg="darkgreen")
   points(1,max(na.omit(data.series[,center.count])),pch=25,lwd=1,col="black",bg="darkgreen")	
   axis(side=1, at=1, labels = paste(names(data.series)[center.count],sep=""))
   dev.off()
   }


#library(stats4)
#mle()
#library(maxLik)
#maxLik()

library(bbmle)

##### Double Pareto 5 parameter
ddoublepareto <- function(x, alpha.par, beta.par, t.par, c.par, m.par, FUN=FALSE) # alpha.par > 0, beta.par > 0, 0 <= c.par < t.par <=  m.par <= inf
  {
  if(FUN==TRUE)
    {
    density.function<-function(x)
      {
      alpha.par=alpha.par
      beta.par=beta.par
      t.par=t.par
      c.par=c.par
      m.par=m.par
      #fun<-(beta.par)/(t.par*(1-((1+(m.par/t.par)^(-alpha.par))/(1+(c.par/t.par)^(-alpha.par)))^(beta.par/alpha.par)))*(((1+(m.par/t.par)^(-alpha.par))^(beta.par/alpha.par))/((1+(x/t.par)^(-alpha.par))^(1+(beta.par/alpha.par))))*((x/t.par)^(-alpha.par-1))

      
      delta.function<-((1+(m.par/t.par)^(-alpha.par))/(1+(c.par/t.par)^(-alpha.par)))^(beta.par/alpha.par)
      eta.function<-(beta.par)/(t.par*(1-delta.function))
      fun<-eta.function*(((1+(m.par/t.par)^(-alpha.par))^(beta.par/alpha.par))/((1+(x/t.par)^(-alpha.par))^(1+(beta.par/alpha.par))))*((x/t.par)^(-alpha.par-1))
      
      return(fun)
      
      }
    return(density.function)
    } else
    {
    d.doublepareto<-(beta.par)/(t.par*(1-((1+(m.par/t.par)^(-alpha.par))/(1+(c.par/t.par)^(-alpha.par)))^(beta.par/alpha.par)))*(((1+(m.par/t.par)^(-alpha.par))^(beta.par/alpha.par))/((1+(x/t.par)^(-alpha.par))^(1+(beta.par/alpha.par))))*((x/t.par)^(-alpha.par-1))
    return(d.doublepareto)
    }
  }

pdoublepareto <- function(x, alpha.par, beta.par, t.par, c.par, m.par, FUN=FALSE) # Integration with integrate
  {
  if(FUN==TRUE)
    {
    cumdensity.function<-function(x)
      {
      density.function<-function(x)
        {
        alpha.par=alpha.par
        beta.par=beta.par
        t.par=t.par
        c.par=c.par
        m.par=m.par
        denfun<-(beta.par)/(t.par*(1-((1+(m.par/t.par)^(-alpha.par))/(1+(c.par/t.par)^(-alpha.par)))^(beta.par/alpha.par)))*(((1+(m.par/t.par)^(-alpha.par))^(beta.par/alpha.par))/((1+(x/t.par)^(-alpha.par))^(1+(beta.par/alpha.par))))*((x/t.par)^(-alpha.par-1))
        return(denfun)
        }
      #fun<-integrate(f=density.function,lower=0.0001,upper=x)
      fun<-integrate(f=density.function,lower=c.par,upper=x)
      
      return(fun$value)
      }
    cumdensity.function<-Vectorize(cumdensity.function)
    return(cumdensity.function)	
    } else
    {
    cumdensity.function<-function(x)
      {
      density.function<-function(x)
        {
        alpha.par=alpha.par
        beta.par=beta.par
        t.par=t.par
        c.par=c.par
        m.par=m.par
        fun<-(beta.par)/(t.par*(1-((1+(m.par/t.par)^(-alpha.par))/(1+(c.par/t.par)^(-alpha.par)))^(beta.par/alpha.par)))*(((1+(m.par/t.par)^(-alpha.par))^(beta.par/alpha.par))/((1+(x/t.par)^(-alpha.par))^(1+(beta.par/alpha.par))))*((x/t.par)^(-alpha.par-1))
        return(fun)
        }
      #p.doublepareto<-integrate(f=density.function,lower=0.0001,upper=x)
      p.doublepareto<-integrate(f=density.function,lower=c.par,upper=x)
      
      return(p.doublepareto$value)
      }
    cumdensity.function<-Vectorize(cumdensity.function)
    return(cumdensity.function(x))	
    }
  }
 
# ## Testing pdoublepareto and ddoublepareto
# c<-pdoublepareto(alpha.par=1.8, beta.par=1.5, t.par=10000, c.par=50, m.par=100000, FUN=TRUE)
# c(seq(50,100000,10))
# pdoublepareto(seq(50,100000,10),alpha.par=1.8, beta.par=1.5, t.par=10000, c.par=50, m.par=100000,FUN=FALSE)
# a<-ddoublepareto(seq(50,100000,10),alpha.par=1.8, beta.par=1.5, t.par=10000, c.par=50, m.par=100000, FUN=TRUE)
# a(seq(50,100000,10))
# 
# values<-seq(50,100000,10)
# den<-ddoublepareto(values,alpha.par=1.8, beta.par=1.5, t.par=10000, c.par=50, m.par=100000)
# plot(values,den,type="l",col="black",log="xy")
# cumden<-pdoublepareto(values,alpha.par=3, beta.par=1.5, t.par=10000, c.par=50, m.par=100000)
# windows()
# plot(values,cumden,type="l",col="black")

uniroot.range<-range(data.series)
# Paramenter value specified by the user, specified_by_user == "YES" in configuration.txt file
if (configuration[30,4]=="YES") (uniroot.range<-as.numeric(configuration[30,2:3]))



rdoublepareto <- function(n, alpha.par, beta.par, t.par, c.par, m.par) 
  {
  cumprobunif<-sort(runif(n))
  fun <- function (x,u)
    {
    alpha=alpha.par
    beta=beta.par
    t=t.par
    c=c.par
    m=m.par
    cumfun<-pdoublepareto(alpha.par=alpha, beta.par=beta, t.par=t, c.par=c, m.par=m,FUN=TRUE)
    return(cumfun(x)-u)
    }
  z <- matrix(NA,nrow=n,ncol=2)
  z[,2]<-cumprobunif
  #colnames(z)<-c("value","cumsumprob")
  for (i in 1:n)
    {
    #print(paste("u=",cumprobunif[i],sep=""))
	if (inherits(try(z[i,1] <- uniroot(fun, lower=uniroot.range[1],upper=uniroot.range[2], tol = 0.0001, u = cumprobunif[i])$root,silent=TRUE),what="try-error"))
		{
		next()
		}	
    #print(paste("value=",z[i,1]))
    #if (z[i,1]<0)
    }
  return(z[,1])
  }

# ## Testing rdoublepareto
# sample<-rdoublepareto(n=1000,alpha.par=1.8, beta.par=1.5, t.par=10000, c.par=50, m.par=100000)
# breaks_type<-nclass.FD(sample)*10
# hist(sample,freq=FALSE,breaks=breaks_type)
# lines(seq(0,100000,10),ddoublepareto(seq(0,100000,10),alpha.par=1.8, beta.par=1.5, t.par=10000, c.par=50, m.par=100000),col="red")
# ## verify with observed data and ks test

# mll.doublepareto <- function(x, alpha.par, beta.par, t.par, c.par, m.par)
#   {
#   misusloglikelihood.doublepareto<--sum(log((beta.par)/(t.par*(1-((1+(m.par/t.par)^(-alpha.par))/(1+(c.par/t.par)^(-alpha.par)))^(beta.par/alpha.par)))*(((1+(m.par/t.par)^(-alpha.par))^(beta.par/alpha.par))/((1+(x/t.par)^(-alpha.par))^(1+(beta.par/alpha.par))))*((x/t.par)^(-alpha.par-1))))
#   return(misusloglikelihood.doublepareto)
#   }

qdoublepareto <- function(x, alpha.par, beta.par, t.par, c.par, m.par, lower.par, upper.par,  FUN=FALSE)
{
  if(FUN==TRUE)
  {
    inverse<-function (f, lower = lower.par, upper = upper.par) # probably better use interval with extendInt=TRUE
    {
      function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
    }
    alpha.par=alpha.par
    beta.par=beta.par
    t.par=t.par
    c.par=c.par
    m.par=m.par
    quantile.function<-inverse(pdoublepareto(alpha.par=alpha.par, beta.par=beta.par,t.par=t.par,c.par=c.par,m.par=m.par,FUN=TRUE), lower=lower.par, upper=upper.par)
    quantile.function<-Vectorize(quantile.function)
    return(quantile.function)  
  } else
  {
    inverse<-function (f, lower = lower.par, upper = upper.par) # probably better use interval with extendInt=TRUE
    {
      function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
    }
    alpha.par=alpha.par
    beta.par=beta.par
    t.par=t.par
    c.par=c.par
    m.par=m.par
    quantile.function<-inverse(pdoublepareto(alpha.par=alpha.par, beta.par=beta.par,t.par=t.par,c.par=c.par,m.par=m.par,FUN=TRUE), lower=lower.par, upper=upper.par)
    quantile.function<-Vectorize(quantile.function)
    quantile.values<-unlist(quantile.function(x))
    names(quantile.values)<-x
    return(quantile.values)	
  }
}

### Testing qdoublepareto
#qdoublepareto(x=c(0.01,0.05,0.25,0.5,0.75,0.95,0.99),alpha.par=1.3,beta.par=1.1, t.par=0.55, c.par=0.01, m.par=1000, lower.par=0.001, upper.par=1000,FUN=FALSE)
#qdoublepareto(x=c(0.01,0.05,0.25,0.5,0.75,0.95,0.99),alpha.par=1.3,beta.par=1.1, t.par=0.55, c.par=0.01, m.par=1000, lower.par=0.001, upper.par=1000,FUN=TRUE)(c(0.01,0.05,0.25,0.5,0.75,0.95,0.99))



mll.doublepareto <- function(x, alpha.par, beta.par, t.par, c.par, m.par)
  {
  delta.function<-((1+(m.par/t.par)^(-alpha.par))/(1+(c.par/t.par)^(-alpha.par)))^(beta.par/alpha.par)
  eta.function<-(beta.par)/(t.par*(1-delta.function))
  density.doublepareto<-eta.function*(((1+(m.par/t.par)^(-alpha.par))^(beta.par/alpha.par))/((1+(x/t.par)^(-alpha.par))^(1+(beta.par/alpha.par))))*((x/t.par)^(-alpha.par-1))
  misusloglikelihood.doublepareto<--sum(log(density.doublepareto))
  return(misusloglikelihood.doublepareto)
  }


##### Simplified Double Pareto (with c.par=0 and m.par=inf)
ddoublepareto.simplified <- function(x, alpha.par, beta.par, t.par,  FUN=FALSE) # alpha.par > 0, beta.par > 0, 0 <= c.par < t.par <=  m.par <= inf
  {
  if(FUN==TRUE)
    {
    density.function<-function(x)
      {
      alpha.par=alpha.par
      beta.par=beta.par
      t.par=t.par
      fun<-(beta.par*(t.par^alpha.par))/((1+((x/t.par)^(-alpha.par)))^(1+(beta.par/alpha.par))*(x)^(alpha.par+1))
      return(fun)
      }
    return(density.function)
    } else
    {
    d.doublepareto.simplified<-(beta.par*(t.par^alpha.par))/((1+((x/t.par)^(-alpha.par)))^(1+(beta.par/alpha.par))*(x)^(alpha.par+1))
    return(d.doublepareto.simplified)
    }
  }

pdoublepareto.simplified <- function(x, alpha.par, beta.par, t.par, FUN=FALSE) # Integration with integrate
  {
  if(FUN==TRUE)
    {
    cumdensity.function<-function(x)
      {
      density.function<-function(x)
        {
        alpha.par=alpha.par
        beta.par=beta.par
        t.par=t.par
        denfun<-(beta.par*(t.par^alpha.par))/((1+((x/t.par)^(-alpha.par)))^(1+(beta.par/alpha.par))*(x)^(alpha.par+1))
        return(denfun)
        }
      fun<-integrate(f=density.function,lower=0.0001,upper=x)
      return(fun$value)
      }
    cumdensity.function<-Vectorize(cumdensity.function)
    return(cumdensity.function)  
    } else
    {
    cumdensity.function<-function(x)
      {
      density.function<-function(x)
        {
        alpha.par=alpha.par
        beta.par=beta.par
        t.par=t.par
        fun<-(beta.par*(t.par^alpha.par))/((1+((x/t.par)^(-alpha.par)))^(1+(beta.par/alpha.par))*(x)^(alpha.par+1))
        return(fun)
        }
      p.doublepareto.simplified<-integrate(f=density.function,lower=0.0001,upper=x)
      return(p.doublepareto.simplified$value)
      }
    cumdensity.function<-Vectorize(cumdensity.function)
    return(cumdensity.function(x))	
    }
  }

# ## Testing pdoublepareto and ddoublepareto
# c<-pdoublepareto.simplified(alpha.par=1.8, beta.par=1.5, t.par=10000, FUN=TRUE)
# c(seq(50,100000,10))
# pdoublepareto.simplified(seq(50,100000,10),alpha.par=1.8, beta.par=1.5, t.par=10000,FUN=FALSE)
# a<-ddoublepareto.simplified(seq(50,100000,10),alpha.par=1.8, beta.par=1.5, t.par=10000, FUN=TRUE)
# a(seq(50,100000,10))
# 
# values<-seq(50,100000,10)
# den<-ddoublepareto.simplified(values,alpha.par=1.8, beta.par=1.5, t.par=10000)
# plot(values,den,type="l",col="black",log="xy")
# cumden<-pdoublepareto.simplified(values,alpha.par=1.8, beta.par=1.5, t.par=10000)
# windows()
# plot(values,cumden,type="l",col="black")

rdoublepareto.simplified <- function(n, alpha.par, beta.par, t.par) 
  {
  cumprobunif<-sort(runif(n))
  fun <- function (x,u)
    {
    alpha=alpha.par
    beta=beta.par
    t=t.par
    cumfun<-pdoublepareto.simplified(alpha.par=alpha, beta.par=beta, t.par=t, FUN=TRUE)
    return(cumfun(x)-u)
    }
  z <- matrix(NA,nrow=n,ncol=2)
  z[,2]<-cumprobunif
  #colnames(z)<-c("value","cumsumprob")
  for (i in 1:n)
    {
    #print(paste("u=",cumprobunif[i],sep=""))
	if (inherits(try(z[i,1] <- uniroot(fun, lower=uniroot.range[1],upper=uniroot.range[2], tol = 0.0001, u = cumprobunif[i])$root,silent=TRUE),what="try-error"))
		{
		next()
		}	
    #print(paste("value=",z[i,1]))
    #if (z[i,1]<0)
    }
  return(z[,1])
  }

# ## Testing rdoublepareto.simplified
# sample<-rdoublepareto.simplified(n=1000,alpha.par=1.8, beta.par=1.5, t.par=10000)
# breaks_type<-nclass.FD(sample)*10
# hist(sample,freq=FALSE,breaks=breaks_type)
# lines(seq(0,100000,10),ddoublepareto.simplified(seq(0,100000,10),alpha.par=1.8, beta.par=1.5, t.par=10000),col="red")
# ## verify with observed data and ks test

qdoublepareto.simplified <- function(x, alpha.par, beta.par, t.par, lower.par, upper.par, FUN=FALSE)
  {
  if(FUN==TRUE)
    {
    inverse<-function (f, lower = lower.par, upper = upper.par) # probably better use interval with extendInt=TRUE
      {
      function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
      }
    alpha.par=alpha.par
    beta.par=beta.par
    t.par=t.par
    quantile.function<-inverse(pdoublepareto.simplified(alpha.par=alpha.par, beta.par=beta.par,t.par=t.par,FUN=TRUE), lower=lower.par, upper=upper.par)
    quantile.function<-Vectorize(quantile.function)
    return(quantile.function)  
    } else
    {
    inverse<-function (f, lower = lower.par, upper = upper.par) # probably better use interval with extendInt=TRUE
      {
      function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
      }
    alpha.par=alpha.par
    beta.par=beta.par
    t.par=t.par
    quantile.function<-inverse(pdoublepareto.simplified(alpha.par=alpha.par, beta.par=beta.par,t.par=t.par,FUN=TRUE), lower=lower.par, upper=upper.par)
    quantile.function<-Vectorize(quantile.function)
    quantile.values<-unlist(quantile.function(x))
    names(quantile.values)<-x
    return(quantile.values)	
    }
  }
  
### Testing qdoublepareto.simplified
#qdoublepareto.simplified(x=c(0.01,0.05,0.25,0.5,0.75,0.95,0.99),alpha.par=1.3,beta.par=1.1, t.par=0.55, lower.par=0.001, upper.par=1000,FUN=FALSE)
#qdoublepareto.simplified(x=c(0.01,0.05,0.25,0.5,0.75,0.95,0.99),alpha.par=1.3,beta.par=1.1, t.par=0.55, lower.par=0.001, upper.par=1000,FUN=TRUE)(c(0.01,0.05,0.25,0.5,0.75,0.95,0.99))


mll.doublepareto.simplified <- function(x, alpha.par, beta.par, t.par)
  {
  misusloglikelihood.doublepareto.simplified<--sum(log((beta.par*(t.par^alpha.par))/((1+((x/t.par)^(-alpha.par)))^(1+(beta.par/alpha.par))*(x)^(alpha.par+1))))
  return(misusloglikelihood.doublepareto.simplified)
  }


##### Inverse gamma (Correponds to Pearson Type 5 distribution also available in PearsonDS)
dinversegamma <- function(x, alpha.par, eta.par, lambda.par,  FUN=FALSE) # lambda.par > 0, alpha.par > 0  (To compare these parameter to those of Malamud et al., 2004: alpha.par=rho, eta.par=sqrt(-s), lambda.par=sqrt(a))
  {
  if(FUN==TRUE)
    {
    density.function<-function(x)
      {
      alpha.par=alpha.par
      eta.par=eta.par
      lambda.par=lambda.par
      fun<-((lambda.par^(2*alpha.par))/gamma(alpha.par))*((1/(x+(eta.par^2)))^(alpha.par+1))*exp(-(lambda.par^2)/(x+(eta.par^2)))
      return(fun)
      }
    return(density.function)
    } else
    {
    d.inversegamma<-((lambda.par^(2*alpha.par))/gamma(alpha.par))*((1/(x+(eta.par^2)))^(alpha.par+1))*exp(-(lambda.par^2)/(x+(eta.par^2)))
    return(d.inversegamma)
    }
  }


pinversegamma <- function(x, alpha.par, eta.par, lambda.par, FUN=FALSE) # Integration with integrate
  {
  if(FUN==TRUE)
    {
    cumdensity.function<-function(x)
      {
      density.function<-function(x)
        {
        alpha.par=alpha.par
        eta.par=eta.par
        lambda.par=lambda.par
        denfun<-((lambda.par^(2*alpha.par))/gamma(alpha.par))*((1/(x+(eta.par^2)))^(alpha.par+1))*exp(-(lambda.par^2)/(x+(eta.par^2)))
        return(denfun)
        }
      fun<-integrate(f=density.function,lower=0.0001,upper=x)
      return(fun$value)
      }
    cumdensity.function<-Vectorize(cumdensity.function)
    return(cumdensity.function)  
    } else
    {
    cumdensity.function<-function(x)
      {
      density.function<-function(x)
        {
        alpha.par=alpha.par
        eta.par=eta.par
        lambda.par=lambda.par
        fun<-((lambda.par^(2*alpha.par))/gamma(alpha.par))*((1/(x+(eta.par^2)))^(alpha.par+1))*exp(-(lambda.par^2)/(x+(eta.par^2)))
        return(fun)
        }
      p.inversegamma<-integrate(f=density.function,lower=0.0001,upper=x)
      return(p.inversegamma$value)
      }
    cumdensity.function<-Vectorize(cumdensity.function)
    return(cumdensity.function(x))  
    }
  }

## Testing pdoublepareto and ddoublepareto
# c<-pinversegamma(alpha.par=1, eta.par=8, lambda.par=30, FUN=TRUE)
# c(seq(50,100000,10))
# pinversegamma(seq(50,100000,10),alpha.par=1, eta.par=8, lambda.par=30,FUN=FALSE)
# a<-dinversegamma(seq(50,100000,10),alpha.par=1, eta.par=8, lambda.par=30,FUN=TRUE)
# a(seq(50,100000,10))


# values<-seq(0.01,100000,10)
# den<-dinversegamma(values,alpha.par=1, eta.par=8, lambda.par=30)
# plot(values,den,type="l",col="black",log="xy")
# cumden<-pinversegamma(values,alpha.par=1, eta.par=8, lambda.par=30)
# windows()
# plot(values,cumden,type="l",col="black")

rinversegamma <- function(n, alpha.par, eta.par, lambda.par) 
	{
	cumprobunif<-sort(runif(n))
	fun <- function (x,u)
		{
		alpha=alpha.par
		eta=eta.par
		lambda=lambda.par
		cumfun<-pinversegamma(alpha.par=alpha, eta.par=eta, lambda.par=lambda, FUN=TRUE)
		return(cumfun(x)-u)
		}
	z <- matrix(NA,nrow=n,ncol=2)
	z[,2]<-cumprobunif
	#colnames(z)<-c("value","cumsumprob")
	for (i in 1:n)
		{
		#print(paste("u=",cumprobunif[i],sep=""))
		if (inherits(try(z[i,1] <- uniroot(fun, lower=uniroot.range[1],upper=uniroot.range[2], tol=0.00001, u = cumprobunif[i])$root,silent=TRUE),what="try-error"))
			{
			next()
			}	
		#print(paste("value=",z[i,1]))
		#if (z[i,1]<0)
		}
	return(z[,1])
	}


# ## Testing rinversegamma
# sample<-rinversegamma(n=1000,alpha.par=1, eta.par=8, lambda.par=30)
# breaks_type<-nclass.FD(sample)*10
# hist(sample,freq=FALSE,breaks=breaks_type)
# lines(seq(0.2,100000,10),dinversegamma(seq(0.2,100000,10),alpha.par=1, eta.par=8, lambda.par=30),col="red")
# ## verify with observed data and ks test

qinversegamma <- function(x, alpha.par=alpha, eta.par=eta, lambda.par=lambda, lower.par, upper.par, FUN=FALSE)
{
  if(FUN==TRUE)
  {
    inverse<-function (f, lower = lower.par, upper = upper.par) # probably better use interval with extendInt=TRUE
    {
      function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
    }
    alpha.par=alpha.par
    eta.par=eta.par
    lambda.par=lambda.par
    quantile.function<-inverse(pinversegamma(alpha.par=alpha.par, eta.par=eta.par,lambda.par=lambda.par,FUN=TRUE), lower=lower.par, upper=upper.par)
    quantile.function<-Vectorize(quantile.function)
    return(quantile.function)  
  } else
  {
    inverse<-function (f, lower = lower.par, upper = upper.par) # probably better use interval with extendInt=TRUE
    {
      function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
    }
    alpha.par=alpha.par
    eta.par=eta.par
    lambda.par=lambda.par
    quantile.function<-inverse(pinversegamma(alpha.par=alpha.par, eta.par=eta.par,lambda.par=lambda.par,FUN=TRUE), lower=lower.par, upper=upper.par)
    quantile.function<-Vectorize(quantile.function)
    quantile.values<-unlist(quantile.function(x))
    names(quantile.values)<-x
    return(quantile.values)	
  }
}

### Testing qinversegamma
#qinversegamma(x=c(0.01,0.05,0.25,0.5,0.75,0.95,0.99),alpha.par=1.3,eta.par=0.1, lambda.par=1, lower.par=0.001, upper.par=1000,FUN=FALSE)
#qinversegamma(x=c(0.01,0.05,0.25,0.5,0.75,0.95,0.99),alpha.par=1.3,eta.par=0.1, lambda.par=1, lower.par=0.001, upper.par=1000,FUN=TRUE)(c(0.01,0.05,0.25,0.5,0.75,0.95,0.99))



mll.inversegamma <- function(x, alpha.par, eta.par, lambda.par)
  {
  misusloglikelihood.inversegamma<--sum(log(((lambda.par^(2*alpha.par))/gamma(alpha.par))*((1/(x+(eta.par^2)))^(alpha.par+1))*exp(-(lambda.par^2)/(x+(eta.par^2)))))
  return(misusloglikelihood.inversegamma)
  }



### Definition of empty matrix for result storing
hde.results<-matrix(data=NA, nrow=14, ncol=4)
colnames(hde.results)<-c("Estimate","Standard_Error","t_value","Pr(>|t|)")
rownames(hde.results)<-c("DoubleParetoSimplified_alpha","DoubleParetoSimplified_beta","DoubleParetoSimplified_t","DoubleParetoSimplified_Rollover","DoublePareto_alpha","DoublePareto_beta","DoublePareto_t","DoublePareto_c","DoublePareto_m","DoublePareto_Rollover","InverseGamma_alpha","InverseGamma_eta","InverseGamma_lambda","InverseGamma_Rollover")

kde.results<-matrix(data=NA, nrow=14, ncol=4)
colnames(kde.results)<-c("Estimate","Standard_Error","t_value","Pr(>|t|)")
rownames(kde.results)<-c("DoubleParetoSimplified_alpha","DoubleParetoSimplified_beta","DoubleParetoSimplified_t","DoubleParetoSimplified_Rollover","DoublePareto_alpha","DoublePareto_beta","DoublePareto_t","DoublePareto_c","DoublePareto_m","DoublePareto_Rollover","InverseGamma_alpha","InverseGamma_eta","InverseGamma_lambda","InverseGamma_Rollover")

mle.results<-matrix(data=NA, nrow=14, ncol=4)
colnames(mle.results)<-c("Estimate","Standard_Error","t_value","Pr(>|t|)")
rownames(mle.results)<-c("DoubleParetoSimplified_alpha","DoubleParetoSimplified_beta","DoubleParetoSimplified_t","DoubleParetoSimplified_Rollover","DoublePareto_alpha","DoublePareto_beta","DoublePareto_t","DoublePareto_c","DoublePareto_m","DoublePareto_Rollover","InverseGamma_alpha","InverseGamma_eta","InverseGamma_lambda","InverseGamma_Rollover")

# Default scale value, specified_by_user == "NO" in configuration.txt file
range.plot.area<-c(10,1000000)
range.plot.density<-c(0.0000001,0.01)

# Scale values specified by the user, specified_by_user == "YES" in configuration.txt file
if (configuration[28,4]=="YES") (range.plot.area<-as.numeric(configuration[28,2:3]))
if (configuration[29,4]=="YES") (range.plot.density<-as.numeric(configuration[29,2:3]))


for (series in 1:dim(data.series)[2])
   {
   # series<-1
   no.events<-length(na.omit(data.series[,series]))
   print(names(data.series)[series])
   log.data.series<-sort(log(na.omit(data.series[,series]),base=10))

   c.par.value<-min(na.omit(data.series[,series]))
   m.par.value<-max(na.omit(data.series[,series]))   
   
# ----------------------------- ECDF estimation  ----------------------------- #
   print(paste("ECDF -",names(data.series)[series]))
   ecdf_function<-ecdf(data.series[,series])
   
   #log_cycles_range<-c(floor(log(min(na.omit(data.series[,series])),10)),ceiling(log(max(na.omit(data.series[,series])),10)))
   #data.series.synthetic<-c(10^log_cycles_range[1],sort((seq(1,10,0.1)%*%t(10^(log_cycles_range[1]:log_cycles_range[2])))[-1,]))
   #data.series.synthetic<-c(min(na.omit(data.series[,series])),data.series.synthetic[which(data.series.synthetic>min(na.omit(data.series[,series])) & data.series.synthetic<max(na.omit(data.series[,series])))],max(na.omit(data.series[,series])))
   
   pdf(file=paste("ECDF_plot.pdf",sep=""))
   #x11()
   plot(data.series[,series],ecdf_function(data.series[,series]),log="x",type="l",main="ECFD", ylab="ECDF", xlab=xlabel,col="blue",xlim=range.plot.area,ylim=c(0,1))
   points(data.series[,series],ecdf_function(data.series[,series]),pch=21,bg="lightskyblue3",col="navyblue",cex=0.9)
   dev.off()
   
   
# ----------------------------- HDE estimation  ----------------------------- #
   print(paste("HDE -",names(data.series)[series]))

   breaks.Sturges<-nclass.Sturges(log.data.series)
   breaks.scott<-nclass.scott(log.data.series)
   breaks.FD<-nclass.FD(log.data.series)
   #breaks.selected<-breaks.scott
    if(bin_method=="Sturges")
      {breaks.method<-breaks.Sturges}
    if(bin_method=="scott")
      {breaks.method<-breaks.scott}
    if(bin_method=="FD")
      {breaks.method<-breaks.FD}

   breaks.selected<-seq(min(log.data.series),max(log.data.series),diff(range(log.data.series))/breaks.method)
   
   ### Generation of synthetic data series
   ## New method
   #log_cycles_range<-c(-1,1)
   log_cycles_range<-c(floor(log(min(na.omit(data.series[,series])),10)),ceiling(log(max(na.omit(data.series[,series])),10)))
   data.series.synthetic<-c(10^log_cycles_range[1],sort((seq(1,10,0.1)%*%t(10^(log_cycles_range[1]:log_cycles_range[2])))[-1,]))
   data.series.synthetic<-c(min(na.omit(data.series[,series])),data.series.synthetic[which(data.series.synthetic>min(na.omit(data.series[,series])) & data.series.synthetic<max(na.omit(data.series[,series])))],max(na.omit(data.series[,series])))

#    ## Old method
#    if (max(na.omit(data.series[,series]))<100000)
# 	{
# 	data.series.synthetic<-round(seq(min(na.omit(data.series[,series])),max(na.omit(data.series[,series])),by=1),0)
# 	} else
# 	{
# 	#data.series.synthetic<-round(seq(min(na.omit(data.series[,series])),max(na.omit(data.series[,series])),by=10),0)
# 	#data.series.synthetic<-round(c(seq(min(na.omit(data.series[,series])),10000,by=10),seq(10000,max(na.omit(data.series[,series])),by=1000)),0)
# 	data.series.synthetic<-round(c(seq(min(na.omit(data.series[,series])),100000,by=10),seq(100000,1000000,by=1000)),seq(1000000,max(na.omit(data.series[,series])),by=100000)),0)
# 	}
# 
#    #data.series.synthetic<-round(seq(min(na.omit(data.series[,series])),max(na.omit(data.series[,series])),by=1),0)
#    # Generating series for volume for which values could be below 1
#    if (round(log(min(na.omit(data.series[,series])),10))<0)
#    		{
#   		if (floor(log(min(na.omit(data.series[,series])),10))==-1) {data.series.synthetic<-c(seq(0.1,0.9,0.1),data.series.synthetic[-1])}   
#   		if (floor(log(min(na.omit(data.series[,series])),10))==-2) {data.series.synthetic<-c(seq(0.01,0.09,0.01),seq(0.1,0.9,0.1),data.series.synthetic[-1])} 
#   		if (floor(log(min(na.omit(data.series[,series])),10))==-3) {data.series.synthetic<-c(seq(0.001,0.009,0.001),seq(0.01,0.09,0.01),seq(0.1,0.9,0.1),data.series.synthetic[-1])}   
#   		if (floor(log(min(na.omit(data.series[,series])),10))==-4) {data.series.synthetic<-c(seq(0.0001,0.0009,0.0001),seq(0.001,0.009,0.001),seq(0.01,0.09,0.01),seq(0.1,0.9,0.1),data.series.synthetic[-1])}   
#   		}
#    #data.series.synthetic<-round(seq(min(na.omit(data.series[,series])),max(na.omit(data.series[,series])),by=10),0)
   
#   if (Sys.info()[1] == "Linux") {x11()}
#   if (Sys.info()[1] == "Windows") {windows()}
#   if (Sys.info()[1] == "Darwin") {quartz()}
   hde.data.series<-hist(log.data.series,breaks=breaks.selected,plot=FALSE)
#   mtext(paste(names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)

	
	
   pdf(file = paste("Histogram_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
   par(cex.axis=0.9)
   hist(log.data.series,breaks=breaks.selected,labels=TRUE)
   mtext(paste(names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
   dev.off()

   # Exlusion of class with 0 counts
   index.exclusion.zero<-which(hde.data.series$counts>0)
   histogram.diff.data.series<-diff(10^hde.data.series$breaks)
   freq.linear.hde.data.series<-hde.data.series$counts[index.exclusion.zero]/histogram.diff.data.series[index.exclusion.zero]
   prob.linear.data.series<-freq.linear.hde.data.series/no.events
   linear.mids.data.series<-10^hde.data.series$mids[index.exclusion.zero]
   

#   ### Fit HDE
#   dps.hde.formula<-formula(prob.linear.data.series ~ (beta.par*(t.par^alpha.par))/((1+((linear.mids.data.series/t.par)^(-alpha.par)))^(1+(beta.par/alpha.par))*(linear.mids.data.series)^(alpha.par+1)))
#   fit.hde.doublepareto.simplified <- nls(dps.hde.formula,algorithm="port",trace=TRUE,lower=c(alpha.par=alpha.range[1], beta.par=beta.range[1],t.par=t.range[1]),upper=c(alpha.par=alpha.range[2], beta.par=beta.range[2],t.par=t.range[2]),start=list(alpha.par=sum(alpha.range)/2, beta.par=sum(beta.range)/2,t.par=sum(t.range)/2/10),control=list(maxiter=1000))
#   summary(fit.hde.doublepareto.simplified)
#   value.fit.hde.doublepareto.simplified<-predict(fit.hde.doublepareto.simplified)
#
#   if (Sys.info()[1] == "Linux") {x11()}
#	if (Sys.info()[1] == "Windows") {windows()}
#	if (Sys.info()[1] == "Darwin") {quartz()}
#   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""), xlab=paste("A [m]",sep=""),col="blue",xlim=range.plot.area,ylim=range.plot.density)
#   lines(linear.mids.data.series,value.fit.hde.doublepareto.simplified,col="red")

	# Default paramenter value, specified_by_user == "NO" in configuration.txt file
	alpha.range<-c(0.1,5)
	beta.range<-c(0.1,5)
	t.range<-c(min(na.omit(data.series[,series]))*2,max(na.omit(data.series[,series]))/10)
	# Paramenter value specified by the user, specified_by_user == "YES" in configuration.txt file
	if (configuration[1,4]=="YES") (alpha.range<-as.numeric(configuration[1,2:3]))
	if (configuration[2,4]=="YES") (beta.range<-as.numeric(configuration[2,2:3]))
	if (configuration[3,4]=="YES") (t.range<-as.numeric(configuration[3,2:3]))

   ### Fit log HDE better  Double Pareto Simplified
   print(paste("HDE (DPS):",names(data.series)[series]))

   
   log.dps.hde.formula<-formula(log10(prob.linear.data.series) ~ log10((beta.par*(t.par^alpha.par))/((1+((linear.mids.data.series/t.par)^(-alpha.par)))^(1+(beta.par/alpha.par))*(linear.mids.data.series)^(alpha.par+1))))
   fit.hde.log.doublepareto.simplified <- nls(log.dps.hde.formula,algorithm="port",trace=TRUE,lower=c(alpha.par=alpha.range[1], beta.par=beta.range[1],t.par=t.range[1]),upper=c(alpha.par=alpha.range[2], beta.par=beta.range[2],t.par=t.range[2]),start=list(alpha.par=sum(alpha.range)/2, beta.par=sum(beta.range)/2,t.par=sum(t.range)/4),control=list(maxiter=1000))

   summary(fit.hde.log.doublepareto.simplified)

   # Storing fitting results
   hde.results[1:3,1:4]<-coef(summary(fit.hde.log.doublepareto.simplified))
   
   value.fit.hde.log.doublepareto.simplified<-predict(fit.hde.log.doublepareto.simplified,newdata=list(linear.mids.data.series=data.series.synthetic))
  
   # Calculating DPS rollover
   hde.index.rollover.dps<-NULL
   hde.index.rollover.dps<-which(10^value.fit.hde.log.doublepareto.simplified==max(10^value.fit.hde.log.doublepareto.simplified))
   if (length(hde.index.rollover.dps)>1) {hde.index.rollover.dps<-hde.index.rollover.dps[1]}
   if (hde.index.rollover.dps>1)
    {
    hde.results[4,1]<-data.series.synthetic[hde.index.rollover.dps]
    }

    # Calculating pdf and cdf for the original data
    data.series.results$hdedps_pdf<-ddoublepareto.simplified(x=data.series[,series],alpha.par=hde.results[1,1],beta.par=hde.results[2,1],t.par=hde.results[3,1])
    data.series.results$hdedps_cdf<-pdoublepareto.simplified(x=data.series[,series],alpha.par=hde.results[1,1],beta.par=hde.results[2,1],t.par=hde.results[3,1])
    #plot(data.series[,series],data.series.results$hdedps_pdf,log="xy")
    #plot(data.series[,series],data.series.results$hdedps_cdf,log="xy")
  
    ### Kolmogorov-Smirnov test
    library(Matching)
    ks_greater_result<-data.frame(matrix(data=NA,nrow=9,ncol=5))
    colnames(ks_greater_result)<-c("Estimation_Distribution","Bootstrap_samples","KS_D","KS_pvalue","KS_pvalue_boot")
    ks_greater_result[,1]<-c("HDE_DPS","HDE_DP","HDE_IG","KDE_DPS","KDE_DP","KDE_IG","MLE_DPS","MLE_DP","MLE_IG")
    
    
    hde.sample.model.dps<-rdoublepareto.simplified(n=no.events, alpha.par=hde.results[1,1], beta.par=hde.results[2,1], t.par=hde.results[3,1]) 
    ks_model_greater<-NULL
    ks_model_greater<-ks.boot(data.series[,series],hde.sample.model.dps,nboots=ks_boot_samples,alternative="two.sided")
    ks_greater_result[1,2:5]<-c(ks_model_greater$nboots,ks_model_greater$ks$statistic,ks_model_greater$ks$p.value,ks_model_greater$ks.boot.pvalue)
    

#   ####### Controllare se sulla versione 2.10.0 se predict.nls integra le funzioni interval e se.fit in particolar modo interval dvrebbe dare i limiti superiore e inferiore e il fit
#   value.fit.hde.log.doublepareto.simplified<-predict(fit.hde.log.doublepareto.simplified,se.fit=TRUE)

#   if (Sys.info()[1] == "Linux") {x11()}
#	if (Sys.info()[1] == "Windows") {windows()}
#	if (Sys.info()[1] == "Darwin") {quartz()}
#   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""), xlab=xlabel,col="blue",xlim=range.plot.area,ylim=range.plot.density)
#   lines(data.series.synthetic,10^value.fit.hde.log.doublepareto.simplified,col="dark green",lty=1)
#   mtext(paste("HDE Double Pareto Simplified -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
#mtext(paste("Boostrapped KS test -> D: ",round(ks_greater_result[1,3],3),"; p-value: ",round(ks_greater_result[1,5],3),sep=""), side=3, col="navyblue", cex=0.8, line=0.5,padj=2.5)


   pdf(file = paste("HDE_Fit_DPS_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
   par(cex.axis=0.9)
   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""), xlab=xlabel,col="blue",xlim=range.plot.area,ylim=range.plot.density)
   lines(data.series.synthetic,10^value.fit.hde.log.doublepareto.simplified,col="dark green",lty=1)
   mtext(paste("HDE Double Pareto Simplified -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
mtext(paste("Boostrapped KS test -> D: ",round(ks_greater_result[1,3],3),"; p-value: ",round(ks_greater_result[1,5],3),sep=""), side=3, col="navyblue", cex=0.8, line=0.5,padj=2.5)
   dev.off()

   ### Fit log HDE better  Double Pareto

   # Default paramenter value, specified_by_user == "NO" in configuration.txt file
   alpha.range<-c(0.1,5)
   beta.range<-c(0.1,5)
   t.range<-c(min(na.omit(data.series[,series]))*2,max(na.omit(data.series[,series]))/10)
   # Paramenter value specified by the user, specified_by_user == "YES" in configuration.txt file
   if (configuration[4,4]=="YES") (alpha.range<-as.numeric(configuration[1,2:3]))
   if (configuration[5,4]=="YES") (beta.range<-as.numeric(configuration[2,2:3]))
   if (configuration[6,4]=="YES") (t.range<-as.numeric(configuration[3,2:3]))
   
   print(paste("HDE (DP):",names(data.series)[series]))
   
   log.dp.hde.formula<-formula(log10(prob.linear.data.series) ~ log10((beta.par)/(t.par*(1-((1+(m.par.value/t.par)^(-alpha.par))/(1+(c.par.value/t.par)^(-alpha.par)))^(beta.par/alpha.par)))*(((1+(m.par.value/t.par)^(-alpha.par))^(beta.par/alpha.par))/((1+(linear.mids.data.series/t.par)^(-alpha.par))^(1+(beta.par/alpha.par))))*((linear.mids.data.series/t.par)^(-alpha.par-1))))
   fit.hde.log.doublepareto <- nls(log.dp.hde.formula,algorithm="port",trace=TRUE,lower=c(alpha.par=alpha.range[1], beta.par=beta.range[1],t.par=t.range[1]),upper=c(alpha.par=alpha.range[2], beta.par=beta.range[2],t.par=t.range[2]),start=list(alpha.par=sum(alpha.range)/2, beta.par=sum(beta.range)/2,t.par=sum(t.range)/4),control=list(maxiter=1000))
   summary(fit.hde.log.doublepareto)
   value.fit.hde.log.doublepareto<-predict(fit.hde.log.doublepareto,newdata=list(linear.mids.data.series=data.series.synthetic))

   # Storing fitting results
   hde.results[5:7,1:4]<-coef(summary(fit.hde.log.doublepareto))
   hde.results[8:9,1]<-rbind(c.par.value,m.par.value)

   # Calculating DP rollover
   hde.index.rollover.dp<-NULL
   hde.index.rollover.dp<-which(10^value.fit.hde.log.doublepareto==max(10^value.fit.hde.log.doublepareto))
   if (length(hde.index.rollover.dp)>1) {hde.index.rollover.dp<-hde.index.rollover.dp[1]}
   if (hde.index.rollover.dp>1)
    {
    hde.results[10,1]<-data.series.synthetic[hde.index.rollover.dp]
    }

    # Calculating pdf and cdf for the original data
    data.series.results$hdedp_pdf<-ddoublepareto(x=data.series[,series],alpha.par=hde.results[5,1],beta.par=hde.results[6,1],t.par=hde.results[7,1],c.par=hde.results[8,1],m.par=hde.results[9,1])
    data.series.results$hdedp_cdf<-pdoublepareto(x=data.series[,series],alpha.par=hde.results[5,1],beta.par=hde.results[6,1],t.par=hde.results[7,1],c.par=hde.results[8,1],m.par=hde.results[9,1])
    #plot(data.series[,series],data.series.results$hdedp_pdf,log="xy")
    #plot(data.series[,series],data.series.results$hdedp_cdf,log="xy")

    ### Kolmogorov-Smirnov test
    hde.sample.model.dp<-rdoublepareto(n=no.events, alpha.par=hde.results[5,1], beta.par=hde.results[6,1], t.par=hde.results[7,1],c.par=hde.results[8,1],m.par=hde.results[9,1]) 
    ks_model_greater<-NULL
    ks_model_greater<-ks.boot(data.series[,series],hde.sample.model.dp,nboots=ks_boot_samples,alternative="two.sided")
    ks_greater_result[2,2:5]<-c(ks_model_greater$nboots,ks_model_greater$ks$statistic,ks_model_greater$ks$p.value,ks_model_greater$ks.boot.pvalue)
    
#	if (Sys.info()[1] == "Linux") {x11()}
#	if (Sys.info()[1] == "Windows") {windows()}
#	if (Sys.info()[1] == "Darwin") {quartz()}
#   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""), xlab=xlabel,col="blue",xlim=range.plot.area,ylim=range.plot.density)
#   lines(data.series.synthetic,10^value.fit.hde.log.doublepareto,col="red",lty=1)
#   mtext(paste("HDE Double Pareto -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
#mtext(paste("Boostrapped KS test -> D: ",round(ks_greater_result[2,3],3),"; p-value: ",round(ks_greater_result[2,5],3),sep=""), side=3, col="navyblue", cex=0.8, line=0.5,padj=2.5)


   pdf(file = paste("HDE_Fit_DP_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
   par(cex.axis=0.9)
   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""), xlab=xlabel,col="blue",xlim=range.plot.area,ylim=range.plot.density)
   lines(data.series.synthetic,10^value.fit.hde.log.doublepareto,col="red",lty=1)
   mtext(paste("HDE Double Pareto -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
  mtext(paste("Boostrapped KS test -> D: ",round(ks_greater_result[2,3],3),"; p-value: ",round(ks_greater_result[2,5],3),sep=""), side=3, col="navyblue", cex=0.8, line=0.5,padj=2.5)
   dev.off()


   ### Fit log HDE better Inverse gamma
   # Default paramenter value, specified_by_user == "NO" in configuration.txt file
   alpha.range<-c(0.1,5)
   eta.range<-c(1,100)
   lambda.range<-c(1,500)
   # Paramenter value specified by the user, specified_by_user == "YES" in configuration.txt file
   if (configuration[7,4]=="YES") (alpha.range<-as.numeric(configuration[7,2:3]))
   if (configuration[8,4]=="YES") (eta.range<-as.numeric(configuration[8,2:3]))
   if (configuration[9,4]=="YES") (lambda.range<-as.numeric(configuration[9,2:3]))
   
   print(paste("HDE (IG):",names(data.series)[series]))

#   ig.hde.formula<-formula(prob.linear.data.series ~ ((lambda.par^(2*alpha.par))/gamma(alpha.par))*((1/(linear.mids.data.series+(eta.par^2)))^(alpha.par+1))*exp(-(lambda.par^2)/(linear.mids.data.series+(eta.par^2))))
#   fit.hde.inversegamma <- nls(ig.hde.formula,algorithm="port",trace=TRUE,lower=c(alpha.par=alpha.range[1],eta.par=eta.range[1],lambda.par=lambda.range[1]),upper=c(alpha.par=alpha.range[2],eta.par=eta.range[2],lambda.par=lambda.range[2]),start=list(alpha.par=sum(alpha.range)/2,eta.par=sum(eta.range)/2,lambda.par=sum(lambda.range)/2),control=list(maxiter=1000))
#   summary(fit.hde.inversegamma)
#   value.fit.hde.inversegamma<-predict(fit.hde.inversegamma,newdata=list(linear.mids.data.series=data.series.synthetic))
#
#   # Storing fitting results
#   hde.results[11:13,1:4]<-coef(summary(fit.hde.inversegamma))


   log.ig.hde.formula<-formula(log10(prob.linear.data.series) ~ log10(((lambda.par^(2*alpha.par))/gamma(alpha.par))*((1/(linear.mids.data.series+(eta.par^2)))^(alpha.par+1))*exp(-(lambda.par^2)/(linear.mids.data.series+(eta.par^2)))))
   fit.hde.log.inversegamma <- nls(log.ig.hde.formula,algorithm="port",trace=TRUE,lower=c(alpha.par=alpha.range[1],eta.par=eta.range[1],lambda.par=lambda.range[1]),upper=c(alpha.par=alpha.range[2],eta.par=eta.range[2],lambda.par=lambda.range[2]),start=list(alpha.par=sum(alpha.range)/2,eta.par=sum(eta.range)/4,lambda.par=sum(lambda.range)/4),control=list(maxiter=1000))

   summary(fit.hde.log.inversegamma)
   value.fit.hde.log.inversegamma<-predict(fit.hde.log.inversegamma,newdata=list(linear.mids.data.series=data.series.synthetic))

   # Storing fitting results
   hde.results[11:13,1:4]<-coef(summary(fit.hde.log.inversegamma))

   # Calculating IG rollover
   hde.index.rollover.ig<-NULL
   hde.index.rollover.ig<-which(10^value.fit.hde.log.inversegamma==max(10^value.fit.hde.log.inversegamma))
   if (length(hde.index.rollover.ig)>1) {hde.index.rollover.ig<-hde.index.rollover.ig[1]}
   if (hde.index.rollover.ig>1)
    {
    hde.results[14,1]<-data.series.synthetic[hde.index.rollover.ig]
    }

    # Calculating pdf and cdf for the original data
    data.series.results$hdeig_pdf<-dinversegamma(x=data.series[,series],alpha.par=hde.results[11,1],eta.par=hde.results[12,1],lambda.par=hde.results[13,1])
    data.series.results$hdeig_cdf<-pinversegamma(x=data.series[,series],alpha.par=hde.results[11,1],eta.par=hde.results[12,1],lambda.par=hde.results[13,1])
    #plot(data.series[,series],data.series.results$hdeig_pdf,log="xy")
    #plot(data.series[,series],data.series.results$hdeig_cdf,log="xy")

    ### Kolmogorov-Smirnov test
    hde.sample.model.ig<-rinversegamma(n=no.events, alpha.par=hde.results[11,1], eta.par=hde.results[12,1], lambda.par=hde.results[13,1])
    ks_model_greater<-NULL
    ks_model_greater<-ks.boot(data.series[,series],hde.sample.model.ig,nboots=ks_boot_samples,alternative="two.sided")
    ks_greater_result[3,2:5]<-c(ks_model_greater$nboots,ks_model_greater$ks$statistic,ks_model_greater$ks$p.value,ks_model_greater$ks.boot.pvalue)
    
#	if (Sys.info()[1] == "Linux") {x11()}
#	if (Sys.info()[1] == "Windows") {windows()}
#	if (Sys.info()[1] == "Darwin") {quartz()}
#   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""), xlab=xlabel,col="blue",xlim=range.plot.area,ylim=range.plot.density)
#   #lines(data.series.synthetic,10^value.fit.hde.log.inversegamma,col="violet",lty=1)
#   lines(data.series.synthetic,value.fit.hde.inversegamma,col="violet",lty=1)
#   mtext(paste("HDE Inverse Gamma -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
#mtext(paste("Boostrapped KS test -> D: ",round(ks_greater_result[3,3],3),"; p-value: ",round(ks_greater_result[3,5],3),sep=""), side=3, col="navyblue", cex=0.8, line=0.5,padj=2.5)



   pdf(file = paste("HDE_Fit_IG_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
   par(cex.axis=0.9)
   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""), xlab=xlabel,col="blue",xlim=range.plot.area,ylim=range.plot.density)
   lines(data.series.synthetic,10^value.fit.hde.log.inversegamma,col="violet",lty=1)
   #lines(data.series.synthetic,value.fit.hde.inversegamma,col="violet",lty=1)
   mtext(paste("HDE Inverse Gamma -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
  mtext(paste("Boostrapped KS test -> D: ",round(ks_greater_result[3,3],3),"; p-value: ",round(ks_greater_result[3,5],3),sep=""), side=3, col="navyblue", cex=0.8, line=0.5,padj=2.5)
   dev.off()



   # Plot comparison HDE different distribution
	if (Sys.info()[1] == "Linux") {x11()}
	if (Sys.info()[1] == "Windows") {windows()}
	if (Sys.info()[1] == "Darwin") {quartz()}
   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""), xlab=xlabel,col="blue",xlim=range.plot.area,ylim=range.plot.density)
   lines(data.series.synthetic,10^value.fit.hde.log.doublepareto.simplified,col="dark green",lty=1)
   lines(data.series.synthetic,10^value.fit.hde.log.doublepareto,col="red",lty=1)
   lines(data.series.synthetic,10^value.fit.hde.log.inversegamma,col="violet",lty=1)
   #lines(data.series.synthetic,value.fit.hde.inversegamma,col="violet",lty=1)
   mtext(paste("HDE Comparison -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
   legend("bottomleft",legend=c("Histogram raw data","Double Pareto simplified","Double Pareto","Inverse Gamma"),pch=c(1,NA_integer_,NA_integer_,NA_integer_),lty=c(FALSE,1,1,1),lwd=1,col=c("blue","dark green","red","violet"),cex=0.8,box.lty=3,box.col="black")

   pdf(file = paste("HDE_Fit_Comparison_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
   par(cex.axis=0.9)
   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""),xlab=xlabel,col="blue",xlim=range.plot.area,ylim=range.plot.density)
   lines(data.series.synthetic,10^value.fit.hde.log.doublepareto.simplified,col="dark green",lty=1)
   lines(data.series.synthetic,10^value.fit.hde.log.doublepareto,col="red",lty=1)
   lines(data.series.synthetic,10^value.fit.hde.log.inversegamma,col="violet",lty=1)
   #lines(data.series.synthetic,value.fit.hde.inversegamma,col="violet",lty=1)
   mtext(paste("HDE Comparison -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
   legend("bottomleft",legend=c("Histogram raw data","Double Pareto simplified","Double Pareto","Inverse Gamma"),pch=c(1,NA_integer_,NA_integer_,NA_integer_),lty=c(FALSE,1,1,1),lwd=1,col=c("blue","dark green","red","violet"),cex=0.8,box.lty=3,box.col="black")
   dev.off()
  
   
  #### -------------- Uncertainty estimation -------------- ####


	### Defining useful funcion and data usefull for the uncertainty evaluation
	#prob_series<-c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99) # They have to be odd and symmetric
	prob_series<-c(0.01,0.05,0.25,0.5,0.75,0.95,0.99) # They have to be odd and symmetric
	#prob_series<-c(0.05,0.25,0.5,0.75,0.95) # They have to be odd and symmetric
	color_vector_pol<-gray(level=1-(1:floor((length(prob_series))/2))/floor(length(prob_series)*1))
	sciNotation <- function(x, digits = 1)
		{ 
		if (length(x) > 1) {return(append(sciNotation(x[1]), sciNotation(x[-1])))} 
		if (!x) return(0) 
		exponent <- floor(log10(x)) 
		base <- round(x / 10^exponent, digits) 
		#as.expression(substitute(base %*% 10^exponent,list(base = base, exponent = exponent))) 
		as.expression(substitute(10^exponent,list(base = base, exponent = exponent))) 
		} 
	quantile_fun<-function (x) quantile(x,prob=prob_series,na.rm=TRUE)
	
	#min.value.x<-min(data.series.synthetic)
	#max.value.x<-max(data.series.synthetic)
	min.value.x<-range.plot.area[1]
	max.value.x<-range.plot.area[2]
	min.value.at.x<-signif(log10(min.value.x),digits=1)
	max.value.at.x<-signif(ceiling(log10(max.value.x)),digits=1)
	at.x<-10^(min.value.at.x:max.value.at.x)
	at.label.x<-rep(10,length(at.x))^log10(at.x)
	at.x.tick<-sort(as.numeric((1:9)%*%t(at.label.x)))
	
	#min.value.y<-min(prob.linear.data.series)
	#max.value.y<-max(prob.linear.data.series)
	min.value.y<-range.plot.density[1]
	max.value.y<-range.plot.density[2]
	min.value.at.y<-signif(log10(min.value.y),digits=1)
	max.value.at.y<-signif(ceiling(log10(max.value.y)),digits=1)
	at.y<-10^(min.value.at.y:max.value.at.y)
	at.label.y<-rep(10,length(at.y))^log10(at.y)
	at.y.tick<-sort(as.numeric((1:9)%*%t(at.label.y)))
	

	# estimating uncentainty
	
  ## New method
  #log_cycles_range<-c(-1,1)
  log_cycles_range<-c(floor(log(min(na.omit(data.series[,series])),10)),ceiling(log(max(na.omit(data.series[,series])),10)))
  data.series.synthetic.sample<-c(10^log_cycles_range[1],sort((seq(1,10,0.1)%*%t(10^(log_cycles_range[1]:log_cycles_range[2])))[-1,]))
  data.series.synthetic.sample<-c(min(na.omit(data.series[,series])),data.series.synthetic.sample[which(data.series.synthetic.sample>min(na.omit(data.series[,series])) & data.series.synthetic.sample<max(na.omit(data.series[,series])))],max(na.omit(data.series[,series])))

  
#   ## Old method
#   #data.series.synthetic.sample<-round(seq(min(na.omit(data.series[,series])),max(na.omit(data.series[,series])),by=10),0) # Based on the observed data
#   #data.series.synthetic.sample<-round(c(seq(min(na.omit(data.series[,series])),10000,by=10),seq(10000,max(na.omit(data.series[,series])),by=1000)),0)
# 
#   # Generation of synthetic data series
#   if (max(na.omit(data.series[,series]))<100000)
#     {
#     data.series.synthetic.sample<-round(seq(min(na.omit(data.series[,series])),max(na.omit(data.series[,series])),by=1),0)
#     } else
#     {
#     #data.series.synthetic<-round(seq(min(na.omit(data.series[,series])),max(na.omit(data.series[,series])),by=10),0)
#     #data.series.synthetic<-round(c(seq(min(na.omit(data.series[,series])),10000,by=10),seq(10000,max(na.omit(data.series[,series])),by=1000)),0)
#     data.series.synthetic.sample<-round(c(seq(min(na.omit(data.series[,series])),100000,by=10),seq(100000,max(na.omit(data.series[,series])),by=100)),0)
#     if (max(na.omit(data.series[,series]))<1000000)
#       {
#       data.series.synthetic.sample<-round(c(seq(min(na.omit(data.series[,series])),100000,by=10),seq(100000,max(na.omit(data.series[,series])),by=100)),0)
#       } else
#       {
#       data.series.synthetic.sample<-round(c(seq(min(na.omit(data.series[,series])),100000,by=10),seq(100000,1000000,by=100),seq(1000000,max(na.omit(data.series[,series])),by=1000)),0)
#       }
#     }
# 
# 
#   if (round(log(min(na.omit(data.series[,series])),10))<0)
#     {
#     if (floor(log(min(na.omit(data.series[,series])),10))==-1) {data.series.synthetic.sample<-c(seq(0.1,0.9,0.1),data.series.synthetic.sample[-1])}   
#     if (floor(log(min(na.omit(data.series[,series])),10))==-2) {data.series.synthetic.sample<-c(seq(0.01,0.09,0.01),seq(0.1,0.9,0.1),data.series.synthetic.sample[-1])} 
#     if (floor(log(min(na.omit(data.series[,series])),10))==-3) {data.series.synthetic.sample<-c(seq(0.001,0.009,0.001),seq(0.01,0.09,0.01),seq(0.1,0.9,0.1),data.series.synthetic.sample[-1])}   
#     if (floor(log(min(na.omit(data.series[,series])),10))==-4) {data.series.synthetic.sample<-c(seq(0.0001,0.0009,0.0001),seq(0.001,0.009,0.001),seq(0.01,0.09,0.01),seq(0.1,0.9,0.1),data.series.synthetic.sample[-1])}   
#     }



	# Number of samples
	samples_number<-100
	
	### Creating object for storing sampling results
	hde_dps_sampling_result<-matrix(NA,nrow=length(data.series.synthetic.sample),ncol=samples_number)
	hde_dps_sampling_result_coef<-matrix(NA,nrow=samples_number,ncol=4)
	colnames(hde_dps_sampling_result_coef)<-c("alpha","beta","t","rollov")
	
	hde_dp_sampling_result<-matrix(NA,nrow=length(data.series.synthetic.sample),ncol=samples_number)
	hde_dp_sampling_result_coef<-matrix(NA,nrow=samples_number,ncol=6)
	colnames(hde_dp_sampling_result_coef)<-c("alpha","beta","t","c","m","rollov")
	
	hde_ig_sampling_result<-matrix(NA,nrow=length(data.series.synthetic.sample),ncol=samples_number)
	hde_ig_sampling_result_coef<-matrix(NA,nrow=samples_number,ncol=4)
	colnames(hde_ig_sampling_result_coef)<-c("alpha","eta","lambda","rollov")
			
	for (countsamples in 1:samples_number)
		{
		#countsamples<-1
		c.par.value.sample<-c.par.value
		m.par.value.sample<-m.par.value
   	print(paste("HDE -",names(data.series)[series],"-> Sample:",countsamples))
		data.series.sample<-sample(data.series[,series],replace=TRUE)
		no.events.sample<-length(data.series.sample)
		log.data.series.sample<-log(data.series.sample,10)
		hde.data.series.sample<-hist(log.data.series.sample,breaks=breaks.selected,plot=FALSE)
		# Exlusion of class with 0 counts
		index.exclusion.zero.sample<-which(hde.data.series.sample$counts>0)
		histogram.diff.data.series<-diff(10^hde.data.series.sample$breaks)
		freq.linear.hde.data.series.sample<-hde.data.series.sample$counts[index.exclusion.zero.sample]/histogram.diff.data.series[index.exclusion.zero.sample]
		prob.linear.data.series.sample<-freq.linear.hde.data.series.sample/no.events.sample
		linear.mids.data.series.sample<-10^hde.data.series.sample$mids[index.exclusion.zero.sample]
		
		### ----------------- Fit Double Pareto Simplified (DPS) ----------------- ###
		# Default paramenter value, specified_by_user == "NO" in configuration.txt file
		alpha.range<-c(0.1,5)
		beta.range<-c(0.1,5)
		t.range<-c(min(na.omit(data.series[,series]))*2,max(na.omit(data.series[,series]))/10)
		# Paramenter value specified by the user, specified_by_user == "YES" in configuration.txt file
		if (configuration[1,4]=="YES") (alpha.range<-as.numeric(configuration[1,2:3]))
		if (configuration[2,4]=="YES") (beta.range<-as.numeric(configuration[2,2:3]))
		if (configuration[3,4]=="YES") (t.range<-as.numeric(configuration[3,2:3]))
		
   		print(paste("HDE (DPS):",names(data.series)[series],"-> Sample:",countsamples))
		log.dps.hde.formula.sample<-formula(log10(prob.linear.data.series.sample) ~ log10((beta.par*(t.par^alpha.par))/((1+((linear.mids.data.series.sample/t.par)^(-alpha.par)))^(1+(beta.par/alpha.par))*(linear.mids.data.series.sample)^(alpha.par+1))))
		if (inherits(try(fit.hde.log.doublepareto.simplified.sample<-nls(log.dps.hde.formula.sample,algorithm="port",trace=FALSE,lower=c(alpha.par=alpha.range[1], beta.par=beta.range[1],t.par=t.range[1]),upper=c(alpha.par=alpha.range[2], beta.par=beta.range[2],t.par=t.range[2]),start=list(alpha.par=sum(alpha.range)/2, beta.par=sum(beta.range)/2,t.par=sum(t.range)/4),control=list(maxiter=1000)),silent=TRUE),what="try-error"))
			{
			next()
			}	
		# Storing fitting results
		hde_dps_sampling_result_coef[countsamples,1:3]<-coef(summary(fit.hde.log.doublepareto.simplified.sample))[,1]
		hde_dps_sampling_result[,countsamples]<-predict(fit.hde.log.doublepareto.simplified.sample,newdata=list(linear.mids.data.series.sample=data.series.synthetic.sample))
		# Calculating DPS rollover
		hde.index.rollover.dps.sample<-NULL
		hde.index.rollover.dps.sample<-which(10^hde_dps_sampling_result[,countsamples]==max(10^hde_dps_sampling_result[,countsamples]))
		if (length(hde.index.rollover.dps.sample)>1) {hde.index.rollover.dps.sample<-hde.index.rollover.dps.sample[1]}
		if (hde.index.rollover.dps.sample>1)
			{
			hde_dps_sampling_result_coef[countsamples,4]<-data.series.synthetic.sample[hde.index.rollover.dps.sample]
			}
		
		### ----------------- Fit Double Pareto (DP) ----------------- ###
		# Default paramenter value, specified_by_user == "NO" in configuration.txt file
		alpha.range<-c(0.1,5)
		beta.range<-c(0.1,5)
		t.range<-c(min(na.omit(data.series[,series]))*2,max(na.omit(data.series[,series]))/10)
		# Paramenter value specified by the user, specified_by_user == "YES" in configuration.txt file
		if (configuration[4,4]=="YES") (alpha.range<-as.numeric(configuration[1,2:3]))
		if (configuration[5,4]=="YES") (beta.range<-as.numeric(configuration[2,2:3]))
		if (configuration[6,4]=="YES") (t.range<-as.numeric(configuration[3,2:3]))
		
		print(paste("HDE (DP):",names(data.series)[series],"-> Sample:",countsamples))
		log.dp.hde.formula.sample<-formula(log10(prob.linear.data.series.sample) ~ log10((beta.par)/(t.par*(1-((1+(m.par.value.sample/t.par)^(-alpha.par))/(1+(c.par.value.sample/t.par)^(-alpha.par)))^(beta.par/alpha.par)))*(((1+(m.par.value.sample/t.par)^(-alpha.par))^(beta.par/alpha.par))/((1+(linear.mids.data.series.sample/t.par)^(-alpha.par))^(1+(beta.par/alpha.par))))*((linear.mids.data.series.sample/t.par)^(-alpha.par-1))))
		if (inherits(try(fit.hde.log.doublepareto.sample<-nls(log.dp.hde.formula.sample,algorithm="port",trace=FALSE,lower=c(alpha.par=alpha.range[1], beta.par=beta.range[1],t.par=t.range[1]),upper=c(alpha.par=alpha.range[2], beta.par=beta.range[2],t.par=t.range[2]),start=list(alpha.par=sum(alpha.range)/2, beta.par=sum(beta.range)/2,t.par=sum(t.range)/4),control=list(maxiter=1000)),silent=TRUE),what="try-error"))
			{
			next()
			}	
		# Storing fitting results
		hde_dp_sampling_result_coef[countsamples,1:3]<-coef(summary(fit.hde.log.doublepareto.sample))[,1]
		hde_dp_sampling_result_coef[countsamples,4:5]<-c(c.par.value.sample,m.par.value.sample)
		hde_dp_sampling_result[,countsamples]<-predict(fit.hde.log.doublepareto.sample,newdata=list(linear.mids.data.series.sample=data.series.synthetic.sample))
		# Calculating DP rollover
		hde.index.rollover.dp.sample<-NULL
		hde.index.rollover.dp.sample<-which(10^hde_dp_sampling_result[,countsamples]==max(10^hde_dp_sampling_result[,countsamples]))
		if (length(hde.index.rollover.dp.sample)>1) {hde.index.rollover.dp.sample<-hde.index.rollover.dp.sample[1]}
		if (hde.index.rollover.dp.sample>1)
			{
			hde_dp_sampling_result_coef[countsamples,6]<-data.series.synthetic.sample[hde.index.rollover.dp.sample]
			}
   
		### ----------------- Fit Inverse Gamma (IG) ----------------- ###
		# Default paramenter value, specified_by_user == "NO" in configuration.txt file
		alpha.range<-c(0.1,5)
		eta.range<-c(1,100)
		lambda.range<-c(1,500)
		# Paramenter value specified by the user, specified_by_user == "YES" in configuration.txt file
		if (configuration[7,4]=="YES") (alpha.range<-as.numeric(configuration[7,2:3]))
		if (configuration[8,4]=="YES") (eta.range<-as.numeric(configuration[8,2:3]))
		if (configuration[9,4]=="YES") (lambda.range<-as.numeric(configuration[9,2:3]))
		
		print(paste("HDE (IG):",names(data.series)[series],"-> Sample:",countsamples))
		log.ig.hde.formula.sample<-formula(log10(prob.linear.data.series.sample) ~ log10(((lambda.par^(2*alpha.par))/gamma(alpha.par))*((1/(linear.mids.data.series.sample+(eta.par^2)))^(alpha.par+1))*exp(-(lambda.par^2)/(linear.mids.data.series.sample+(eta.par^2)))))
		if (inherits(try(fit.hde.log.inversegamma.sample<-nls(log.ig.hde.formula.sample,algorithm="port",trace=FALSE,lower=c(alpha.par=alpha.range[1],eta.par=eta.range[1],lambda.par=lambda.range[1]),upper=c(alpha.par=alpha.range[2],eta.par=eta.range[2],lambda.par=lambda.range[2]),start=list(alpha.par=sum(alpha.range)/2,eta.par=sum(eta.range)/4,lambda.par=sum(lambda.range)/4),control=list(maxiter=1000)),silent=TRUE),what="try-error"))
			{
			next()
			}	
		# Storing fitting results
		hde_ig_sampling_result_coef[countsamples,1:3]<-coef(summary(fit.hde.log.inversegamma.sample))[,1]
		hde_ig_sampling_result[,countsamples]<-predict(fit.hde.log.inversegamma.sample,newdata=list(linear.mids.data.series.sample=data.series.synthetic.sample))
		# Calculating IG rollover
		hde.index.rollover.ig.sample<-NULL
		hde.index.rollover.ig.sample<-which(10^hde_ig_sampling_result[,countsamples]==max(10^hde_ig_sampling_result[,countsamples]))
		if (length(hde.index.rollover.ig.sample)>1) {hde.index.rollover.ig.sample<-hde.index.rollover.ig.sample[1]}
		if (hde.index.rollover.ig.sample>1)
			{
			hde_ig_sampling_result_coef[countsamples,4]<-data.series.synthetic.sample[hde.index.rollover.ig.sample]
			}
		print(paste("HDE -",names(data.series)[series],"-> Sample:",countsamples,"Progress:",round(countsamples/samples_number*100,1),"%"))
		}

		quantiles_hde_dps_sampling_result_coef<-apply(hde_dps_sampling_result_coef,MARGIN=2,FUN=quantile_fun)
		quantiles_hde_dp_sampling_result_coef<-apply(hde_dp_sampling_result_coef,MARGIN=2,FUN=quantile_fun)
		quantiles_hde_ig_sampling_result_coef<-apply(hde_ig_sampling_result_coef,MARGIN=2,FUN=quantile_fun)
		
		write.table(data.frame(percentile=rownames(quantiles_hde_dps_sampling_result_coef),round(quantiles_hde_dps_sampling_result_coef,2)),paste("HDE_Fit_DPS_parameter_uncertainty_",names(data.series)[series],".txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t")
		write.table(data.frame(percentile=rownames(quantiles_hde_dp_sampling_result_coef),round(quantiles_hde_dp_sampling_result_coef,2)),paste("HDE_Fit_DP_parameter_uncertainty_",names(data.series)[series],".txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t")
		write.table(data.frame(percentile=rownames(quantiles_hde_ig_sampling_result_coef),round(quantiles_hde_ig_sampling_result_coef,2)),paste("HDE_Fit_IG_parameter_uncertainty_",names(data.series)[series],".txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t")
		
		quantiles_hde_dps_sampling_result<-apply(hde_dps_sampling_result,MARGIN=1,FUN=quantile_fun)
		quantiles_hde_dp_sampling_result<-apply(hde_dp_sampling_result,MARGIN=1,FUN=quantile_fun)
		quantiles_hde_ig_sampling_result<-apply(hde_ig_sampling_result,MARGIN=1,FUN=quantile_fun)

		
		### Plot HDE estimation with uncertainty
		pdf(file = paste("HDE_Fit_DPS_uncertainty_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
		par(cex.axis=0.9)
		#if (Sys.info()[1] == "Linux") {x11()}
		#if (Sys.info()[1] == "Windows") {windows()}
		#if (Sys.info()[1] == "Darwin") {quartz()}
		#xlabel<-expression(A (m^{2}))
		ylabel<-paste("Probability density",sep="")
		plot(NA, NA, log="xy", main="Probability densities", ylab=ylabel, xlab=xlabel,xlim=range.plot.area,ylim=range.plot.density,xaxt="n",yaxt="n")
		axis(side=1,at=at.x,label=sciNotation(at.x,1))
		axis(side=2,at=at.y,label=sciNotation(at.y,1))
		axis(1, at.x.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		axis(2, at.y.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		axis(3, at.x.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		axis(4, at.y.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		for (count_q in 1:floor((length(prob_series))/2))
			{
			#count_q<-1	
			polygon(x=c(sort(data.series.synthetic.sample),rev(sort(data.series.synthetic.sample))),y=10^c(quantiles_hde_dps_sampling_result[length(prob_series)+1-count_q,],rev(quantiles_hde_dps_sampling_result[count_q,])),border=NA,col=color_vector_pol[count_q])
			}
		lines(data.series.synthetic,10^value.fit.hde.log.doublepareto.simplified,col="darkred",lty=1,lwd=2)
		lines(data.series.synthetic.sample,10^quantiles_hde_dps_sampling_result["50%",],col="black",lty="dashed",lwd=1)
		points(linear.mids.data.series, prob.linear.data.series,pch=20,cex=1,col="navyblue")
		mtext(paste("HDE Double Pareto Simplified -> ",names(data.series)[series],sep=""), side=3, col="darkred", cex=0.8, line=0.5)
    mtext(paste("Boostrapped KS test -> D: ",round(ks_greater_result[1,3],3),"; p-value: ",round(ks_greater_result[1,5],3),sep=""), side=3, col="navyblue", cex=0.8, line=0.5,padj=2.5)
  legend("bottomleft",c("Histogram empirical data","Distribution fit","Median values",paste("Quantiles: ",prob_series[1:floor((length(prob_series))/2)],"-",prob_series[(length(prob_series)+1)-1:floor((length(prob_series))/2)],sep="")),pt.cex=1,pch=c(20,NA,NA,rep(NA,floor((length(prob_series))/2))),lty=c(NA,1,2,rep(1,floor((length(prob_series))/2))),lwd=c(NA,2,1,rep(6,floor((length(prob_series))/2))),col=c("navyblue","darkred","black",color_vector_pol),bty="o",bg="transparent",box.col="transparent",cex=0.9)
		dev.off()
		
		
		pdf(file = paste("HDE_Fit_DP_uncertainty_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
		par(cex.axis=0.9)
		#if (Sys.info()[1] == "Linux") {x11()}
		#if (Sys.info()[1] == "Windows") {windows()}
		#if (Sys.info()[1] == "Darwin") {quartz()}
		#xlabel<-expression(A (m^{2}))
		ylabel<-paste("Probability density",sep="")
		plot(NA, NA, log="xy", main="Probability densities", ylab=ylabel, xlab=xlabel,xlim=range.plot.area,ylim=range.plot.density,xaxt="n",yaxt="n")
		axis(side=1,at=at.x,label=sciNotation(at.x,1))
		axis(side=2,at=at.y,label=sciNotation(at.y,1))
		axis(1, at.x.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		axis(2, at.y.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		axis(3, at.x.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		axis(4, at.y.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		for (count_q in 1:floor((length(prob_series))/2))
			{
			#count_q<-1	
			polygon(x=c(sort(data.series.synthetic.sample),rev(sort(data.series.synthetic.sample))),y=10^c(quantiles_hde_dp_sampling_result[length(prob_series)+1-count_q,],rev(quantiles_hde_dp_sampling_result[count_q,])),border=NA,col=color_vector_pol[count_q])
			}
		lines(data.series.synthetic,10^value.fit.hde.log.doublepareto,col="darkred",lty=1,lwd=2)
		lines(data.series.synthetic.sample,10^quantiles_hde_dp_sampling_result["50%",],col="black",lty="dashed",lwd=1)
		points(linear.mids.data.series, prob.linear.data.series,pch=20,cex=1,col="navyblue")
		mtext(paste("HDE Double Pareto -> ",names(data.series)[series],sep=""), side=3, col="darkred", cex=0.8, line=0.5)
    mtext(paste("Boostrapped KS test -> D: ",round(ks_greater_result[2,3],3),"; p-value: ",round(ks_greater_result[2,5],3),sep=""), side=3, col="navyblue", cex=0.8, line=0.5,padj=2.5)
		legend("bottomleft",c("Histogram empirical data","Distribution fit","Median values",paste("Quantiles: ",prob_series[1:floor((length(prob_series))/2)],"-",prob_series[(length(prob_series)+1)-1:floor((length(prob_series))/2)],sep="")),pt.cex=1,pch=c(20,NA,NA,rep(NA,floor((length(prob_series))/2))),lty=c(NA,1,2,rep(1,floor((length(prob_series))/2))),lwd=c(NA,2,1,rep(6,floor((length(prob_series))/2))),col=c("navyblue","darkred","black",color_vector_pol),bty="o",bg="transparent",box.col="transparent",cex=0.9)
		dev.off()
		
		pdf(file = paste("HDE_Fit_IG_uncertainty_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
		par(cex.axis=0.9)
		#if (Sys.info()[1] == "Linux") {x11()}
		#if (Sys.info()[1] == "Windows") {windows()}
		#if (Sys.info()[1] == "Darwin") {quartz()}
		#xlabel<-expression(A (m^{2}))
		ylabel<-paste("Probability density",sep="")
		plot(NA, NA, log="xy", main="Probability densities", ylab=ylabel, xlab=xlabel,xlim=range.plot.area,ylim=range.plot.density,xaxt="n",yaxt="n")
		axis(side=1,at=at.x,label=sciNotation(at.x,1))
		axis(side=2,at=at.y,label=sciNotation(at.y,1))
		axis(1, at.x.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		axis(2, at.y.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		axis(3, at.x.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		axis(4, at.y.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		for (count_q in 1:floor((length(prob_series))/2))
			{
			#count_q<-1	
			polygon(x=c(sort(data.series.synthetic.sample),rev(sort(data.series.synthetic.sample))),y=10^c(quantiles_hde_ig_sampling_result[length(prob_series)+1-count_q,],rev(quantiles_hde_ig_sampling_result[count_q,])),border=NA,col=color_vector_pol[count_q])
			}
		lines(data.series.synthetic,10^value.fit.hde.log.inversegamma,col="darkred",lty=1,lwd=2)
		lines(data.series.synthetic.sample,10^quantiles_hde_ig_sampling_result["50%",],col="black",lty="dashed",lwd=1)
		points(linear.mids.data.series, prob.linear.data.series,pch=20,cex=1,col="navyblue")
		mtext(paste("HDE Inverse Gamma -> ",names(data.series)[series],sep=""), side=3, col="darkred", cex=0.8, line=0.5)
    mtext(paste("Boostrapped KS test -> D: ",round(ks_greater_result[3,3],3),"; p-value: ",round(ks_greater_result[3,5],3),sep=""), side=3, col="navyblue", cex=0.8, line=0.5,padj=2.5)
		legend("bottomleft",c("Histogram empirical data","Distribution fit","Median values",paste("Quantiles: ",prob_series[1:floor((length(prob_series))/2)],"-",prob_series[(length(prob_series)+1)-1:floor((length(prob_series))/2)],sep="")),pt.cex=1,pch=c(20,NA,NA,rep(NA,floor((length(prob_series))/2))),lty=c(NA,1,2,rep(1,floor((length(prob_series))/2))),lwd=c(NA,2,1,rep(6,floor((length(prob_series))/2))),col=c("navyblue","darkred","black",color_vector_pol),bty="o",bg="transparent",box.col="transparent",cex=0.9)
		dev.off()
	
# ----------------------------- KDE estimation  ----------------------------- #
   print(paste("KDE -",names(data.series)[series]))

   num.kde<-1000
   bandwidth.data.series<-bw.nrd0(log.data.series)
   #bandwidth.data.series<-bw.ucv(log.data.series)
   #bandwidth.data.series<-bw.nrd(log.data.series)
   #bandwidth.data.series<-bandwidth.data.series*1.5

   kde.data.series<-density(log.data.series,bw=bandwidth.data.series,from=min(log.data.series),to=max(log.data.series),cut=range(log.data.series),adjust=1, kernel=c("g"), n=num.kde, give.Rkern=FALSE)
   range(log.data.series)
   range(kde.data.series$x)
   kernel.diff.data.series<-diff(10^kde.data.series$x)
   kde.data.series.linear<-kde.data.series$y[1:(length(kde.data.series$y)-1)]*((kde.data.series$x[2:length(kde.data.series$x)])-(kde.data.series$x[1:(length(kde.data.series$x)-1)]))/kernel.diff.data.series

   freq.linear.data.series<-kde.data.series$y[1:(length(kde.data.series$y)-1)]*no.events*((kde.data.series$x[2:length(kde.data.series$x)])-(kde.data.series$x[1:(length(kde.data.series$x)-1)]))/kernel.diff.data.series
   linear.x.data.series<-10^(kde.data.series$x[1:(length(kde.data.series$x)-1)])

   kde.data.series.linear.log<-log(kde.data.series.linear,base=10)
   freq.linear.log.kde.data.series<-log(freq.linear.data.series,base=10)
   linear.x.log.kde.data.series<-log(linear.x.data.series,base=10)

#   if (Sys.info()[1] == "Linux") {x11()}
#   if (Sys.info()[1] == "Windows") {windows()}
#   if (Sys.info()[1] == "Darwin") {quartz()}
#   plot(kde.data.series,type="l",main=paste("Kernel density of log(A)", sep=""), xlab=paste("Log(A) [m]",sep=""), ylab="Kernel density")
#   mtext(paste(names(data.series)[series]), side=3, col="red", cex=1, line=0.5)
#   rug(log.data.series)

   pdf(file = paste("KDE_Density_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
   par(cex.axis=0.9)
   plot(kde.data.series,type="l",main=paste("Kernel density of log(A)", sep=""), xlab=paste("Log(A) [m]",sep=""), ylab="Kernel density")
   mtext(paste(names(data.series)[series]), side=3, col="red", cex=1, line=0.5)
   rug(log.data.series)
   dev.off()

#   if (Sys.info()[1] == "Linux") {x11()}
#   if (Sys.info()[1] == "Windows") {windows()}
#   if (Sys.info()[1] == "Darwin") {quartz()}
#   plot(kde.data.series$x,kde.data.series$y*no.events,type="l",main=paste("Frequency densities of log(A)", sep=""), xlab=paste("Log (A) [m]",sep=""),ylab="Frequency density")
#   mtext(paste(names(data.series)[series]), side=3, col="red", cex=1, line=0.5)
#   rug(log.data.series)

   pdf(file = paste("KDE_FrequencyDensity_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
   par(cex.axis=0.9)
   plot(kde.data.series$x,kde.data.series$y*no.events,type="l",main=paste("Frequency densities of log(A)", sep=""), xlab=paste("Log (A) [m]",sep=""),ylab="Frequency density")
   mtext(paste(names(data.series)[series]), side=3, col="red", cex=1, line=0.5)
   rug(log.data.series)
   dev.off()

   
   ### Fit log KDE better  Double Pareto Simplified

   # Default paramenter value, specified_by_user == "NO" in configuration.txt file
   alpha.range<-c(0.1,5)
   beta.range<-c(0.1,5)
   t.range<-c(min(na.omit(data.series.sample))*2,max(na.omit(data.series.sample))/10)
   # Paramenter value specified by the user, specified_by_user == "YES" in configuration.txt file
   if (configuration[10,4]=="YES") (alpha.range<-as.numeric(configuration[10,2:3]))
   if (configuration[11,4]=="YES") (beta.range<-as.numeric(configuration[11,2:3]))
   if (configuration[12,4]=="YES") (t.range<-as.numeric(configuration[12,2:3]))
   
   
   print(paste("KDE (DPS):",names(data.series)[series]))

   log.dps.kde.formula<-formula(log10(kde.data.series.linear) ~ log10((beta.par*(t.par^alpha.par))/((1+((linear.x.data.series/t.par)^(-alpha.par)))^(1+(beta.par/alpha.par))*(linear.x.data.series)^(alpha.par+1))))
   fit.kde.log.doublepareto.simplified <- nls(log.dps.kde.formula,algorithm="port",trace=TRUE,lower=c(alpha.par=alpha.range[1], beta.par=beta.range[1],t.par=t.range[1]),upper=c(alpha.par=alpha.range[2], beta.par=beta.range[2],t.par=t.range[2]),start=list(alpha.par=sum(alpha.range)/2, beta.par=sum(beta.range)/2,t.par=sum(t.range)/4),control=list(maxiter=1000))
   summary(fit.kde.log.doublepareto.simplified)

   #value.fit.kde.log.doublepareto.simplified<-predict(fit.kde.log.doublepareto.simplified)
   value.fit.kde.log.doublepareto.simplified<-predict(fit.kde.log.doublepareto.simplified,newdata=list(linear.x.data.series=data.series.synthetic))


   # Storing fitting results
   kde.results[1:3,1:4]<-coef(summary(fit.kde.log.doublepareto.simplified))

   # Calculating DPS rollover
   kde.index.rollover.dps<-NULL
   kde.index.rollover.dps<-which(10^value.fit.kde.log.doublepareto.simplified==max(10^value.fit.kde.log.doublepareto.simplified))
   if (length(kde.index.rollover.dps)>1) {kde.index.rollover.dps<-kde.index.rollover.dps[1]}
   if (kde.index.rollover.dps>1)
    {
    kde.results[4,1]<-data.series.synthetic[kde.index.rollover.dps]
    }

    # Calculating pdf and cdf for the original data
    data.series.results$kdedps_pdf<-ddoublepareto.simplified(x=data.series[,series],alpha.par=kde.results[1,1],beta.par=kde.results[2,1],t.par=kde.results[3,1])
    data.series.results$kdedps_cdf<-pdoublepareto.simplified(x=data.series[,series],alpha.par=kde.results[1,1],beta.par=kde.results[2,1],t.par=kde.results[3,1])
    #plot(data.series[,series],data.series.results$kdedps_pdf,log="xy")
    #plot(data.series[,series],data.series.results$kdedps_cdf,log="xy")

    ### Kolmogorov-Smirnov test
   
    kde.sample.model.dps<-rdoublepareto.simplified(n=no.events, alpha.par=kde.results[1,1], beta.par=kde.results[2,1], t.par=kde.results[3,1]) 
    ks_model_greater<-NULL
    ks_model_greater<-ks.boot(data.series[,series],kde.sample.model.dps,nboots=ks_boot_samples,alternative="two.sided")
    ks_greater_result[4,2:5]<-c(ks_model_greater$nboots,ks_model_greater$ks$statistic,ks_model_greater$ks$p.value,ks_model_greater$ks.boot.pvalue)
    
#	if (Sys.info()[1] == "Linux") {x11()}
#	if (Sys.info()[1] == "Windows") {windows()}
#	if (Sys.info()[1] == "Darwin") {quartz()}
#   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""), xlab=xlabel,col="blue",xlim=range.plot.area,ylim=range.plot.density)
#   lines(linear.x.data.series, kde.data.series.linear, col="orange",lty=1,lwd=1)
#   lines(linear.x.data.series,10^value.fit.kde.log.doublepareto.simplified,col="dark green",lty=2)
#   mtext(paste("KDE Double Pareto Simplified-> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
#mtext(paste("Boostrapped KS test -> D: ",round(ks_greater_result[4,3],3),"; p-value: ",round(ks_greater_result[4,5],3),sep=""), side=3, col="navyblue", cex=0.8, line=0.5,padj=2.5)


   pdf(file = paste("KDE_Fit_DPS_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
   par(cex.axis=0.9)
   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""), xlab=xlabel,col="blue",xlim=range.plot.area,ylim=range.plot.density)
   lines(linear.x.data.series, kde.data.series.linear, col="orange",lty=1,lwd=1)
   lines(data.series.synthetic,10^value.fit.kde.log.doublepareto.simplified,col="dark green",lty=2)
   mtext(paste("KDE Double Pareto Simplified-> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
    mtext(paste("Boostrapped KS test -> D: ",round(ks_greater_result[4,3],3),"; p-value: ",round(ks_greater_result[4,5],3),sep=""), side=3, col="navyblue", cex=0.8, line=0.5,padj=2.5)
   dev.off()


   ### Fit log KDE better  Double Pareto

   # Default paramenter value, specified_by_user == "NO" in configuration.txt file
   alpha.range<-c(0.1,5)
   beta.range<-c(0.1,5)
   t.range<-c(min(na.omit(data.series.sample))*2,max(na.omit(data.series.sample))/10)
   # Paramenter value specified by the user, specified_by_user == "YES" in configuration.txt file
   if (configuration[13,4]=="YES") (alpha.range<-as.numeric(configuration[13,2:3]))
   if (configuration[14,4]=="YES") (beta.range<-as.numeric(configuration[14,2:3]))
   if (configuration[15,4]=="YES") (t.range<-as.numeric(configuration[15,2:3]))
   
   print(paste("KDE (DP):",names(data.series)[series]))

   log.dp.kde.formula<-formula(log10(kde.data.series.linear) ~ log10((beta.par)/(t.par*(1-((1+(m.par.value/t.par)^(-alpha.par))/(1+(c.par.value/t.par)^(-alpha.par)))^(beta.par/alpha.par)))*(((1+(m.par.value/t.par)^(-alpha.par))^(beta.par/alpha.par))/((1+(linear.x.data.series/t.par)^(-alpha.par))^(1+(beta.par/alpha.par))))*((linear.x.data.series/t.par)^(-alpha.par-1))))
   fit.kde.log.doublepareto <- nls(log.dp.kde.formula,algorithm="port",trace=TRUE,lower=c(alpha.par=alpha.range[1], beta.par=beta.range[1],t.par=t.range[1]),upper=c(alpha.par=alpha.range[2], beta.par=beta.range[2],t.par=t.range[2]),start=list(alpha.par=sum(alpha.range)/2, beta.par=sum(beta.range)/2,t.par=sum(t.range)/4),control=list(maxiter=1000))
   summary(fit.kde.log.doublepareto)

   #value.fit.kde.log.doublepareto<-predict(fit.kde.log.doublepareto)
  value.fit.kde.log.doublepareto<-predict(fit.kde.log.doublepareto,newdata=list(linear.x.data.series=data.series.synthetic))


   # Storing fitting results
   kde.results[5:7,1:4]<-coef(summary(fit.kde.log.doublepareto))
   kde.results[8:9,1]<-rbind(c.par.value,m.par.value)

   # Calculating DP rollover
   kde.index.rollover.dp<-NULL
   kde.index.rollover.dp<-which(10^value.fit.kde.log.doublepareto==max(10^value.fit.kde.log.doublepareto))
   if (length(kde.index.rollover.dp)>1) {kde.index.rollover.dp<-kde.index.rollover.dp[1]}
   if (kde.index.rollover.dp>1)
    {
    kde.results[10,1]<-data.series.synthetic[kde.index.rollover.dp]
    }

    # Calculating pdf and cdf for the original data
    data.series.results$kdedp_pdf<-ddoublepareto(x=data.series[,series],alpha.par=kde.results[5,1],beta.par=kde.results[6,1],t.par=kde.results[7,1],c.par=kde.results[8,1],m.par=kde.results[9,1])
    data.series.results$kdedp_cdf<-pdoublepareto(x=data.series[,series],alpha.par=kde.results[5,1],beta.par=kde.results[6,1],t.par=kde.results[7,1],c.par=kde.results[8,1],m.par=kde.results[9,1])
    #plot(data.series[,series],data.series.results$kdedp_pdf,log="xy")
    #plot(data.series[,series],data.series.results$kdedp_cdf,log="xy")

    ### Kolmogorov-Smirnov test
    
    kde.sample.model.dp<-rdoublepareto(n=no.events, alpha.par=kde.results[5,1], beta.par=kde.results[6,1], t.par=kde.results[7,1],c.par=kde.results[8,1],m.par=kde.results[9,1]) 
    ks_model_greater<-NULL
    ks_model_greater<-ks.boot(data.series[,series],kde.sample.model.dp,nboots=ks_boot_samples,alternative="two.sided")
    ks_greater_result[5,2:5]<-c(ks_model_greater$nboots,ks_model_greater$ks$statistic,ks_model_greater$ks$p.value,ks_model_greater$ks.boot.pvalue)
    
#	if (Sys.info()[1] == "Linux") {x11()}
#	if (Sys.info()[1] == "Windows") {windows()}
#	if (Sys.info()[1] == "Darwin") {quartz()}
#   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""), xlab=xlabel,col="blue",xlim=range.plot.area,ylim=range.plot.density)
#   lines(linear.x.data.series, kde.data.series.linear, col="orange",lty=1,lwd=1)
#   lines(linear.x.data.series,10^value.fit.kde.log.doublepareto,col="red",lty=2)
#   mtext(paste("KDE Double Pareto -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
#mtext(paste("Boostrapped KS test -> D: ",round(ks_greater_result[5,3],3),"; p-value: ",round(ks_greater_result[5,5],3),sep=""), side=3, col="navyblue", cex=0.8, line=0.5,padj=2.5)


   pdf(file = paste("KDE_Fit_DP_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
   par(cex.axis=0.9)
   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""), xlab=xlabel,col="blue",xlim=range.plot.area,ylim=range.plot.density)
   lines(linear.x.data.series, kde.data.series.linear, col="orange",lty=1,lwd=1)
   lines(data.series.synthetic,10^value.fit.kde.log.doublepareto,col="red",lty=2)
   mtext(paste("KDE Double Pareto -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
    mtext(paste("Boostrapped KS test -> D: ",round(ks_greater_result[5,3],3),"; p-value: ",round(ks_greater_result[5,5],3),sep=""), side=3, col="navyblue", cex=0.8, line=0.5,padj=2.5)
   dev.off()


   ### Fit log KDE better Inverse gamma

   # Default paramenter value, specified_by_user == "NO" in configuration.txt file
   alpha.range<-c(0.1,5)
   eta.range<-c(1,100)
   lambda.range<-c(1,500)
   # Paramenter value specified by the user, specified_by_user == "YES" in configuration.txt file
   if (configuration[16,4]=="YES") (alpha.range<-as.numeric(configuration[16,2:3]))
   if (configuration[17,4]=="YES") (eta.range<-as.numeric(configuration[17,2:3]))
   if (configuration[18,4]=="YES") (lambda.range<-as.numeric(configuration[18,2:3]))
   

   print(paste("KDE (IG):",names(data.series)[series]))


#   ig.kde.formula<-formula(kde.data.series.linear ~ ((lambda.par^(2*alpha.par))/gamma(alpha.par))*((1/(linear.x.data.series+(eta.par^2)))^(alpha.par+1))*exp(-(lambda.par^2)/(linear.x.data.series+(eta.par^2))))
#
#   fit.kde.inversegamma <- nls(ig.kde.formula,algorithm="port",trace=TRUE,lower=c(alpha.par=alpha.range[1],eta.par=eta.range[1],lambda.par=lambda.range[1]),upper=c(alpha.par=alpha.range[2],eta.par=eta.range[2],lambda.par=lambda.range[2]),start=list(alpha.par=sum(alpha.range)/2,eta.par=sum(eta.range)/2,lambda.par=sum(lambda.range)/2),control=list(maxiter=1000))
#
#   summary(fit.kde.inversegamma)
#   value.fit.kde.inversegamma<-predict(fit.kde.inversegamma)
#
#   # Storing fitting results
#   kde.results[9:11,1:4]<-coef(summary(fit.kde.inversegamma))




   log.ig.kde.formula<-formula(log10(kde.data.series.linear) ~ log10(((lambda.par^(2*alpha.par))/gamma(alpha.par))*((1/(linear.x.data.series+(eta.par^2)))^(alpha.par+1))*exp(-(lambda.par^2)/(linear.x.data.series+(eta.par^2)))))

   fit.kde.log.inversegamma <- nls(log.ig.kde.formula,algorithm="port",trace=TRUE,lower=c(alpha.par=alpha.range[1],eta.par=eta.range[1],lambda.par=lambda.range[1]),upper=c(alpha.par=alpha.range[2],eta.par=eta.range[2],lambda.par=lambda.range[2]),start=list(alpha.par=sum(alpha.range)/2,eta.par=sum(eta.range)/2,lambda.par=sum(lambda.range)/2),control=list(maxiter=1000))

   summary(fit.kde.log.inversegamma)

   #value.fit.kde.log.inversegamma<-predict(fit.kde.log.inversegamma)
  value.fit.kde.log.inversegamma<-predict(fit.kde.log.inversegamma,newdata=list(linear.x.data.series=data.series.synthetic))



   # Storing fitting results
   kde.results[11:13,1:4]<-coef(summary(fit.kde.log.inversegamma))

   # Calculating IG rollover
   kde.index.rollover.ig<-NULL
   kde.index.rollover.ig<-which(10^value.fit.kde.log.inversegamma==max(10^value.fit.kde.log.inversegamma))
   if (length(kde.index.rollover.ig)>1) {kde.index.rollover.ig<-kde.index.rollover.ig[1]}
   if (kde.index.rollover.ig>1)
    {
    kde.results[14,1]<-data.series.synthetic[kde.index.rollover.ig]
    }

    # Calculating pdf and cdf for the original data
    data.series.results$kdeig_pdf<-dinversegamma(x=data.series[,series],alpha.par=kde.results[11,1],eta.par=kde.results[12,1],lambda.par=kde.results[13,1])
    data.series.results$kdeig_cdf<-pinversegamma(x=data.series[,series],alpha.par=kde.results[11,1],eta.par=kde.results[12,1],lambda.par=kde.results[13,1])
    #plot(data.series[,series],data.series.results$kdeig_pdf,log="xy")
    #plot(data.series[,series],data.series.results$kdeig_cdf,log="xy")

    ### Kolmogorov-Smirnov test
    kde.sample.model.ig<-rinversegamma(n=no.events, alpha.par=kde.results[11,1], eta.par=kde.results[12,1], lambda.par=kde.results[13,1])
    ks_model_greater<-NULL
    ks_model_greater<-ks.boot(data.series[,series],kde.sample.model.ig,nboots=ks_boot_samples,alternative="two.sided")
    ks_greater_result[6,2:5]<-c(ks_model_greater$nboots,ks_model_greater$ks$statistic,ks_model_greater$ks$p.value,ks_model_greater$ks.boot.pvalue)
    
#	if (Sys.info()[1] == "Linux") {x11()}
#	if (Sys.info()[1] == "Windows") {windows()}
#	if (Sys.info()[1] == "Darwin") {quartz()}
#   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""), xlab=xlabel,col="blue",xlim=range.plot.area,ylim=range.plot.density)
#   lines(linear.x.data.series, kde.data.series.linear, col="orange",lty=1,lwd=1)
#   lines(linear.x.data.series,10^value.fit.kde.log.inversegamma,col="violet",lty=2)
#   #lines(linear.x.data.series,value.fit.kde.inversegamma,col="violet",lty=2)
#   mtext(paste("KDE Inverse Gamma -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
#mtext(paste("Boostrapped KS test -> D: ",round(ks_greater_result[6,3],3),"; p-value: ",round(ks_greater_result[6,5],3),sep=""), side=3, col="navyblue", cex=0.8, line=0.5,padj=2.5)


   pdf(file = paste("KDE_Fit_IG_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
   par(cex.axis=0.9)
   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""), xlab=xlabel,col="blue",xlim=range.plot.area,ylim=range.plot.density)
   lines(linear.x.data.series, kde.data.series.linear, col="orange",lty=1,lwd=1)
   lines(data.series.synthetic,10^value.fit.kde.log.inversegamma,col="violet",lty=2)
   #lines(linear.x.data.series,value.fit.kde.inversegamma,col="violet",lty=2)
   mtext(paste("KDE Inverse Gamma -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
    mtext(paste("Boostrapped KS test -> D: ",round(ks_greater_result[6,3],3),"; p-value: ",round(ks_greater_result[6,5],3),sep=""), side=3, col="navyblue", cex=0.8, line=0.5,padj=2.5)
   dev.off()


   # Plot comparison KDE different distribution
	if (Sys.info()[1] == "Linux") {x11()}
	if (Sys.info()[1] == "Windows") {windows()}
	if (Sys.info()[1] == "Darwin") {quartz()}
   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""), xlab=xlabel,col="blue",xlim=range.plot.area,ylim=range.plot.density)
   lines(linear.x.data.series, kde.data.series.linear, col="orange",lty=1,lwd=1)
   lines(data.series.synthetic,10^value.fit.kde.log.doublepareto.simplified,col="dark green",lty=2)
   lines(data.series.synthetic,10^value.fit.kde.log.doublepareto,col="red",lty=2)
   lines(data.series.synthetic,10^value.fit.kde.log.inversegamma,col="violet",lty=2)
   #lines(linear.x.data.series,value.fit.kde.inversegamma,col="violet",lty=2)
   mtext(paste("KDE Comparison -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
   legend("bottomleft",legend=c("Histogram raw data","Kernel density raw data","Double Pareto simplified","Double Pareto","Inverse Gamma"),pch=c(1,NA_integer_,NA_integer_,NA_integer_,NA_integer_),lty=c(FALSE,1,2,2,2),lwd=1,col=c("blue","orange","dark green","red","violet"),cex=0.8,box.lty=3,box.col="black")

   pdf(file = paste("KDE_Fit_Comparison_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
   par(cex.axis=0.9)
   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""), xlab=xlabel,col="blue",xlim=range.plot.area,ylim=range.plot.density)
   lines(linear.x.data.series, kde.data.series.linear, col="orange",lty=1,lwd=1)
   lines(data.series.synthetic,10^value.fit.kde.log.doublepareto.simplified,col="dark green",lty=2)
   lines(data.series.synthetic,10^value.fit.kde.log.doublepareto,col="red",lty=2)
   lines(data.series.synthetic,10^value.fit.kde.log.inversegamma,col="violet",lty=2)
   #lines(linear.x.data.series,value.fit.kde.inversegamma,col="violet",lty=2)
   mtext(paste("KDE Comparison -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
   legend("bottomleft",legend=c("Histogram raw data","Kernel density raw data","Double Pareto simplified","Double Pareto","Inverse Gamma"),pch=c(1,NA_integer_,NA_integer_,NA_integer_,NA_integer_),lty=c(FALSE,1,2,2,2),lwd=1,col=c("blue","orange","dark green","red","violet"),cex=0.8,box.lty=3,box.col="black")
   dev.off()
   
   
  #### -------------- Uncertainty estimation -------------- ####
	   
	### Creating object for storing sampling results
	kde_dps_sampling_result<-matrix(NA,nrow=length(data.series.synthetic.sample),ncol=samples_number)
	kde_dps_sampling_result_coef<-matrix(NA,nrow=samples_number,ncol=4)
	colnames(kde_dps_sampling_result_coef)<-c("alpha","beta","t","rollov")
	  
	kde_dp_sampling_result<-matrix(NA,nrow=length(data.series.synthetic.sample),ncol=samples_number)
	kde_dp_sampling_result_coef<-matrix(NA,nrow=samples_number,ncol=6)
	colnames(kde_dp_sampling_result_coef)<-c("alpha","beta","t","c","m","rollov")
	   
	kde_ig_sampling_result<-matrix(NA,nrow=length(data.series.synthetic.sample),ncol=samples_number)
	kde_ig_sampling_result_coef<-matrix(NA,nrow=samples_number,ncol=4)
	colnames(kde_ig_sampling_result_coef)<-c("alpha","eta","lambda","rollov")
   
	for (countsamples in 1:samples_number)
		{
		#countsamples<-1
		c.par.value.sample<-c.par.value
		m.par.value.sample<-m.par.value
		print(paste("KDE -",names(data.series)[series],"-> Sample:",countsamples))
		data.series.sample<-sample(data.series[,series],replace=TRUE)
		no.events.sample<-length(data.series.sample)
		log.data.series.sample<-log(data.series.sample,10)
		
		#bandwidth.data.series.sample<-bw.nrd0(log.data.series.sample)
		#kde.data.series.sample<-density(log.data.series.sample,bw=bandwidth.data.series.sample,from=min(log.data.series.sample),to=max(log.data.series.sample),cut=range(log.data.series.sample),adjust=1, kernel=c("g"), n=num.kde, give.Rkern=FALSE)
		
		kde.data.series.sample<-density(log.data.series.sample,bw=bandwidth.data.series,from=min(log.data.series.sample),to=max(log.data.series.sample),cut=range(log.data.series.sample),adjust=1, kernel=c("g"), n=num.kde, give.Rkern=FALSE)
		kernel.diff.data.series.sample<-diff(10^kde.data.series.sample$x)
		kde.data.series.linear.sample<-kde.data.series.sample$y[1:(length(kde.data.series.sample$y)-1)]*((kde.data.series.sample$x[2:length(kde.data.series.sample$x)])-(kde.data.series.sample$x[1:(length(kde.data.series.sample$x)-1)]))/kernel.diff.data.series
		
		freq.linear.data.series.sample<-kde.data.series.sample$y[1:(length(kde.data.series.sample$y)-1)]*no.events*((kde.data.series.sample$x[2:length(kde.data.series.sample$x)])-(kde.data.series.sample$x[1:(length(kde.data.series.sample$x)-1)]))/kernel.diff.data.series
		linear.x.data.series.sample<-10^(kde.data.series.sample$x[1:(length(kde.data.series.sample$x)-1)])
		
		kde.data.series.linear.log.sample<-log(kde.data.series.linear.sample,base=10)
		freq.linear.log.kde.data.series.sample<-log(freq.linear.data.series.sample,base=10)
		linear.x.log.kde.data.series.sample<-log(linear.x.data.series.sample,base=10)

		### ----------------- Fit Double Pareto Simplified (DPS) ----------------- ###
		# Default paramenter value, specified_by_user == "NO" in configuration.txt file
		alpha.range<-c(0.1,5)
		beta.range<-c(0.1,5)
		t.range<-c(min(na.omit(data.series.sample))*2,max(na.omit(data.series.sample))/10)
		# Paramenter value specified by the user, specified_by_user == "YES" in configuration.txt file
		if (configuration[13,4]=="YES") (alpha.range<-as.numeric(configuration[13,2:3]))
		if (configuration[14,4]=="YES") (beta.range<-as.numeric(configuration[14,2:3]))
		if (configuration[15,4]=="YES") (t.range<-as.numeric(configuration[15,2:3]))
		
		print(paste("KDE (DPS):",names(data.series)[series],"-> Sample:",countsamples))
		log.dps.kde.formula.sample<-formula(log10(kde.data.series.linear.sample) ~ log10((beta.par*(t.par^alpha.par))/((1+((linear.x.data.series.sample/t.par)^(-alpha.par)))^(1+(beta.par/alpha.par))*(linear.x.data.series.sample)^(alpha.par+1))))
		if (inherits(try(fit.kde.log.doublepareto.simplified.sample<-nls(log.dps.kde.formula.sample,algorithm="port",trace=FALSE,lower=c(alpha.par=alpha.range[1], beta.par=beta.range[1],t.par=t.range[1]),upper=c(alpha.par=alpha.range[2], beta.par=beta.range[2],t.par=t.range[2]),start=list(alpha.par=sum(alpha.range)/2, beta.par=sum(beta.range)/2,t.par=sum(t.range)/4),control=list(maxiter=1000)),silent=TRUE),what="try-error"))
			{
			next()
			}	
		# Storing fitting results
		kde_dps_sampling_result_coef[countsamples,1:3]<-coef(summary(fit.kde.log.doublepareto.simplified.sample))[,1]
		kde_dps_sampling_result[,countsamples]<-predict(fit.kde.log.doublepareto.simplified.sample,newdata=list(linear.x.data.series.sample=data.series.synthetic.sample))
		# Calculating DPS rollover
		kde.index.rollover.dps.sample<-NULL
		kde.index.rollover.dps.sample<-which(10^kde_dps_sampling_result[,countsamples]==max(10^kde_dps_sampling_result[,countsamples]))
		if (length(kde.index.rollover.dps.sample)>1) {kde.index.rollover.dps.sample<-kde.index.rollover.dps.sample[1]}
		if (kde.index.rollover.dps.sample>1)
			{
			kde_dps_sampling_result_coef[countsamples,4]<-data.series.synthetic.sample[kde.index.rollover.dps.sample]
			}

		### ----------------- Fit Double Pareto (DP) ----------------- ###
		# Default paramenter value, specified_by_user == "NO" in configuration.txt file
		alpha.range<-c(0.1,5)
		beta.range<-c(0.1,5)
		t.range<-c(min(na.omit(data.series.sample))*2,max(na.omit(data.series.sample))/10)
		# Paramenter value specified by the user, specified_by_user == "YES" in configuration.txt file
		if (configuration[13,4]=="YES") (alpha.range<-as.numeric(configuration[13,2:3]))
		if (configuration[14,4]=="YES") (beta.range<-as.numeric(configuration[14,2:3]))
		if (configuration[15,4]=="YES") (t.range<-as.numeric(configuration[15,2:3]))
		
		print(paste("KDE (DP):",names(data.series)[series],"-> Sample:",countsamples))
		log.dp.kde.formula.sample<-formula(log10(kde.data.series.linear.sample) ~ log10((beta.par*(t.par^alpha.par))/((1+((linear.x.data.series.sample/t.par)^(-alpha.par)))^(1+(beta.par/alpha.par))*(linear.x.data.series.sample)^(alpha.par+1))))
		if (inherits(try(fit.kde.log.doublepareto.sample<-nls(log.dp.kde.formula.sample,algorithm="port",trace=FALSE,lower=c(alpha.par=alpha.range[1], beta.par=beta.range[1],t.par=t.range[1]),upper=c(alpha.par=alpha.range[2], beta.par=beta.range[2],t.par=t.range[2]),start=list(alpha.par=sum(alpha.range)/2, beta.par=sum(beta.range)/2,t.par=sum(t.range)/4),control=list(maxiter=1000)),silent=TRUE),what="try-error"))
			{
			next()
			}	
		# Storing fitting results
		kde_dp_sampling_result_coef[countsamples,1:3]<-coef(summary(fit.kde.log.doublepareto.sample))[,1]
		kde_dp_sampling_result_coef[countsamples,4:5]<-c(c.par.value.sample,m.par.value.sample) 
		kde_dp_sampling_result[,countsamples]<-predict(fit.kde.log.doublepareto.sample,newdata=list(linear.x.data.series.sample=data.series.synthetic.sample))
		# Calculating DPS rollover
		kde.index.rollover.dp.sample<-NULL
		kde.index.rollover.dp.sample<-which(10^kde_dp_sampling_result[,countsamples]==max(10^kde_dp_sampling_result[,countsamples]))
		if (length(kde.index.rollover.dp.sample)>1) {kde.index.rollover.dp.sample<-kde.index.rollover.dp.sample[1]}
		if (kde.index.rollover.dp.sample>1)
			{
			kde_dp_sampling_result_coef[countsamples,6]<-data.series.synthetic.sample[kde.index.rollover.dp.sample]
			}
		
		### ----------------- Fit Inverse Gamma (IG) ----------------- ###
		# Default paramenter value, specified_by_user == "NO" in configuration.txt file
		alpha.range<-c(0.1,5)
		eta.range<-c(1,100)
		lambda.range<-c(1,500)
		# Paramenter value specified by the user, specified_by_user == "YES" in configuration.txt file
		if (configuration[16,4]=="YES") (alpha.range<-as.numeric(configuration[16,2:3]))
		if (configuration[17,4]=="YES") (eta.range<-as.numeric(configuration[17,2:3]))
		if (configuration[18,4]=="YES") (lambda.range<-as.numeric(configuration[18,2:3]))
		
		print(paste("KDE (IG):",names(data.series)[series],"-> Sample:",countsamples))
		log.ig.kde.formula.sample<-formula(log10(kde.data.series.linear.sample) ~ log10(((lambda.par^(2*alpha.par))/gamma(alpha.par))*((1/(linear.x.data.series.sample+(eta.par^2)))^(alpha.par+1))*exp(-(lambda.par^2)/(linear.x.data.series.sample+(eta.par^2)))))
		if (inherits(try(fit.kde.log.inversegamma.sample<-nls(log.ig.kde.formula.sample,algorithm="port",trace=FALSE,lower=c(alpha.par=alpha.range[1],eta.par=eta.range[1],lambda.par=lambda.range[1]),upper=c(alpha.par=alpha.range[2],eta.par=eta.range[2],lambda.par=lambda.range[2]),start=list(alpha.par=sum(alpha.range)/2,eta.par=sum(eta.range)/2,lambda.par=sum(lambda.range)/2),control=list(maxiter=1000)),silent=TRUE),what="try-error"))
			{
			next()
			}	
		# Storing fitting results
		kde_ig_sampling_result_coef[countsamples,1:3]<-coef(summary(fit.kde.log.inversegamma.sample))[,1]
		kde_ig_sampling_result[,countsamples]<-predict(fit.kde.log.inversegamma.sample,newdata=list(linear.x.data.series.sample=data.series.synthetic.sample))
		# Calculating DPS rollover
		kde.index.rollover.ig.sample<-NULL
		kde.index.rollover.ig.sample<-which(10^kde_ig_sampling_result[,countsamples]==max(10^kde_ig_sampling_result[,countsamples]))
		if (length(kde.index.rollover.ig.sample)>1) {kde.index.rollover.ig.sample<-kde.index.rollover.ig.sample[1]}
		if (kde.index.rollover.ig.sample>1)
			{
			kde_ig_sampling_result_coef[countsamples,4]<-data.series.synthetic.sample[kde.index.rollover.ig.sample]
			}
		}	
		
		quantiles_kde_dps_sampling_result_coef<-apply(kde_dps_sampling_result_coef,MARGIN=2,FUN=quantile_fun)
		quantiles_kde_dp_sampling_result_coef<-apply(kde_dp_sampling_result_coef,MARGIN=2,FUN=quantile_fun)
		quantiles_kde_ig_sampling_result_coef<-apply(kde_ig_sampling_result_coef,MARGIN=2,FUN=quantile_fun)
					
		write.table(data.frame(percentile=rownames(quantiles_kde_dps_sampling_result_coef),round(quantiles_kde_dps_sampling_result_coef,2)),paste("KDE_Fit_DPS_parameter_uncertainty_",names(data.series)[series],".txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t")
		write.table(data.frame(percentile=rownames(quantiles_kde_dp_sampling_result_coef),round(quantiles_kde_dp_sampling_result_coef,2)),paste("KDE_Fit_DP_parameter_uncertainty_",names(data.series)[series],".txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t")
		write.table(data.frame(percentile=rownames(quantiles_kde_ig_sampling_result_coef),round(quantiles_kde_ig_sampling_result_coef,2)),paste("KDE_Fit_IG_parameter_uncertainty_",names(data.series)[series],".txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t")
		
		quantiles_kde_dps_sampling_result<-apply(kde_dps_sampling_result,MARGIN=1,FUN=quantile_fun)
		quantiles_kde_dp_sampling_result<-apply(kde_dp_sampling_result,MARGIN=1,FUN=quantile_fun)
		quantiles_kde_ig_sampling_result<-apply(kde_ig_sampling_result,MARGIN=1,FUN=quantile_fun)
 
				
		### Plot KDE estimation with uncertainty
		 
		pdf(file = paste("KDE_Fit_DPS_uncertainty_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
		par(cex.axis=0.9)
		#if (Sys.info()[1] == "Linux") {x11()}
		#if (Sys.info()[1] == "Windows") {windows()}
		#if (Sys.info()[1] == "Darwin") {quartz()}
		#xlabel<-expression(A (m^{2}))
		ylabel<-paste("Probability density",sep="")
		plot(NA, NA, log="xy", main="Probability densities", ylab=ylabel, xlab=xlabel,xlim=range.plot.area,ylim=range.plot.density,xaxt="n",yaxt="n")
		axis(side=1,at=at.x,label=sciNotation(at.x,1))
		axis(side=2,at=at.y,label=sciNotation(at.y,1))
		axis(1, at.x.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		axis(2, at.y.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		axis(3, at.x.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		axis(4, at.y.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		for (count_q in 1:floor((length(prob_series))/2))
			{
			#count_q<-1	
			polygon(x=c(sort(data.series.synthetic.sample),rev(sort(data.series.synthetic.sample))),y=10^c(quantiles_kde_dps_sampling_result[length(prob_series)+1-count_q,],rev(quantiles_kde_dps_sampling_result[count_q,])),border=NA,col=color_vector_pol[count_q])
			}
		lines(linear.x.data.series, kde.data.series.linear, col="orange",lty=1,lwd=1)
		lines(data.series.synthetic,10^value.fit.kde.log.doublepareto.simplified,col="darkred",lty=1,lwd=2)
		lines(data.series.synthetic.sample,10^quantiles_kde_dps_sampling_result["50%",],col="black",lty="dashed",lwd=1)
		points(linear.mids.data.series, prob.linear.data.series,pch=20,cex=1,col="navyblue")
		mtext(paste("KDE Double Pareto Simplified -> ",names(data.series)[series],sep=""), side=3, col="darkred", cex=0.8, line=0.5)
    mtext(paste("Boostrapped KS test -> D: ",round(ks_greater_result[4,3],3),"; p-value: ",round(ks_greater_result[4,5],3),sep=""), side=3, col="navyblue", cex=0.8, line=0.5,padj=2.5)
		legend("bottomleft",c("Histogram empirical data","Raw kernel density data","Distribution fit","Median values",paste("Quantiles: ",prob_series[1:floor((length(prob_series))/2)],"-",prob_series[(length(prob_series)+1)-1:floor((length(prob_series))/2)],sep="")),pt.cex=1,pch=c(20,NA,NA,NA,rep(NA,floor((length(prob_series))/2))),lty=c(NA,1,1,2,rep(1,floor((length(prob_series))/2))),lwd=c(NA,1,2,1,rep(6,floor((length(prob_series))/2))),col=c("navyblue","orange","darkred","black",color_vector_pol),bty="o",bg="transparent",box.col="transparent",cex=0.9)
		dev.off()
 

		pdf(file = paste("KDE_Fit_DP_uncertainty_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
		par(cex.axis=0.9)
		#if (Sys.info()[1] == "Linux") {x11()}
		#if (Sys.info()[1] == "Windows") {windows()}
		#if (Sys.info()[1] == "Darwin") {quartz()}
		#xlabel<-expression(A (m^{2}))
		ylabel<-paste("Probability density",sep="")
		plot(NA, NA, log="xy", main="Probability densities", ylab=ylabel, xlab=xlabel,xlim=range.plot.area,ylim=range.plot.density,xaxt="n",yaxt="n")
		axis(side=1,at=at.x,label=sciNotation(at.x,1))
		axis(side=2,at=at.y,label=sciNotation(at.y,1))
		axis(1, at.x.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		axis(2, at.y.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		axis(3, at.x.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		axis(4, at.y.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		for (count_q in 1:floor((length(prob_series))/2))
			{
			#count_q<-1	
			polygon(x=c(sort(data.series.synthetic.sample),rev(sort(data.series.synthetic.sample))),y=10^c(quantiles_kde_dp_sampling_result[length(prob_series)+1-count_q,],rev(quantiles_kde_dp_sampling_result[count_q,])),border=NA,col=color_vector_pol[count_q])
			}
		lines(linear.x.data.series, kde.data.series.linear, col="orange",lty=1,lwd=1)
		lines(data.series.synthetic,10^value.fit.kde.log.doublepareto,col="darkred",lty=1,lwd=2)
		lines(data.series.synthetic.sample,10^quantiles_kde_dp_sampling_result["50%",],col="black",lty="dashed",lwd=1)
		points(linear.mids.data.series, prob.linear.data.series,pch=20,cex=1,col="navyblue")
		mtext(paste("KDE Double Pareto -> ",names(data.series)[series],sep=""), side=3, col="darkred", cex=0.8, line=0.5)
    mtext(paste("Boostrapped KS test -> D: ",round(ks_greater_result[5,3],3),"; p-value: ",round(ks_greater_result[5,5],3),sep=""), side=3, col="navyblue", cex=0.8, line=0.5,padj=2.5)
		legend("bottomleft",c("Histogram empirical data","Raw kernel density data","Distribution fit","Median values",paste("Quantiles: ",prob_series[1:floor((length(prob_series))/2)],"-",prob_series[(length(prob_series)+1)-1:floor((length(prob_series))/2)],sep="")),pt.cex=1,pch=c(20,NA,NA,NA,rep(NA,floor((length(prob_series))/2))),lty=c(NA,1,1,2,rep(1,floor((length(prob_series))/2))),lwd=c(NA,1,2,1,rep(6,floor((length(prob_series))/2))),col=c("navyblue","orange","darkred","black",color_vector_pol),bty="o",bg="transparent",box.col="transparent",cex=0.9)
		dev.off()

		pdf(file = paste("KDE_Fit_IG_uncertainty_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
		par(cex.axis=0.9)
		#if (Sys.info()[1] == "Linux") {x11()}
		#if (Sys.info()[1] == "Windows") {windows()}
		#if (Sys.info()[1] == "Darwin") {quartz()}
		#xlabel<-expression(A (m^{2}))
		ylabel<-paste("Probability density",sep="")
		plot(NA, NA, log="xy", main="Probability densities", ylab=ylabel, xlab=xlabel,xlim=range.plot.area,ylim=range.plot.density,xaxt="n",yaxt="n")
		axis(side=1,at=at.x,label=sciNotation(at.x,1))
		axis(side=2,at=at.y,label=sciNotation(at.y,1))
		axis(1, at.x.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		axis(2, at.y.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		axis(3, at.x.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		axis(4, at.y.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		for (count_q in 1:floor((length(prob_series))/2))
			{
			#count_q<-1	
			polygon(x=c(sort(data.series.synthetic.sample),rev(sort(data.series.synthetic.sample))),y=10^c(quantiles_kde_ig_sampling_result[length(prob_series)+1-count_q,],rev(quantiles_kde_ig_sampling_result[count_q,])),border=NA,col=color_vector_pol[count_q])
			}
		lines(linear.x.data.series, kde.data.series.linear, col="orange",lty=1,lwd=1)
		
		lines(data.series.synthetic,10^value.fit.kde.log.inversegamma,col="darkred",lty=1,lwd=2)
		lines(data.series.synthetic.sample,10^quantiles_kde_ig_sampling_result["50%",],col="black",lty="dashed",lwd=1)
		points(linear.mids.data.series, prob.linear.data.series,pch=20,cex=1,col="navyblue")
		mtext(paste("KDE Inverse Gamma -> ",names(data.series)[series],sep=""), side=3, col="darkred", cex=0.8, line=0.5)
    mtext(paste("Boostrapped KS test -> D: ",round(ks_greater_result[6,3],3),"; p-value: ",round(ks_greater_result[6,5],3),sep=""), side=3, col="navyblue", cex=0.8, line=0.5,padj=2.5)
		legend("bottomleft",c("Histogram empirical data","Raw kernel density data","Distribution fit","Median values",paste("Quantiles: ",prob_series[1:floor((length(prob_series))/2)],"-",prob_series[(length(prob_series)+1)-1:floor((length(prob_series))/2)],sep="")),pt.cex=1,pch=c(20,NA,NA,NA,rep(NA,floor((length(prob_series))/2))),lty=c(NA,1,1,2,rep(1,floor((length(prob_series))/2))),lwd=c(NA,1,2,1,rep(6,floor((length(prob_series))/2))),col=c("navyblue","orange","darkred","black",color_vector_pol),bty="o",bg="transparent",box.col="transparent",cex=0.9)
		dev.off()


# ----------------------------- MLE estimation  ----------------------------- #
   print(paste("MLE -",names(data.series)[series]))


   #### MLE using Double Pareto Simplified

   # Default paramenter value, specified_by_user == "NO" in configuration.txt file
   alpha.range<-c(min(c(hde.results[1,1],kde.results[1,1]))/1.5,max(c(hde.results[1,1],kde.results[1,1]))*1.5)
   beta.range<-c(min(c(hde.results[2,1],kde.results[2,1]))/1.5,max(c(hde.results[2,1],kde.results[2,1]))*1.5)
   t.range<-c(min(c(hde.results[3,1],kde.results[2,1]))/1.5,max(c(hde.results[3,1],kde.results[3,1]))*1.5)
   # Paramenter value specified by the user, specified_by_user == "YES" in configuration.txt file
   if (configuration[19,4]=="YES") (alpha.range<-as.numeric(configuration[19,2:3]))
   if (configuration[20,4]=="YES") (beta.range<-as.numeric(configuration[20,2:3]))
   if (configuration[21,4]=="YES") (t.range<-as.numeric(configuration[21,2:3]))
   
   print(paste("MLE (DPS):",names(data.series)[series]))   
   
   mle.data.series.dps<-mle2(mll.doublepareto.simplified,method="L-BFGS-B",lower=c(alpha.par=alpha.range[1], beta.par=beta.range[1],t.par=t.range[1]),upper=c(alpha.par=alpha.range[2], beta.par=beta.range[2],t.par=t.range[2]),start=list(alpha.par=sum(alpha.range)/2, beta.par=sum(beta.range)/2,t.par=sum(t.range)/2),data=list(x=na.omit(data.series[,series])),control=list(maxit=1000000,lmm=5,factr=10^-8,pgtol=0,trace=6,parscale=c(alpha.par=1, beta.par=1,t.par=1)))
   summary(mle.data.series.dps)

  # PDF generated using parameters estimated by MLE using a Simplified Double Pareto distribution (dps)
   pdf.mle.data.series.dps<-ddoublepareto.simplified(x=data.series.synthetic, alpha.par=mle.data.series.dps@coef[1], beta.par=mle.data.series.dps@coef[2], t.par=mle.data.series.dps@coef[3])
   cdf.mle.data.series.dps<-pdoublepareto.simplified(x=data.series.synthetic, alpha.par=mle.data.series.dps@coef[1], beta.par=mle.data.series.dps@coef[2], t.par=mle.data.series.dps@coef[3])
#lines(data.series.synthetic,cdf.mle.data.series.dps,col="green")

   # Storing fitting results
   mle.results[1:3,1:4]<-coef(summary(mle.data.series.dps))

   # Calculating DPS rollover
   mle.index.rollover.dps<-NULL
   mle.index.rollover.dps<-which(pdf.mle.data.series.dps==max(pdf.mle.data.series.dps))
   if (length(mle.index.rollover.dps)>1) {mle.index.rollover.dps<-mle.index.rollover.dps[1]}
   if (mle.index.rollover.dps>1)
    {
    mle.results[4,1]<-data.series.synthetic[mle.index.rollover.dps]
    }

    # Calculating pdf and cdf for the original data
    data.series.results$mledps_pdf<-ddoublepareto.simplified(x=data.series[,series],alpha.par=mle.results[1,1],beta.par=mle.results[2,1],t.par=mle.results[3,1])
    data.series.results$mledps_cdf<-pdoublepareto.simplified(x=data.series[,series],alpha.par=mle.results[1,1],beta.par=mle.results[2,1],t.par=mle.results[3,1])
    #plot(data.series[,series],data.series.results$mledps_pdf,log="xy")
    #plot(data.series[,series],data.series.results$mledps_cdf,log="xy")

    ### Kolmogorov-Smirnov test
    mle.sample.model.dps<-rdoublepareto.simplified(n=no.events, alpha.par=mle.results[1,1], beta.par=mle.results[2,1], t.par=mle.results[3,1]) 
    ks_model_greater<-NULL
    ks_model_greater<-ks.boot(data.series[,series],mle.sample.model.dps,nboots=ks_boot_samples,alternative="two.sided")
    ks_greater_result[7,2:5]<-c(ks_model_greater$nboots,ks_model_greater$ks$statistic,ks_model_greater$ks$p.value,ks_model_greater$ks.boot.pvalue)
    
#   if (Sys.info()[1] == "Linux") {x11()}
#	if (Sys.info()[1] == "Windows") {windows()}
#	if (Sys.info()[1] == "Darwin") {quartz()}
#   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""), xlab=xlabel,col="blue",xlim=range.plot.area,ylim=range.plot.density)  # Histogram raw
#   lines(data.series.synthetic,pdf.mle.data.series.dps,col="dark green",lty=3)
#   mtext(paste("MLE Double Pareto Simplified -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
#mtext(paste("Boostrapped KS test -> D: ",round(ks_greater_result[7,3],3),"; p-value: ",round(ks_greater_result[7,5],3),sep=""), side=3, col="navyblue", cex=0.8, line=0.5,padj=2.5)

   pdf(file = paste("MLE_Fit_DPS_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
   par(cex.axis=0.9)
   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""), xlab=xlabel,col="blue",xlim=range.plot.area,ylim=range.plot.density)  # Histogram raw
   lines(data.series.synthetic,pdf.mle.data.series.dps,col="dark green",lty=3)
   mtext(paste("MLE Double Pareto Simplified -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
    mtext(paste("Boostrapped KS test -> D: ",round(ks_greater_result[7,3],3),"; p-value: ",round(ks_greater_result[7,5],3),sep=""), side=3, col="navyblue", cex=0.8, line=0.5,padj=2.5)
   dev.off()

#   ### Plot of DPS cdf
#	if (Sys.info()[1] == "Linux") {x11()}
#	if (Sys.info()[1] == "Windows") {windows()}
#	if (Sys.info()[1] == "Darwin") {quartz()}
#   plot(data.series.synthetic,cdf.mle.data.series.dps, log="x", main="Cumulative distribution function", ylab=paste("Cumulative probability density",sep=""), xlab=xlabel,type="l",col="dark green",xlim=range.plot.area,ylim=c(0,1))
#   mtext(paste("CDF MLE Double Pareto Simplified -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)

   pdf(file = paste("CDF_MLE_DPS_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
   plot(data.series.synthetic,cdf.mle.data.series.dps, log="x", main="Cumulative distribution function", ylab=paste("Cumulative probability density",sep=""), xlab=xlabel,type="l",col="dark green",xlim=range.plot.area,ylim=c(0,1))
   mtext(paste("CDF MLE Double Pareto Simplified -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
   dev.off()
   
   pdf(file = paste("CDF_MLE_DPS_EntireRange_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
   data.series.synthetic.cdf<-c(10^log_cycles_range[1],sort((seq(1,10,0.1)%*%t(10^(log_cycles_range[1]:log_cycles_range[2])))[-1,]))
   cdf.mle.data.series.dps.entire.range<-pdoublepareto.simplified(x=data.series.synthetic.cdf, alpha.par=mle.data.series.dps@coef[1], beta.par=mle.data.series.dps@coef[2], t.par=mle.data.series.dps@coef[3])

   plot(data.series.synthetic,cdf.mle.data.series.dps, log="x", main="Cumulative distribution function", ylab=paste("Cumulative probability density",sep=""), xlab=xlabel,type="l",xlim=range.plot.area,ylim=c(0,1),xaxt="n",yaxt="n",col="dark green")
   index.cdf.mle.data.fill.min<-which(data.series.synthetic.cdf < min(data.series.synthetic))
   lines(data.series.synthetic.cdf[index.cdf.mle.data.fill.min],cdf.mle.data.series.dps.entire.range[index.cdf.mle.data.fill.min],lty="dashed",col="dark green")
   index.cdf.mle.data.fill.max<-which(data.series.synthetic.cdf > max(data.series.synthetic))
   lines(data.series.synthetic.cdf[index.cdf.mle.data.fill.max],cdf.mle.data.series.dps.entire.range[index.cdf.mle.data.fill.max],lty="dashed",col="dark green")
   cdf.percentiles.sizes.qdoublepareto.simplified<-qdoublepareto.simplified(x=cdf_percentiles,alpha.par=mle.data.series.dps@coef[1], beta.par=mle.data.series.dps@coef[2], t.par=mle.data.series.dps@coef[3],lower.par=uniroot.range[1]/10, upper.par=uniroot.range[2]*10,FUN=FALSE)
   points(cdf.percentiles.sizes.qdoublepareto.simplified,cdf_percentiles,pch=21,bg="forestgreen",col="black")
   lines(cdf.percentiles.sizes.qdoublepareto.simplified,cdf_percentiles,type="h",lty="dotted",col="forestgreen")
   text(cdf.percentiles.sizes.qdoublepareto.simplified, cdf_percentiles, labels=paste("(",signif(cdf.percentiles.sizes.qdoublepareto.simplified,2),",",cdf_percentiles,")",sep=""), cex= 0.6,pos=3)
   axis(side=1,at=at.x,label=sciNotation(at.x,1))
   at.y.cdf<-seq(0,1,0.1)
   at.y.tick.cdf<-seq(0,1,0.01)
   axis(side=2,at=at.y.cdf,label=at.y.cdf)
   axis(1, at.x.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
   axis(2, at.y.tick.cdf, labels = NA, lty = 1, lwd = 1, tck = -0.01)
   axis(3, at.x.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
   axis(4, at.y.tick.cdf, labels = NA, lty = 1, lwd = 1, tck = -0.01)
   mtext(paste("CDF MLE Double Pareto Simplified -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
   dev.off()
   
   
   #### MLE using Double Pareto (Fixing c.par and m.par)

   # Default paramenter value, specified_by_user == "NO" in configuration.txt file
   alpha.range<-c(min(c(hde.results[5,1],kde.results[5,1]))/2,max(c(hde.results[5,1],kde.results[5,1]))*2)
   beta.range<-c(min(c(hde.results[6,1],kde.results[6,1]))/4,max(c(hde.results[6,1],kde.results[6,1]))*4)
   t.range<-c(min(c(hde.results[7,1],kde.results[7,1]))/2,max(c(hde.results[7,1],kde.results[7,1]))*2)
   # Paramenter value specified by the user, specified_by_user == "YES" in configuration.txt file
   if (configuration[22,4]=="YES") (alpha.range<-as.numeric(configuration[22,2:3]))
   if (configuration[23,4]=="YES") (beta.range<-as.numeric(configuration[23,2:3]))
   if (configuration[24,4]=="YES") (t.range<-as.numeric(configuration[24,2:3]))
      
   print(paste("MLE (DP):",names(data.series)[series]))

   mle.data.series.dp<-mle2(mll.doublepareto,method="L-BFGS-B",fixed=list(c.par=c.par.value, m.par=m.par.value), lower=c(alpha.par=alpha.range[1], beta.par=beta.range[1],t.par=t.range[1]),upper=c(alpha.par=alpha.range[2], beta.par=beta.range[2],t.par=t.range[2]),start=list(alpha.par=sum(alpha.range)/2, beta.par=sum(beta.range)/2,t.par=sum(t.range)/2),data=list(x=na.omit(data.series[,series])),control=list(maxit=1000000,lmm=5,factr=10^-8,pgtol=0,trace=6,parscale=c(alpha.par=1, beta.par=1,t.par=1)))
   summary(mle.data.series.dp)

   # PDF generated using parameters estimated by MLE using a Double Pareto distribution (dp)
   pdf.mle.data.series.dp<-ddoublepareto(x=data.series.synthetic, alpha.par=mle.data.series.dp@coef[1], beta.par=mle.data.series.dp@coef[2], t.par=mle.data.series.dp@coef[3], c.par=c.par.value, m.par=m.par.value)
   cdf.mle.data.series.dp<-pdoublepareto(x=data.series.synthetic, alpha.par=mle.data.series.dp@coef[1], beta.par=mle.data.series.dp@coef[2], t.par=mle.data.series.dp@coef[3], c.par=c.par.value, m.par=m.par.value)
#plot(data.series.synthetic,cdf.mle.data.series.dp,col="red")

   
   mle.results[5:7,1:4]<-coef(summary(mle.data.series.dp))
   mle.results[8:9,1]<-rbind(c.par.value,m.par.value)

   # Calculating DP rollover
   mle.index.rollover.dp<-NULL
   mle.index.rollover.dp<-which(pdf.mle.data.series.dp==max(pdf.mle.data.series.dp))
   if (length(mle.index.rollover.dp)>1) {mle.index.rollover.dp<-mle.index.rollover.dp[1]}
   if (mle.index.rollover.dp>1)
    {
    mle.results[10,1]<-data.series.synthetic[mle.index.rollover.dp]
    }

    # Calculating pdf and cdf for the original data
    data.series.results$mledp_pdf<-ddoublepareto(x=data.series[,series],alpha.par=mle.results[5,1],beta.par=mle.results[6,1],t.par=mle.results[7,1],c.par=mle.results[8,1],m.par=mle.results[9,1])
    data.series.results$mledp_cdf<-pdoublepareto(x=data.series[,series],alpha.par=mle.results[5,1],beta.par=mle.results[6,1],t.par=mle.results[7,1],c.par=mle.results[8,1],m.par=mle.results[9,1])
    #plot(data.series[,series],data.series.results$mledp_pdf,log="xy")
    #plot(data.series[,series],data.series.results$mledp_cdf,log="xy")

    ### Kolmogorov-Smirnov test
    mle.sample.model.dp<-rdoublepareto(n=no.events, alpha.par=mle.results[5,1], beta.par=mle.results[6,1], t.par=mle.results[7,1],c.par=mle.results[8,1],m.par=mle.results[9,1]) 
    ks_model_greater<-NULL
    ks_model_greater<-ks.boot(data.series[,series],mle.sample.model.dp,nboots=ks_boot_samples,alternative="two.sided")
    ks_greater_result[8,2:5]<-c(ks_model_greater$nboots,ks_model_greater$ks$statistic,ks_model_greater$ks$p.value,ks_model_greater$ks.boot.pvalue)
    
#	if (Sys.info()[1] == "Linux") {x11()}
#	if (Sys.info()[1] == "Windows") {windows()}
#	if (Sys.info()[1] == "Darwin") {quartz()}
#   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""), xlab=xlabel,col="blue",xlim=range.plot.area,ylim=range.plot.density)  # Histogram raw
#   lines(data.series.synthetic,pdf.mle.data.series.dp,col="red",lty=3)
#   mtext(paste("MLE Double Pareto -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
#mtext(paste("Boostrapped KS test -> D: ",round(ks_greater_result[8,3],3),"; p-value: ",round(ks_greater_result[8,5],3),sep=""), side=3, col="navyblue", cex=0.8, line=0.5,padj=2.5)


   pdf(file = paste("MLE_Fit_DP_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
   par(cex.axis=0.9)
   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""), xlab=xlabel,col="blue",xlim=range.plot.area,ylim=range.plot.density)  # Histogram raw
   lines(data.series.synthetic,pdf.mle.data.series.dp,col="red",lty=3)
   mtext(paste("MLE Double Pareto -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
    mtext(paste("Boostrapped KS test -> D: ",round(ks_greater_result[8,3],3),"; p-value: ",round(ks_greater_result[8,5],3),sep=""), side=3, col="navyblue", cex=0.8, line=0.5,padj=2.5)
   dev.off()

#   ### Plot of DP cdf
#	if (Sys.info()[1] == "Linux") {x11()}
#	if (Sys.info()[1] == "Windows") {windows()}
#	if (Sys.info()[1] == "Darwin") {quartz()}
#   plot(data.series.synthetic,cdf.mle.data.series.dp, log="x", main="Cumulative distribution function", ylab=paste("Cumulative probability density",sep=""), xlab=xlabel,type="l",col="red",xlim=range.plot.area,ylim=c(0,1))
#   mtext(paste("CDF MLE Double Pareto -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)

   pdf(file = paste("CDF_MLE_DP_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
   plot(data.series.synthetic,cdf.mle.data.series.dp, log="x", main="Cumulative distribution function", ylab=paste("Cumulative probability density",sep=""), xlab=xlabel,type="l",col="red",xlim=range.plot.area,ylim=c(0,1))
   mtext(paste("CDF MLE Double Pareto -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
   dev.off()

   pdf(file = paste("CDF_MLE_DP_EntireRange_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
   data.series.synthetic.cdf<-c(10^log_cycles_range[1],sort((seq(1,10,0.1)%*%t(10^(log_cycles_range[1]:log_cycles_range[2])))[-1,]))
   cdf.mle.data.series.dp.entire.range<-pdoublepareto(x=data.series.synthetic.cdf, alpha.par=mle.results[5,1],beta.par=mle.results[6,1],t.par=mle.results[7,1],c.par=mle.results[8,1],m.par=mle.results[9,1])
   plot(data.series.synthetic,cdf.mle.data.series.dp, log="x", main="Cumulative distribution function", ylab=paste("Cumulative probability density",sep=""), xlab=xlabel,type="l",xlim=range.plot.area,ylim=c(0,1),xaxt="n",yaxt="n",col="red")
   #index.cdf.mle.data.fill.min<-which(data.series.synthetic.cdf < min(data.series.synthetic))
   index.cdf.mle.data.fill.min<-which(data.series.synthetic.cdf < min(data.series.synthetic) & data.series.synthetic.cdf >= mle.results[8,1]) # only for DP
   lines(data.series.synthetic.cdf[index.cdf.mle.data.fill.min],cdf.mle.data.series.dp.entire.range[index.cdf.mle.data.fill.min],lty="dashed",col="red")
   #index.cdf.mle.data.fill.max<-which(data.series.synthetic.cdf > max(data.series.synthetic))
   index.cdf.mle.data.fill.max<-which(data.series.synthetic.cdf > max(data.series.synthetic) & data.series.synthetic.cdf <= mle.results[9,1]) # only for DP
   lines(data.series.synthetic.cdf[index.cdf.mle.data.fill.max],cdf.mle.data.series.dp.entire.range[index.cdf.mle.data.fill.max],lty="dashed",col="red")
   cdf.percentiles.sizes.qdoublepareto<-qdoublepareto(x=cdf_percentiles,alpha.par=mle.results[5,1],beta.par=mle.results[6,1],t.par=mle.results[7,1],c.par=mle.results[8,1],m.par=mle.results[9,1],lower.par=uniroot.range[1]/100, upper.par=uniroot.range[2]*100,FUN=FALSE)
   points(cdf.percentiles.sizes.qdoublepareto,cdf_percentiles,pch=21,bg="red",col="black")
   lines(cdf.percentiles.sizes.qdoublepareto,cdf_percentiles,type="h",lty="dotted",col="red")
   text(cdf.percentiles.sizes.qdoublepareto, cdf_percentiles, labels=paste("(",signif(cdf.percentiles.sizes.qdoublepareto,2),",",cdf_percentiles,")",sep=""), cex= 0.6,pos=3)
   axis(side=1,at=at.x,label=sciNotation(at.x,1))
   at.y.cdf<-seq(0,1,0.1)
   at.y.tick.cdf<-seq(0,1,0.01)
   axis(side=2,at=at.y.cdf,label=at.y.cdf)
   axis(1, at.x.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
   axis(2, at.y.tick.cdf, labels = NA, lty = 1, lwd = 1, tck = -0.01)
   axis(3, at.x.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
   axis(4, at.y.tick.cdf, labels = NA, lty = 1, lwd = 1, tck = -0.01)
   mtext(paste("CDF MLE Double Pareto -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
   dev.off()
   
   
   #### MLE using Inverse gamma
   # Default paramenter value, specified_by_user == "NO" in configuration.txt file
   alpha.range<-c(min(c(hde.results[11,1],kde.results[11,1]))/1.5,max(c(hde.results[11,1],kde.results[11,1]))*1.5)
   eta.range<-c(min(c(hde.results[12,1],kde.results[12,1]))/3,max(c(hde.results[12,1],kde.results[12,1]))*3)
   lambda.range<-c(min(c(hde.results[13,1],kde.results[13,1]))/3,max(c(hde.results[13,1],kde.results[13,1]))*3)
   # Paramenter value specified by the user, specified_by_user == "YES" in configuration.txt file
   if (configuration[25,4]=="YES") (alpha.range<-as.numeric(configuration[25,2:3]))
   if (configuration[26,4]=="YES") (eta.range<-as.numeric(configuration[26,2:3]))
   if (configuration[27,4]=="YES") (lambda.range<-as.numeric(configuration[27,2:3]))
   
   print(paste("MLE (IG):",names(data.series)[series]))

   mle.data.series.ig<-mle2(mll.inversegamma,method="L-BFGS-B",lower=c(alpha.par=alpha.range[1], eta.par=eta.range[1], lambda.par=lambda.range[1]),upper=c(alpha.par=alpha.range[2], eta.par=eta.range[2], lambda.par=lambda.range[2]),start=list(alpha.par=sum(alpha.range)/2, eta.par=sum(eta.range)/2, lambda.par=sum(lambda.range)/2),data=list(x=na.omit(data.series[,series])),control=list(maxit=1000000,lmm=5,factr=10^-8,ndeps=c(0.00001,0.00001,0.00001),pgtol=0,trace=6,parscale=c(alpha.par=1, eta.par=1,lambda.par=1)))
   summary(mle.data.series.ig)

#   # SANN optimization method
#   mle.data.series.ig<-mle2(mll.inversegamma,method="SANN",start=list(alpha.par=sum(alpha.range)/2, eta.par=sum(eta.range)/2, lambda.par=sum(lambda.range)/2),data=list(x=na.omit(data.series.sample)),control=list(maxit=50000,trace=6))
#   summary(mle.data.series.ig)
#
#   # Nelder-Mead optimization method
#   mle.data.series.ig<-mle2(mll.inversegamma,method="Nelder-Mead",start=list(alpha.par=1, eta.par=100, lambda.par=100),data=list(x=na.omit(data.series.sample)),control=list(maxit=100000,trace=6))
#   summary(mle.data.series.ig)
#

   # PDF generated using parameters estimated by MLE using a Simplified Double Pareto distribution (dps)
   pdf.mle.data.series.ig<-dinversegamma(x=data.series.synthetic, alpha.par=mle.data.series.ig@coef[1], eta.par=mle.data.series.ig@coef[2], lambda.par=mle.data.series.ig@coef[3])
   cdf.mle.data.series.ig<-pinversegamma(x=data.series.synthetic, alpha.par=mle.data.series.ig@coef[1], eta.par=mle.data.series.ig@coef[2], lambda.par=mle.data.series.ig@coef[3])
#plot(data.series.synthetic,cdf.mle.data.series.ig,type="l")

   mle.results[11:13,1:4]<-coef(summary(mle.data.series.ig))

   # Calculating IG rollover
   mle.index.rollover.ig<-NULL
   mle.index.rollover.ig<-which(pdf.mle.data.series.ig==max(pdf.mle.data.series.ig))
   if (length(mle.index.rollover.ig)>1) {mle.index.rollover.ig<-mle.index.rollover.ig[1]}
   if (mle.index.rollover.ig>1)
    {
    mle.results[14,1]<-data.series.synthetic[mle.index.rollover.ig]
    }

    # Calculating pdf and cdf for the original data
    data.series.results$mleig_pdf<-dinversegamma(x=data.series[,series],alpha.par=mle.results[11,1],eta.par=mle.results[12,1],lambda.par=mle.results[13,1])
    data.series.results$mleig_cdf<-pinversegamma(x=data.series[,series],alpha.par=mle.results[11,1],eta.par=mle.results[12,1],lambda.par=mle.results[13,1])
    #plot(data.series[,series],data.series.results$mleig_pdf,log="xy")
    #plot(data.series[,series],data.series.results$mleig_cdf,log="xy")

    ### Kolmogorov-Smirnov test
    mle.sample.model.ig<-rinversegamma(n=no.events, alpha.par=mle.results[11,1], eta.par=mle.results[12,1], lambda.par=mle.results[13,1])
    ks_model_greater<-NULL
    ks_model_greater<-ks.boot(data.series[,series],mle.sample.model.ig,nboots=ks_boot_samples,alternative="two.sided")
    ks_greater_result[9,2:5]<-c(ks_model_greater$nboots,ks_model_greater$ks$statistic,ks_model_greater$ks$p.value,ks_model_greater$ks.boot.pvalue)

#	if (Sys.info()[1] == "Linux") {x11()}
#	if (Sys.info()[1] == "Windows") {windows()}
#	if (Sys.info()[1] == "Darwin") {quartz()}
#   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""), xlab=xlabel,col="blue",xlim=range.plot.area,ylim=range.plot.density)  # Histogram raw
#   lines(data.series.synthetic,pdf.mle.data.series.ig,col="violet")
#   mtext(paste("MLE Inverse Gamma -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
#mtext(paste("Boostrapped KS test -> D: ",round(ks_greater_result[9,3],3),"; p-value: ",round(ks_greater_result[9,5],3),sep=""), side=3, col="navyblue", cex=0.8, line=0.5,padj=2.5)


   pdf(file = paste("MLE_Fit_IG_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""), xlab=xlabel,col="blue",xlim=range.plot.area,ylim=range.plot.density)  # Histogram raw
   lines(data.series.synthetic,pdf.mle.data.series.ig,col="violet")
   mtext(paste("MLE Inverse Gamma -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
    mtext(paste("Boostrapped KS test -> D: ",round(ks_greater_result[9,3],3),"; p-value: ",round(ks_greater_result[9,5],3),sep=""), side=3, col="navyblue", cex=0.8, line=0.5,padj=2.5)
   dev.off()

#   ### Plot of IG cdf
#	if (Sys.info()[1] == "Linux") {x11()}
#	if (Sys.info()[1] == "Windows") {windows()}
#	if (Sys.info()[1] == "Darwin") {quartz()}
#   plot(data.series.synthetic,cdf.mle.data.series.ig, log="x", main="Cumulative distribution function", ylab=paste("Cumulative probability density",sep=""), xlab=xlabel,type="l",col="violet",xlim=range.plot.area,ylim=c(0,1))
#   mtext(paste("CDF MLE Inverse Gamma -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)

   pdf(file = paste("CDF_MLE_IG_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
   plot(data.series.synthetic,cdf.mle.data.series.ig, log="x", main="Cumulative distribution function", ylab=paste("Cumulative probability density",sep=""), xlab=xlabel,type="l",col="violet",xlim=range.plot.area,ylim=c(0,1))
   mtext(paste("CDF MLE Inverse Gamma -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
   dev.off()
   
   
   pdf(file = paste("CDF_MLE_IG_EntireRange_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
   data.series.synthetic.cdf<-c(10^log_cycles_range[1],sort((seq(1,10,0.1)%*%t(10^(log_cycles_range[1]:log_cycles_range[2])))[-1,]))
   cdf.mle.data.series.ig.entire.range<-pinversegamma(x=data.series.synthetic.cdf, alpha.par=mle.results[11,1], eta.par=mle.results[12,1], lambda.par=mle.results[13,1])
   
   plot(data.series.synthetic,cdf.mle.data.series.ig, log="x", main="Cumulative distribution function", ylab=paste("Cumulative probability density",sep=""), xlab=xlabel,type="l",xlim=range.plot.area,ylim=c(0,1),xaxt="n",yaxt="n",col="violet")
   index.cdf.mle.data.fill.min<-which(data.series.synthetic.cdf < min(data.series.synthetic))
   lines(data.series.synthetic.cdf[index.cdf.mle.data.fill.min],cdf.mle.data.series.ig.entire.range[index.cdf.mle.data.fill.min],lty="dashed",col="violet")
   index.cdf.mle.data.fill.max<-which(data.series.synthetic.cdf > max(data.series.synthetic))
   lines(data.series.synthetic.cdf[index.cdf.mle.data.fill.max],cdf.mle.data.series.ig.entire.range[index.cdf.mle.data.fill.max],lty="dashed",col="violet")
   cdf.percentiles.sizes.qinversegamma<-qinversegamma(x=cdf_percentiles, alpha.par=mle.results[11,1], eta.par=mle.results[12,1], lambda.par=mle.results[13,1],lower.par=uniroot.range[1]/10, upper.par=uniroot.range[2]*10,FUN=FALSE)
   points(cdf.percentiles.sizes.qinversegamma,cdf_percentiles,pch=21,bg="violet",col="black")
   lines(cdf.percentiles.sizes.qinversegamma,cdf_percentiles,type="h",lty="dotted",col="violet")
   text(cdf.percentiles.sizes.qinversegamma, cdf_percentiles, labels=paste("(",signif(cdf.percentiles.sizes.qinversegamma,2),",",cdf_percentiles,")",sep=""), cex= 0.6,pos=3)
   axis(side=1,at=at.x,label=sciNotation(at.x,1))
   at.y.cdf<-seq(0,1,0.1)
   at.y.tick.cdf<-seq(0,1,0.01)
   axis(side=2,at=at.y.cdf,label=at.y.cdf)
   axis(1, at.x.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
   axis(2, at.y.tick.cdf, labels = NA, lty = 1, lwd = 1, tck = -0.01)
   axis(3, at.x.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
   axis(4, at.y.tick.cdf, labels = NA, lty = 1, lwd = 1, tck = -0.01)
   mtext(paste("CDF MLE Inverse Gamma -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
   dev.off()
   

   # Plot comparison MLE different distribution
	if (Sys.info()[1] == "Linux") {x11()}
	if (Sys.info()[1] == "Windows") {windows()}
	if (Sys.info()[1] == "Darwin") {quartz()}
   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""), xlab=xlabel,col="blue",xlim=range.plot.area,ylim=range.plot.density)  # Histogram raw
   lines(data.series.synthetic,pdf.mle.data.series.dps,col="dark green",lty=3)
   lines(data.series.synthetic,pdf.mle.data.series.dp,col="red",lty=3)
   lines(data.series.synthetic,pdf.mle.data.series.ig,col="violet",lty=3)
   mtext(paste("MLE Comparison -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
   legend("bottomleft",legend=c("Histogram raw data","Double Pareto simplified","Double Pareto","Inverse Gamma"),pch=c(1,NA_integer_,NA_integer_,NA_integer_),lty=c(FALSE,3,3,3),lwd=1,col=c("blue","dark green","red","violet"),cex=0.8,box.lty=3,box.col="black")

   pdf(file = paste("MLE_Fit_Comparison_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
   par(cex.axis=0.9)
   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""), xlab=xlabel,col="blue",xlim=range.plot.area,ylim=range.plot.density)  # Histogram raw
   lines(data.series.synthetic,pdf.mle.data.series.dps,col="dark green",lty=3)
   lines(data.series.synthetic,pdf.mle.data.series.dp,col="red",lty=3)
   lines(data.series.synthetic,pdf.mle.data.series.ig,col="violet",lty=3)
   mtext(paste("MLE Comparison -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
   legend("bottomleft",legend=c("Histogram raw data","Double Pareto simplified","Double Pareto","Inverse Gamma"),pch=c(1,NA_integer_,NA_integer_,NA_integer_),lty=c(FALSE,3,3,3),lwd=1,col=c("blue","dark green","red","violet"),cex=0.8,box.lty=3,box.col="black")
   dev.off()


   # CDF plot comparison MLE different distribution
	if (Sys.info()[1] == "Linux") {x11()}
	if (Sys.info()[1] == "Windows") {windows()}
	if (Sys.info()[1] == "Darwin") {quartz()}
   plot(data.series.synthetic,cdf.mle.data.series.dps, log="x", main="Cumulative distribution function", ylab=paste("Cumulative probability density",sep=""), xlab=xlabel,type="l",col="dark green",xlim=range.plot.area,ylim=c(0,1))
   lines(data.series.synthetic,cdf.mle.data.series.dp,col="red")
   lines(data.series.synthetic,cdf.mle.data.series.ig,col="violet")
   mtext(paste("CDF MLE Comparison -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
   legend("topleft",legend=c("Double Pareto simplified","Double Pareto","Inverse Gamma"),pch=c(NA_integer_,NA_integer_,NA_integer_),lty=c(1,1,1),lwd=1,col=c("dark green","red","violet"),cex=0.8,box.lty=3,box.col="black")


   pdf(file = paste("CDF_MLE_Fit_Comparison_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
   plot(data.series.synthetic,cdf.mle.data.series.dps, log="x", main="Cumulative distribution function", ylab=paste("Cumulative probability density",sep=""), xlab=xlabel,type="l",col="dark green",xlim=range.plot.area,ylim=c(0,1))
   lines(data.series.synthetic,cdf.mle.data.series.dp,col="red")
   lines(data.series.synthetic,cdf.mle.data.series.ig,col="violet")
   mtext(paste("CDF MLE Comparison -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
   legend("topleft",legend=c("Double Pareto simplified","Double Pareto","Inverse Gamma"),pch=c(NA_integer_,NA_integer_,NA_integer_),lty=c(1,1,1),lwd=1,col=c("dark green","red","violet"),cex=0.8,box.lty=3,box.col="black")
   dev.off()

  #### -------------- Uncertainty estimation -------------- ####

	### Creating object for storing sampling results
	mle_dps_sampling_result<-matrix(NA,nrow=length(data.series.synthetic.sample),ncol=samples_number)
	mle_dps_sampling_result_coef<-matrix(NA,nrow=samples_number,ncol=4)
	colnames(mle_dps_sampling_result_coef)<-c("alpha","beta","t","rollov")
  
	mle_dp_sampling_result<-matrix(NA,nrow=length(data.series.synthetic.sample),ncol=samples_number)
	mle_dp_sampling_result_coef<-matrix(NA,nrow=samples_number,ncol=6)
	colnames(mle_dp_sampling_result_coef)<-c("alpha","beta","t","c","m","rollov")
 
	mle_ig_sampling_result<-matrix(NA,nrow=length(data.series.synthetic.sample),ncol=samples_number)
	mle_ig_sampling_result_coef<-matrix(NA,nrow=samples_number,ncol=4)
	colnames(mle_ig_sampling_result_coef)<-c("alpha","eta","lambda","rollov")

	for (countsamples in 1:samples_number)
		{
		#countsamples<-1
		c.par.value.sample<-c.par.value
		m.par.value.sample<-m.par.value
		print(paste("MLE -",names(data.series)[series],"-> Sample:",countsamples))
		data.series.sample<-sample(data.series[,series],replace=TRUE)
		no.events.sample<-length(data.series.sample)
	   
		### ----------------- Fit Double Pareto Simplified (DPS) ----------------- ###
		# Default paramenter value, specified_by_user == "NO" in configuration.txt file
		alpha.range<-c(min(c(hde.results[1,1],kde.results[1,1]))/1.5,max(c(hde.results[1,1],kde.results[1,1]))*1.5)
		beta.range<-c(min(c(hde.results[2,1],kde.results[2,1]))/1.5,max(c(hde.results[2,1],kde.results[2,1]))*1.5)
		t.range<-c(min(c(hde.results[3,1],kde.results[2,1]))/1.5,max(c(hde.results[3,1],kde.results[3,1]))*1.5)
		# Paramenter value specified by the user, specified_by_user == "YES" in configuration.txt file
		if (configuration[19,4]=="YES") (alpha.range<-as.numeric(configuration[19,2:3]))
		if (configuration[20,4]=="YES") (beta.range<-as.numeric(configuration[20,2:3]))
		if (configuration[21,4]=="YES") (t.range<-as.numeric(configuration[21,2:3]))
   
		print(paste("MLE (DPS):",names(data.series)[series],"-> Sample:",countsamples))
		if (inherits(try(mle.data.series.dps.sample<-mle2(mll.doublepareto.simplified,method="L-BFGS-B",trace=FALSE,lower=c(alpha.par=alpha.range[1], beta.par=beta.range[1],t.par=t.range[1]),upper=c(alpha.par=alpha.range[2], beta.par=beta.range[2],t.par=t.range[2]),start=list(alpha.par=sum(alpha.range)/2, beta.par=sum(beta.range)/2,t.par=sum(t.range)/2),data=list(x=na.omit(data.series.sample)),control=list(maxit=1000000,lmm=5,factr=10^-8,pgtol=0,trace=6,parscale=c(alpha.par=1, beta.par=1,t.par=1))),silent=TRUE),what="try-error"))
			{
			next()
			}	
		# Storing fitting results
		mle_dps_sampling_result_coef[countsamples,1:3]<-coef(summary(mle.data.series.dps.sample))[,1]
		# PDF generated using parameters estimated by MLE using a Simplified Double Pareto distribution (dps)
		mle_dps_sampling_result[,countsamples]<-ddoublepareto.simplified(x=data.series.synthetic.sample, alpha.par=mle.data.series.dps.sample@coef[1], beta.par=mle.data.series.dps.sample@coef[2], t.par=mle.data.series.dps.sample@coef[3])	
		
		#Calculating DPS rollover
		mle.index.rollover.dps.sample<-NULL
		mle.index.rollover.dps.sample<-which(10^mle_dps_sampling_result[,countsamples]==max(10^mle_dps_sampling_result[,countsamples]))
		if (length(mle.index.rollover.dps.sample)>1) {mle.index.rollover.dps.sample<-mle.index.rollover.dps.sample[1]}
		if (mle.index.rollover.dps.sample>1)
			{
			mle_dps_sampling_result_coef[countsamples,4]<-data.series.synthetic.sample[mle.index.rollover.dps.sample]
			}

		### ----------------- Fit Double Pareto (DP) ----------------- ###
		# Default paramenter value, specified_by_user == "NO" in configuration.txt file
		alpha.range<-c(min(c(hde.results[5,1],kde.results[5,1]))/2,max(c(hde.results[5,1],kde.results[5,1]))*2)
		beta.range<-c(min(c(hde.results[6,1],kde.results[6,1]))/4,max(c(hde.results[6,1],kde.results[6,1]))*4)
		t.range<-c(min(c(hde.results[7,1],kde.results[7,1]))/2,max(c(hde.results[7,1],kde.results[7,1]))*2)
		# Paramenter value specified by the user, specified_by_user == "YES" in configuration.txt file
		if (configuration[22,4]=="YES") (alpha.range<-as.numeric(configuration[22,2:3]))
		if (configuration[23,4]=="YES") (beta.range<-as.numeric(configuration[23,2:3]))
		if (configuration[24,4]=="YES") (t.range<-as.numeric(configuration[24,2:3]))

		print(paste("MLE (DP):",names(data.series)[series],"-> Sample:",countsamples))
		if (inherits(try(mle.data.series.dp.sample<-mle2(mll.doublepareto,method="L-BFGS-B",trace=FALSE,fixed=list(c.par=c.par.value, m.par=m.par.value), lower=c(alpha.par=alpha.range[1], beta.par=beta.range[1],t.par=t.range[1]),upper=c(alpha.par=alpha.range[2], beta.par=beta.range[2],t.par=t.range[2]),start=list(alpha.par=sum(alpha.range)/2, beta.par=sum(beta.range)/2,t.par=sum(t.range)/2),data=list(x=na.omit(data.series.sample)),control=list(maxit=1000000,lmm=5,factr=10^-8,pgtol=0,trace=6,parscale=c(alpha.par=1, beta.par=1,t.par=1))),silent=TRUE),what="try-error"))
			{
			next()
			}	
		# Storing fitting results
		mle_dp_sampling_result_coef[countsamples,1:3]<-coef(summary(mle.data.series.dp.sample))[,1]
		# PDF generated using parameters estimated by MLE using a Double Pareto distribution (dp)
		mle_dp_sampling_result[,countsamples]<-ddoublepareto(x=data.series.synthetic.sample, alpha.par=mle.data.series.dp.sample@coef[1], beta.par=mle.data.series.dp.sample@coef[2], t.par=mle.data.series.dp.sample@coef[3], c.par=c.par.value.sample, m.par=m.par.value.sample)
		mle_dp_sampling_result_coef[countsamples,4:5]<-c(c.par.value.sample,m.par.value.sample)
		
		#Calculating DPS rollover
		mle.index.rollover.dp.sample<-NULL
		mle.index.rollover.dp.sample<-which(10^mle_dp_sampling_result[,countsamples]==max(10^mle_dp_sampling_result[,countsamples]))
		if (length(mle.index.rollover.dp.sample)>1) {mle.index.rollover.dp.sample<-mle.index.rollover.dp.sample[1]}
		if (mle.index.rollover.dp.sample>1)
			{
			mle_dp_sampling_result_coef[countsamples,6]<-data.series.synthetic.sample[mle.index.rollover.dp.sample]
			}

		### ----------------- Fit Inverse Gamma (IG) ----------------- ###
		# Default paramenter value, specified_by_user == "NO" in configuration.txt file
		alpha.range<-c(min(c(hde.results[11,1],kde.results[11,1]))/1.5,max(c(hde.results[11,1],kde.results[11,1]))*1.5)
		eta.range<-c(min(c(hde.results[12,1],kde.results[12,1]))/3,max(c(hde.results[12,1],kde.results[12,1]))*3)
		lambda.range<-c(min(c(hde.results[13,1],kde.results[13,1]))/3,max(c(hde.results[13,1],kde.results[13,1]))*3)
		# Paramenter value specified by the user, specified_by_user == "YES" in configuration.txt file
		if (configuration[25,4]=="YES") (alpha.range<-as.numeric(configuration[25,2:3]))
		if (configuration[26,4]=="YES") (eta.range<-as.numeric(configuration[26,2:3]))
		if (configuration[27,4]=="YES") (lambda.range<-as.numeric(configuration[27,2:3]))

		print(paste("MLE (IG):",names(data.series)[series],"-> Sample:",countsamples))
		if (inherits(try(mle.data.series.ig.sample<-mle2(mll.inversegamma,method="L-BFGS-B",trace=FALSE,lower=c(alpha.par=alpha.range[1], eta.par=eta.range[1], lambda.par=lambda.range[1]),upper=c(alpha.par=alpha.range[2], eta.par=eta.range[2], lambda.par=lambda.range[2]),start=list(alpha.par=sum(alpha.range)/2, eta.par=sum(eta.range)/2, lambda.par=sum(lambda.range)/2),data=list(x=na.omit(data.series.sample)),control=list(maxit=1000000,lmm=5,factr=10^-8,ndeps=c(0.00001,0.00001,0.00001),pgtol=0,trace=6,parscale=c(alpha.par=1, eta.par=1,lambda.par=1))),silent=TRUE),what="try-error"))
			{
			next()
			}	
	
		# Storing fitting results
		mle_ig_sampling_result_coef[countsamples,1:3]<-coef(summary(mle.data.series.ig.sample))[,1]
		# PDF generated using parameters estimated by MLE using a Double Pareto distribution (ig)
		mle_ig_sampling_result[,countsamples]<-dinversegamma(x=data.series.synthetic.sample, alpha.par=mle.data.series.ig.sample@coef[1], eta.par=mle.data.series.ig.sample@coef[2], lambda.par=mle.data.series.ig.sample@coef[3])
		#Calculating DPS rollover
		mle.index.rollover.ig.sample<-NULL
		mle.index.rollover.ig.sample<-which(10^mle_ig_sampling_result[,countsamples]==max(10^mle_ig_sampling_result[,countsamples]))
		if (length(mle.index.rollover.ig.sample)>1) {mle.index.rollover.ig.sample<-mle.index.rollover.ig.sample[1]}
		if (mle.index.rollover.ig.sample>1)
			{
			mle_ig_sampling_result_coef[countsamples,4]<-data.series.synthetic.sample[mle.index.rollover.ig.sample]
			}
		}	

		quantiles_mle_dps_sampling_result_coef<-apply(mle_dps_sampling_result_coef,MARGIN=2,FUN=quantile_fun)
		quantiles_mle_dp_sampling_result_coef<-apply(mle_dp_sampling_result_coef,MARGIN=2,FUN=quantile_fun)
		quantiles_mle_ig_sampling_result_coef<-apply(mle_ig_sampling_result_coef,MARGIN=2,FUN=quantile_fun)
		 
		write.table(data.frame(percentile=rownames(quantiles_mle_dps_sampling_result_coef),round(quantiles_mle_dps_sampling_result_coef,2)),paste("MLE_Fit_DPS_parameter_uncertainty_",names(data.series)[series],".txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t")
		write.table(data.frame(percentile=rownames(quantiles_mle_dp_sampling_result_coef),round(quantiles_mle_dp_sampling_result_coef,2)),paste("MLE_Fit_DP_parameter_uncertainty_",names(data.series)[series],".txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t")
		write.table(data.frame(percentile=rownames(quantiles_mle_ig_sampling_result_coef),round(quantiles_mle_ig_sampling_result_coef,2)),paste("MLE_Fit_IG_parameter_uncertainty_",names(data.series)[series],".txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t")
		 
		quantiles_mle_dps_sampling_result<-apply(mle_dps_sampling_result,MARGIN=1,FUN=quantile_fun)
		quantiles_mle_dp_sampling_result<-apply(mle_dp_sampling_result,MARGIN=1,FUN=quantile_fun)
		quantiles_mle_ig_sampling_result<-apply(mle_ig_sampling_result,MARGIN=1,FUN=quantile_fun)
   
 
		### Plot MLE estimation with uncertainty
		 
		pdf(file = paste("MLE_Fit_DPS_uncertainty_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
		par(cex.axis=0.9)
		#if (Sys.info()[1] == "Linux") {x11()}
		#if (Sys.info()[1] == "Windows") {windows()}
		#if (Sys.info()[1] == "Darwin") {quartz()}
		#xlabel<-expression(A (m^{2}))
		ylabel<-paste("Probability density",sep="")
		plot(NA, NA, log="xy", main="Probability densities", ylab=ylabel, xlab=xlabel,xlim=range.plot.area,ylim=range.plot.density,xaxt="n",yaxt="n")
		axis(side=1,at=at.x,label=sciNotation(at.x,1))
		axis(side=2,at=at.y,label=sciNotation(at.y,1))
		axis(1, at.x.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		axis(2, at.y.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		axis(3, at.x.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		axis(4, at.y.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		for (count_q in 1:floor((length(prob_series))/2))
			{
			#count_q<-1	
		polygon(x=c(sort(data.series.synthetic.sample),rev(sort(data.series.synthetic.sample))),y=c(quantiles_mle_dps_sampling_result[length(prob_series)+1-count_q,],rev(quantiles_mle_dps_sampling_result[count_q,])),border=NA,col=color_vector_pol[count_q])
		}
		lines(linear.x.data.series, kde.data.series.linear, col="orange",lty=1,lwd=1)
		lines(data.series.synthetic,pdf.mle.data.series.dps,col="darkred",lty=1,lwd=2)
		#lines(data.series.synthetic.sample,10^quantiles_mle_dps_sampling_result["50%",],col="black",lty="dashed",lwd=1)
    lines(data.series.synthetic.sample,quantiles_mle_dps_sampling_result["50%",],col="black",lty="dashed",lwd=1)
		points(linear.mids.data.series, prob.linear.data.series,pch=20,cex=1,col="navyblue")
		mtext(paste("MLE Double Pareto Simplified -> ",names(data.series)[series],sep=""), side=3, col="darkred", cex=0.8, line=0.5)
    mtext(paste("Boostrapped KS test -> D: ",round(ks_greater_result[7,3],3),"; p-value: ",round(ks_greater_result[7,5],3),sep=""), side=3, col="navyblue", cex=0.8, line=0.5,padj=2.5)
		legend("bottomleft",c("Histogram empirical data","Raw kernel density data","Distribution fit","Median values",paste("Quantiles: ",prob_series[1:floor((length(prob_series))/2)],"-",prob_series[(length(prob_series)+1)-1:floor((length(prob_series))/2)],sep="")),pt.cex=1,pch=c(20,NA,NA,NA,rep(NA,floor((length(prob_series))/2))),lty=c(NA,1,1,2,rep(1,floor((length(prob_series))/2))),lwd=c(NA,1,2,1,rep(6,floor((length(prob_series))/2))),col=c("navyblue","orange","darkred","black",color_vector_pol),bty="o",bg="transparent",box.col="transparent",cex=0.9)
		dev.off()
   		  
		pdf(file = paste("MLE_Fit_DP_uncertainty_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
		par(cex.axis=0.9)
		#if (Sys.info()[1] == "Linux") {x11()}
		#if (Sys.info()[1] == "Windows") {windows()}
		#if (Sys.info()[1] == "Darwin") {quartz()}
		#xlabel<-expression(A (m^{2}))
		ylabel<-paste("Probability density",sep="")
		plot(NA, NA, log="xy", main="Probability densities", ylab=ylabel, xlab=xlabel,xlim=range.plot.area,ylim=range.plot.density,xaxt="n",yaxt="n")
		axis(side=1,at=at.x,label=sciNotation(at.x,1))
		axis(side=2,at=at.y,label=sciNotation(at.y,1))
		axis(1, at.x.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		axis(2, at.y.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		axis(3, at.x.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		axis(4, at.y.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		for (count_q in 1:floor((length(prob_series))/2))
		{
		#count_q<-1	
		polygon(x=c(sort(data.series.synthetic.sample),rev(sort(data.series.synthetic.sample))),y=c(quantiles_mle_dp_sampling_result[length(prob_series)+1-count_q,],rev(quantiles_mle_dp_sampling_result[count_q,])),border=NA,col=color_vector_pol[count_q])
		}
		lines(linear.x.data.series, kde.data.series.linear, col="orange",lty=1,lwd=1)
		lines(data.series.synthetic,pdf.mle.data.series.dp,col="darkred",lty=1,lwd=2)
		#lines(data.series.synthetic.sample,10^quantiles_mle_dp_sampling_result["50%",],col="black",lty="dashed",lwd=1)
    lines(data.series.synthetic.sample,quantiles_mle_dp_sampling_result["50%",],col="black",lty="dashed",lwd=1)
		points(linear.mids.data.series, prob.linear.data.series,pch=20,cex=1,col="navyblue")
		mtext(paste("MLE Double Pareto -> ",names(data.series)[series],sep=""), side=3, col="darkred", cex=0.8, line=0.5)
    mtext(paste("Boostrapped KS test -> D: ",round(ks_greater_result[8,3],3),"; p-value: ",round(ks_greater_result[8,5],3),sep=""), side=3, col="navyblue", cex=0.8, line=0.5,padj=2.5)
		legend("bottomleft",c("Histogram empirical data","Raw kernel density data","Distribution fit","Median values",paste("Quantiles: ",prob_series[1:floor((length(prob_series))/2)],"-",prob_series[(length(prob_series)+1)-1:floor((length(prob_series))/2)],sep="")),pt.cex=1,pch=c(20,NA,NA,NA,rep(NA,floor((length(prob_series))/2))),lty=c(NA,1,1,2,rep(1,floor((length(prob_series))/2))),lwd=c(NA,1,2,1,rep(6,floor((length(prob_series))/2))),col=c("navyblue","orange","darkred","black",color_vector_pol),bty="o",bg="transparent",box.col="transparent",cex=0.9)
		dev.off()
		  
		pdf(file = paste("MLE_Fit_IG_uncertainty_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
		par(cex.axis=0.9)
		#if (Sys.info()[1] == "Linux") {x11()}
		#if (Sys.info()[1] == "Windows") {windows()}
		#if (Sys.info()[1] == "Darwin") {quartz()}
		#xlabel<-expression(A (m^{2}))
		ylabel<-paste("Probability density",sep="")
		plot(NA, NA, log="xy", main="Probability densities", ylab=ylabel, xlab=xlabel,xlim=range.plot.area,ylim=range.plot.density,xaxt="n",yaxt="n")
		axis(side=1,at=at.x,label=sciNotation(at.x,1))
		axis(side=2,at=at.y,label=sciNotation(at.y,1))
		axis(1, at.x.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		axis(2, at.y.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		axis(3, at.x.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		axis(4, at.y.tick, labels = NA, lty = 1, lwd = 1, tck = -0.01)
		for (count_q in 1:floor((length(prob_series))/2))
		{
		#count_q<-1	
		polygon(x=c(sort(data.series.synthetic.sample),rev(sort(data.series.synthetic.sample))),y=c(quantiles_mle_ig_sampling_result[length(prob_series)+1-count_q,],rev(quantiles_mle_ig_sampling_result[count_q,])),border=NA,col=color_vector_pol[count_q])
		}
		lines(linear.x.data.series, kde.data.series.linear, col="orange",lty=1,lwd=1)
		lines(data.series.synthetic,pdf.mle.data.series.ig,col="darkred",lty=1,lwd=2)
		#lines(data.series.synthetic.sample,10^quantiles_mle_ig_sampling_result["50%",],col="black",lty="dashed",lwd=1)
    lines(data.series.synthetic.sample,quantiles_mle_ig_sampling_result["50%",],col="black",lty="dashed",lwd=1)
		points(linear.mids.data.series, prob.linear.data.series,pch=20,cex=1,col="navyblue")
		mtext(paste("MLE Inverse Gamma -> ",names(data.series)[series],sep=""), side=3, col="darkred", cex=0.8, line=0.5)
    mtext(paste("Boostrapped KS test -> D: ",round(ks_greater_result[9,3],3),"; p-value: ",round(ks_greater_result[9,5],3),sep=""), side=3, col="navyblue", cex=0.8, line=0.5,padj=2.5)
		legend("bottomleft",c("Histogram empirical data","Raw kernel density data","Distribution fit","Median values",paste("Quantiles: ",prob_series[1:floor((length(prob_series))/2)],"-",prob_series[(length(prob_series)+1)-1:floor((length(prob_series))/2)],sep="")),pt.cex=1,pch=c(20,NA,NA,NA,rep(NA,floor((length(prob_series))/2))),lty=c(NA,1,1,2,rep(1,floor((length(prob_series))/2))),lwd=c(NA,1,2,1,rep(6,floor((length(prob_series))/2))),col=c("navyblue","orange","darkred","black",color_vector_pol),bty="o",bg="transparent",box.col="transparent",cex=0.9)
		dev.off()
   
# ------------------------ HDE, KDE, MLE Comparison  ------------------------ #

   # Plot comparison DPS
	if (Sys.info()[1] == "Linux") {x11()}
	if (Sys.info()[1] == "Windows") {windows()}
	if (Sys.info()[1] == "Darwin") {quartz()}
   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""), xlab=xlabel,col="blue",xlim=range.plot.area,ylim=range.plot.density)  # Histogram raw
   lines(linear.x.data.series, kde.data.series.linear, col="orange",lty=1,lwd=1)
   lines(data.series.synthetic,10^value.fit.hde.log.doublepareto.simplified,col="dark green",lty=1)
   lines(data.series.synthetic,10^value.fit.kde.log.doublepareto.simplified,col="dark green",lty=2)
   lines(data.series.synthetic,pdf.mle.data.series.dps,col="dark green",lty=3)
   mtext(paste("Double Pareto Simplified Comparison -> ",names(data.series)[series],sep=""), side=3, col="dark green", cex=0.8, line=0.5)
   legend("bottomleft",legend=c("Histogram raw data","Kernel density raw data","HDE","KDE","MLE"),pch=c(1,NA_integer_,NA_integer_,NA_integer_,NA_integer_),lty=c(FALSE,1,1,2,3),lwd=1,col=c("blue","orange","dark green","dark green","dark green"),cex=0.8,box.lty=3,box.col="black")

   pdf(file = paste("Model_Comparison_DPS_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
   par(cex.axis=0.9)
   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""), xlab=xlabel,col="blue",xlim=range.plot.area,ylim=range.plot.density)  # Histogram raw
   lines(linear.x.data.series, kde.data.series.linear, col="orange",lty=1,lwd=1)
   lines(data.series.synthetic,10^value.fit.hde.log.doublepareto.simplified,col="dark green",lty=1)
   lines(data.series.synthetic,10^value.fit.kde.log.doublepareto.simplified,col="dark green",lty=2)
   lines(data.series.synthetic,pdf.mle.data.series.dps,col="dark green",lty=3)
   mtext(paste("Double Pareto Simplified Comparison -> ",names(data.series)[series],sep=""), side=3, col="dark green", cex=0.8, line=0.5)
   legend("bottomleft",legend=c("Histogram raw data","Kernel density raw data","HDE","KDE","MLE"),pch=c(1,NA_integer_,NA_integer_,NA_integer_,NA_integer_),lty=c(FALSE,1,1,2,3),lwd=1,col=c("blue","orange","dark green","dark green","dark green"),cex=0.8,box.lty=3,box.col="black")
   dev.off()


   # Plot comparison DP
	if (Sys.info()[1] == "Linux") {x11()}
	if (Sys.info()[1] == "Windows") {windows()}
	if (Sys.info()[1] == "Darwin") {quartz()}
   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""), xlab=xlabel,col="blue",xlim=range.plot.area,ylim=range.plot.density)  # Histogram raw
   lines(linear.x.data.series, kde.data.series.linear, col="orange",lty=1,lwd=1)
   lines(data.series.synthetic,10^value.fit.hde.log.doublepareto,col="red",lty=1)
   lines(data.series.synthetic,10^value.fit.kde.log.doublepareto,col="red",lty=2)
   lines(data.series.synthetic,pdf.mle.data.series.dp,col="red",lty=3)
   mtext(paste("Double Pareto Comparison -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
   legend("bottomleft",legend=c("Histogram raw data","Kernel density raw data","HDE","KDE","MLE"),pch=c(1,NA_integer_,NA_integer_,NA_integer_,NA_integer_),lty=c(FALSE,1,1,2,3),lwd=1,col=c("blue","orange","red","red","red"),cex=0.8,box.lty=3,box.col="black")

   pdf(file = paste("Model_Comparison_DP_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
   par(cex.axis=0.9)
   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""), xlab=xlabel,col="blue",xlim=range.plot.area,ylim=range.plot.density)  # Histogram raw
   lines(linear.x.data.series, kde.data.series.linear, col="orange",lty=1,lwd=1)
   lines(data.series.synthetic,10^value.fit.hde.log.doublepareto,col="red",lty=1)
   lines(data.series.synthetic,10^value.fit.kde.log.doublepareto,col="red",lty=2)
   lines(data.series.synthetic,pdf.mle.data.series.dp,col="red",lty=3)
   mtext(paste("Double Pareto Comparison -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
   legend("bottomleft",legend=c("Histogram raw data","Kernel density raw data","HDE","KDE","MLE"),pch=c(1,NA_integer_,NA_integer_,NA_integer_,NA_integer_),lty=c(FALSE,1,1,2,3),lwd=1,col=c("blue","orange","red","red","red"),cex=0.8,box.lty=3,box.col="black")
   dev.off()


   # Plot comparison IG
	if (Sys.info()[1] == "Linux") {x11()}
	if (Sys.info()[1] == "Windows") {windows()}
	if (Sys.info()[1] == "Darwin") {quartz()}
   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""), xlab=xlabel,col="blue",xlim=range.plot.area,ylim=range.plot.density)
   lines(linear.x.data.series, kde.data.series.linear, col="orange",lty=1,lwd=1)
   lines(data.series.synthetic,10^value.fit.hde.log.inversegamma,col="violet",lty=1)
   #lines(data.series.synthetic,value.fit.hde.inversegamma,col="violet",lty=1)
   lines(data.series.synthetic,10^value.fit.kde.log.inversegamma,col="violet",lty=2)
   #lines(linear.x.data.series,value.fit.kde.inversegamma,col="violet",lty=2)
   lines(data.series.synthetic,pdf.mle.data.series.ig,col="violet",lty=3)
   mtext(paste("Inverse Gamma Comparison -> ",names(data.series)[series],sep=""), side=3, col="red", cex=0.8, line=0.5)
   legend("bottomleft",legend=c("Histogram raw data","Kernel density raw data","HDE","KDE","MLE"),pch=c(1,NA_integer_,NA_integer_,NA_integer_,NA_integer_),lty=c(FALSE,1,1,2,3),lwd=1,col=c("blue","orange","violet","violet","violet"),cex=0.8,box.lty=3,box.col="black")

   pdf(file = paste("Model_Comparison_IG_",names(data.series)[series],".pdf",sep=""), width = 6, height = 6, onefile = TRUE, family = "Helvetica", fonts = NULL, paper = "special", pagecentre=TRUE)
   par(cex.axis=0.9)
   plot(linear.mids.data.series, prob.linear.data.series, log="xy", main="Probability densities", ylab=paste("Probability density",sep=""), xlab=xlabel,col="blue",xlim=range.plot.area,ylim=range.plot.density)
   lines(linear.x.data.series, kde.data.series.linear, col="orange",lty=1,lwd=1)
   lines(data.series.synthetic,10^value.fit.hde.log.inversegamma,col="violet",lty=1)
   #lines(data.series.synthetic,value.fit.hde.inversegamma,col="violet",lty=1)
   lines(data.series.synthetic,10^value.fit.kde.log.inversegamma,col="violet",lty=2)
   #lines(linear.x.data.series,value.fit.kde.inversegamma,col="violet",lty=2)
   lines(data.series.synthetic,pdf.mle.data.series.ig,col="violet",lty=3)
   mtext(paste("Inverse Gamma Comparison -> ",names(data.series)[series],sep=""), side=3, col="violet", cex=0.8, line=0.5)
   legend("bottomleft",legend=c("Histogram raw data","Kernel density raw data","HDE","KDE","MLE"),pch=c(1,NA_integer_,NA_integer_,NA_integer_,NA_integer_),lty=c(FALSE,1,1,2,3),lwd=1,col=c("blue","orange","violet","violet","violet"),cex=0.8,box.lty=3,box.col="black")
   dev.off()


   #### Exporting results
   write.table(hde.results,file=paste("HDE_Results_",names(data.series)[series],".txt",sep=''), quote=FALSE, sep = "\t", row.names=TRUE, col.names=NA)
   write.table(kde.results,file=paste("KDE_Results_",names(data.series)[series],".txt",sep=''), quote=FALSE, sep = "\t", row.names=TRUE, col.names=NA)
   write.table(mle.results,file=paste("MLE_Results_",names(data.series)[series],".txt",sep=''), quote=FALSE, sep = "\t", row.names=TRUE, col.names=NA)

  #### Exporting KS results 
    ks_greater_result[,3:5]<-round(ks_greater_result[,3:5],3)
    write.table(ks_greater_result,file=paste("Bootstrapped_KS_Test_Results_",names(data.series)[series],".txt",sep=''), quote=FALSE, sep = "\t", row.names=FALSE, col.names=TRUE)

    #### Percentile CDF size values estimated using parameters obtained with MLE 
    cdf.percentiles.sizes.results<-rbind(cdf.percentiles.sizes.qdoublepareto.simplified,cdf.percentiles.sizes.qdoublepareto,cdf.percentiles.sizes.qinversegamma)
    rownames(cdf.percentiles.sizes.results)<-NULL
    cdf.percentiles.sizes.results<-data.frame(estimation_function=c("MLE_CDF_DPS_sizes","MLE_CDF_DP_sizes","MLE_CDF_IG_sizes"),cdf.percentiles.sizes.results)
    colnames(cdf.percentiles.sizes.results)<-c("estimation_function",paste(as.numeric(cdf_percentiles)*100,"th",sep=""))
    write.table(cdf.percentiles.sizes.results,file=paste("CDF_MLE_PercentileSizes_",names(data.series)[series],".txt",sep=''), quote=FALSE, sep = "\t", row.names=FALSE, col.names=TRUE)
    

  ### Export table probability values
  if (use_shape==FALSE)
    {
    write.table(data.series.results,file=paste("ProbabilityResults_",names(data.series)[series],".txt",sep=''), quote=FALSE, sep = "\t", row.names=FALSE, col.names=TRUE)
    } else 
    {
      write.table(data.series.results,file=paste("ProbabilityResults_",names(data.series)[series],".txt",sep=''), quote=FALSE, sep = "\t", row.names=FALSE, col.names=TRUE)
      data.shape.results<-data.shape
      data.shape.results@data<-data.series.results
      writeOGR(data.shape.results,dsn=getwd(),layer=paste("ProbabilityResults_",names(data.series)[series],sep=''),driver="ESRI Shapefile",overwrite_layer=TRUE)
    }
   }


### ---------------- CDF sensitivity analysis ---------------- ###

if(executing_CDF_sensitivity_analysis==TRUE)
  {
  
  sample.size.cdf.analysis<-dim(data.series[series])[1]*10
  data.series.cdf.analysis<-rdoublepareto.simplified(n=sample.size.cdf.analysis, alpha.par=mle.results[1,1], beta.par=mle.results[2,1], t.par=mle.results[3,1])	
  data.series.synthetic.cdf<-c(10^log_cycles_range[1],sort((seq(1,10,0.1)%*%t(10^(log_cycles_range[1]:log_cycles_range[2])))[-1,]))
  data.filtering.synthetic.cdf<-c(10^log_cycles_range[1],sort((seq(1,10,1)%*%t(10^(log_cycles_range[1]:log_cycles_range[2])))[-1,]))
  
  ### Small size - High-Pass filtering
  #size.filtering.vector<-data.series.synthetic.cdf[which(data.series.synthetic.cdf<mle.results[4,1])]
  size.filtering.vector<-data.filtering.synthetic.cdf[which(data.filtering.synthetic.cdf<mle.results[4,1])]
  results.synthetic.cdf.matrix.small<-matrix(as.numeric(NA),nrow=length(data.series.synthetic.cdf),ncol=length(size.filtering.vector))
  colnames(results.synthetic.cdf.matrix.small)<-size.filtering.vector
  
  results.synthetic.cdf.coefficient.matrix.small<-matrix(as.numeric(NA),nrow=length(size.filtering.vector),ncol=10)
  colnames(results.synthetic.cdf.coefficient.matrix.small)<-c("alpha","beta","t","alpha_pvalue","beta_pvalue","t_pvalue","nboot_samples","KS_D","KS_pvalue","KS_pvalue_boot")

  data.series.cdf.analysis.filtered.old<-NULL
  for(count in 1:length(size.filtering.vector))
    {
    # count<-1
    size.filter.value<-size.filtering.vector[count]
    data.series.cdf.analysis.filtered<-data.series.cdf.analysis[which(data.series.cdf.analysis>size.filter.value)]
    print(paste("CDF Size filter value: ",size.filter.value,"; Done: ",round(count/length(size.filtering.vector)*100,1),"%",sep=""))
  
    if(identical(data.series.cdf.analysis.filtered,data.series.cdf.analysis.filtered.old) | length(data.series.cdf.analysis.filtered)==0) {print("Skipping value");next()}
    data.series.cdf.analysis.filtered.old<-data.series.cdf.analysis.filtered
  
    alpha.range<-c(min(c(hde.results[1,1],kde.results[1,1]))/1.5,max(c(hde.results[1,1],kde.results[1,1]))*1.5)
    beta.range<-c(min(c(hde.results[2,1],kde.results[2,1]))/1.5,max(c(hde.results[2,1],kde.results[2,1]))*1.5)
    t.range<-c(min(c(hde.results[3,1],kde.results[2,1]))/1.5,max(c(hde.results[3,1],kde.results[3,1]))*1.5)
    # Paramenter value specified by the user, specified_by_user == "YES" in configuration.txt file
    if (configuration[19,4]=="YES") (alpha.range<-as.numeric(configuration[19,2:3]))
    if (configuration[20,4]=="YES") (beta.range<-as.numeric(configuration[20,2:3]))
    if (configuration[21,4]=="YES") (t.range<-as.numeric(configuration[21,2:3]))
    
    if (inherits(try(mle.cdf.analysis<-mle2(mll.doublepareto.simplified,method="L-BFGS-B",trace=FALSE,lower=c(alpha.par=alpha.range[1], beta.par=beta.range[1],t.par=t.range[1]),upper=c(alpha.par=alpha.range[2], beta.par=beta.range[2],t.par=t.range[2]),start=list(alpha.par=sum(alpha.range)/2, beta.par=sum(beta.range)/2,t.par=sum(t.range)/2),data=list(x=na.omit(data.series.cdf.analysis.filtered)),control=list(maxit=1000000,lmm=5,factr=10^-8,pgtol=0,trace=6,parscale=c(alpha.par=1, beta.par=1,t.par=1))),silent=TRUE),what="try-error"))
      {
      print(paste("Warning: MLE not coverged",sep=""))
      next()
      }	
    mle.cdf.analysis.coefficients<-coef(summary(mle.cdf.analysis))[,1]
    results.synthetic.cdf.coefficient.matrix.small[count,1:3]<-mle.cdf.analysis.coefficients
    results.synthetic.cdf.coefficient.matrix.small[count,4:6]<-coef(summary(mle.cdf.analysis))[,4]
    results.synthetic.cdf.matrix.small[,count]<-pdoublepareto.simplified(x=data.series.synthetic.cdf, alpha.par=mle.cdf.analysis.coefficients["alpha.par"], beta.par=mle.cdf.analysis.coefficients["beta.par"], t.par=mle.cdf.analysis.coefficients["t.par"], FUN=FALSE)
    
    require(Matching)
    ks_boot_samples_cdf<-100
    data.series.cdf.analysis.filtered.random<-rdoublepareto.simplified(n=length(data.series.cdf.analysis.filtered), alpha.par=mle.cdf.analysis.coefficients["alpha.par"], beta.par=mle.cdf.analysis.coefficients["beta.par"], t.par=mle.cdf.analysis.coefficients["t.par"]) 
    ks.cdf.analysis<-NULL
    ks.cdf.analysis<-ks.boot(data.series.cdf.analysis.filtered,data.series.cdf.analysis.filtered.random,nboots=ks_boot_samples_cdf,alternative="two.sided")
    results.synthetic.cdf.coefficient.matrix.small[count,7:10]<-c(ks.cdf.analysis$nboots,ks.cdf.analysis$ks$statistic,ks.cdf.analysis$ks$p.value,ks.cdf.analysis$ks.boot.pvalue)
    }
  
  results.synthetic.cdf.matrix.small<-data.frame(cbind(x=data.series.synthetic.cdf,results.synthetic.cdf.matrix.small))
  names(results.synthetic.cdf.matrix.small)<-c("x",size.filtering.vector)
  results.synthetic.cdf.matrix.small<-results.synthetic.cdf.matrix.small[,is.finite(colSums(results.synthetic.cdf.matrix.small))]
  write.table(results.synthetic.cdf.matrix.small,"CDF_MLE_Sensitivity_SmallSizes_HighPass_Filtering_series.txt",sep="\t",row.names=FALSE,col.names=TRUE)
  sel.row.finite<-is.finite(rowSums(results.synthetic.cdf.coefficient.matrix.small))
  results.synthetic.cdf.coefficient.matrix.small<-results.synthetic.cdf.coefficient.matrix.small[sel.row.finite,]
  results.synthetic.cdf.coefficient.matrix.small<-data.frame(filter=paste("Size>",size.filtering.vector[sel.row.finite],sep=""),results.synthetic.cdf.coefficient.matrix.small)
  write.table(data.frame(results.synthetic.cdf.coefficient.matrix.small),"CDF_MLE_Sensitivity_SmallSizes_HighPass_Filtering_coefficients.txt",sep="\t",row.names=FALSE,col.names=TRUE)
 
  
  str(results.synthetic.cdf.matrix.small)
  library(ggplot2)
  library(reshape2)
  d.small<-melt(results.synthetic.cdf.matrix.small, id.vars="x")
  str(d.small)
  # Everything on the same plot
  pdf("CDF_MLE_Sensitivity_SmallSizes_HighPass_Filtering.pdf")
  min.value.x<-range(results.synthetic.cdf.matrix.small$x)[1]
  max.value.x<-range(results.synthetic.cdf.matrix.small$x)[2]
  min.value.at.x<-signif(log10(min.value.x),digits=1)
  max.value.at.x<-signif(ceiling(log10(max.value.x)),digits=1)
  at.x<-10^(min.value.at.x:max.value.at.x)
  at.label.x<-rep(10,length(at.x))^log10(at.x)
  at.x.tick<-sort(as.numeric((1:9)%*%t(at.label.x)))
  at.y.cdf<-seq(0,1,0.1)
  at.y.tick.cdf<-seq(0,1,0.01)
  ggplot(d.small, aes(x,value,col=variable)) + scale_x_log10(limits=c(min.value.x,max.value.x),breaks=at.label.x,labels=sciNotation(at.x,1)) + scale_y_continuous(breaks=at.y.cdf,labels=at.y.cdf) + geom_line() + labs(title="CDF sensitivity: high-pass filtering",x =xlabel, y = "Cumulative Distribution Function") + scale_colour_discrete(name = "Threshold size") + guides(col=guide_legend(ncol=2))
  dev.off()
  
  ### Large size - Low Pass filtering
  #size.filtering.vector<-data.series.synthetic.cdf[which(data.series.synthetic.cdf>mle.results[4,1])]
  size.filtering.vector<-data.filtering.synthetic.cdf[which(data.filtering.synthetic.cdf>mle.results[4,1])]
  results.synthetic.cdf.matrix.large<-matrix(as.numeric(NA),nrow=length(data.series.synthetic.cdf),ncol=length(size.filtering.vector))
  colnames(results.synthetic.cdf.matrix.large)<-size.filtering.vector
  
  results.synthetic.cdf.coefficient.matrix.large<-matrix(as.numeric(NA),nrow=length(size.filtering.vector),ncol=10)
  colnames(results.synthetic.cdf.coefficient.matrix.large)<-c("alpha","beta","t","alpha_pvalue","beta_pvalue","t_pvalue","nboot_samples","KS_D","KS_pvalue","KS_pvalue_boot")
  
  
  data.series.cdf.analysis.filtered.old<-NULL
  for(count in 1:length(size.filtering.vector))
  {
    # count<-1
    size.filter.value<-size.filtering.vector[count]
    data.series.cdf.analysis.filtered<-data.series.cdf.analysis[which(data.series.cdf.analysis>size.filter.value)]
    print(paste("CDF Size filter value: ",size.filter.value,"; Done: ",round(count/length(size.filtering.vector)*100,1),"%",sep=""))
    
    if(identical(data.series.cdf.analysis.filtered,data.series.cdf.analysis.filtered.old) | length(data.series.cdf.analysis.filtered)==0) {print("Skipping value");next()}
    data.series.cdf.analysis.filtered.old<-data.series.cdf.analysis.filtered
    
    alpha.range<-c(min(c(hde.results[1,1],kde.results[1,1]))/1.5,max(c(hde.results[1,1],kde.results[1,1]))*1.5)
    beta.range<-c(min(c(hde.results[2,1],kde.results[2,1]))/1.5,max(c(hde.results[2,1],kde.results[2,1]))*1.5)
    t.range<-c(min(c(hde.results[3,1],kde.results[2,1]))/1.5,max(c(hde.results[3,1],kde.results[3,1]))*1.5)
    # Paramenter value specified by the user, specified_by_user == "YES" in configuration.txt file
    if (configuration[19,4]=="YES") (alpha.range<-as.numeric(configuration[19,2:3]))
    if (configuration[20,4]=="YES") (beta.range<-as.numeric(configuration[20,2:3]))
    if (configuration[21,4]=="YES") (t.range<-as.numeric(configuration[21,2:3]))
    
    t.range[2]<-t.range[2]*100
    
    if (inherits(try(mle.cdf.analysis<-mle2(mll.doublepareto.simplified,method="L-BFGS-B",trace=FALSE,lower=c(alpha.par=alpha.range[1], beta.par=beta.range[1],t.par=t.range[1]),upper=c(alpha.par=alpha.range[2], beta.par=beta.range[2],t.par=t.range[2]),start=list(alpha.par=sum(alpha.range)/2, beta.par=sum(beta.range)/2,t.par=sum(t.range)/2),data=list(x=na.omit(data.series.cdf.analysis.filtered)),control=list(maxit=1000000,lmm=5,factr=10^-8,pgtol=0,trace=6,parscale=c(alpha.par=1, beta.par=1,t.par=1))),silent=TRUE),what="try-error"))
    {
      print(paste("Warning: MLE not coverged",sep=""))
      next()
    }	
    mle.cdf.analysis.coefficients<-coef(summary(mle.cdf.analysis))[,1]
    results.synthetic.cdf.coefficient.matrix.large[count,1:3]<-mle.cdf.analysis.coefficients
    results.synthetic.cdf.coefficient.matrix.large[count,4:6]<-coef(summary(mle.cdf.analysis))[,4]
    results.synthetic.cdf.matrix.large[,count]<-pdoublepareto.simplified(x=data.series.synthetic.cdf, alpha.par=mle.cdf.analysis.coefficients["alpha.par"], beta.par=mle.cdf.analysis.coefficients["beta.par"], t.par=mle.cdf.analysis.coefficients["t.par"], FUN=FALSE)
    
    require(Matching)
    ks_boot_samples_cdf<-100
    data.series.cdf.analysis.filtered.random<-rdoublepareto.simplified(n=length(data.series.cdf.analysis.filtered), alpha.par=mle.cdf.analysis.coefficients["alpha.par"], beta.par=mle.cdf.analysis.coefficients["beta.par"], t.par=mle.cdf.analysis.coefficients["t.par"]) 
    ks.cdf.analysis<-NULL
    ks.cdf.analysis<-ks.boot(data.series.cdf.analysis.filtered,data.series.cdf.analysis.filtered.random,nboots=ks_boot_samples_cdf,alternative="two.sided")
    results.synthetic.cdf.coefficient.matrix.large[count,7:10]<-c(ks.cdf.analysis$nboots,ks.cdf.analysis$ks$statistic,ks.cdf.analysis$ks$p.value,ks.cdf.analysis$ks.boot.pvalue)
    }
  results.synthetic.cdf.matrix.large<-data.frame(cbind(x=data.series.synthetic.cdf,results.synthetic.cdf.matrix.large))
  names(results.synthetic.cdf.matrix.large)<-c("x",size.filtering.vector)
  results.synthetic.cdf.matrix.large<-results.synthetic.cdf.matrix.large[,is.finite(colSums(results.synthetic.cdf.matrix.large))]
  write.table(results.synthetic.cdf.matrix.large,"CDF_MLE_Sensitivity_LargeSizes_LowPass_Filtering_series.txt",sep="\t",row.names=FALSE,col.names=TRUE)
  sel.row.finite<-is.finite(rowSums(results.synthetic.cdf.coefficient.matrix.large))
  results.synthetic.cdf.coefficient.matrix.large<-results.synthetic.cdf.coefficient.matrix.large[sel.row.finite,]
  results.synthetic.cdf.coefficient.matrix.large<-data.frame(filter=paste("Size<",size.filtering.vector[sel.row.finite],sep=""),results.synthetic.cdf.coefficient.matrix.large)
  write.table(results.synthetic.cdf.coefficient.matrix.large,"CDF_MLE_Sensitivity_LargeSizes_LowPass_Filtering_coefficients.txt",sep="\t",row.names=FALSE,col.names=TRUE)
  
  
  library(ggplot2)
  library(reshape2)
  d.large<-melt(results.synthetic.cdf.matrix.large, id.vars="x")
  str(d.large)
  # Everything on the same plot
  scientific_10 <- function(x) {parse(text=gsub("e", " %*% 10^", scientific_format()(x)))}
  
  pdf("CDF_MLE_Sensitivity_LargeSizes_LowPass_Filtering.pdf")
  min.value.x<-range(results.synthetic.cdf.matrix.large$x)[1]
  max.value.x<-range(results.synthetic.cdf.matrix.large$x)[2]
  min.value.at.x<-signif(log10(min.value.x),digits=1)
  max.value.at.x<-signif(ceiling(log10(max.value.x)),digits=1)
  at.x<-10^(min.value.at.x:max.value.at.x)
  at.label.x<-rep(10,length(at.x))^log10(at.x)
  at.x.tick<-sort(as.numeric((1:9)%*%t(at.label.x)))
  at.y.cdf<-seq(0,1,0.1)
  at.y.tick.cdf<-seq(0,1,0.01)
  ggplot(d.large, aes(x,value,col=variable)) + scale_x_log10(limits=c(min.value.x,max.value.x),breaks=at.label.x,labels=sciNotation(at.x,1)) + scale_y_continuous(breaks=at.y.cdf,labels=at.y.cdf) + geom_line() + labs(title="CDF sensitivity: low-pass filtering",x =xlabel, y = "Cumulative Distribution Function") + scale_colour_discrete(name = "Threshold size") + guides(col=guide_legend(ncol=2))
  dev.off()
  }
