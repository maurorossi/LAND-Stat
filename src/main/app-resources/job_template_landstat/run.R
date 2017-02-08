#!/opt/anaconda/bin/Rscript --vanilla --slave --quiet

# load rciop library to access the developer cloud sandbox functions
library("rciop")

# load any other R library required
library("rgeos")

# load the parametee values with rciop.getparam() function
myparam1 <- rciop.getparam("param1")
myparam1 <- as.numeric(myparam1)

myparam2 <- rciop.getparam("param2")
myparam2 <- as.numeric(myparam2)

param_file_name <- rciop.getparam("file_name")
####
res <- rciop.copy(param_file_name, TMPDIR, uncompress=TRUE)
if (res$exit.code==0) local.url <- res$output
tmp.df <- read.table(local.url,sep="\t",dec=".",header=TRUE,stringsAsFactors=FALSE)
#str(tmp.df)

output_dir<-paste(TMPDIR,"output",sep="/")
dir.create(output_dir)

pdf(paste(output_dir,"histogram.pdf",sep="/"))
hist(tmp.df[,1])
dev.off()

histogram<-hist(tmp.df[,1])
write.table(data.frame(histogram$density),paste(output_dir,"histogram.txt",sep="/"),sep="\t")

#input <- "http://landsat.usgs.gov/documents/L7_60m20090422.tgz"
#res <- rciop.copy(input, TMPDIR, uncompress=TRUE)
#if (res$exit.code==0) local.url <- res$output
#local.prms <- apply(as.data.frame(prm), 1, function(url) { rciop.copy(url, TMPDIR)$output })
#tmp.df <- read.table(textConnection(local.prms))
####

# add a log message
rciop.log("DEBUG", paste("I'm running a job with parameter values:", myparam1, myparam2,param_file_name, sep=" "))

# # read the inputs coming from stdin
# f <- file("stdin")
# open(f)

# while(length(input <- readLines(f, n=1)) > 0) {
  
#   rciop.log("INFO", paste("processing input:", input, sep=" "))
  
#   # copy the input to the process temporary folder TMPDIR
#   res <- rciop.copy(input, TMPDIR, uncompress=TRUE)
  
#   if (res$exit.code==0) local.url <- res$output
  
#   mycsv <- read.csv(local.url)
  
#   # do something with the downloaded csv here in TMPDIR/output
  
#   # publish the any results done 
#   rciop.publish(paste(TMPDIR,"output", sep="/"), recursive=TRUE, metalink=FALSE)
 
# }
data_random<-rnorm(100,mean=myparam2,sd=myparam1)
data_random2<-rnorm(100,mean=myparam2,sd=myparam1)
a<-hist(data_random)

dt<-data.frame(x=data_random,y=data_random2)
rciop.publish(output_dir, recursive=TRUE, metalink=FALSE)



