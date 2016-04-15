rm(list=ls())

setwd("/Users/syedrahman/Documents/Spring2015/Hipergator/07062015")

code.generate = function(par.vals, pat.nams, path="", path.opt=F,
                            templatename, out_name){ 
  for(i in 1:length(par.vals)){
    SourceCode=readLines(templatename) #read desired code

    #replace "AA" with desired value
    SourceCode=gsub(pattern=pat.nams[i],replacement=par.vals[i],x=SourceCode) 
  
    #output changed code
    path=ifelse(path==""|path.opt==F,"",paste(path,"/",sep=""))
    cat(SourceCode,file=paste(path, out_name,
                            paste(par.vals[i],collapse=""),
                            '.R',sep=""),sep='\n') 
  }
}  

######################################################################

numvals = 1:100
nams = rep("AA",times=length(numvals))

code.generate(numvals,nams,templatename='tdist07062015p100n200.R', out_name='tdist07062015p100n200_')
