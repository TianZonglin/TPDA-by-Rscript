
###################################################
#                                                 #
#                   统 计 画 图                   #
#                                                 #
###################################################
#install.packages('gcookbook')#R数据可视化手册书中的数据集
#install.packages('gmodels')  #需要gmodels包 
##############################################################################



sourceFile<-read.csv("D:\\RST\\Sample.csv",header = F,stringsAsFactors = FALSE)[-1,]
ComplexZhengShu<-read.csv("D:\\RST\\SampleAfterNumeric.csv",header = F,stringsAsFactors = FALSE)[-1,]
QuJian<-read.csv("D:\\RST\\SampleAfterLisan.csv",header = F,stringsAsFactors = FALSE)[-1,]
ZhengShu<-read.csv("D:\\RST\\SampleAfterChangeNB.csv",header = F,stringsAsFactors = FALSE)[-1,]

rownames(sourceFile)<-sourceFile[,1]
sourceFile<-sourceFile[,-1]


#########################################
line<-as.factor(z)
fqc<-data.frame(1,1)
colnames(fqc)<-c('name','value')
numfqc<-1


for(p in 1:length(line)){
  print(p) 
  
  for(k in 1:length(fqc$name)){
    if(fqc$name[k]==line[p]){
      fqc$value[k]=fqc$value[k]+1
    }else{
      #numfqc<-numfqc+1
      #fqc[numfqc,1]<-line[p]
      #fqc[numfqc,2]<-1
      fqc<-c(fqc,line[p],1)
    }
  }
}



library(ggplot2)
library(gcookbook) 
#ggplot(z, aes(x=group, y=weight)) + geom_bar(stat="identity")
ggplot(z, aes(x=V1,y=V2)) + geom_bar(stat="identity")
z<-as.matrix(z)
length(z[z[,1]==1])
