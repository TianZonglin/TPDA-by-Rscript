

library(igraph)
library(bnlearn)
library(Rgraphviz) 

####################################################################################

mydata<- read.csv("D:\\RST\\8基因实验\\Gene8.csv", header=FALSE)

mydata<-data.frame(mydata)
 
mydata<-t(mydata)
colnames(mydata)<-mydata[1,]
mydata<-mydata[-1,]
mydata<-data.frame(mydata)

t <- data.frame(matrix(as.numeric(unlist(mydata)),ncol = length(mydata[1,])))
rownames(t) <- rownames(mydata)
colnames(t) <- colnames(mydata)
mydata <- t
mydata <- discretize(t(mydata),method = 'interval',breaks= 7 ) #离散化

#print(mydata)
##################################################################################
OPinfo<-function(bb,aa){
  a<-data.frame(aa)
  NineInfo<-data.frame(1,2,3)
  for(i in 1:length(a[,1])){
      NineInfo[i,1]<-bb
      NineInfo[i,2]<-as.character(a[i,1])
      NineInfo[i,3]<-as.character(a[i,2])
  }
  return(NineInfo)
}


OPinfo('gs',gs(mydata)$arcs) 
OPinfo('hc',hc(mydata)$arcs) 
OPinfo('iamb',iamb(mydata)$arcs) 
OPinfo('mmpc',mmpc(mydata)$arcs) 
OPinfo('rsmax',rsmax2(mydata)$arcs) 
OPinfo('tabu',tabu(mydata)$arcs)
OPinfo('fastiamb',fast.iamb(mydata)$arcs)
OPinfo('interiamb',inter.iamb(mydata)$arcs) 
OPinfo('mmhc',mmhc(mydata)$arcs) 


##################################################################################
