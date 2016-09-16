

###################################################
#                                                 #
#                                                 #
#                   NO.33333333                   #
#                                                 #
#                                                 #
###################################################

#########################################################################################################
#
#   源数据是G100.csv
#   通过排序后的数据是G100_after.csv
#   两者应该是不一样的
#   此R程序即完成排序功能
#
#########################################################################################################

ruleBack<-read.csv("D:\\RST\\Gene100.csv")
ruleBack<-data.frame(ruleBack)
row.names(ruleBack)<-ruleBack[,1]
ruleBack<-t(ruleBack)
ruleBack<-ruleBack[-1,]

 

gen<-ruleBack

gen<-data.frame(t(gen))
rownum<-nrow(gen)#行数
gen<-data.frame(c('s'),gen,stringsAsFactors = F)#加一行
gen <- data.frame(gen[,1:2])




for(s in 1:rownum){
  gen[s,1] <- paste(as.character(rownames(gen)[s]),sep="")#赴名字
}
gen[,1]<-as.character(gen[,1])
order<-c(nrow(gen))
gen<-data.frame(order,gen)
gen <- data.frame(gen[,1:2])


#读入关联规则的结果
hc <- read.csv("D:\\RST\\rules\\sum_data.csv",header = T,stringsAsFactors = F)
hc <- hc[,-2] 
hc <- hc[,-3] 
hc <- hc[,-3]   

order<-c(nrow(hc))
from=hc[,1]

to=hc[,2]
row_names<-c(row.names(hc))
hc40ft<-data.frame(row_names,from,to,order,stringsAsFactors = F)




#################
count=0;
gencount=0;
tt=0;
mycount<-c();
mygencount<-c();
###################
hc40ft<-data.frame(hc40ft,tag=c(0))
mylist<-c()
for(ix in 1:nrow(hc40ft)){
  for(iy in 1:nrow(hc40ft)){
    
    if(ix==79&&iy==98){
      print(hc40ft[79,2])
      print(hc40ft[98,3])
      if(hc40ft[ix,2]==hc40ft[iy,3]){
        print("kao le")
        print(hc40ft[ix,5])
      }
      
    }
    if(hc40ft[ix,2]==hc40ft[iy,3]&&hc40ft[ix,5]==0){
      hc40ft[ix,5]=1
      mylist<-c(mylist,hc40ft[ix,1])
    }
  }
}

mylist2<-c()
for(xi in 1:nrow(hc40ft)){
  t=0
  for(yi in 1:length(mylist)){
    if(hc40ft[xi,1]==mylist[yi]){
      t=1
      break
    }
  }
  if(t==0){
    mylist2<-c(mylist2,hc40ft[xi,1])
  }
}
if(length(mylist2)==0){
  print("---wu tou jie dian----you huan-1-------")
  stop();
}
hc40ft_1<-hc40ft[mylist2,]  #头结点(1 jie dian)
#####################################
##hc40ft 头节点编号
for(xi in 1:nrow(hc40ft_1)){
  for(yi in 1:nrow(hc40ft)){
    if(hc40ft_1[xi,1]==hc40ft[yi,1]&&hc40ft[yi,4]>count){
      count=count+1;
      hc40ft_1[xi,4]=count;
      hc40ft[yi,4]=count;
      break;
    }
  }
}
###############################
## gen  头节点编号
for(xi in 1:nrow(hc40ft_1)){
  for(yi in 1:nrow(gen)){
    if(hc40ft_1[xi,2]==gen[yi,2]&&gen[yi,1]>gencount){
      gencount=gencount+1;
      gen[yi,1]=gencount;
      break;
    }
  }
}
##gencount
#######
mycount<-c(mycount,count);
mygencount<-c(mygencount,gencount);
######
for(zz in 1:1000){
  #####################################
  print(zz)
  #再对头结点的后继编号 (di er ci)
  for(xi in 1:nrow(hc40ft_1)){
    for(yi in 1:nrow(gen)){
      if(hc40ft_1[xi,3]==gen[yi,2]&&gen[yi,1]>gencount){
        gencount=gencount+1;
        gen[yi,1]=gencount;
        break;
      }
    }
  }
  #######################
  mygencount<-c(mygencount,gencount);
  tt=count;
  ##选出头结点的后继的节点
  mylist3<-c()
  mylist4<-c()
  for(xi in 1:nrow(hc40ft_1)){
    for(yi in 1:nrow(hc40ft)){
      if(hc40ft_1[xi,3]==hc40ft[yi,2]&&hc40ft[yi,4]>count){
        mylist3<-c(mylist3,hc40ft_1[xi,4])
        mylist4<-c(mylist4,hc40ft[yi,1])
      }
    }
  }
  if(length(mylist3)==0&&length(mylist4)==0){ break;}
  mydata_f_t<-data.frame(mylist3,mylist4)
  #hc40ft_2<-hc40ft[mylist4,]  #(2 jie dian)
  ######################################
  ###########hc40ft 第二次编号
  for(xi in 1:nrow(mydata_f_t)){
    for(yi in 1:nrow(hc40ft)){
      if(hc40ft[yi,1]==mydata_f_t[xi,2]&&hc40ft[yi,4]>count){
        count=count+1;
        hc40ft[yi,4]=count;
        break;
      }
    }
  }
  mycount<-c(mycount,count);
  if(tt==count||gencount==(nrow(gen)-1)){break;}
  ##############################
  ##消除重复
  hc40ft_2<-hc40ft[mylist4,]
  
  p=0
  for(xi in 1:nrow(hc40ft_2)){
    q=0
    ll<-c(strsplit(row.names(hc40ft_2)[xi-p],""))          
    ####消除重复
    #if(ll!=null){
    for(yi in 1:length(ll[[1]])){     
      if(ll[[1]][yi]=="."){q=1}
    }
    if(q==1){
      hc40ft_2<-hc40ft_2[-(xi-p),]
      p=p+1
    }
    if((xi-p)>=nrow(hc40ft_2)){break;}
  }
  hc40ft_1<-data.frame();
  mydata_f_t<-data.frame();
  hc40ft_1<-hc40ft_2;
  
}
######################################################################################################
cc=0
for(ix in 1:nrow(gen)){
  if(gen[ix,1]==nrow(gen))
    print(paste(gen[ix,1]," ",nrow(gen)))
  cc=cc+1;
}



names(gen)<-c('order','gene')
gen <- cbind(gen$gene, gen$order)
#现在的gene表结构为两列，左为节点名右为节点序，随即写入文件



write.csv(gen,"D:\\RST\\rules\\genNode.csv",row.names = F)
if(length(mylist3)==0&&length(mylist4)==0){
  print("-------关系找完，有不存在的基因-------")
}

######################################################################################################

sx <- read.csv("D:\\RST\\rules\\genNode.csv",header = T)

zero <- ruleBack

zero<-data.frame(t(zero))

one <- data.frame(sx,zero)
two <- one[order(one[,2],decreasing=T),] #顺序【降序 100...1 】
two <- two[,-2]

write.csv(two,"D:\\RST\\rules\\express100_after.csv",row.names = F)


#####################################【2016-6-24 结果改为调整源数据】####################################
#########################################################################################################






