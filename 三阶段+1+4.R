###################################################
#                                                 #
#                                                 #
#                   NO.44444444                   #
#                                                 #
#                                                 #
################################################### 如果不进行排序，那么直接运行本文件

library(igraph)
library(bnlearn)
library(Rgraphviz) 

##############################################################
Sort1_Or_NotSort2 <- FALSE    #是否运行了关联规则算法排序
            LSdim <- 7        #离散值 { 使 LSdim<-0 为 不离散 }
         weight_1 <- 0.01     #互信息阀值
         weight_2 <- 0.01     #互信息阀值
         weight_3 <- 0.01     #互信息阀值
##############################################################

         
#预处理
####################################################################################################################

Vitamin9<-read.csv("D:\\RST\\weishengsu_add_1_and_4.csv",header = F,stringsAsFactors = FALSE)[-1,]

if(Sort1_Or_NotSort2==TRUE){
  Data100<-read.csv("D:\\RST\\rules\\express100_after.csv",header = F,stringsAsFactors = FALSE)[-1,]
}else{
  Data100<-read.csv("D:\\RST\\Gene100.csv",header = F,stringsAsFactors = FALSE)[-1,]
}
SampleA<-rbind(Vitamin9,Data100) #100+9
write.csv(SampleA,"D:\\RST\\Sample.csv",row.names = F)
####################################################################################################################

#正式开始

z<-read.csv("D:\\RST\\Sample.csv",header = F,stringsAsFactors = FALSE)[-1,]
row.names(z)<-z[,1]
z<-z[-1]

#将源数据z离散 
#并按区间置换其结果
#即：n值离散后，z的每一列变为n种区间，将n种区间按大小置换为整数1到n
#Ex：[2:4],[0:2],[4:9].共3个区间（无序），应对应将[0:2]换为1，[2:4]换为2，[7:9]换为3
####################################################################################

ExchangeSplit<-function(complexx){
  countPt<-length(complexx)#离散数
  sttr<-complexx[,1]
  r<-data.frame()
  for(index in 1:countPt){
    ta<-strsplit(sttr[index],split=",")
    ta<-unlist(ta)
    r[index,1]<-as.numeric(chartr(old = '(',new=' ',chartr(old='[',new = ' ',ta[1])))
    r[index,2]<-as.numeric(chartr(old = ')',new=' ',chartr(old=']',new = ' ',ta[2])))
    r[index,3]<-0.5*(r[index,1]+r[index,2]) #用区间 中值 做比较量
    r[index,4]<-0
    #将区间端点映射为数值
  }
  for(ValueChge in 1:length(r[,1])){
    minn<-255 #假设最小值
    for(pp in 1:length(r[,1])){
      if(r[pp,4]==0&&r[pp,3]<minn){
        minn<-r[pp,3]
      }
    }
    if(minn!=255){
      for(kk in 1:length(r[,1])){
        if(r[kk,3]==minn)r[kk,4]=ValueChge
      }
    }else{}
  }
   return (r[,4])#返回第四列value列
}
 
ChangeNumeric<-function(tempz,dls){
  tempz<-data.frame(tempz)
  t <- data.frame(matrix(as.numeric(unlist(tempz)),ncol = length(tempz[1,])))
  rownames(t) <- rownames(tempz)
  colnames(t) <- colnames(tempz)
  tempz <- t
  
  #统计输出(涉及统计画图)...
  
  write.csv(tempz,'D:\\RST\\SampleAfterNumeric.csv',row.names = F) 
  tempz<-data.frame(tempz)
  #注意保留原来的列名（即基因名称）
  tempName <- discretize(tempz,method = 'interval',breaks=dls) #离散化
  rownames(tempName) <- rownames(tempz)
  colnames(tempName) <- colnames(tempz)
  tempz <- tempName
  write.csv(tempz,'D:\\RST\\SampleAfterLisan.csv',row.names = F) 
  tempz<-as.matrix(tempz)
  for(p in 1:length(tempz[1,])){#列数    
    print(paste('row:',p,' complete'))    
    splitt<-data.frame(unique(tempz[,p]))
    tp<-as.matrix(splitt)
    
    ch<-ExchangeSplit(tp)#置换数据标准    
    
    for(k in 1:length(tempz[,1])){#行数    
      for(kp in 1:length(tp[,1])){
        if(tempz[,p][k]==tp[,1][kp]){ 
          tempz[,p][k]=ch[kp]   
        }else{}  
      }    
    }  
  }  
  return(tempz)
}
##############################################################################


if(LSdim!=0){
  z<-data.frame(ChangeNumeric(z,LSdim)) #二值离散，是or否
  write.csv(z,'D:\\RST\\SampleAfterChangeNB.csv',row.names = F) 
}

##############################################################################

# 计算列的 协方差矩阵的 行列式的值
getDoubleHang_value<-function(mylist1,mylist2){
  Dx1=var(mylist1)
  Dx2=var(mylist2)
  double_h_value=Dx1*Dx2-cov(mylist1,mylist2)*cov(mylist1,mylist2)
}
## 计算初始互信息，最后的结果relation即三元组[from,to,value]
mylist_a<-c()
mylist_b<-c()
mylist_info<-c()
for(xi in 1:nrow(z)){
  mylist_a<-as.double(c(z[xi,]))
  Dxa=var(mylist_a)
  for(yi in 1:nrow(z)){
    #print(yi)
    mylist_b<-as.double(c(z[yi,]))
    Dxb=var(mylist_b)
    cov_ab=getDoubleHang_value(mylist_a,mylist_b)
    mi=0.5*log(Dxa*Dxb/cov_ab)
    mylist_info<-c(mylist_info,mi)
  }
}
info<-matrix(mylist_info,nrow(z),nrow(z),byrow = TRUE)

rname<-c(row.names(z))
row.names(info)<-row.names(z)
colnames(info)<-row.names(z)
# zhao guan xi
relation<-data.frame()
for(xi in 1:(nrow(info)-1)){
  yi=xi
  while(TRUE){
    yi=yi+1
    if(nrow(relation)==0){
      from<-c(rname[xi])
      to<-c(rname[yi])
      value<-c(info[xi,yi])
      relation<-data.frame(from,to,value,stringsAsFactors = FALSE)
      relation<-data.frame(xi=c(t(relation)),stringsAsFactors = FALSE) 
    }else
      relation<-data.frame(relation,xi=c(rname[xi],rname[yi],info[xi,yi]),stringsAsFactors = FALSE)
    if(yi==nrow(info)) break
  }
}
row.names(relation)<-c("from","to","value")
relation<-data.frame(t(relation),stringsAsFactors = FALSE)


######################################################################################################

AllResult<-function(){
  from<-relation[,1][order(as.numeric(relation$value),decreasing = T)] ##降序排序
  to<-relation[,2][order(as.numeric(relation$value),decreasing = T)]
  value<-relation[,3][order(as.numeric(relation$value),decreasing = T)]
  relation2<-data.frame(from,to,value)
  mark1 <- 1
  for(wx in 1:(nrow(relation2))){
    if( as.numeric(as.character(relation2[wx,3])) > weight_1 ){
      print(paste(as.character(relation2[wx,3]),weight_1))
      mark1 = mark1 +1 
    }
  }
  result<-data.frame()
  result<-relation2[1:mark1,]
  return(result)
} 
# 计算过程中直接丢弃“基因-基因”关系
PartlyResult<-function(){
  #step1. 去掉原来relation中的“基因-基因”关系，
  #     留下的关系中必有一方含有5个特殊表型的一个，存在rlt中
  Gx<-c('LUT','ZEA','Bcry','AC','BC')
  Gx<-data.frame((Gx))
  rlt <- data.frame(1,1,1)
  cp <- 0
  for(i in 1:length(relation[,1])){
    cp = cp+1
    for(j in 1:length(Gx[,1])){
      
      if( (as.character(Gx[j,])==as.character(relation$to[i])) || (as.character(Gx[j,])==as.character(relation$from[i])) ){
        rlt<-data.frame(rlt,i=c((relation$from[i]),(relation$to[i]),(relation$value[i])),stringsAsFactors = FALSE) 
      }else{
        print(paste(as.character(Gx[j,]),"  ",as.character(relation$from[i])))
        
        cp = cp-1
      }
      
    }
  }
  #step2.  与处理relation的方法一样，得到处理后的rlt（水平结构）
  row.names(rlt)<-c("from","to","value")
  rlt<-data.frame(t(rlt),stringsAsFactors = FALSE)
  rlt <- data.frame((rlt))
  rlt <- rlt[-(1:3),]
  #step3. 获得排序后的rlt2（竖直）
  from<-rlt[,1][order(as.numeric(rlt$value),decreasing = T)] ##pai xu (jiang xu)
  to<-rlt[,2][order(as.numeric(rlt$value),decreasing = T)]
  value<-rlt[,3][order(as.numeric(rlt$value),decreasing = T)]
  rlt2<-data.frame(from,to,value)
  
  mark1 <- 1
  for(wx in 1:(nrow(rlt2))){
    if( as.numeric(as.character(rlt2[wx,3])) > weight_1 ){
      print(paste(as.character(rlt2[wx,3]),weight_1))
      mark1 = mark1 +1 
    }
  }
  result<-data.frame()
  result<-rlt2[1:mark1,]
  return(result)
}
#######################################################################################################

# 上述的【relation】和【rlt2】是相同的结构，前者是全部，后者是过滤后的结果
# 按需要调用响应函数
# PartlyResult()得到的是全部
# AllResult()得到的是部分

tempRst<-AllResult() #函数调用
write.csv(tempRst,"D:\\RST\\c2links.csv",row.names = F) 


##########################################################################
    ###############################   ###############################
     ########################### ##   ## ##########################
          ###################### ##   ## ###################
          ###########       #### ##   ## ####       ########
          ###################### ##   ## ###################
    ###############################   ###############################
##########################################################################



#weight_2 <- 0.025   #互信息阀值
#weight_3 <- 0.025  #互信息阀值

#gene是原始基因，build_data是根据互信息计算的排名
gene<-read.csv("D:\\RST\\Sample.csv",header = FALSE,stringsAsFactors = FALSE)
build_data<-read.csv("D:\\RST\\c2links.csv",header = T,stringsAsFactors = F)
ck <- union(build_data$from,build_data$to)

lname<-gene[,1]
gene<-gene[-1,]

labe<-intersect(lname,ck)

#实现将一个n.x转化成GRMZM2G151227
Translat <- function(data_nx){
  retdata <- data_nx
  for(i in 1:length(data_nx)){
    retdata[i] <- labe[which( trans_labe == data_nx[i])]
  }
  return(retdata)
}

#名字转化: GRMZM2G162755类 转为 n.x类
trans_labe<-c()
for(i in 1:length(labe)){ trans_labe<-c(trans_labe,paste("n.",i,sep = "")) }

ct<-c()

snum<-length(build_data$from)
for(i in 1:snum){ ct<-c(ct,build_data[i,1],build_data[i,2]) }

for(i in 1:length(ct)){ 
  #print(paste(i,ct[i]))
  ct[i]<-paste("n.",which(labe == ct[i]),sep = "") 
}

############
# 第一阶段 #
############

graE<-c()
graR<-c()

for(i in seq(1,length(ct),2)){
  if( length(union(graE,graE)) !=length(union(graE,c(ct[i],ct[i+1])))  ){
    graE<-c(graE,ct[i],ct[i+1])
    g<- graph(graE, directed=T)
  }else{
    #print(paste(ct[i],"+",ct[i+1]))
    if(edge_connectivity(g, source = ct[i], target = ct[i+1], checks = TRUE) == 0){
      graE<-c(graE,ct[i],ct[i+1])
      g<- graph(graE, directed=T)
    }
    else{
      graR<-c(graR,ct[i],ct[i+1])
    }
  }
}

#plot(g,layout=layout.fruchterman.reingold, vertex.size=10,vertex.color="green")

#######################################################################################
#计算互信息，源点，目标，割集，所有数据
Info <- function(gxi,gyi,cutset,gen){
  cutset<-data.frame(cutset,stringsAsFactors = F)
  x<-c()
  y<-c()
  cutset_table<-data.frame()
  for(xi in 1:nrow(gen)){
    if(gxi==gen[xi,1]){
      x<-as.numeric(c(gen[xi,2:length(gen)]))
    }
    if(gyi==gen[xi,]){
      y<-as.numeric(gen[xi,2:length(gen)])
    }
    for(yi in 1:nrow(cutset)){
      if(gen[xi,1]==cutset[yi,1]){
        if(nrow(cutset_table)==0)
          cutset_table<-data.frame(t(gen[xi,]),stringsAsFactors = F)
        else
          cutset_table<-data.frame(cutset_table,t(gen[xi,]),stringsAsFactors = F)
      }
    }
  }
  cutset_table<-data.frame(t(cutset_table),stringsAsFactors = F)
  cutset_table<-cutset_table[,-1]
  to_list<-c()
  cut<-data.frame()
  for(xi in 1:nrow(cutset_table)){
    to_list<-as.numeric(c(cutset_table[xi,]))
    if(nrow(cut)==0)
      cut<-data.frame(to_list)
    else
      cut<-data.frame(cut,to_list)
  }
  
  #代码翻译计算公式
  #实际计算中有零项，此处区别对待之
  
  if(length(x)==0){
    cut_x<-data.frame(cut)
  }else{
    cut_x<-data.frame(cut,x)
  }
  if(length(y)==0){
    cut_y<-data.frame(cut)
  }else{
    cut_y<-data.frame(cut,y)
  }
  if(length(x)==0 && length(y)==0 ){
    cut_x_y<-data.frame(cut)
  }else if(length(x)==0){
    cut_x_y<-data.frame(cut,y)
  }else if(length(y)==0){
    cut_x_y<-data.frame(cut,x)
  }else{
    cut_x_y<-data.frame(cut,x,y)
  }
  t=det(cov(cut_x))*det(cov(cut_y))/(det(cov(cut))*det(cov(cut_x_y)))
  cmi=0.5*log(t)
}
#######################################################################################


############
# 第二阶段 #
############

g<- graph(graE, directed=F)
for(i in  seq(1,length(graR),2)){
  #找到路径，并存储为n.x到one_path中
  shortpa<-shortest_paths(g, from = graR[i], to = graR[i+1], mode = c("all"))$vpath
  one_path<-names(V(g))[as.integer(shortpa[[1]])]
  #输入n.x输出GRMZM2G162755  for(j in 1:length(one_path)){ one_path[j] <- labe[which( trans_labe == one_path[j])]}
  #计算割点以及互信息
  brek <- one_path[-c(1,length(one_path))]
  if(length(brek) == 0){
    #print(" brek=0! ")          ################      
    info = 1                     #     出错     #
  }else{                         #     位置     #
    #print(" do! ")              ################
    info <- Info(Translat(graR[i]),Translat(graR[i+1]),Translat(brek),gene)
    #print(info) 
  }
  print(info)  #2016.8.8
  if(info > weight_2){ graE <- c(graE,graR[i],graR[i+1]) } #如果互信息大于weight_2则添加
}

rm("one_path","i","shortpa","info","brek","graR")


############
# 第三阶段 #
############

g<- graph(graE, directed=F)

#起点和终点，其中每一个相邻的两个节点都代表着一条边，除的关系是(i+1)/2,其中所有的边都没有重复的边
for(i in  seq(1,length(graE),2)){
  g_d <- g - edge(paste(graE[i],"|",graE[i+1],sep = ""))  #删除当前边并存放在g_d中
  #在 g_d 中查找是不是还存在路径,存在的话就进行计算互信息，不存在则跳过
  if( edge_connectivity(g_d, source = graE[i], target = graE[i+1], checks = TRUE) > 0){
    #找到路径，并存储为n.x到one_path中
    shortpa<-shortest_paths(g_d, from = graE[i], to = graE[i+1], mode = c("all"))$vpath
    one_path<-names(V(g))[as.integer(shortpa[[1]])]
    #计算割点以及互信息
    brek <- one_path[-c(1,length(one_path))]
    if(length(brek) > 0){
      print(brek)
      info <- Info(Translat(graE[i]),Translat(graE[i+1]),Translat(brek),gene)
      #print(info)
      if(info < weight_3){ g <- g_d } #如果互信息小于weight_3，那么久确定删除
    }else{
      print(i)
    }
  }
}
#最终结果在result中
result <- as_edgelist(g, names = TRUE)
for(i in 1:length(result[,1])){
  
  result[i,] <- Translat(result[i,])
}
write.csv(result,"D:\\RST\\XResultX.csv",row.names = F)

##########################################################################################
   ##################################################################################
         ######################################################################
         #########        “基因表型”实验中   最后的比较计算部分        ######## 
   ##################################################################################
##########################################################################################



TDBX5<-c('BC')#特定表型
TDJY4<-c('GRMZM2G012966','GRMZM2G152135','GRMZM2G300348','GRMZM2G108457')#特定基因
CheckInArr<-function(nm,arr){
  cnt<-length(arr)
  for(ci in 1:cnt){
    if(nm[,1]==arr[ci]||nm[,2]==arr[ci]){
      return (TRUE)
    }
  }
  return (FALSE)
}

#筛选满足特定表型的关系
tempBX<-data.frame()
countYes<-0
for(k in 1:length(resultw[,1])){
  if(CheckInArr(resultw[k,],TDBX5)==TRUE){
    countYes <- countYes+1
    #print(resultw[k,])
    tempBX[countYes,1]<-resultw[k,1]
    tempBX[countYes,2]<-resultw[k,2]
  }
}

#在上述结果基础上 筛选满足特定基因的关系
tempG<-data.frame()
countYes<-0
for(k in 1:length(tempBX[,1])){
  if(CheckInArr(tempBX[k,],TDJY4)==TRUE){
    countYes <- countYes+1
    tempG[countYes,1]<-tempBX[k,1]
    tempG[countYes,2]<-tempBX[k,2]
  }
  
}
tempG<-tempG[order(tempG[,1],decreasing=T),] #排个序
#print(tempBX)#含表型BC的关系数(M)


cat(paste('\n流程中是否涉及排序\t',Sort1_Or_NotSort2,
          '\n源抽样数据基因数为\t',(100+4),
          '\n其中含有特定基因数\t X =',4,
          '\n其中不含特定基因数\t Y =',100,
          '\n本次全部的关系数为\t',length(resultw[,1]),
          '\n三阶段的【阈值】为\t',weight_1,' ',weight_2,' ',weight_3,' ',
          '\n此阶段【离散度】为\t',LSdim,
          '\n含表型BC的关系数有\t M =',length(tempBX[,1]),
          '\n其中含有特定基因数\t P =',countYes,
          '\n  无关基因的比重为\t',(length(tempBX[,1])-countYes)/100,
          '\n【特定基因】比重为\t',countYes/4
))
print(tempG)#M中含4个特定基因的关系如下



