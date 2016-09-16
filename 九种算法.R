






NineAlgorithmTest<-function(methodd,lsd,sourcee){
  library(igraph)
  library(bnlearn)
  library(Rgraphviz) 
  methodd<-methodd
  #lsd<-2
  sourcee<-sourcee
  
  ####################################################################################
  Sort1_Or_NotSort2 <- FALSE    #是否运行了关联规则算法排序
  LSdim <-lsd       #离散值 { 使 LSdim<-0 为 不离散 }
  NotSort_File <-sourcee                     #不排序的·输入
  Sort_File <- 'D:\\RST\\rules\\express100_after.csv'       #排序后的·输入
  ####################################################################################
  
  
  #预处理 100 + 4 + 5
  ####################################################################################################################
  
  Vitamin9<-read.csv("D:\\RST\\weishengsu_add_5_and_4.csv",header = F,stringsAsFactors = FALSE)[-1,]
  
  if(Sort1_Or_NotSort2==TRUE){
    Data100<-read.csv( Sort_File ,header = F,stringsAsFactors = FALSE)[-1,]
    #排序的话·每次抽样的结果为express100_after.csv
  }else{
    Data100<-read.csv( NotSort_File ,header = F,stringsAsFactors = FALSE)[-1,]
    #不排序的话·每次抽样的结果为Gene100.csv
  }
  SampleA<-rbind(Vitamin9,Data100) #100+9
  write.csv(SampleA,"D:\\RST\\SampleNine.csv",row.names = F)
  ####################################################################################################################
  
  #正式开始
  zp<-read.csv("D:\\RST\\SampleNine.csv",header = F,stringsAsFactors = FALSE)[-1,]
  row.names(zp)<-zp[,1]
  zp<-zp[-1]
  mydata<-data.frame(zp)  
  tq <- data.frame(matrix(as.numeric(unlist(mydata)),ncol = length(mydata[1,])))
  rownames(tq) <- rownames(mydata)
  colnames(tq) <- colnames(mydata)
  mydata <- tq
  mydata<-data.frame(t(mydata))
  
  tpN <- discretize(mydata,method =methodd,breaks= lsd ) #离散化
  rownames(tpN) <- rownames(mydata)
  colnames(tpN) <- colnames(mydata)
  mydata <- data.frame((tpN))
  
  s<-gs(mydata)
  #sink("D:\\RST\\任意取名.txt",append=TRUE,split=TRUE)
  
  
  #与之前的比较方式不同，此处直接对各算法的结果进行比较查找
  #找出满足条件的关系的数量，用函数CompareG()完成
  
  
  TDBX5<-c('LUT','ZEA','Bcry','AC','BC')#特定表型
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
  
  ToArrayStr<-function(df){
    sstr<-c()
    longstr<-c()
    for(ssi in 1:length(df[,1])){
      sstr<-c(sstr,paste('[',df[ssi,1],',',df[ssi,2],']',sep=''))
      longstr<-paste(longstr,sstr[ssi])
    }
    return (longstr)
  }
  
  
  
  
  OPinfo<-function(mtd,mdata){
    dataa<-data.frame()
    timestartNine<-Sys.time()
    if(mtd=='gs'){
      dataa<-gs(mdata)$arcs
    }else if(mtd=='hc'){
      dataa<-hc(mdata)$arcs 
    }else if(mtd=='iamb'){
      dataa<-iamb(mdata)$arcs 
    }else if(mtd=='mmpc'){
      dataa<-mmpc(mdata)$arcs 
    }else if(mtd=='rsmax2'){
      dataa<-rsmax2(mdata)$arcs 
    }else if(mtd=='tabu'){
      dataa<-tabu(mdata)$arcs
    }else if(mtd=='fastiamb'){
      dataa<-fast.iamb(mdata)$arcs
    }else if(mtd=='interiamb'){
      dataa<-inter.iamb(mdata)$arcs 
    }else if(mtd=='mmhc'){
      dataa<-mmhc(mdata)$arcs 
    }
    print(paste(mtd,'调用完成..'))
    timeendNine<-Sys.time()
    timeNine<-c(timeendNine-timestartNine)#timeNine
    
    resultw<-data.frame(dataa)
    
    tempBX<-data.frame()
    count<-0
    for(k in 1:length(resultw[,1])){
      if(CheckInArr(resultw[k,],TDBX5)==TRUE){
        count <- count+1
        tempBX[count,1]<-resultw[k,1]
        tempBX[count,2]<-resultw[k,2]
      }
    }
    
    tempG<-data.frame()
    countYes<-0
    for(k in 1:length(tempBX[,1])){
      if(CheckInArr(tempBX[k,],TDJY4)==TRUE){
        countYes <- countYes+1
        
        tempG[countYes,1]<-tempBX[k,1]
        tempG[countYes,2]<-tempBX[k,2]
      }
    }
    print(paste(mtd,'比较完成..'))
    #如果tempG为空，那么直接打印或者看其长度都会出错
    
    #tempG<-tempG[order(tempG[,1],decreasing=T),] #排个序
    #print(tempBX)#含表型BC的关系数(M)
    
    
    #cat(paste('\n本程序段使用算法为\t',mtd,
    #          '\n抽样中【表型数】为\t Z =',length(TDBX5),
    #          '\n其中含有特定基因数\t Y =',4,
    #          '\n其中不含特定基因数\t X =',100,
    #          '\n本次全部的关系数为\t',length(resultw[,1]),
    #          '\n此阶段【离散度】为\t',LSdim,
    #          '\n含表型BC的关系数有\t M =',length(tempBX[,1]),
    #          '\n其中含有特定基因数\t N =',countYes,
    #          '\n  无关基因的比重为\t p`=',(length(tempBX[,1])-countYes)/100,
    #          '\n【特定基因】比重为\t P =',countYes/4
    #))
    #cat('\n- - - - - - - - - - - - - - - - - - - - - \n')
    
    NineInfo<-data.frame(1,2,3,4,5,6,7,8,9,10,11,12,13)
    NineInfo[,1]<-mtd
    NineInfo[,6]<-length(resultw[,1])
    NineInfo[,7]<-timeNine
    NineInfo[,8]<-length(tempBX[,1])
    NineInfo[,9]<-countYes
    NineInfo[,10]<-paste('p`=',(length(tempBX[,1])-countYes)/100)
    NineInfo[,11]<-paste('P =',countYes/4)
    NineInfo[,12]<-as.character(Sys.time())
    NineInfo[,13]<-ToArrayStr(tempG)
    
    
    return(NineInfo)
  }
  
  outputN<-data.frame(1,2,3,4,5,6,7,8,9,10,11,12,13)
  
  outputN[1,]<-OPinfo('gs',mydata)
  outputN[2,]<-OPinfo('hc',mydata)
  outputN[3,]<-OPinfo('iamb',mydata)
  outputN[4,]<-OPinfo('mmpc',mydata)
  outputN[5,]<-OPinfo('rsmax2',mydata)
  outputN[6,]<-OPinfo('tabu',mydata)
  outputN[7,]<-OPinfo('fastiamb',mydata)
  outputN[8,]<-OPinfo('interiamb',mydata)
  outputN[9,]<-OPinfo('mmhc',mydata)
  
  outputN[,2]<-methodd
  outputN[,3]<-LSdim
  outputN[,4]<-sourcee
  outputN[,5]<-as.character(Sort1_Or_NotSort2)   
  
  write.table(outputN,"D:\\RST\\分阶段存储数据ssssssssssssss.txt",row.names = F,append=TRUE,col.names = F,sep=' ')
  
  
}



NineAlgorithmTest('hartemink',2,'D:\\RST\\Gene100_4.csv')
#op[1+NineCnt]<-NineAlgorithmTest('interval',2,'D:\\RST\\Gene100_4.csv')



#op<-NineAlgorithmTest('hartemink',2,'D:\\RST\\Gene100_4.csv')















