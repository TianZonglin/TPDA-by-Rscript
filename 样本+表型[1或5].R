Sort1_Or_NotSort2=TRUE
Vitamin9<-read.csv("D:\\RST\\weishengsu_add_5_and_4.csv",header = F,stringsAsFactors = FALSE)[-1,]
if(Sort1_Or_NotSort2==TRUE){
 print("1111111111")
}else{
  Data100<-read.csv("D:\\RST\\Gene100.csv",header = F,stringsAsFactors = FALSE)[-1,]
}
SampleA<-rbind(Vitamin9,Data100) #100+9
write.csv(SampleA,"D:\\RST\\Sample.csv",row.names = F)


#抽样集构造完成

