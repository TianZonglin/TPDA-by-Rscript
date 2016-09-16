

###################################################
#                                                 #
#                                                 #
#                   NO.11111111                   #
#                                                 #
#                                                 #
###################################################
                #NO_1_FIRSTSTEP.JAVA

AllG<-read.csv("D:\\RST\\weishengsu_expression_change.csv",header = F,stringsAsFactors = FALSE)#D:\\RST\\weishengsu_add_1_and_4.csv
AllG<-data.frame(t(AllG))
AllG<-AllG[,-1]
G100<-sample(AllG,100)
G100<-data.frame(t(G100))
write.csv(G100,"D:\\RST\\Gene100_6.csv",row.names = F)
###############################################

