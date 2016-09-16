

###################################################
#                                                 #
#                                                 #
#                   NO.22222222                   #
#                                                 #
#                                                 #
###################################################
#NO_2_SECONDSTEP.JAVA

#上一步使用java中的vector来遍历csv中的每一个元素，并统一转成Vi*X的格式
#为什么不用R：R很适合处理表的行列变换，但是需要遍历每个元素时效率很低，要避免之

#此步
#处理之后的原数据已经可以作为apriori函数的输入（此输入有固定格式）
#此简短代码就是使用已有封装apriori函数来完成关联规则的分析，并把结果转存，作为后续使用

setwd("D:/RST/rules/")
library("arules")
tr <-  read.transactions("Mult_data_100.txt", format="basket", sep=",")
rules=apriori(tr,parameter=list(support=0.15,confidence=0.45,minlen=2,maxlen=2)) #原来是 0.3  0.6
summary(rules)    #察看求得的关联规则之摘要
outdata <- as(rules,"data.frame")  #转换数据格式
write.table(outdata,"Result.csv",sep=",")  #输出结果






