require('ggplot2')
require('reshape2')
setwd('/Users/erfan/Main/Library/polytomytest/')
d <- read.csv('data/comapre-s_tree.polytomy.csv',sep=" ",header=F)
d2 <- d[d$V4 != 10,]
#d2$V4<-as.factor(d2$V4)
d2$V8<-as.numeric(as.character(d2$V8))
d2$V6<-as.factor(d2$V6)
ggplot(data=d2,aes(x=V4,y=V8,fill=V6))+geom_boxplot()
d2<-d2[!is.na(d2$V8),]


ggplot(data=d2[d2$V3=="true",],aes(x=as.factor(V4),y=V8,fill=V6))+
  geom_jitter(alpha=0.333,width = 0.1)+
  facet_wrap(~V6,nrow=2)+theme_bw()+
  geom_hline(yintercept = 0.05,color="red",linetype=2)+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=13),
        legend.text = element_text(size=12),
        legend.position = "",
        strip.text=element_text(size=12))+
        xlab("number of genes")+ylab("p-values")
ggsave('polytomy-model1and2-jitter-true.pdf',width = 12.5, height= 7.01)

ggplot(data=d2[d2$V3=="estimated",],aes(x=as.factor(V4),y=V8,fill=V6))+
  geom_jitter(alpha=0.333,width = 0.1)+
  facet_wrap(~V6,nrow=2)+theme_bw()+
  geom_hline(yintercept = 0.05,color="red",linetype=2)+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=13),
        legend.text = element_text(size=12),
        legend.position = "",
        strip.text=element_text(size=12))+
  xlab("number of genes")+ylab("p-values")
ggsave('polytomy-model1and2-jitter-estimated.pdf',width = 12.5, height= 7..1)








ggplot(data=d2[d2$V3 == "true",],aes(x=V8,color=as.factor(V4)))+stat_ecdf()+
  facet_wrap(~V6,scales = "free",nrow = 2)+theme_bw()+
  geom_vline(xintercept = 0.05,color="#ff0000",linetype=2)+
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=12),
        legend.text = element_text(size=11),
        legend.position = "bottom",
        strip.text=element_text(size=11))+
  scale_color_manual(name="",values = 
    c('#b2182b','#ef8a62','#ffb9b7','#b1e5f0','#67a9cf','#2166ac'))+
  xlab("p-values")+ylab("ecdf")

ggsave('polytomy-model1and2-true.pdf',width = 12.5, height= 7.01)




ggplot(data=d2[d2$V3 == "true",],aes(x=V8,color=as.factor(V4)))+stat_ecdf()+
  facet_wrap(~V6,nrow = 2)+theme_bw()+
  geom_vline(xintercept = 0.05,color="#ff0000",linetype=2)+
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=12),
        legend.text = element_text(size=11),
        legend.position = "bottom",
        strip.text=element_text(size=11))+
  scale_color_manual(name="",values = 
                       c('#b2182b','#ef8a62','#ffb9b7','#b1e5f0','#67a9cf','#2166ac'))+
  xlab("p-values")+ylab("ecdf")

ggsave('polytomy-model1and2-true.pdf',width = 12.5, height= 7.01)



ggplot(data=d2[d2$V3 != "true",],aes(x=V8,color=as.factor(V4)))+stat_ecdf()+
  facet_wrap(~V6,scales = "free",nrow = 2)+theme_bw()+
  geom_vline(xintercept = 0.05,color="#ff0000",linetype=2)+
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=12),
        legend.text = element_text(size=11),
        legend.position = "bottom",
        strip.text=element_text(size=11))+
  scale_color_manual(name="",values = 
                       c('#b2182b','#ef8a62','#ffb9b7','#b1e5f0','#67a9cf','#2166ac'))+
  xlab("p-values")+ylab("ecdf")

ggsave('polytomy-model1and2-est.pdf',width = 12.5, height= 7.01)




a<-dcast(data=d2,V1+V3+V4+V6~.,fun.aggregate = function(x){length(x[x<0.05])},value.var="V8")
head(d2)




