require(ggplot2)
require(reshape2)
require(plyr)
setwd('/Users/erfan/Main/Library/polytomytest')

f<-read.csv('data/scores.csv',sep=" ",header=F)
f<-f[! (f$V1 == "model.200.500000.0.000001" &f$V2 %in% c("8","15","49")),]
f1<-f[f$V1 %in% c("model.200.500000.0.000001","model.200.500000.0.0000001"), ]
f2<-f[f$V1 %in% c("model.200.2000000.0.000001","model.200.2000000.0.0000001"), ]
f3<-f[f$V1 %in% c("model.200.10000000.0.000001","model.200.10000000.0.0000001"), ]
f1$V1<-"0.5M"
f2$V1<-"2M"
f3$V1<-"10M"
f<-rbind(f1,f2)
f<-rbind(f,f3)
af<-dcast(f,V1~.,fun.aggregate=mean,value.var="V3")



d<-read.csv('data/pvalues-polytomytest.csv',header=F,sep=" ")

d$V5<-d$V5/200000

dtarget1<-d[d$V1 %in% c("model.200.500000.0.0000001"),]
dtarget2<-d[d$V1 %in% c("model.200.500000.0.000001"),]
dtarget<-d[grepl("model.200.[1-2]",d$V1),]
dtarget2<-dtarget2[!dtarget2$V2 %in% c("8","15","49"),]
dtarget<-rbind(dtarget,dtarget1)
dtarget<-rbind(dtarget,dtarget2)
dtarget$V8<-(as.numeric(as.character(dtarget$V7))<0.05)
dtarget<-dtarget[!is.na(dtarget$V7),]
dtarget$V9<-">=0.05"
dtarget[dtarget$V8,]$V9<-"<0.05"
dtarget$V4<-as.factor(dtarget$V4)
dtarget$V3<-as.factor(dtarget$V3)
levels(dtarget$V4)<-list("50gt"="50","200gt"="200","1000gt"="1000")
levels(dtarget$V3)<-list("true gene trees"="truegt","estimated gene trees"="estimatedgt")
ggplot(data=dtarget[dtarget$V1 %in% c("model.200.500000.0.0000001","model.200.500000.0.000001"),],
       aes(fill=V9,x=(as.numeric(as.character(V5)))))+
  stat_bin(binwidth = 0.1,position = "stack" )+facet_grid(V3~V4)+
  theme_bw()+scale_fill_brewer(name="",palette = "Set1")+
  xlab('branch lengths (cu.)')+ylab('count')+
  theme(axis.text.x = element_text(size=10,hjust=1,angle = 40),
        axis.text.y = element_text(size=10),
        title = element_text(size=10),
        legend.position = c(0.1,0.8))+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,0.5,2,15))+geom_vline(xintercept = 0.1)
ggsave('figure/results-simulation-0.5M.pdf')



ggplot(data=dtarget[dtarget$V1 %in% c("model.200.2000000.0.0000001","model.200.2000000.0.000001"),],
       aes(fill=V9,x=(as.numeric(as.character(V5)))))+
  stat_bin(binwidth = 0.1,position = "stack" )+facet_grid(V3~V4)+
  theme_bw()+scale_fill_brewer(name="",palette = "Set1")+
  xlab('branch lengths (cu.)')+ylab('count')+
  theme(axis.text.x = element_text(size=10,hjust=1,angle = 40),
        axis.text.y = element_text(size=10),
        title = element_text(size=10),
        legend.position = c(0.1,0.8))+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,0.5,2,15))+geom_vline(xintercept = 0.1)
ggsave('figure/results-simulation-2M.pdf')



ggplot(data=dtarget[dtarget$V1 %in% c("model.200.10000000.0.000001","model.200.10000000.0.0000001"),],
       aes(fill=V9,x=(as.numeric(as.character(V5)))))+
  stat_bin(binwidth = 0.1,position = "stack" )+facet_grid(V3~V4)+
  theme_bw()+scale_fill_brewer(name="",palette = "Set1")+
  xlab('branch lengths (cu.)')+ylab('count')+
  theme(axis.text.x = element_text(size=10,hjust=1,angle = 40),
        axis.text.y = element_text(size=10),
        title = element_text(size=10),
        legend.position = c(0.07,0.8))+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,0.5,2,15))+geom_vline(xintercept = 0.1)
ggsave('figure/results-simulation-10M.pdf')



ggplot(data=dtarget[dtarget$V1 %in% c("model.200.10000000.0.000001","model.200.10000000.0.0000001"),],
       aes(fill=V9,x=(as.numeric(as.character(V5)))))+
  stat_bin(binwidth = 0.1,position = "stack" )+facet_grid(V3~V4)+
  theme_bw()+scale_fill_brewer(name="",palette = "Set1")+
  xlab('branch lengths (cu.)')+ylab('count')+
  theme(axis.text.x = element_text(size=10,hjust=1,angle = 40),
        axis.text.y = element_text(size=10),
        title = element_text(size=10),
        legend.position = c(0.07,0.8))+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,0.5,2,15))+geom_vline(xintercept = 0.1)


dt1<-dtarget[dtarget$V1 %in% c("model.200.500000.0.0000001","model.200.500000.0.000001"),]
a=20
dt1$V10<-cut(dt1$V5,quantile(round(dt1$V5*10^8)/10^8,0:a*1/a))

cls<-c("#7bbee6","#4191c5","#08519c")
ggplot(data=dt1,
       aes(color=V4,linetype=V3,group=interaction(V3,V4),x=V10,y=(as.numeric(as.character(V7)))))+
  stat_summary(fun.y=function(x){length(x[x<0.05])/length(x)*100},geom="line",size=0.6)+
#  stat_bin(binwidth = 0.1,position = "stack" )+facet_grid(V3~V4)+
  theme_bw()+scale_fill_brewer(name="",palette = "Set1")+
  xlab('branch lengths (cu.)')+ylab('percent of branches getting rejected for each bin')+
  theme(axis.text.x = element_text(size=10,hjust=1,angle = 40),
        axis.text.y = element_text(size=10),
        title = element_text(size=10),
        legend.position = c(0.83,0.4))+
  scale_color_grey(start=0.8, end=0.2,name="")+
  scale_linetype_manual(name="",values=c(2,1))
ggsave('figure/results-simulation-0.5M-line-grey.pdf')


dt2<-dtarget[dtarget$V1 %in% c("model.200.2000000.0.0000001","model.200.2000000.0.000001"),]
a=20
dt2$V10<-cut(dt2$V5,quantile(round(dt2$V5*10^6)/10^6,0:a*1/a))

ggplot(data=dt2,
       aes(color=V4,linetype=V3,group=interaction(V3,V4),x=V10,y=(as.numeric(as.character(V7)))))+
  stat_summary(fun.y=function(x){length(x[x<0.05])/length(x)*100},geom="line",size=0.6)+
  #  stat_bin(binwidth = 0.1,position = "stack" )+facet_grid(V3~V4)+
  theme_bw()+scale_fill_brewer(name="",palette = "Set1")+
  xlab('branch lengths (cu.)')+ylab('percent of branches getting rejected for each bin')+
  theme(axis.text.x = element_text(size=10,hjust=1,angle = 40),
        axis.text.y = element_text(size=10),
        title = element_text(size=10),
        legend.position = c(0.83,0.4))+
  geom_vline(xintercept = 31,color="red",size=0.6)+scale_color_grey(start=0.8, end=0.2,name="")+
  scale_linetype_manual(name="",values=c(2,1))
ggsave('figure/results-simulation-2M-line-grey.pdf')



dt3<-dtarget[dtarget$V1 %in% c("model.200.10000000.0.0000001","model.200.10000000.0.000001"),]
a=20
l<-unname(quantile(round(dt3$V5*10^5)/10^5,0:a*1/a))
l[1]<-min(dt3$V5)-min(dt3$V5)/20
dt3$V10<-cut(dt3$V5,l)

ggplot(data=dt3,
       aes(color=V4,linetype=V3,group=interaction(V3,V4),x=V10,y=(as.numeric(as.character(V7)))))+
  stat_summary(fun.y=function(x){length(x[x<0.05])/length(x)*100},geom="line",size=0.6)+
  #  stat_bin(binwidth = 0.1,position = "stack" )+facet_grid(V3~V4)+
  theme_bw()+scale_fill_brewer(name="",palette = "Set1")+
  xlab('branch lengths (cu.)')+ylab('percent of branches getting rejected for each bin')+
  theme(axis.text.x = element_text(size=10,hjust=1,angle = 40),
        axis.text.y = element_text(size=10),
        title = element_text(size=10),
        legend.position = c(0.83,0.4))+
  geom_vline(xintercept = 31,color="red",size=0.6)+scale_color_grey(start=0.8, end=0.2,name="")+
  scale_linetype_manual(name="",values=c(2,1))
ggsave('figure/results-simulation-10M-line-grey.pdf')

dcast(data=dtarget[dtarget$V5>=0.1 & dtarget$V5<=0.2  & dtarget$V4=="1000gt",],V3~.,
      fun.aggregate = function(x){length(x[x<0.05])/length(x)},value.var="V7")

h<-read.csv('data/allCompare.csv',header=F,sep=" ")

h1<-h[h$V1=="1KP",]
h2<-h[h$V1=="Avian_biological",]
h3<-h[h$V1=="Insects",]
h4<-h[h$V1=="Xenoturbella_Cannon",]
h5<-h[h$V1=="Xenoturbella_Greg",]
h1$V8 <- 844
h2$V8 <- 2022
h3$V8 <- 1478
h4$V8 <- 212
h5$V8 <- 393
ht<-rbind(h1,h2)
ht<-rbind(ht,h3)
ht<-rbind(ht,h4)
ht<-rbind(ht,h5)

ht$V9<-round(ht$V8*ht$V2/100)

ht$V5<-as.numeric(as.character(ht$V5))
ht<-ht[ht$V5<5000,]

ht$V10<-">=0.05"
ht[ht$V5<0.05,]$V10<-"<0.05"
ht2<-dcast(ht,V1+V2+V6+V7+V8~.,fun.aggregate=mean,value.var="V5")
ht3<-dcast(ht,V1+V2+V6+V7+V8~.,fun.aggregate=mean,value.var="V4")
names(ht2)<-c("V1","V2","V6","V7","V8","V5")
names(ht3)<-c("V1","V2","V6","V7","V8","V4")

htf<-join(ht3,ht2)
htf$V10<-">=0.05"
htf[htf$V5<0.05,]$V10<-"<0.05"
htf$V4<-htf$V4+min(htf[htf$V4>0,]$V4)/20
htf$V4<-htf$V4-min(htf$V4)
ggplot(data=htf[htf$V2=="100",],
       aes(fill=V10,x=(as.numeric(as.character(V4)))))+
  stat_bin(binwidth = 0.1,position = "stack" )+facet_wrap(~V1)+
  theme_bw()+scale_fill_brewer(name="",palette = "Set1")+
  xlab('branch lengths (cu.)')+ylab('count')+
  theme(axis.text.x = element_text(size=10,hjust=1,angle = 40),
        axis.text.y = element_text(size=10),
        title = element_text(size=10),
        legend.position = c(0.1,0.8))+
  scale_x_log10(breaks=c(0.000001, 0.00001,0.0001,0.001,0.01,0.1,0.2,0.5,2,15))+geom_vline(xintercept = 0.1)

htf$V11<-paste(htf$V1,htf$V8,sep=" (")
htf$V11<-paste(htf$V11,")",sep="")
htf$V11<-as.factor(htf$V11)
levels(htf$V11)<-list("Plants (844gt)"="1KP (844)","Avian (2022gt)"="Avian_biological (2022)","Insects (1478gt)"="Insects (1478)","Xenoturbella (Cannon et. al. 212gt)"="Xenoturbella_Cannon (212)","(Rouse et. al. 393gt)"="Xenoturbella_Greg (393)"  )
ggplot(data=htf[htf$V2=="100"& htf$V4>0,],
       aes(color=V10,x=(as.numeric(as.character(V4))),y=V5))+
  geom_point(size=1.5)+facet_wrap(~V11,nrow=1)+
  theme_bw()+scale_fill_brewer(name="",palette = "Set1")+
  xlab('branch lengths (cu.)')+ylab('p-value')+
  theme_classic()+theme(axis.text.x = element_text(size=12,hjust=1,angle = 50),
                        axis.text.y = element_text(size=12),
                        title = element_text(size=12),
                        legend.position = c(0.98,0.8))+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,0.5,2,15))+
  scale_color_brewer(name="",palette = "Set1")
ggsave('figure/biological-point.pdf',width=16, height=3.7 )

mp<-read.csv('data/map',sep="\t",header=F)
ht2<-merge(x=htf,y=mp,by.x=c("V1","V6","V7"),by.y=c("V1","V2","V3"))
levels(ht2$V1)<-list("Plants"="1KP","Avian"="Avian_biological","Insects"="Insects","Xenoturbella (Cannon et. al.)"="Xenoturbella_Cannon","Xenoturbella (Rouse et. al.)"="Xenoturbella_Greg")

ggplot(data=ht2[ht2$V2>4 & ht2$V1=="Plants",],aes(x=V8*V2/100,y=V5,color=as.factor(V4.y)))+geom_point()+#facet_wrap(~V4.y)+
  theme_bw()+facet_wrap(~V1)+
  stat_smooth(se=F)+theme(legend.position = "bottom")+xlab("Number of genes")+ylab("p-value")+
  geom_hline(yintercept = 0.05,size=0.6)+scale_color_brewer(name="",palette = "Paired")+xlab("Number of genes")+ylab("p-value")+
  guides(colour=guide_legend(nrow=3))
ggsave('figure/plants-biological.pdf',width=7.26  ,height=6.72)


ggplot(data=ht2[ht2$V2>4 & ht2$V1=="Insects" & ht2$V4.y %in% c("G","E","Ingroup","P","J","L"),],aes(x=V8*V2/100,y=V5,color=V4.y))+geom_point()+theme_bw()+
  facet_wrap(~V1)+
  stat_smooth(se=F)+theme(legend.position = "bottom")+guides(colour=guide_legend(nrow=2))+xlab("Number of genes")+ylab("p-value")+
  geom_hline(yintercept = 0.05,size=0.6)+facet_wrap(~V1)+
  scale_color_brewer(name="",palette = "Paired",labels=c("Holometabola","Acercaria+Hymenoptera","Hexapoda","Orthopteroidea","Pterygota","Psocodea+Holometabola"))+
  xlab("Number of genes")+ylab("p-value")
ggsave('figure/insects-biological.pdf',width=7.26  ,height=6.67)

# ggplot(data=ht2[ht2$V1=="Insects" & ht2$V4.y %in% c("G","E","Ingroup","P","J","L"),],aes(x=V8*V2/100,y=log10(V5),color=V4.y))+geom_point()+theme_bw()+
#   scale_y_continuous(labels = function(x){print(x);10^(x)})+stat_smooth(se=F)+
#   geom_hline(yintercept = log10(0.05),size=0.6)
  
# ggplot(data=ht2[ht2$V1%in% c("Xenoturbella_Cannon","Xenoturbella_Greg"),],aes(x=V8*V2/100,y=log10(V5),color=V4.y))+
#   geom_point()+facet_wrap(~V1,scales = "free_x")+theme_bw()+scale_y_continuous(labels = function(x){print(x);10^(x)})+
#   geom_hline(yintercept =log10(0.05) )+stat_smooth(se=F)
ggplot(data=ht2[ht2$V2>4 & ht2$V1%in% c("Xenoturbella (Cannon et. al.)","Xenoturbella (Rouse et. al.)"),],aes(x=V8*V2/100,y=V5,color=V4.y))+
  geom_point()+facet_wrap(~V1,scales = "free_x")+theme_bw()+
  #scale_y_continuous(labels = function(x){print(x);1-10^(x)},breaks = c(0,log10(0.95),log10(0.5),log10(0.9)))+
  stat_smooth(se=F)+theme(legend.position = "bottom")+xlab("Number of genes")+ylab("p-value")+
  geom_hline(yintercept = 0.05,size=0.6)+scale_color_brewer(name="",palette = "Paired")+xlab("Number of genes")+ylab("p-value")
ggsave('figure/xenoturbella-biological.pdf',width=7.26  ,height=6.2)

ggplot(data=ht2[ht2$V2>4 & ht2$V1=="Avian" & 
                  ! ht2$V4.y %in% c("Core Waterbirds+Phaethontimorphae","Neoaves-Columbea-Otidimorphae","Columbimorphae+Otidimorphae"),],
       aes(x=V8*V2/100,y=V5,color=V4.y))+geom_point()+theme_bw()+facet_wrap(~V1)+
  #scale_y_continuous(labels = function(x){print(x);1-10^(x)},breaks = c(0,log10(0.95),log10(0.5),log10(0.9)))+
  stat_smooth(se=F)+theme(legend.position = "bottom")+
  geom_hline(yintercept = 0.05,size=0.6)+scale_color_brewer(name="",palette = "Paired")+xlab("Number of genes")+ylab("p-value")+
  guides(colour=guide_legend(nrow=3))
ggsave('figure/avian-biological.pdf',width=7.26 ,height= 6.67)
