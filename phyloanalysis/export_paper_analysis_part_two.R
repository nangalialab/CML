### Part Two
## Generate tree and related main figures for the paper

source("load_and_annotate_tree.R")
library("TeachingDemos")
source("paper_plot_functions.R")
source("revised_sigs.R")
source("cml_sigs.R")
VERSION=readLines("../export/VERSION.txt")
library("tidyverse")
library("ggpubr")
library("scales")
library("ggrepel")

PDD=readRDS("../cache/PDDxs.RDS")
## 4c and d
scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))
}
scientific2 <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("1e", " 10^", scientific_format()(x)))))
}
## Pull in data per donor
dfcomb=read.table(sprintf("../export/per_donor_summary_v%s.txt",VERSION),head=TRUE,stringsAsFactors = FALSE,comment.char = "")
dfcomb=dfcomb[order(dfcomb$t_diagnosis),]
dfcomb=dfcomb %>% mutate(Donor=factor(Donor,levels=Donor)) %>% mutate(label2=sprintf("%s(%syrs)",Donor,age_at_diagnosis)) %>% mutate(label2=factor(label2,levels=label2))
dfcomb=dfcomb %>% mutate(label=sprintf("%d",age_at_diagnosis)) %>% mutate(label=factor(label,levels=label))

#selres=read.table(sprintf("../export/all_growth_rates_v%s.txt",VERSION),head=TRUE,stringsAsFactors = FALSE,comment.char = "")
donor.exclude="PD57333"
dfcomb=dfcomb %>% mutate(doubling_time=log(2)/(log(1+S_ltt_median)))
model=lm(data=dfcomb,log(log(1+S_ltt_median))~log(latency_median))
model.ticks=seq(3,16,0.01)
df=data.frame(latency_median=model.ticks)
model.predict=predict(model,newdata=df,interval="confidence")
model.predict=cbind(df,model.predict) %>% mutate(S_lb=exp(exp(lwr))-1,S=exp(exp(fit))-1,S_ub=exp(exp(upr))-1) %>%
  mutate(td_lb=log(2)/log(1+S_ub),td=log(2)/log(1+S),td_ub=log(2)/log(1+S_lb))
print(summary(model))
fig4c=ggplot(dfcomb %>% filter(!is.na(S_ltt_median)),aes(x=latency_median,y=100*S_ltt_median,col=Donor))+
  scale_y_log10(label=scientific2,breaks=10**(seq(0,8,1)),limits=c(50,5e8))+
  theme_bw()+annotation_logticks(sides = "l")+ylab("Growth Rate(%)")+xlab("Latency (Years)")+
  geom_ribbon(data=model.predict,aes(x=latency_median,ymin=100*S_lb,ymax=100*S_ub,y=S),alpha=0.5,fill="lightgrey",color=NA)+
  geom_line(data=model.predict,aes(x=latency_median,y=100*S),color="black")+geom_point(size=3)+
  geom_text(aes(label=label2),hjust="left",nudge_x =0.15,show.legend=FALSE,size=3.88)+
  scale_color_manual(values=get_colors(dfcomb %>% filter(!(Donor %in% donor.exclude )),"Donor"))+scale_x_continuous(breaks=c(4,8,12,16),minor_breaks=seq(3,16,1),limits=c(3,16))
figure4d_doubling=ggplot(dfcomb %>% filter(!is.na(S_ltt_median)),aes(x=latency_median,y=doubling_time*12,col=Donor))+
  geom_ribbon(data=model.predict,aes(x=latency_median,ymin=12*td_lb,ymax=12*td_ub,y=12*td),alpha=0.5,fill="lightgrey",color=NA)+
  geom_line(data=model.predict,aes(x=latency_median,y=12*td),color="black")+geom_point(size=3)+
  geom_text_repel(aes(label=label2),show.legend=FALSE)+
  scale_color_manual(values=get_colors(dfcomb %>% filter(!(Donor %in% donor.exclude )),"Donor"))+theme_bw()+
  xlab("Latency (Years)")+ylab("Doubling Time (Months)")+scale_y_continuous(breaks=seq(0,15,3),minor_breaks = seq(1,15))+
  scale_x_continuous(breaks=c(4,8,12,16),minor_breaks=seq(3,16,1))+coord_cartesian(ylim = c(-2,15), xlim = c(3,16))
pdf(sprintf("../figures/figure_4c_growth_vs_latency_v%s.pdf",VERSION),w=8,height=6)
print(fig4c)
dev.off()
pdf(sprintf("../figures/figure_4d_doubling_vs_latency_v%s.pdf",VERSION),w=8,height=6)
print(figure4d_doubling)
dev.off()

pa=ggplot(dfcomb %>% filter(!(Donor %in% donor.exclude )),aes(x=Donor,y=latency_median,fill=label2))+geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=latency_lb,ymax=latency_ub))+ylab("Latency to Diagnosis (Years)")+xlab("")+theme_bw()+ggtitle("Latency to CML Diagnosis")+
  scale_fill_manual(values=get_colors(dfcomb %>% filter(!(Donor %in% donor.exclude )),field="label2"))+guides(fill=guide_legend(title="Donor(Age @ Diagnosis)"))

pdf(sprintf("../figures/figure_4a_latency_%s.pdf",VERSION),w=8,h=6)
print(pa)
dev.off()

BRK=c(10,100,1000,10000,100000,1e6,1e7,1e8)
#BRKL=c("10","100","1,000","10,000","100,000","1,000,000","10,000,000","100,000,000")

s2b=ggplot(dfcomb %>% mutate(S=100*S_ltt_median) %>% filter(!is.na(S)),
           aes(x=label2,y=S,fill=label2))+geom_bar(stat="identity")+geom_errorbar(aes(ymin=100*S_ltt_lb,ymax=100*S_ltt_ub))+theme_bw()+ylab("Growth Rate (%)")+
  scale_y_log10(breaks=BRK,labels=scientific2)+ggtitle("CML Clones: Estimated Annual Growth Rates")+xlab("")+
  scale_fill_manual(values=get_colors(dfcomb  %>% filter(!is.na(S_ltt_median)),field = "label2"))+
  guides(fill=guide_legend(title="Donor(Age at Diagnosis)"))

pdf(sprintf("../figures/figure_4b_v%s.pdf",VERSION),w=8,h=6)
print(s2b)
dev.off()
### Tree plots (Figure 2)

# Present on appropriate order..
DONORS=c("PD57334", "PD56961","PD51633", "PD57335" ,"PD51634" , "PD51632" ,"PD51635" ,"PD57333", "PD57332" )
PDD2=lapply(PDD,function(PD){PD$pdx$tree_ml=PD$fit$poisson_tree$altmodel$ultratree;PD})
plot_all_trees_scaled_v2(PDD[DONORS],bw=0.01,mar = c(1,2,1,3),expand.top.tree = "mt",figure.path=sprintf("../figures/figure_2_v%s.pdf",VERSION))
plot_all_trees_scaled_v2(PDD2[DONORS],bw=0.01,mar = c(1,2,1,3),m=c(55,60,85),expand.top.tree = "ultra",figure.path=sprintf("../figures/figure_2_ultra_v%s.pdf",VERSION))

pdf(sprintf("../figures/figure_3b_c_embedded_and_extended_figure_3b_sigs_v%s.pdf",VERSION),w=16,h=18,pointsize=11)
par(mfrow=c(5,2),oma=c(2,2,2,2))
for(PD in PDD[1:9]){
  plot_cml_sig_tree_posterior(PD,b.ultra = FALSE)
}
dev.off()

pdf(sprintf("../figures/figure_3d_embedded_and_extended_figure_4_ctcpg_v%s.pdf",VERSION),w=16,h=18,pointsize=12)
par(mfrow=c(5,2),oma=c(2,2,2,2))
for(PD in PDD[1:9]){
  plot_cml_sig_tree_posterior(PD,b.add.sigs = FALSE,b.ultra = TRUE)
}
dev.off()

## Supplementary note 3 figures.
