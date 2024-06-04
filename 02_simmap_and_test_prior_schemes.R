library(phytools)
library(geiger)
setwd('/Users/rs155/Dropbox/my_research/mating-systems/REVISION/unrandomized_analyses')

tr <- read.tree("/Users/rs155/Dropbox/my_research/mating-systems/MCC_of_100_hackett_trees.tre")

dat <- read.csv("/Users/rs155/Dropbox/my_research/mating-systems/REVISION/unrandomized_analyses/0randdata.csv")
#dat$Mating_system <- gsub(x = dat$Mating_system, pattern = " ", replacement = "")
rownames(dat) <- dat$Species
rownames(dat) <- gsub(x = rownames(dat), pattern = " ", replacement = "_")
dat=dat[,c(2:4)]
#colnames(dat)[2]='Lit_mating'

prior_schemes=read.csv('/Users/rs155/Dropbox/my_research/mating-systems/REVISION/prior_schemes.csv')
row.names(prior_schemes)=NULL
colnames(prior_schemes)=prior_schemes[1,]
ps1=prior_schemes[2:13,]
row.names(ps1)=NULL
ps2=prior_schemes[16:27,]
row.names(ps2)=NULL
ps3=prior_schemes[30:41,]
row.names(ps3)=NULL


#######PRIOR SCHEME 1########

####assign priors#####
priormatrix1 <- as.data.frame(dat)
priormatrix1 <- priormatrix1[c("Lit_mating", "Data_quality")]
priormatrix1$M=rep(NA,nrow(priormatrix1))
priormatrix1$P=rep(NA,nrow(priormatrix1))
priormatrix1$L=rep(NA,nrow(priormatrix1))


for (sp in rownames(priormatrix1)){
  #cat(sp, '\n')
  if (priormatrix1[sp,'Lit_mating']=='U'){priormatrix1[sp,3:5]=ps1[11,3:5]}
  if (priormatrix1[sp,'Lit_mating']=='PU'){priormatrix1[sp,3:5]=ps1[10,3:5]}
  
  if (priormatrix1[sp,'Lit_mating']=='M' & priormatrix1[sp,'Data_quality']==3){priormatrix1[sp,3:5]=ps1[1,3:5]}
  if (priormatrix1[sp,'Lit_mating']=='P' & priormatrix1[sp,'Data_quality']==3){priormatrix1[sp,3:5]=ps1[2,3:5]}
  if (priormatrix1[sp,'Lit_mating']=='L' & priormatrix1[sp,'Data_quality']==3){priormatrix1[sp,3:5]=ps1[3,3:5]}
  
  if (priormatrix1[sp,'Lit_mating']=='M' & priormatrix1[sp,'Data_quality']==2){priormatrix1[sp,3:5]=ps1[4,3:5]}
  if (priormatrix1[sp,'Lit_mating']=='P' & priormatrix1[sp,'Data_quality']==2){priormatrix1[sp,3:5]=ps1[5,3:5]}
  if (priormatrix1[sp,'Lit_mating']=='L' & priormatrix1[sp,'Data_quality']==2){priormatrix1[sp,3:5]=ps1[6,3:5]}
  
  if (priormatrix1[sp,'Lit_mating']=='M' & priormatrix1[sp,'Data_quality']==1){priormatrix1[sp,3:5]=ps1[7,3:5]}
  if (priormatrix1[sp,'Lit_mating']=='P' & priormatrix1[sp,'Data_quality']==1){priormatrix1[sp,3:5]=ps1[8,3:5]}
  if (priormatrix1[sp,'Lit_mating']=='L' & priormatrix1[sp,'Data_quality']==1){priormatrix1[sp,3:5]=ps1[9,3:5]}
}

#####simmap with PS1#####

priormatrix1=priormatrix1[,3:5]
priormatrix1[,1]=as.numeric(priormatrix1[,1])
priormatrix1[,2]=as.numeric(priormatrix1[,2])
priormatrix1[,3]=as.numeric(priormatrix1[,3])
td <- treedata(phy = tr, data = priormatrix1)


Sys.time()
m.ard.ps1 <- make.simmap(tree = td$phy, x = td$data, model = "ARD", nsim = 1000, pi = "estimated")
Sys.time()
md.ard.ps1 <- describe.simmap(m.ard.ps1)
Sys.time()


#######PRIOR SCHEME 2########

####assign priors#####
priormatrix2 <- as.data.frame(dat)
priormatrix2 <- priormatrix2[c("Lit_mating", "Data_quality")]
priormatrix2$M=rep(NA,nrow(priormatrix2))
priormatrix2$P=rep(NA,nrow(priormatrix2))
priormatrix2$L=rep(NA,nrow(priormatrix2))


for (sp in rownames(priormatrix2)){
  #cat(sp, '\n')
  if (priormatrix2[sp,'Lit_mating']=='U'){priormatrix2[sp,3:5]=ps2[11,3:5]}
  if (priormatrix2[sp,'Lit_mating']=='PU'){priormatrix2[sp,3:5]=ps2[10,3:5]}
  
  if (priormatrix2[sp,'Lit_mating']=='M' & priormatrix2[sp,'Data_quality']==3){priormatrix2[sp,3:5]=ps2[1,3:5]}
  if (priormatrix2[sp,'Lit_mating']=='P' & priormatrix2[sp,'Data_quality']==3){priormatrix2[sp,3:5]=ps2[2,3:5]}
  if (priormatrix2[sp,'Lit_mating']=='L' & priormatrix2[sp,'Data_quality']==3){priormatrix2[sp,3:5]=ps2[3,3:5]}
  
  if (priormatrix2[sp,'Lit_mating']=='M' & priormatrix2[sp,'Data_quality']==2){priormatrix2[sp,3:5]=ps2[4,3:5]}
  if (priormatrix2[sp,'Lit_mating']=='P' & priormatrix2[sp,'Data_quality']==2){priormatrix2[sp,3:5]=ps2[5,3:5]}
  if (priormatrix2[sp,'Lit_mating']=='L' & priormatrix2[sp,'Data_quality']==2){priormatrix2[sp,3:5]=ps2[6,3:5]}
  
  if (priormatrix2[sp,'Lit_mating']=='M' & priormatrix2[sp,'Data_quality']==1){priormatrix2[sp,3:5]=ps2[7,3:5]}
  if (priormatrix2[sp,'Lit_mating']=='P' & priormatrix2[sp,'Data_quality']==1){priormatrix2[sp,3:5]=ps2[8,3:5]}
  if (priormatrix2[sp,'Lit_mating']=='L' & priormatrix2[sp,'Data_quality']==1){priormatrix2[sp,3:5]=ps2[9,3:5]}
}

#####simmap with PS2#####

priormatrix2=priormatrix2[,3:5]
priormatrix2[,1]=as.numeric(priormatrix2[,1])
priormatrix2[,2]=as.numeric(priormatrix2[,2])
priormatrix2[,3]=as.numeric(priormatrix2[,3])
td <- treedata(phy = tr, data = priormatrix2)


Sys.time()
m.ard.ps2 <- make.simmap(tree = td$phy, x = td$data, model = "ARD", nsim = 1000, pi = "estimated")
Sys.time()
md.ard.ps2 <- describe.simmap(m.ard.ps2)
Sys.time()


#######PRIOR SCHEME 3########

####assign priors#####
priormatrix3 <- as.data.frame(dat)
priormatrix3 <- priormatrix3[c("Lit_mating", "Data_quality")]
priormatrix3$M=rep(NA,nrow(priormatrix3))
priormatrix3$P=rep(NA,nrow(priormatrix3))
priormatrix3$L=rep(NA,nrow(priormatrix3))


for (sp in rownames(priormatrix3)){
  #cat(sp, '\n')
  if (priormatrix3[sp,'Lit_mating']=='U'){priormatrix3[sp,3:5]=ps3[11,3:5]}
  if (priormatrix3[sp,'Lit_mating']=='PU'){priormatrix3[sp,3:5]=ps3[10,3:5]}
  
  if (priormatrix3[sp,'Lit_mating']=='M' & priormatrix3[sp,'Data_quality']==3){priormatrix3[sp,3:5]=ps3[1,3:5]}
  if (priormatrix3[sp,'Lit_mating']=='P' & priormatrix3[sp,'Data_quality']==3){priormatrix3[sp,3:5]=ps3[2,3:5]}
  if (priormatrix3[sp,'Lit_mating']=='L' & priormatrix3[sp,'Data_quality']==3){priormatrix3[sp,3:5]=ps3[3,3:5]}
  
  if (priormatrix3[sp,'Lit_mating']=='M' & priormatrix3[sp,'Data_quality']==2){priormatrix3[sp,3:5]=ps3[4,3:5]}
  if (priormatrix3[sp,'Lit_mating']=='P' & priormatrix3[sp,'Data_quality']==2){priormatrix3[sp,3:5]=ps3[5,3:5]}
  if (priormatrix3[sp,'Lit_mating']=='L' & priormatrix3[sp,'Data_quality']==2){priormatrix3[sp,3:5]=ps3[6,3:5]}
  
  if (priormatrix3[sp,'Lit_mating']=='M' & priormatrix3[sp,'Data_quality']==1){priormatrix3[sp,3:5]=ps3[7,3:5]}
  if (priormatrix3[sp,'Lit_mating']=='P' & priormatrix3[sp,'Data_quality']==1){priormatrix3[sp,3:5]=ps3[8,3:5]}
  if (priormatrix3[sp,'Lit_mating']=='L' & priormatrix3[sp,'Data_quality']==1){priormatrix3[sp,3:5]=ps3[9,3:5]}
}

#####simmap with ps3#####

priormatrix3=priormatrix3[,3:5]
priormatrix3[,1]=as.numeric(priormatrix3[,1])
priormatrix3[,2]=as.numeric(priormatrix3[,2])
priormatrix3[,3]=as.numeric(priormatrix3[,3])
td <- treedata(phy = tr, data = priormatrix3)


Sys.time()
m.ard.ps3 <- make.simmap(tree = td$phy, x = td$data, model = "ARD", nsim = 1000, pi = "estimated")
Sys.time()
md.ard.ps3 <- describe.simmap(m.ard.ps3)
Sys.time()

save.image(file="/Users/rs155/Dropbox/my_research/mating-systems/REVISION/02_simmap_and_test_prior_schemes_REZ.Rdata")
#load(file="/Users/rs155/Dropbox/my_research/mating-systems/REVISION/02_simmap_and_test_prior_schemes_REZ.Rdata")


####set tips to most likely state in results####

ps1.tips=as.data.frame(md.ard.ps1$tips)
ps1.tips$state=rep(x=NA)
for (t in 1:nrow(ps1.tips)) {
  s <- which.max(ps1.tips[t, ])
  ps1.tips[t,4] <- colnames(ps1.tips)[s]
}
colnames(ps1.tips)=paste('ps1',colnames(ps1.tips),sep='.')


ps2.tips=as.data.frame(md.ard.ps2$tips)
ps2.tips$state=rep(x=NA)
for (t in 1:nrow(ps2.tips)) {
  s <- which.max(ps2.tips[t, ])
  ps2.tips[t,4] <- colnames(ps2.tips)[s]
}
colnames(ps2.tips)=paste('ps2',colnames(ps2.tips),sep='.')


ps3.tips=as.data.frame(md.ard.ps3$tips)
ps3.tips$state=rep(x=NA)
for (t in 1:nrow(ps3.tips)) {
  s <- which.max(ps3.tips[t, ])
  ps3.tips[t,4] <- colnames(ps3.tips)[s]
}
colnames(ps3.tips)=paste('ps3',colnames(ps3.tips),sep='.')

####merge everything####

rez=merge(dat,ps1.tips,all=F,by=0)
rez$Species=NULL
rez=merge(x=rez,y=ps2.tips,by.x=1,by.y=0)
rez=merge(x=rez,y=ps3.tips,by.x=1,by.y=0)
rownames(rez)=rez$Row.names



#######find instances where PSs disagree#######

#lit X ps1
litVps1=as.data.frame(matrix(nrow=0,ncol=2))
  
for (sp in rez$Row.names){
  if (rez[sp,'Lit_mating']!='U' && rez[sp,'Lit_mating']!='PU'){
  if (rez[sp,'Lit_mating']!=rez[sp,'ps1.state']){
  litVps1[nrow(litVps1)+1,]=rez[sp,c(1,8)]

  }}}
colnames(litVps1)=c('sp','ps1')


#lit X ps2
litVps2=as.data.frame(matrix(nrow=0,ncol=2))
  
for (sp in rez$Row.names){
  if (rez[sp,'Lit_mating']!='U' && rez[sp,'Lit_mating']!='PU'){
  if (rez[sp,'Lit_mating']!=rez[sp,'ps2.state']){
  litVps2[nrow(litVps2)+1,]=rez[sp,c(1,12)]

  }}}
colnames(litVps2)=c('sp','ps2')


#lit X ps3
litVps3=as.data.frame(matrix(nrow=0,ncol=2))
  
for (sp in rez$Row.names){
  if (rez[sp,'Lit_mating']!='U' && rez[sp,'Lit_mating']!='PU'){
  if (rez[sp,'Lit_mating']!=rez[sp,'ps3.state']){
  litVps3[nrow(litVps3)+1,]=rez[sp,c(1,16)]

  }}}
colnames(litVps3)=c('sp','ps3')


litVps=merge(x=dat[,c(1:2,4)],y=litVps3,by=0,all=T)
litVps=merge(x=litVps,y=litVps2,by=1,all=T)
litVps=merge(x=litVps,y=litVps1,by=1,all=T)
litVps=litVps[,is.na(litVps$ps3)==T]
rows_with_values <- !is.na(litVps$ps3) | !is.na(litVps$ps2) | !is.na(litVps$ps1)
litVps=litVps[rows_with_values,]
litVps=litVps[,c(2,3,4,6,7,8)]


             
  




