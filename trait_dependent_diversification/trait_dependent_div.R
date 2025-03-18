##############################
#Manuscript: "How to be a big genus? Ficus L. as an emerging model"
#Authors: Nicole Mitidieri-Rivera, Elliot Gardner, Julianne Phipps, Nyree Zerega, Leandro C. Pederneiras, Alexander Damián-Parizaca, Kenneth J. Sytsma
#For correspondence. mitidieririv@wisc.edu
##############################

#This script was modified from Revell & Harmon (2022) to run trait-dependent diversification analyses in Mitidieri et al. (2025)

###################################################################
#Initialization
rm(list = ls())

#getwd() 
#setwd()

#install and call the appropriate R packages for phylogenetic analyses.

required.libraries <- c("phytools","geiger","ape","phangorn","corHMM", "lmtest",
                        "diversitree", "BAMMtools", "hisse", "devtools","RColorBrewer")
needed.libraries <- required.libraries[!(required.libraries %in% installed.packages()[,"Package"])]
if(length(needed.libraries)) install.packages(needed.libraries)

#load all required libraries at once
lapply(required.libraries, require, character.only = TRUE)

###################################################################
# Trait-dependent analyses on Moraceae and the outgroup using the phylogeny from Gardner et al. (2023)

#read the tree
mor.tree<-read.tree(file="mor_and_outgroup_520.tre")
mor.tree

#coerse a phylogenetic tree to be ultrametric
mor.tree <- read.tree(file="mor508.tre")
mor.tree<-force.ultrametric(mor.tree)

invol<-read.csv(file="~/Desktop/big_genera_ficus/trait_dependent_diversification/data/raw/bisse_hisse/involucre.csv", row.names = 1)
invol1<-setNames(invol[,1],row.names(invol))

chk<-name.check(mor.tree,invol)
chk

#set a sampling fraction based on the number of species within each state
#represented in the true phylogeny, not in the reconstructed phylogeny
rho.bisse.inv<-setNames(c(0.293,0.70), #1 = 29.3%; 2 = 70%
                        1:2) 

#make BiSSE likelihood function
bisse.model.inv <- make.bisse(mor.tree,invol1,sampling.f = rho.bisse.inv)

#find reasonable parameter values for optimization
p <- starting.point.bisse(mor.tree)
p

#optimize BiSSE model
bisse.mle.inv <- find.mle(bisse.model.inv,p)
bisse.mle.inv
coef(bisse.mle.inv)

#create constrained null model 
bissenull.model.inv <- constrain(bisse.model.inv,
                                 lambda1~lambda0,mu1~mu0)

#optimize null model
bissenull.mle.inv <- find.mle(bissenull.model.inv,
                              p[c(-2,-4)])

coef(bissenull.mle.inv)
logLik(bissenull.mle.inv)

#run likelihood-ratio test
bisseAnova.inv <- anova(null=bissenull.mle.inv, bisse.mle.inv)
bisseAnova.inv # The likelihood for the full BiSSE model is sufficiently improved compared to the null, and the P-value for our likelihood ratio test is significant (0.008449)

#pull Akaike Information Criterion (AIC) values from this table and compute Akaike weights
aicw(setNames(bisseAnova.inv$AIC,
              rownames(bisseAnova.inv))) #the AIC weight of the BiSSE model is larger than the null, indicating relatively high weight of evidence in support of the state-dependent model.

#Analyzing a BiSSE model using Bayesian MCMC
prior.i <- make.prior.exponential(1/(2*0.4))
prior.i

#run Bayesian MCMC 
#(produces a posterior sample for each of the six parameters in our model)
bisse.mcmc.i <- mcmc(bisse.model.inv,bisse.mle.inv$par,
                     nsteps=10000, prior=prior.i,w=0.1,
                     print.every=100)

#Visualize the posterior distribution of speciation & extinction in our model
png("bisse_mcmc_invol.png", width = 8, height = 5, units = 'in', res = 700)
#subdivide plot and set margins
par(mfrow=c(1,2),mar=c(5.1,4.1,3.1,2.1))
#set colors for plotting
col <- setNames(c("blue","red"),
                c("no involucre","involucre"))
#create graph of posterior sample for lambda
profiles.plot(bisse.mcmc.i[,c("lambda0","lambda1")],
              col.line = col,las=1,bty="n",
              xlab=expression(lambda),cex.axis=0.7)
#add legend & panel label
legend("topright",names(col),pch=15,col=col,
       pt.cex=1.5,bty="n",cex=0.7)
mtext("a)",line=0.5,adj=0)
#create graph of posterior sample for mu
profiles.plot(bisse.mcmc.i[,c("mu0","mu1")],
              col.line=col,las=1,bty="n",
              xlab=expression(mu),cex.axis=0.7)
#add legend & panel label
legend("topright",names(col),pch=15,col=col,
       pt.cex=1.5,bty="n",cex=0.7)
mtext("b)",line=0.5,adj=0)
dev.off()

#Visualize the posterior distribution of net diversification in our model
png("bisse_mcmc_net_diversification_invol.png", width = 5, height = 5, units = 'in', res = 700)
net.div<-bisse.mcmc.i[,grep("lambda",colnames(bisse.mcmc.i))]-
  bisse.mcmc.i[,grep("mu",colnames(bisse.mcmc.i))]
colnames(net.div)<-paste("lambda-mu(",1:2,")",sep="")
profiles.plot(net.div, 
              xlab="Net diversification rate", ylab="Probability density",
              legend.pos="topright",col.line=setNames(col,colnames(net.div)),
              lty=1)
dev.off()

#Calculate the posterior probability that lambda1 > lambda0 or mu1 > mu0
#lambda
sum(bisse.mcmc.i$lambda1>bisse.mcmc.i$lambda0)/length(bisse.mcmc.i$lambda1) #0.0044

#mu
sum(bisse.mcmc.i$mu1>bisse.mcmc.i$mu0)/length(bisse.mcmc.i$mu1) #0.0019

#HiSSE (Beaulieu & O'Meara, 2016)
#create input data frame for hisse
hs<-read.csv(file="~/Desktop/big_genera_ficus/trait_dependent_diversification/data/raw/bisse_hisse/involucre.csv", row.names=1)
hs1<-data.frame(Genus.species=rownames(hs),
                x=hs[,"involucre"])
head(hs1)
#hs1<-setNames(invol[,1],rownames(invol))

#create HiSSE design matrix
rates.hisse<-TransMatMakerHiSSE(hidden.traits=1)
rates.hisse

#create hisse design matrix for BiSSE model
rates.bisse<-TransMatMakerHiSSE(hidden.traits = 0)
#fit BiSSE model using HiSSE
bisse.hmle<-hisse(mor.tree,hs1,turnover=c(1,2),
                  eps=c(1,2),hidden.states = FALSE,
                  trans.rate = rates.bisse)
bisse.hmle

#custom function to back-transform turnover and extinction-fraction 
#to lambda & mu
repar.bd<-function(object,k=2){
  pars<-object$solution
  tt<-pars[grep("turnover",names(pars))] [1:k]
  ee<-pars[grep("eps",names(pars))] [1:k]
  lambda<-tt/(1+ee)
  mu<-tt-lambda
  nn<-sapply(strsplit(names(tt),"turnover"),
             function(x) x[2])
  matrix(c(lambda,mu),k,2,dimnames=list(nn,
                                        c("lambda","mu")))
}
repar.bd(bisse.hmle)

#fit CID (character-independent diversification) using hisse; this is the same as 
#the constant-rate birth-death model.
cid.mle<-hisse(mor.tree,hs1,turnover = c(1,1),
               eps=c(1,1),hidden.states = FALSE,
               trans.rate = rates.bisse)
cid.mle #this fitted model parameters match what we obtained for our null model using diversitree
#and those from fit.bd in phytools (see fit.bd(cid.mle,1) below to confirm)

#fit birth-death model using phytools
fit.bd(mor.tree) #lambda = 0.07389 

#reparameterize CID model in terms of lambda and mu
repar.bd(cid.mle,1) #lambda = 0.07389121

#create CID-2 design matrix; in this model, there's no influence of our measured trait
#on diversification but there is an influence of a hidden character
rates.cid2<-rates.hisse
rates.cid2[!is.na(rates.cid2)]<-1
rates.cid2

#fit CID-2 model using hisse
cid2.mle<-hisse(mor.tree,hs1,f=c(1,1),
                turnover=c(1,1,2,2),eps=c(1,1,2,2),
                hidden.states = TRUE,trans.rate = rates.cid2)
cid2.mle

#reparameterize to lambda & mu
repar.bd(cid2.mle,k=4)

#let's create a four-level hidden-rate model, 
#independent of the value of any discrete character state
#create design matrix for CID-4 model 
rates.cid4<-TransMatMakerHiSSE(hidden.traits = 3)
rates.cid4[!is.na(rates.cid4)]<-1
rates.cid4

#fit CID-4 model
cid4.mle<-hisse(mor.tree,hs1,f=c(1,1),
                turnover = c(1,1,2,2,3,3,4,4),
                eps=c(1,1,2,2,3,3,4,4),
                hidden.states = TRUE,
                trans.rate = rates.cid4)
cid4.mle

#reparameterize to lambda & mu
repar.bd(cid4.mle,8)

#fit full HiSSE model
hisse.mle<-hisse(mor.tree,hs1,f=c(1,1),
                 hidden.states = TRUE,
                 turnover = c(1,2,3,4,5,6,7,8),
                 eps=c(1,2,3,4,5,6,7,8),
                 trans.rate = rates.cid4)
hisse.mle
repar.bd(hisse.mle,8)

#our logLik methods
logLik.hisse.fit<-function(x,...){
  lik<-x$loglik
  attr(lik,"df")<-(x$AIC+2*lik)/2
  lik
}
#print a table of results
data.frame(
  model=c("CID","BiSSE","HiSSE CID-2","HiSSE CID-4","HiSSE full"),
  logL=sapply(list(cid.mle,bisse.hmle,
                   cid2.mle,cid4.mle,hisse.mle),
              logLik),
  k=sapply(list(cid.mle,bisse.hmle,
                cid2.mle,cid4.mle,hisse.mle),
           function(x) attr(logLik(x),"df")),
  AIC=aic<-sapply(list(cid.mle,bisse.hmle,
                       cid2.mle,cid4.mle,hisse.mle),
                  AIC),
  Akaike.weight=unclass(aic.w(aic))
)

