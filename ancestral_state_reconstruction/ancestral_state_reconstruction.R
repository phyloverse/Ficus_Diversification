##############################
#Manuscript: "How to be a big genus? Ficus L. as an emerging model"
#Authors: Nicole Mitidieri-Rivera, Elliot Gardner, Julianne Phipps, Nyree Zerega, Leandro C. Pederneiras, Alexander Damián-Parizaca, Kenneth J. Sytsma
#For correspondence. mitidieririv@wisc.edu
##############################

#This script was modified from Revell & Harmon (2022) to run Ancestral State Reconstructions in Mitidieri et al. (2025)

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
# ASR analyses on Moraceae and the outgroup using the phylogeny from Gardner et al. (2023)

#read and (optionally) plot the tree
mor.tree<-read.tree(file="mor_and_outgroup_520.tre")
mor.tree
#png("moraceae_tree.png", width = 4, height = 4, units = 'in', res = 700)
#plotTree(mor.tree,ftype="i", fsize=0.3)
#dev.off()

###################################################################
# ANCESTRAL STATE RECONSTRUCTION 

#Replace "trait" with the character of interest (e.g., phyllotaxis) and its actual file name

trait<-read.csv(file="trait.csv", row.names = 1)
trait1<-setNames(infl[,1],row.names(trait))
trait2 <- as.factor(trait1)
cols<-setNames(c("green","purple","orange"),levels(trait2))

chk<-name.check(mor.tree,trait)
chk

#choose a character model of evolution

#fit ER model
#png("fitER_trait.png", width = 4, height = 4, units = 'in', res = 700)
fitER_trait<-fitMk(mor.tree,trait2,model="ER")
fitER_trait
plot(fitER_trait)
#dev.off()

#fit ARD model
#png("fitARD_infl.png", width = 4, height = 4, units = 'in', res = 700)
fitARD_trait<-fitMk(mor.tree,trait2,model="ARD")
fitARD_trait
plot(fitARD_trait)
#dev.off()

#fit SYM model
#png("fitSYM_trait.png", width = 4, height = 4, units = 'in', res = 700)
fitSYM_trait<-fitMk(mor.tree,trait2,model="SYM")
fitSYM_trait
plot(fitSYM_trait)
#dev.off()

#likelihood-ratio test comparing ER & ARD
lrtest(fitER_trait,fitARD_trait)

#likelihood-ratio test comparing ER & SYM
lrtest(fitER_trait,fitSYM_trait)

#likelihood-ratio test comparing SYM & ARD
lrtest(fitSYM_trait,fitARD_trait)

#extract AIC values for each model
aic<-c(AIC(fitER_trait),AIC(fitARD_trait),AIC(fitSYM_trait))
aicw <- round(aic.w(aic),digits=2)
data.frame(model=c("ER","ARD","SYM"),
           logL=c(logLik(fitER_trait),logLik(fitARD_trait),logLik(fitSYM_trait)),
           AIC=aic,delta.AIC=aic-min(aic))

#ASR using marginal state reconstruction

#estimate marginal ancestral states under a ARD model
fit.marginal<-ancr(fitARD_trait,type="marginal")
fit.marginal

cols<-setNames(viridisLite::viridis(n=3),levels(infl2))
png("marginal_asr_trait.png", width = 20, height = 20, units = 'in', res = 700)
plotTree.datamatrix(mor.tree,as.data.frame(trait2),
                    colors=list(cols),header=FALSE, fsize=0.25)
legend("topright",legend=levels(trait2),pch=22,
       pt.cex=1.5,pt.bg=cols,bty="n",cex=0.7)
nodelabels(pie=fit.marginal$ace,piecol=cols,cex=0.25)
dev.off()

#ASR using a MCMC approach (Huelsenbeck et al. 2003)

#generate one stochastic character history
mtree<-make.simmap(mor.tree,trait2,model="ARD")

#plot single stochastic map
plot(mtree,cols,fsize=0.4,ftype="i",lwd=2,offset=0.4,
    ylim=c(-1,Ntip(mor.tree)))

#generate 1,000 stochastic character maps
mtrees<-make.simmap(mor.tree,trait2,model="ARD",nsim=10000) 
mtrees.trait<-make.simmap(mor.tree,trait2,model="ARD",nsim=10000)
pd.trait<-summary(mtrees.trait)

#compute posterior probabilities at nodes
pd<-summary(mtrees)
pd

#create a plot showing PP at all nodes of the tree
#png("stochastic_trait.png", width = 20, height = 20, units = 'in', res = 700)
plot(pd,colors=cols,fsize=0.1,ftype="i",lwd=1,
     offset=0.1,ylim=c(-1,Ntip(mor.tree)),
     cex=c(0.3,0.1))

#add a legend
legend("bottomleft",legend=levels(trait2),pch=22,
       pt.cex=1.7,pt.bg=cols,bty="n",cex=0.8)
#dev.off()

#map a multi-state discrete character onto the edges of a tree using stochastic mapping & transparent colors
#create a visualization of the posterior sample

#png("densityMap_trait.png", width = 20, height = 20, units = 'in', res = 700)
plotTree(mor.tree,ftype="i",fsize=0.1,offset=0.1,
         lwd=1)
par(fg="transparent",lend=1)
plotTree(mor.tree,ftype="i",fsize=0.1,offset=0.1,
         lwd=1,color="white",add=TRUE)
#now plot our 100 stochastic map trees with 99% transparency
for(i in 1:length(mtrees)) plot(mtrees[[i]],
                                   colors=sapply(cols,make.transparent,alpha=0.01),
                                   add=TRUE,lwd=1,ftype="i",fsize=0.1,offset=0.1)
par(fg="black")
nodelabels(pie=summary(mtrees)$ace,piecol=cols,
           cex=0.2)
legend(x="bottomleft",levels(trait2),pch=22,
       pt.bg=cols,pt.cex=1.5,bty="n",cex=0.7)
dev.off()