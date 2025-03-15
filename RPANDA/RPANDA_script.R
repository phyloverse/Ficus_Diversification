##############################
#Manuscript: "How to be a big genus? Ficus L. as an emerging model"
#Authors: Nicole Mitidieri-Rivera, Elliot Gardner, Julianne Phipps, Nyree Zerega, Leandro C. Pederneiras, Alexander Damián-Parizaca, Kenneth J. Sytsma
#For correspondence. mitidieririv@wisc.edu
##############################

#This script was modified from Condamine et al. (2013) to run time- and temperature-dependent diversification analyses in Mitidieri et al. (2025)

#install.packages("phytools")
#install_github("JClavel/mvMORPH", build_vignettes = TRUE)
#install_github("hmorlon/PANDA")
library(phytools)
library(devtools)
library(mvMORPH)
library(RPANDA)
library(picante)
library(pspline)
library(qpcR)
library("DDD")
library(TreeSim)
library(geiger)
#library(TreePar)

tree <- read.tree(file="tree.tre")
tot_time<-max(node.age(tree)$ages)

# 1) BCST (Pure birth)
print("BCST")
f.lamb<-function(x,y){y[1]}
f.mu<-function(x,y){0}
lamb_par<-c(0.1)
mu_par<-c()
cst.lamb=T; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T

treei_BCST<-fit_bd(tree,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=508/1285,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")

print(treei_BCST)
plot_fit_bd(treei_BCST,tot_time)

# 2) BCST DCST (constant Birth-death)
print("BCST DCST")
f.lamb<-function(x,y){y[1]}
f.mu<-function(x,y){y[1]}
lamb_par<-c(treei_BCST$lamb_par[1])
mu_par<-c(0.01)
cst.lamb=T; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F

treei_BCSTDCST<-fit_bd(tree,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=508/1285,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")

print(treei_BCSTDCST)
plot_fit_bd(treei_BCSTDCST,tot_time)

###############################################
###### Time Dependence (exponential variation) ######
###############################################

# 3) BTimeVar EXPO
print("BTimeVar EXPO")
f.lamb<-function(x,y){y[1]*exp(y[2]*x)}
f.mu<-function(x,y){0}
lamb_par<-c(treei_BCSTDCST$lamb_par[1],0.01)
mu_par<-c()
cst.lamb=F; cst.mu=T; expo.lamb=T; expo.mu=F; fix.mu=T

treei_BTimeVar_EXPO<-fit_bd(tree,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=508/1285,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")

print(treei_BTimeVar_EXPO)
plot_fit_bd(treei_BTimeVar_EXPO,tot_time)

# 4) BTimeVar DCST EXPO
print("BTimeVar DCST EXPO")
f.lamb<-function(x,y){y[1]*exp(y[2]*x)}
f.mu<-function(x,y){y[1]}
lamb_par<-c(treei_BTimeVar_EXPO$lamb_par[1],treei_BTimeVar_EXPO$lamb_par[2])
mu_par<-c(treei_BCSTDCST$mu_par[1])
mu_par<-c(0.01)
cst.lamb=F; cst.mu=T; expo.lamb=T; expo.mu=F; fix.mu=F

treei_BTimeVarDCST_EXPO<-fit_bd(tree,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=508/1285,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")

print(treei_BTimeVarDCST_EXPO)
plot_fit_bd(treei_BTimeVarDCST_EXPO,tot_time)

# 5) BCST DTimeVar EXPO
print("BCST DTimeVar EXPO")
f.lamb<-function(x,y){y[1]}
f.mu<-function(x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(treei_BCSTDCST$lamb_par[1])
mu_par<-c(0.01,0.001)
cst.lamb=T; cst.mu=F; expo.lamb=F; expo.mu=T; fix.mu=F

treei_BCSTDTimeVar_EXPO<-fit_bd(tree,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=508/1285,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")

print(treei_BCSTDTimeVar_EXPO)
plot_fit_bd(treei_BCSTDTimeVar_EXPO,tot_time)

# 6) BTimeVar DTimeVar EXPO
print("BTimeVar DTimeVar EXPO")
f.lamb<-function(x,y){y[1]*exp(y[2]*x)}
f.mu<-function(x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(treei_BTimeVarDCST_EXPO$lamb_par[1],0.001)
mu_par<-c(0.05,0.001)
cst.lamb=F; cst.mu=F; expo.lamb=T; expo.mu=T; fix.mu=F

treei_BTimeVarDTimeVar_EXPO<-fit_bd(tree,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=508/1285,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")

print(treei_BTimeVarDTimeVar_EXPO)
plot_fit_bd(treei_BTimeVarDTimeVar_EXPO,tot_time)

##########################################
###### Time Dependence (linear variation) ######
##########################################

#Not included in the analyses

# 7) BTimeVar LIN
#print("BTimeVar LIN")
#f.lamb<-function(x,y){y[1]+(y[2]*x)}
#f.mu<-function(x,y){0}
#lamb_par<-c(0.01,0.001)
#mu_par<-c()
#cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T

#treei_BTimeVar_LIN<-fit_bd(tree,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=508/1285,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")

#print(treei_BTimeVar_LIN)
#plot_fit_bd(treei_BTimeVar_LIN,tot_time)

# 8) BTimeVar DCST LIN
#print("BTimeVar DCST LIN")
#f.lamb<-function(x,y){y[1]+(y[2]*x)}
#f.mu<-function(x,y){y[1]}
#lamb_par<-c(abs(treei_BTimeVar_LIN$lamb_par[1]),treei_BTimeVar_LIN$lamb_par[2])
#mu_par<-c(0.01)
#cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F

#treei_BTimeVarDCST_LIN<-fit_bd(tree,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=508/1285,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")

#print(treei_BTimeVarDCST_LIN)
#plot_fit_bd(treei_BTimeVarDCST_LIN,tot_time)

# 9) BCST DTimeVar LIN
#print("BCST DTimeVar LIN")

#f.lamb<-function(x,y){y[1]}
#f.mu<-function(x,y){y[1]+(y[2]*x)}
#lamb_par<-c(treei_BCSTDCST$lamb_par[1])
#mu_par<-c(0.001,0.001)
#cst.lamb=T; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F

#treei_BCSTDTimeVar_LIN<-fit_bd(tree,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=508/1285,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")

#print(treei_BCSTDTimeVar_LIN)
#plot_fit_bd(treei_BCSTDTimeVar_LIN,tot_time)

# 10) BTimeVar DTimeVar LIN
#print("BTimeVar DTimeVar LIN")
#f.lamb<-function(x,y){y[1]+(y[2]*x)}
#f.mu<-function(x,y){y[1]+(y[2]*x)}
#lamb_par<-c(abs(treei_BTimeVarDCST_LIN$lamb_par[1]),0.001)
#mu_par<-c(0.05,-0.001)
#cst.lamb=F; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F

#treei_BTimeVarDTimeVar_LIN<-fit_bd(tree,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=508/1285,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")

#print(treei_BTimeVarDTimeVar_LIN)
#plot_fit_bd(treei_BTimeVarDTimeVar_LIN,tot_time)

############# RESULTS ###########################################

results<-matrix(NA,6,8)
colnames(results)<-c("Models","NP","logL","AICc","Lambda","AlphaTime","Mu","BetaTime")

#Models
results[,1]<-c("BCST","BCSTDCST","BTimeVar_EXPO","BTimeVarDCST_EXPO","BCSTDTimeVar_EXPO","BTimeVarDTimeVar_EXPO")

#NP
results[1,2]<-1
results[2,2]<-2
results[3,2]<-2
results[4,2]<-3
results[5,2]<-3
results[6,2]<-4

#logL
results[1,3]<-round(treei_BCST$LH,3)
results[2,3]<-round(treei_BCSTDCST$LH,3)
results[3,3]<-round(treei_BTimeVar_EXPO$LH,3)
results[4,3]<-round(treei_BTimeVarDCST_EXPO$LH,3)
results[5,3]<-round(treei_BCSTDTimeVar_EXPO$LH,3)
results[6,3]<-round(treei_BTimeVarDTimeVar_EXPO$LH,3)

#AICc
results[1,4]<-round(treei_BCST$aicc,3)
results[2,4]<-round(treei_BCSTDCST$aicc,3)
results[3,4]<-round(treei_BTimeVar_EXPO$aicc,3)
results[4,4]<-round(treei_BTimeVarDCST_EXPO$aicc,3)
results[5,4]<-round(treei_BCSTDTimeVar_EXPO$aicc,3)
results[6,4]<-round(treei_BTimeVarDTimeVar_EXPO$aicc,3)

#Lambda0
results[1,5]<-round(abs(treei_BCST$lamb_par[1]),4)
results[2,5]<-round(abs(treei_BCSTDCST$lamb_par[1]),4)
results[3,5]<-round(abs(treei_BTimeVar_EXPO$lamb_par[1]),4)
results[4,5]<-round(abs(treei_BTimeVarDCST_EXPO$lamb_par[1]),4)
results[5,5]<-round(abs(treei_BCSTDTimeVar_EXPO$lamb_par[1]),4)
results[6,5]<-round(abs(treei_BTimeVarDTimeVar_EXPO$lamb_par[1]),4)

#Alpha Time
results[3,6]<-round(treei_BTimeVar_EXPO$lamb_par[2],5)
results[4,6]<-round(treei_BTimeVarDCST_EXPO$lamb_par[2],5)
results[6,6]<-round(treei_BTimeVarDTimeVar_EXPO$lamb_par[2],5)

#Mu0
results[2,7]<-round(abs(treei_BCSTDCST$mu_par[1]),5)
results[4,7]<-round(abs(treei_BTimeVarDCST_EXPO$mu_par[1]),5)
results[5,7]<-round(abs(treei_BCSTDTimeVar_EXPO$mu_par[1]),5)
results[6,7]<-round(abs(treei_BTimeVarDTimeVar_EXPO$mu_par[1]),5)

#Beta Time
results[5,8]<-round(treei_BCSTDTimeVar_EXPO$mu_par[2],5)
results[6,8]<-round(treei_BTimeVarDTimeVar_EXPO$mu_par[2],5)

#Akaike weights
all_AICc<-c(results[,4])
results[1,9]<-round(akaike.weights(as.numeric(all_AICc))$weights[1],3)
results[2,9]<-round(akaike.weights(as.numeric(all_AICc))$weights[2],3)
results[3,9]<-round(akaike.weights(as.numeric(all_AICc))$weights[3],3)
results[4,9]<-round(akaike.weights(as.numeric(all_AICc))$weights[4],3)
results[5,9]<-round(akaike.weights(as.numeric(all_AICc))$weights[5],3)
results[6,9]<-round(akaike.weights(as.numeric(all_AICc))$weights[6],3)

final_time<-results
print(results)

write.table(final_time,paste("Results_Moraceae_time",".txt"),quote=FALSE,sep="\t",row.names=FALSE)

print(final_time)

## Condamine's models: paleoenvironment-dependent diversification models with continuous rates through time
# If using this approach, please cite:
# Condamine F.L., Rolland J., Morlon H. 2013. Macroevolutionary perspectives to environmental change. Ecol. Lett. 16: 72–85.

# BCST DCST (constant Birth-death)
print("BCST DCST")
f.lamb<-function(x,y){y[1]}
f.mu<-function(x,y){y[1]}
lamb_par<-c(0.1)
mu_par<-c(0.01)
cst.lamb=T; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F

treei_BCSTDCST<-fit_bd(tree,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=508/1285,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")

print(treei_BCSTDCST)

###############################################
###### Temp Dependence (exponential variation) ######
###############################################

env_data <- read.csv("tempdata.csv") #d18o data
dof<-smooth.spline(env_data[,1], env_data[,2])$df

# 11) BEnv.Var EXPO
print("BEnv.Var EXPO")
f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){0}
lamb_par<-c(0.1, 0.0)
mu_par<-c()
cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T

treei_BEnv.Var_EXPO<-fit_env(tree,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=508/1285,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")

print(treei_BEnv.Var_EXPO)
plot_fit_env(treei_BEnv.Var_EXPO,env_data,tot_time)

# 12) BEnv.Var DCST EXPO
print("BEnv.Var DCST EXPO")
f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y[1]}
lamb_par<-c(abs(treei_BEnv.Var_EXPO$lamb_par[1]), 0.0)
mu_par<-c(0.01)
cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F

treei_BEnv.VarDCST_EXPO<-fit_env(tree,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=508/1285,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")

print(treei_BEnv.VarDCST_EXPO)
plot_fit_env(treei_BEnv.VarDCST_EXPO,env_data,tot_time)

# 13) BCST DEnv.Var EXPO
print("BCST DEnv.Var EXPO")
f.lamb<-function(t,x,y){y[1]}
f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(treei_BCSTDCST$lamb_par[1])
mu_par<-c(0.01,0.001)
cst.lamb=T; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F

treei_BCSTDEnv.Var_EXPO<-fit_env(tree,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=508/1285,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")

print(treei_BCSTDEnv.Var_EXPO)
plot_fit_env(treei_BCSTDEnv.Var_EXPO,env_data,tot_time)

# 14) BEnv.Var DEnv.Var EXPO
print("BEnv.Var DEnv.Var EXPO")
f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(abs(treei_BEnv.VarDCST_EXPO$lamb_par[1]), 0.0)
mu_par<-c(0.01, 0.0)
cst.lamb=F; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F

treei_BEnv.VarDEnv.Var_EXPO<-fit_env(tree,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=508/1285,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")

print(treei_BEnv.VarDEnv.Var_EXPO)
plot_fit_env(treei_BEnv.VarDEnv.Var_EXPO,env_data,tot_time)

##########################################
###### Temp Dependence (linear variation) ######
##########################################

#Not included in the analyses

# 15) BEnv.Var LIN
#print("BEnv.Var LIN")
#f.lamb<-function(t,x,y){y[1]+y[2]*x}
#f.mu<-function(t,x,y){0}
#lamb_par<-c(abs(treei_BEnv.Var_EXPO$lamb_par[1]), 0.0)
#mu_par<-c()
#cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T

#treei_BEnv.Var_LIN<-fit_env(tree,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=508/1285,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")

#print(treei_BEnv.Var_LIN)
#plot_fit_env(treei_BEnv.Var_LIN,env_data,tot_time)

# 16) BEnv.Var DCST LIN
#print("BEnv.Var DCST LIN")
#f.lamb<-function(t,x,y){y[1]+y[2]*x}
#f.mu<-function(t,x,y){y[1]}
#lamb_par<-c(abs(treei_BEnv.Var_LIN$lamb_par[1]), 0.0)
#mu_par<-c(0.01)
#cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F

#treei_BEnv.VarDCST_LIN<-fit_env(tree,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=508/1285,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")

#print(treei_BEnv.VarDCST_LIN)
#plot_fit_env(treei_BEnv.VarDCST_LIN,env_data,tot_time)

# 17) BCST DEnv.Var LIN
#print("BCST DEnv.Var LIN")
#f.lamb<-function(t,x,y){y[1]}
#f.mu<-function(t,x,y){y[1]+y[2]*x}
#lamb_par<-c(treei_BCSTDCST$lamb_par[1])
#mu_par<-c(0.02, 0.0)
#cst.lamb=T; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F

#treei_BCSTDEnv.Var_LIN<-fit_env(tree,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=508/1285,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")

#print(treei_BCSTDEnv.Var_LIN)
#plot_fit_env(treei_BCSTDEnv.Var_LIN,env_data,tot_time)

# 18) BEnv.Var DEnv.Var LIN
#print("BEnv.Var DEnv.Var LIN")
#f.lamb<-function(t,x,y){y[1]+y[2]*x}
#f.mu<-function(t,x,y){y[1]+y[2]*x}
#lamb_par<-c(abs(treei_BEnv.VarDCST_LIN$lamb_par[1]), 0.0)
#mu_par<-c(0.02, 0.0)
#cst.lamb=F; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F

#treei_BEnv.VarDEnv.Var_LIN<-fit_env(tree,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=508/1285,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")

#print(treei_BEnv.VarDEnv.Var_LIN)
#plot_fit_env(treei_BEnv.VarDEnv.Var_LIN,env_data,tot_time)

############# RESULTS ###########################################

results<-matrix(NA,4,9)
colnames(results)<-c("Models","NP","logL","AICc","Lambda","AlphaEnv.","Mu","BetaEnv.","Akaike_w")

#Models
results[,1]<-c("BEnv.Var_EXPO","BEnv.VarDCST_EXPO","BCSTDEnv.Var_EXPO","BEnv.VarDEnv.Var_EXPO")

#NP
results[1,2]<-2
results[2,2]<-3
results[3,2]<-3
results[4,2]<-4

#logL
results[1,3]<-treei_BEnv.Var_EXPO$LH
results[2,3]<-treei_BEnv.VarDCST_EXPO$LH
results[3,3]<-treei_BCSTDEnv.Var_EXPO$LH
results[4,3]<-treei_BEnv.VarDEnv.Var_EXPO$LH

#AICc
results[1,4]<-treei_BEnv.Var_EXPO$aicc
results[2,4]<-treei_BEnv.VarDCST_EXPO$aicc
results[3,4]<-treei_BCSTDEnv.Var_EXPO$aicc
results[4,4]<-treei_BEnv.VarDEnv.Var_EXPO$aicc

#Lambda0
results[1,5]<-abs(treei_BEnv.Var_EXPO$lamb_par[1])
results[2,5]<-abs(treei_BEnv.VarDCST_EXPO$lamb_par[1])
results[3,5]<-abs(treei_BCSTDEnv.Var_EXPO$lamb_par[1])
results[4,5]<-abs(treei_BEnv.VarDEnv.Var_EXPO$lamb_par[1])

#Alpha Env.
results[1,6]<-treei_BEnv.Var_EXPO$lamb_par[2]
results[2,6]<-treei_BEnv.VarDCST_EXPO$lamb_par[2]
results[4,6]<-treei_BEnv.VarDEnv.Var_EXPO$lamb_par[2]

#Mu0
results[2,7]<-abs(treei_BEnv.VarDCST_EXPO$mu_par[1])
results[3,7]<-abs(treei_BCSTDEnv.Var_EXPO$mu_par[1])
results[4,7]<-abs(treei_BEnv.VarDEnv.Var_EXPO$mu_par[1])

#Beta Env.
results[3,8]<-treei_BCSTDEnv.Var_EXPO$mu_par[2]
results[4,8]<-treei_BEnv.VarDEnv.Var_EXPO$mu_par[2]

#Akaike weights
all_AICc<-c(results[,4])
results[1,9]<-round(akaike.weights(as.numeric(all_AICc))$weights[1],3)
results[2,9]<-round(akaike.weights(as.numeric(all_AICc))$weights[2],3)
results[3,9]<-round(akaike.weights(as.numeric(all_AICc))$weights[3],3)
results[4,9]<-round(akaike.weights(as.numeric(all_AICc))$weights[4],3)

final_env<-results
print(results)

write.table(final_env,paste("Results_Moraceae_env",".txt"),quote=FALSE,sep="\t",row.names=FALSE)
