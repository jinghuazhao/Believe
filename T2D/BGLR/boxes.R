#Box 1: Loading BGLR and wheat data
library(BGLR)
load("wheat599.RData")
ls()  #list objects

#Boxplot for grain yield
boxplot(Yield~Env,data=Pheno,xlab="Environment",
        ylab="Grain yield (ton/ha)")


#Box 2: Fitting Bayesian Ridge Regression (BRR)

#Specify linear predictor
EtaR<-list(markers=list(X=X,model="BRR"))

#Grain yield in environment 4
y<-as.vector(subset(Pheno,Env==4)$Yield)

#Set random seed
set.seed(456)

#Fit the model
fmR<-BGLR(y=y,ETA=EtaR, nIter=10000,burnIn=5000,thin=10,verbose=FALSE)

#Estimated marker effects
betaHat<-fmR$ETA$markers$b

#Plot estimated marker effects
plot(betaHat,xlab="Marker",ylab="Estimated marker effect")
abline(h=0,col="red",lwd=2)

#Box 3: Fitting GBLUP

#A genomic relationship matrix
Z<-scale(X,center=TRUE,scale=TRUE)
G<-tcrossprod(Z)/ncol(Z)

#Specify linear predictor
EtaG<-list(markers=list(K=G,model="RKHS"))

#Grain yield in environment 4
y<-as.vector(subset(Pheno,Env==4)$Yield)

#Set random seed
set.seed(789)

#Fit the model
fmG<-BGLR(y=y,ETA=EtaG, nIter=10000,burnIn=5000,thin=10,verbose=FALSE)

#BLUPs
uHat<-fmG$ETA$markers$u

#Variance parameters
fmG$varE  					#Residual
fmG$ETA[[1]]$markers$varU	#Genotypes


#Box 4 Reaction norm model 
library(BGLR)
load("wheat599.RData")

#incidence matrix for main eff. of environments.
Pheno$Env<-as.factor(Pheno$Env)
ZE<-model.matrix(~Pheno$Env-1)  

#incidence matrix for main eff. of lines.
Pheno$Var<-as.factor(Pheno$Var)
ZVar<-model.matrix(~Pheno$Var-1)

Z<-scale(X,center=TRUE,scale=TRUE)
G<-tcrossprod(Z)/ncol(Z)

K1<-ZVar%*%G%*%t(ZVar)

ZEZE<-tcrossprod(ZE)
K2<-K1*ZEZE
  
EtaRN<-list(ENV=list(X=ZE,model="BRR"),
            Grm=list(K=K1,model="RKHS"),
            EGrm=list(K=K2,model="RKHS"))

fmRN<-BGLR(y=Pheno$Yield,ETA=EtaRN, nIter=10000,
           burnIn=5000,thin=10,verbose=FALSE)
           
plot(fmRN$y,fmRN$yHat,
     pch=as.integer(Pheno$Env),
     col=Pheno$Env,
     xlab="Grain yield observed (ton/ha)",
     ylab="Grain yield predicted (ton/ha)")

legend("topleft",legend=c(1,2,3,4),title="Environment",
       bty="n",col=1:4,pch=1:4)
       
#Variance parameters
fmRN$varE			#Residual
fmRN$ETA$ENV$varB	#Environment
fmRN$ETA$Grm$varU	#Genotypes
fmRN$ETA$EGrm$varU	#Interaction

#Box 5 Multi trait model (MTM)

#install MTM package
install.packages("remotes")             #install MTM
remotes::install_github("QuantGen/MTM") #install MTM
library(MTM)
load("wheat599.RData")

#Phenotypes
y1<-subset(Pheno,Env==1)$Yield
y2<-subset(Pheno,Env==2)$Yield
y3<-subset(Pheno,Env==3)$Yield
y4<-subset(Pheno,Env==4)$Yield

Y<-cbind(y1,y2,y3,y4)

#Genotypes
Z<-scale(X,center=TRUE,scale=TRUE)
G<-tcrossprod(Z)/ncol(Z)

#Linear predictor
EtaM<-list(
           list(K = G, COV = list(type = "UN",df0 = 4,S0 = diag(4)))
          )

#Residual
residual<-list(type = "UN",S0 = diag(4),df0 = 4)


fmM <- MTM(Y = Y, K=EtaM,resCov=residual,
		  nIter=10000,burnIn=5000,thin=10)

#Predictions of phenotypical values
fmM$YHat 

#Predictions of random effects
fmM$K[[1]]$U

#Residual covariance matrix
fmM$resCov$R 

#Genetic covariance matrix
fmM$K[[1]]$G


#Box 6 Multi trait model (MTM) with factor analytic (FA) model

library(MTM)
load("wheat599.RData")

#Phenotypes
y1<-subset(Pheno,Env==1)$Yield
y2<-subset(Pheno,Env==2)$Yield
y3<-subset(Pheno,Env==3)$Yield
y4<-subset(Pheno,Env==4)$Yield

Y<-cbind(y1,y2,y3,y4)

#Genotypes
Z<-scale(X,center=TRUE,scale=TRUE)
G<-tcrossprod(Z)/ncol(Z)

#Linear predictor using factor analytic for G strucure
EtaFA<-list(
  list(K = G, COV = list(type = 'FA', nF=1,M=matrix(nrow = 4,ncol = 1, TRUE),df0 = rep(1,4),S0 = rep(1,4),var=100)))
#Residual
residual<-list(type = "UN",S0 = diag(4),df0 = 4)


fmFA <- MTM(Y = Y, K=EtaFA,resCov=residual,
		  nIter=10000,burnIn=5000,thin=10)

#Predictions of phenotypical values
fmFA$YHat 

#Predictions of random effects
fmFA$K[[1]]$U

#Residual covariance matrix
fmFA$resCov$R 

#Genetic covariance matrix
fmFA$K[[1]]$G

