
#Y = vector of TP and FN test results, all in cases, length = 2^n
#S = vector of sensitivities
#C = vector of specificities

#have a matrix of possible test results (0000,0001, etc), ncol = n, nrows = 2^n
#loop over each row, which corresponds to Y
#for test i:
#if Y[i]=1, use S[i] and 1-C[i]
#if Y[i]=0, use 1-S[i] and C[i]
#can take ((S[i])**T[i])*(1-S[i])**(1-T[i]) #where T is the test result
#when T[i]==1, it will be S[i], when T[i]==0, will be (1-S[i])

#Sens = beta(truepos,falseneg)
#this will be the sum of Y[i] correspondign to the rows in the matrix where test i = 1
#Spec = beta(trueneg,falseneg)
library(tidyverse)

load('data_lcm.RData')

n=5

number2binary = function(number, noBits) {
  binary_vector = rev(as.numeric(intToBits(number)))
  if(missing(noBits)) {
    return(binary_vector)
  } else {
    binary_vector[-(1:(length(binary_vector) - noBits))]
  }
}

test.matrix = matrix((NA),nrow = 2^n, ncol=n)
for (i in 0:(2^n -1)){
  test.matrix[i+1,] = number2binary(i,n)
}  

## need to generate the result vector
trim.data = data.lcm 

result.vector= rep(0,dim(test.matrix)[1])
for (i in 1:dim(test.matrix)[1]){
  for (j in 1:dim(trim.data)[1]){
    if (sum(as.numeric(as.numeric(trim.data[j,])== as.numeric(test.matrix[i,])))==dim(test.matrix)[2]){
      result.vector[i] = result.vector[i] + 1
    }
  }
}


N= sum(result.vector)

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

#culture, widal, testit,alter_etiol_flip, DPP

S.alpha= c(57,
           as.numeric(estBetaParams(0.486,((0.545-0.425)/4)**2)[1]),
           as.numeric(estBetaParams(0.62,((0.675-0.56)/4)**2)[1]),
           1,
           1)
S.beta= c(38,
          as.numeric(estBetaParams(0.486,((0.545-0.425)/4)**2)[2]),
          as.numeric(estBetaParams(0.62,((0.675-0.56)/4)**2)[2]),
          1,
          1)
C.alpha=c(99999,
          as.numeric(estBetaParams(0.837,((0.784-0.881)/4)**2)[1]),
          as.numeric(estBetaParams(0.997,((1-0.986)/4)**2)[1]),
          as.numeric(estBetaParams(0.95,((1-0.9)/4)**2)[1]),
          1)
C.beta=c(1,
         as.numeric(estBetaParams(0.837,((0.784-0.881)/4)**2)[2]),
         as.numeric(estBetaParams(0.997,((1-0.986)/4)**2)[2]),
         as.numeric(estBetaParams(0.95,((1-0.9)/4)**2)[2]),
         1)
pi.alpha=1
pi.beta=1


S.vec = rep(NA,n)
C.vec= rep(NA,n)
Y.vec = rep(NA,2^n)
for (i in 1:n){
  S.vec[i]=rbeta(1,S.alpha[i],S.beta[i])
  C.vec[i]=rbeta(1,C.alpha[i],C.beta[i])
}
pi=0.5

sims= 200000

pi.vec= rep(NA,sims)
S.mat= matrix((NA),nrow=n,ncol=sims)
C.mat= matrix((NA),nrow=n,ncol=sims)


for (i in 1:sims){ #loop over simulations
  for (k in 1:(2^n)){ #loop over possible test result combinations
    tp = pi
    fn= 1-pi
    for (m in 1:n){ #loop over tests
      tp = tp*(S.vec[m]**test.matrix[k,m])*((1-S.vec[m])**(1-test.matrix[k,m]))
      fn= fn*((1-C.vec[m])**(test.matrix[k,m]))*(C.vec[m]**(1-test.matrix[k,m]))
    }
    Y.vec[k] = rbinom(1,result.vector[k],tp/(tp+fn))
  }
  pi= rbeta(1,sum(Y.vec)+pi.alpha,N-sum(Y.vec)+pi.beta)
  for (j in 1:n){
    S.vec[j]= rbeta(1,sum(test.matrix[,j]*Y.vec) + S.alpha[j],sum((1-test.matrix[,j])*Y.vec) + S.beta[j]) #multiple vector of test results by the test matrix
    C.vec[j]= rbeta(1,sum((1-test.matrix[,j])*result.vector) - sum((1-test.matrix[,j])*Y.vec) + C.alpha[j], 
                    sum(test.matrix[,j]*result.vector) - sum(test.matrix[,j]*Y.vec) + C.beta[j]) #first number is all negatives minus false negatives; second number is all positives minus true positives 
  }
  pi.vec[i]=pi
  S.mat[,i]=S.vec
  C.mat[,i]=C.vec
  }


#blood.culture,SD.IgM,Typhidot.IgM,Widal,Enterocheck,TestIt,CTK.IgG,Tubex,Spectrum.IgM
subsample = seq(sims/2 + 1,sims,100)

df = data.frame(
  pi=pi.vec,
  S.bloodculture= S.mat[1,],
  S.Widal= S.mat[2,],
  S.TestIt = S.mat[3,],
  S.AltEtiol = S.mat[4,],
  S.DPP = S.mat[5,],
  C.bloodculture= C.mat[1,],
  C.Widal= C.mat[2,],
  C.TestIt = C.mat[3,],
  C.AltEtiol = C.mat[4,],
  C.DPP = C.mat[5,],
  A.bloodculture = (S.mat[1,] + C.mat[1,])/2,
  A.Widal = (S.mat[2,] + C.mat[2,])/2,
  A.TestIt = (S.mat[3,] + C.mat[3,])/2,
  A.AltEtiol = (S.mat[4,] + C.mat[4,])/2,
  A.DPP = (S.mat[5,] + C.mat[5,])/2,
  sim= 1:sims
)

df.final = tibble(df[subsample,])

df.long = pivot_longer(df,cols=1:16)

results2=df.long %>%
  group_by(name) %>%
  summarise(med= median(value),q2.5= quantile(value,0.025),q97.5=quantile(value,0.975))

results2$param = c(rep("Balanced Accuracy",5),rep("Specificity",5),rep("Sensitivity",5),"Prevalence")
results2$name = gsub("C.","",results2$name)
results2$name = gsub("S.","",results2$name)
results2$name = gsub("A.","",results2$name)

#### Convergence diagnostics
library(coda)
geweke.diag(df[,1:16],frac1=0.5,frac2 = 0.5)


